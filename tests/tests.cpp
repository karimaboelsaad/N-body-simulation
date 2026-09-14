#include "gravity.hpp"
#include "integrators.hpp"
#include "quadtree.hpp"
#include "scenarios.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message){
    if(!condition){
        throw std::runtime_error(message);
    }
}

void require_near(
    double actual,
    double expected,
    double tolerance,
    const std::string& message){
    if(std::abs(actual - expected) > tolerance){
        throw std::runtime_error(message);
    }
}

void test_square_contains_points(){
    const Square square{{0.0, 0.0}, 2.0};

    require(square.contains({0.0, 0.0}), "centre should be inside");
    require(square.contains({2.0, -2.0}), "boundary should be inside");
    require(!square.contains({2.01, 0.0}), "point beyond right edge is outside");
    require(!square.contains({0.0, -2.01}), "point below bottom edge is outside");
}

void test_direct_force_symmetry(){
    const std::vector<Particle> particles{
        Particle{1.0, {-0.5, 0.0}, {0.0, 0.0}},
        Particle{1.0, {0.5, 0.0}, {0.0, 0.0}}
    };
    const std::vector<Vec2> accelerations = direct_accelerations(particles);

    require(accelerations[0].x > 0.0, "left particle should accelerate right");
    require(accelerations[1].x < 0.0, "right particle should accelerate left");
    require_near(
        accelerations[0].x,
        -accelerations[1].x,
        1.0e-12,
        "equal masses should have opposite x accelerations"
    );
    require_near(accelerations[0].y, 0.0, 1.0e-12, "y acceleration should be zero");
    require_near(accelerations[1].y, 0.0, 1.0e-12, "y acceleration should be zero");
}

void test_quadtree_mass_and_centre_of_mass(){
    const std::vector<Particle> particles{
        Particle{2.0, {-1.0, 0.0}, {0.0, 0.0}},
        Particle{1.0, {2.0, 0.0}, {0.0, 0.0}},
        Particle{3.0, {0.0, 3.0}, {0.0, 0.0}}
    };
    const std::unique_ptr<QuadTreeNode> root = build_quadtree(particles);

    require(!root->is_leaf(), "three separated particles should subdivide the root");
    require_near(root->total_mass, 6.0, 1.0e-12, "root mass should equal total mass");
    require_near(root->centre_of_mass.x, 0.0, 1.0e-12, "incorrect centre-of-mass x");
    require_near(root->centre_of_mass.y, 1.5, 1.0e-12, "incorrect centre-of-mass y");

    for(const Particle& particle : particles){
        require(
            root->region.contains(particle.position),
            "root square should contain every particle"
        );
    }
}

void test_barnes_hut_matches_direct_with_small_theta(){
    const std::vector<Particle> particles = make_uniform_cloud(200, 42);
    const std::vector<Vec2> direct = direct_accelerations(particles);
    const std::vector<Vec2> barnes_hut =
        barnes_hut_accelerations(particles, 1.0e-6);
    const double error = relative_acceleration_error(direct, barnes_hut);

    require(error < 1.0e-12, "small theta should closely match direct forces");
}

void test_default_barnes_hut_accuracy(){
    const std::vector<Particle> particles = make_uniform_cloud(1000, 42);
    const std::vector<Vec2> direct = direct_accelerations(particles);
    const std::vector<Vec2> barnes_hut = barnes_hut_accelerations(particles);
    const double error = relative_acceleration_error(direct, barnes_hut);

    require(error < 0.01, "default Barnes-Hut error should stay below one percent");
}

double separation(const Particle& first, const Particle& second){
    const double dx = second.position.x - first.position.x;
    const double dy = second.position.y - first.position.y;
    return std::sqrt(dx * dx + dy * dy);
}

void test_velocity_verlet_preserves_two_body_orbit(){
    constexpr double dt = 0.001;
    constexpr std::size_t number_of_steps = 5000;
    const std::vector<Particle> initial_particles = make_two_body_orbit();
    std::vector<Particle> euler_particles = initial_particles;
    std::vector<Particle> verlet_particles = initial_particles;
    const double initial_energy = total_energy(initial_particles);
    const double initial_separation =
        separation(initial_particles[0], initial_particles[1]);

    for(std::size_t step = 0; step < number_of_steps; ++step){
        euler_step(euler_particles, dt, ForceMethod::direct);
        velocity_verlet_step(verlet_particles, dt, ForceMethod::direct);
    }

    const double euler_drift =
        relative_energy_drift(total_energy(euler_particles), initial_energy);
    const double verlet_drift =
        relative_energy_drift(total_energy(verlet_particles), initial_energy);
    const double verlet_separation_drift = std::abs(
        separation(verlet_particles[0], verlet_particles[1]) - initial_separation
    ) / initial_separation;

    require(euler_drift > 1.0e-3, "Euler should show measurable energy drift");
    require(verlet_drift < 1.0e-8, "Velocity Verlet energy drift is too large");
    require(verlet_drift < euler_drift, "Velocity Verlet should outperform Euler");
    require(
        verlet_separation_drift < 1.0e-4,
        "Velocity Verlet should preserve orbital separation"
    );
}

void run_test(
    const std::string& name,
    void (*test_function)(),
    int& failure_count){
    try{
        test_function();
        std::cout << "[PASS] " << name << '\n';
    } catch(const std::exception& error){
        ++failure_count;
        std::cerr << "[FAIL] " << name << ": " << error.what() << '\n';
    }
}

} // namespace

int main(){
    int failure_count = 0;

    run_test("square containment", test_square_contains_points, failure_count);
    run_test("direct force symmetry", test_direct_force_symmetry, failure_count);
    run_test(
        "quadtree mass and centre of mass",
        test_quadtree_mass_and_centre_of_mass,
        failure_count
    );
    run_test(
        "Barnes-Hut small-theta agreement",
        test_barnes_hut_matches_direct_with_small_theta,
        failure_count
    );
    run_test(
        "Barnes-Hut default accuracy",
        test_default_barnes_hut_accuracy,
        failure_count
    );
    run_test(
        "Velocity Verlet orbital stability",
        test_velocity_verlet_preserves_two_body_orbit,
        failure_count
    );

    if(failure_count != 0){
        std::cerr << failure_count << " test(s) failed\n";
        return 1;
    }

    std::cout << "All tests passed\n";
    return 0;
}
