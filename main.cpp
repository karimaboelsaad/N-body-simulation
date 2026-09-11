#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

constexpr double G = 1.0;
constexpr double softening = 1.0e-3;


struct Vec2{
    double x;
    double y;
};

struct Particle{
    double mass;
    Vec2 position;
    Vec2 velocity;
};


Vec2 gravitational_acceleration(const Particle& target, const Particle& source){
    Vec2 acceleration{};
    double dx = source.position.x - target.position.x;
    double dy = source.position.y - target.position.y;
    double distance_squared =dx * dx + dy * dy + softening * softening;
    double distance = std::sqrt(distance_squared);
    double scale = G*source.mass/(distance_squared*distance);
    acceleration.x = dx*scale;
    acceleration.y = dy*scale;

    return acceleration;
}

std::vector<Vec2> direct_accelerations(const std::vector<Particle>& particles){
    std::vector<Vec2> accelerations(
        particles.size(),
        Vec2{0.0, 0.0}
    );

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        for (std::size_t j = 0; j < particles.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }

            Vec2 contribution =
                gravitational_acceleration(particles[i], particles[j]);

            accelerations[i].x += contribution.x;
            accelerations[i].y += contribution.y;
        }
    }

    return accelerations;
}



void euler_step(std::vector<Particle>& particles, double dt){
    const std::vector<Vec2> accelerations =
        direct_accelerations(particles);

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        particles[i].position.x+=particles[i].velocity.x*dt;
        particles[i].position.y+=particles[i].velocity.y*dt;
        particles[i].velocity.x+=accelerations[i].x*dt;
        particles[i].velocity.y+=accelerations[i].y*dt;
    }

}

void velocity_verlet_step(std::vector<Particle>& particles, double dt){
    const std::vector<Vec2> old_accelerations =
        direct_accelerations(particles);
    const double half_dt_squared = 0.5 * dt * dt;

    for(std::size_t i = 0; i < particles.size(); ++i){
        particles[i].position.x +=
            particles[i].velocity.x * dt +
            old_accelerations[i].x * half_dt_squared;
        particles[i].position.y +=
            particles[i].velocity.y * dt +
            old_accelerations[i].y * half_dt_squared;
    }

    const std::vector<Vec2> new_accelerations =
        direct_accelerations(particles);

    for(std::size_t i = 0; i < particles.size(); ++i){
        particles[i].velocity.x +=
            0.5 * (old_accelerations[i].x + new_accelerations[i].x) * dt;
        particles[i].velocity.y +=
            0.5 * (old_accelerations[i].y + new_accelerations[i].y) * dt;
    }
}

double kinetic_energy(const std::vector<Particle>& particles){
    double energy{0};
    for(const Particle& particle:particles){
        energy+=0.5*particle.mass*(particle.velocity.x*particle.velocity.x+particle.velocity.y*particle.velocity.y);
    }

    return energy;
}

double potential_energy(const std::vector<Particle>& particles){
    double energy{0};
    for(std::size_t i = 0; i < particles.size(); ++i){
        for(std::size_t j = i + 1; j < particles.size(); ++j){
            double dx = particles[j].position.x - particles[i].position.x;
            double dy = particles[j].position.y - particles[i].position.y;
            double distance_squared = dx * dx + dy * dy + softening * softening;
            energy -= G * particles[i].mass * particles[j].mass /
                      std::sqrt(distance_squared);
        }
    }

    return energy;
}

double total_energy(const std::vector<Particle>& particles){
    return kinetic_energy(particles)+potential_energy(particles);
}

std::vector<Particle> make_two_body_orbit(){
    const double orbital_speed = std::sqrt(0.5);

    return {
        Particle{1.0, {-0.5, 0.0}, {0.0, -orbital_speed}},
        Particle{1.0, {0.5, 0.0}, {0.0, orbital_speed}}
    };
}

std::vector<Particle> make_figure_eight(){
    
    return {
        Particle{1.0,{-0.97000436, 0.24308753},{0.466203685, 0.432365730}},
        Particle{1.0,{0.97000436, -0.24308753},{0.466203685, 0.432365730}},
        Particle{1.0,{0.0, 0.0},{-0.932407370, -0.864731460}}
    };
}

std::vector<Particle> make_random_system(std::size_t particle_count,unsigned int seed){
    if(particle_count == 0){
        throw std::invalid_argument("particle_count must be greater than zero");
    }

    std::vector<Particle> particles;
    particles.reserve(particle_count);

    const double central_mass = 1.0;
    particles.push_back(Particle{central_mass, {0.0, 0.0}, {0.0, 0.0}});

    if(particle_count == 1){
        return particles;
    }

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> angle_distribution(
        0.0,
        2.0 * std::acos(-1.0)
    );
    std::uniform_real_distribution<double> radius_squared_distribution(
        0.09,
        1.0
    );

    // The orbiters are deliberately much lighter than the central body.
    // They still attract one another, but close encounters cannot easily
    // fling them out of the system during this visual demo.
    const double total_orbiting_mass = 0.001;
    const double orbiting_particle_mass =
        total_orbiting_mass / static_cast<double>(particle_count - 1);
    Vec2 orbiting_momentum{0.0, 0.0};

    for(std::size_t i = 1; i < particle_count; ++i){
        const double angle = angle_distribution(generator);
        const double radius = std::sqrt(radius_squared_distribution(generator));
        const double orbital_speed = std::sqrt(G * central_mass / radius);

        const Vec2 position{
            radius * std::cos(angle),
            radius * std::sin(angle)
        };
        const Vec2 velocity{
            -orbital_speed * std::sin(angle),
            orbital_speed * std::cos(angle)
        };

        particles.push_back(Particle{
            orbiting_particle_mass,
            position,
            velocity
        });

        orbiting_momentum.x += orbiting_particle_mass * velocity.x;
        orbiting_momentum.y += orbiting_particle_mass * velocity.y;
    }

    // Give the central body the opposite momentum so the whole system's
    // centre of mass does not drift across the screen.
    particles[0].velocity.x = -orbiting_momentum.x / central_mass;
    particles[0].velocity.y = -orbiting_momentum.y / central_mass;

    Vec2 centre_of_mass{0.0, 0.0};
    const double total_mass = central_mass + total_orbiting_mass;
    for(const Particle& particle : particles){
        centre_of_mass.x += particle.mass * particle.position.x;
        centre_of_mass.y += particle.mass * particle.position.y;
    }

    centre_of_mass.x /= total_mass;
    centre_of_mass.y /= total_mass;
    for(Particle& particle : particles){
        particle.position.x -= centre_of_mass.x;
        particle.position.y -= centre_of_mass.y;
    }

    return particles;
}

double relative_energy_drift(double energy, double initial_energy){
    if(std::abs(initial_energy) < 1.0e-15){
        return 0.0;
    }

    return std::abs(energy - initial_energy) / std::abs(initial_energy);
}

void write_state(
    std::ostream& output,
    double time,
    const std::vector<Particle>& particles){
    for(std::size_t id = 0; id < particles.size(); ++id){
        output << time << ','
               << id << ','
               << particles[id].position.x << ','
               << particles[id].position.y << '\n';
    }
}

int main(int argc, char* argv[]){
    const std::string scenario = argc > 1 ? argv[1] : "random";
    std::size_t particle_count = 0;
    unsigned int seed = std::random_device{}();

    if(scenario == "random"){
        if(argc > 2){
            std::cerr << "Usage: " << argv[0]
                      << " [random | figure-eight]\n";
            return 1;
        }

        std::cerr << "Choose the number of bodies (1-100): " << std::flush;
        if(!(std::cin >> particle_count) ||
           particle_count < 1 || particle_count > 100){
            std::cerr << "\nPlease enter a whole number between 1 and 100.\n";
            return 1;
        }
    } else if(scenario != "figure-eight"){
        std::cerr << "Usage: " << argv[0]
                  << " [random | figure-eight]\n";
        return 1;
    } else if(argc > 2){
        std::cerr << "Usage: " << argv[0]
                  << " [random | figure-eight]\n";
        return 1;
    }

    const double dt = 0.001;
    const std::size_t number_of_steps =
        scenario == "random" ? 10000 : 20000;

    const std::vector<Particle> initial_particles =
        scenario == "random"
            ? make_random_system(particle_count, seed)
            : make_figure_eight();

    if(scenario == "random"){
        std::cerr << "Random system: " << particle_count
                  << " bodies (seed " << seed << ")\n";
    }

    std::vector<Particle> euler_particles = initial_particles;
    std::vector<Particle> verlet_particles = initial_particles;
    const double initial_energy = total_energy(initial_particles);

    std::ofstream trajectory_file{"trajectory.csv"};
    if(!trajectory_file){
        std::cerr << "Could not open trajectory.csv\n";
        return 1;
    }

    trajectory_file << "time,particle_id,x,y\n";
    write_state(trajectory_file, 0.0, verlet_particles);

    std::cout << "time,euler_energy,euler_drift,verlet_energy,verlet_drift\n";
    std::cout << 0.0 << ','
              << initial_energy << ','
              << 0.0 << ','
              << initial_energy << ','
              << 0.0 << '\n';

    for (std::size_t step = 0; step < number_of_steps; ++step)
    {
        euler_step(euler_particles, dt);
        velocity_verlet_step(verlet_particles, dt);

        if((step + 1) % 10 == 0){
            const double trajectory_time =
                static_cast<double>(step + 1) * dt;
            write_state(trajectory_file, trajectory_time, verlet_particles);
        }

        if ((step + 1) % 100 == 0)
        {
            const double time = static_cast<double>(step + 1) * dt;
            const double euler_energy = total_energy(euler_particles);
            const double verlet_energy = total_energy(verlet_particles);
            const double euler_drift =
                relative_energy_drift(euler_energy, initial_energy);
            const double verlet_drift =
                relative_energy_drift(verlet_energy, initial_energy);

            std::cout << time << ','
                      << euler_energy << ','
                      << euler_drift << ','
                      << verlet_energy << ','
                      << verlet_drift << '\n';
        }
    }

    return 0;
}
