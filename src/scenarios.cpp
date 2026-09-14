#include "scenarios.hpp"

#include "gravity.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

std::vector<Particle> make_two_body_orbit(){
    const double orbital_speed = std::sqrt(0.5);

    return {
        Particle{1.0, {-0.5, 0.0}, {0.0, -orbital_speed}},
        Particle{1.0, {0.5, 0.0}, {0.0, orbital_speed}}
    };
}

std::vector<Particle> make_figure_eight(){
    return {
        Particle{1.0, {-0.97000436, 0.24308753},
                      {0.466203685, 0.432365730}},
        Particle{1.0, {0.97000436, -0.24308753},
                      {0.466203685, 0.432365730}},
        Particle{1.0, {0.0, 0.0},
                      {-0.932407370, -0.864731460}}
    };
}

std::vector<Particle> make_random_system(
    std::size_t particle_count,
    unsigned int seed){
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

    const double total_orbiting_mass = 0.001;
    const double orbiting_particle_mass =
        total_orbiting_mass / static_cast<double>(particle_count - 1);
    Vec2 orbiting_momentum{0.0, 0.0};

    for(std::size_t index = 1; index < particle_count; ++index){
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

std::vector<Particle> make_uniform_cloud(
    std::size_t particle_count,
    unsigned int seed){
    if(particle_count == 0){
        throw std::invalid_argument("particle_count must be greater than zero");
    }

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> position_distribution(-1.0, 1.0);
    std::vector<Particle> particles;
    particles.reserve(particle_count);

    const double particle_mass =
        1.0 / static_cast<double>(particle_count);

    for(std::size_t index = 0; index < particle_count; ++index){
        particles.push_back(Particle{
            particle_mass,
            Vec2{
                position_distribution(generator),
                position_distribution(generator)
            },
            Vec2{0.0, 0.0}
        });
    }

    return particles;
}
