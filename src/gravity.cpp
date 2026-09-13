#include "gravity.hpp"

#include "quadtree.hpp"

#include <cmath>

Vec2 gravitational_acceleration(
    const Particle& target,
    const Particle& source){
    const double dx = source.position.x - target.position.x;
    const double dy = source.position.y - target.position.y;
    const double distance_squared =
        dx * dx + dy * dy + softening * softening;
    const double distance = std::sqrt(distance_squared);
    const double scale = G * source.mass / (distance_squared * distance);

    return Vec2{dx * scale, dy * scale};
}

std::vector<Vec2> direct_accelerations(
    const std::vector<Particle>& particles){
    std::vector<Vec2> accelerations(
        particles.size(),
        Vec2{0.0, 0.0}
    );

    for(std::size_t i = 0; i < particles.size(); ++i){
        for(std::size_t j = 0; j < particles.size(); ++j){
            if(i == j){
                continue;
            }

            const Vec2 contribution =
                gravitational_acceleration(particles[i], particles[j]);
            accelerations[i].x += contribution.x;
            accelerations[i].y += contribution.y;
        }
    }

    return accelerations;
}

std::vector<Vec2> calculate_accelerations(
    const std::vector<Particle>& particles,
    ForceMethod method){
    if(method == ForceMethod::direct){
        return direct_accelerations(particles);
    }

    return barnes_hut_accelerations(particles);
}

double relative_acceleration_error(
    const std::vector<Vec2>& exact,
    const std::vector<Vec2>& approximate){
    double squared_error = 0.0;
    double squared_exact = 0.0;

    for(std::size_t index = 0; index < exact.size(); ++index){
        const double dx = approximate[index].x - exact[index].x;
        const double dy = approximate[index].y - exact[index].y;
        squared_error += dx * dx + dy * dy;
        squared_exact +=
            exact[index].x * exact[index].x +
            exact[index].y * exact[index].y;
    }

    if(squared_exact == 0.0){
        return 0.0;
    }

    return std::sqrt(squared_error / squared_exact);
}

double kinetic_energy(const std::vector<Particle>& particles){
    double energy = 0.0;

    for(const Particle& particle : particles){
        const double speed_squared =
            particle.velocity.x * particle.velocity.x +
            particle.velocity.y * particle.velocity.y;
        energy += 0.5 * particle.mass * speed_squared;
    }

    return energy;
}

double potential_energy(const std::vector<Particle>& particles){
    double energy = 0.0;

    for(std::size_t i = 0; i < particles.size(); ++i){
        for(std::size_t j = i + 1; j < particles.size(); ++j){
            const double dx =
                particles[j].position.x - particles[i].position.x;
            const double dy =
                particles[j].position.y - particles[i].position.y;
            const double distance_squared =
                dx * dx + dy * dy + softening * softening;

            energy -= G * particles[i].mass * particles[j].mass /
                      std::sqrt(distance_squared);
        }
    }

    return energy;
}

double total_energy(const std::vector<Particle>& particles){
    return kinetic_energy(particles) + potential_energy(particles);
}

double relative_energy_drift(double energy, double initial_energy){
    if(std::abs(initial_energy) < 1.0e-15){
        return 0.0;
    }

    return std::abs(energy - initial_energy) / std::abs(initial_energy);
}
