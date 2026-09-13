#pragma once

#include "types.hpp"

#include <vector>

inline constexpr double G = 1.0;
inline constexpr double softening = 1.0e-3;

enum class ForceMethod {
    direct,
    barnes_hut
};

Vec2 gravitational_acceleration(
    const Particle& target,
    const Particle& source
);

std::vector<Vec2> direct_accelerations(
    const std::vector<Particle>& particles
);

std::vector<Vec2> calculate_accelerations(
    const std::vector<Particle>& particles,
    ForceMethod method
);

double relative_acceleration_error(
    const std::vector<Vec2>& exact,
    const std::vector<Vec2>& approximate
);

double kinetic_energy(const std::vector<Particle>& particles);
double potential_energy(const std::vector<Particle>& particles);
double total_energy(const std::vector<Particle>& particles);
double relative_energy_drift(double energy, double initial_energy);
