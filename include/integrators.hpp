#pragma once

#include "gravity.hpp"
#include "types.hpp"

#include <vector>

void euler_step(
    std::vector<Particle>& particles,
    double dt,
    ForceMethod force_method
);

void velocity_verlet_step(
    std::vector<Particle>& particles,
    double dt,
    ForceMethod force_method
);
