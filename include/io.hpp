#pragma once

#include "types.hpp"

#include <ostream>
#include <vector>

void write_state(
    std::ostream& output,
    double time,
    const std::vector<Particle>& particles
);
