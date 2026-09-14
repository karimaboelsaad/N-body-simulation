#pragma once

#include "types.hpp"

#include <cstddef>
#include <vector>

std::vector<Particle> make_two_body_orbit();
std::vector<Particle> make_figure_eight();
std::vector<Particle> make_random_system(
    std::size_t particle_count,
    unsigned int seed
);

std::vector<Particle> make_uniform_cloud(
    std::size_t particle_count,
    unsigned int seed
);
