#pragma once

#include "types.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

inline constexpr double barnes_hut_theta = 0.5;

struct Square {
    Vec2 centre;
    double half_size;

    bool contains(const Vec2& point) const;
};

struct QuadTreeNode {
    Square region;
    double total_mass{0.0};
    Vec2 centre_of_mass{0.0, 0.0};
    int particle_index{-1};
    std::array<std::unique_ptr<QuadTreeNode>, 4> children{};

    explicit QuadTreeNode(const Square& node_region);

    bool is_leaf() const;
    std::size_t child_index_for(const Vec2& point) const;
    void subdivide();
    void add_to_mass(const Particle& particle);
    bool insert(
        std::size_t new_particle_index,
        const std::vector<Particle>& particles
    );
};

Square bounding_square(const std::vector<Particle>& particles);

std::unique_ptr<QuadTreeNode> build_quadtree(
    const std::vector<Particle>& particles
);

Vec2 acceleration_from_node(
    std::size_t target_index,
    const std::vector<Particle>& particles,
    const QuadTreeNode& node,
    double theta
);

std::vector<Vec2> barnes_hut_accelerations(
    const std::vector<Particle>& particles,
    double theta = barnes_hut_theta
);
