#include "quadtree.hpp"

#include "gravity.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

bool Square::contains(const Vec2& point) const {
    const double left = centre.x - half_size;
    const double right = centre.x + half_size;
    const double bottom = centre.y - half_size;
    const double top = centre.y + half_size;

    return point.x >= left &&
           point.x <= right &&
           point.y >= bottom &&
           point.y <= top;
}

QuadTreeNode::QuadTreeNode(const Square& node_region)
    : region(node_region) {
}

bool QuadTreeNode::is_leaf() const {
    return children[0] == nullptr;
}

std::size_t QuadTreeNode::child_index_for(const Vec2& point) const {
    const bool east = point.x >= region.centre.x;
    const bool north = point.y >= region.centre.y;

    if(north){
        return east ? 1 : 0;
    }

    return east ? 3 : 2;
}

void QuadTreeNode::subdivide() {
    const double child_half_size = region.half_size / 2.0;
    const double x = region.centre.x;
    const double y = region.centre.y;

    children[0] = std::make_unique<QuadTreeNode>(
        Square{{x - child_half_size, y + child_half_size}, child_half_size}
    );
    children[1] = std::make_unique<QuadTreeNode>(
        Square{{x + child_half_size, y + child_half_size}, child_half_size}
    );
    children[2] = std::make_unique<QuadTreeNode>(
        Square{{x - child_half_size, y - child_half_size}, child_half_size}
    );
    children[3] = std::make_unique<QuadTreeNode>(
        Square{{x + child_half_size, y - child_half_size}, child_half_size}
    );
}

void QuadTreeNode::add_to_mass(const Particle& particle) {
    const double new_total_mass = total_mass + particle.mass;

    centre_of_mass.x =
        (centre_of_mass.x * total_mass +
         particle.position.x * particle.mass) /
        new_total_mass;
    centre_of_mass.y =
        (centre_of_mass.y * total_mass +
         particle.position.y * particle.mass) /
        new_total_mass;
    total_mass = new_total_mass;
}

bool QuadTreeNode::insert(
    std::size_t new_particle_index,
    const std::vector<Particle>& particles){
    const Particle& new_particle = particles[new_particle_index];

    if(!region.contains(new_particle.position)){
        return false;
    }

    add_to_mass(new_particle);

    if(is_leaf() && particle_index == -1){
        particle_index = static_cast<int>(new_particle_index);
        return true;
    }

    if(is_leaf()){
        if(region.half_size < 1.0e-12){
            throw std::runtime_error(
                "Particles are too close to separate in the quadtree"
            );
        }

        const std::size_t old_particle_index =
            static_cast<std::size_t>(particle_index);
        particle_index = -1;
        subdivide();

        const std::size_t old_child_index =
            child_index_for(particles[old_particle_index].position);
        children[old_child_index]->insert(old_particle_index, particles);
    }

    const std::size_t new_child_index =
        child_index_for(new_particle.position);
    return children[new_child_index]->insert(new_particle_index, particles);
}

Square bounding_square(const std::vector<Particle>& particles){
    if(particles.empty()){
        throw std::invalid_argument("Cannot build a quadtree with no particles");
    }

    double minimum_x = particles[0].position.x;
    double maximum_x = particles[0].position.x;
    double minimum_y = particles[0].position.y;
    double maximum_y = particles[0].position.y;

    for(const Particle& particle : particles){
        minimum_x = std::min(minimum_x, particle.position.x);
        maximum_x = std::max(maximum_x, particle.position.x);
        minimum_y = std::min(minimum_y, particle.position.y);
        maximum_y = std::max(maximum_y, particle.position.y);
    }

    const Vec2 centre{
        0.5 * (minimum_x + maximum_x),
        0.5 * (minimum_y + maximum_y)
    };
    const double width = std::max(
        maximum_x - minimum_x,
        maximum_y - minimum_y
    );
    const double half_size = 0.5 * width + softening;

    return Square{centre, half_size};
}

std::unique_ptr<QuadTreeNode> build_quadtree(
    const std::vector<Particle>& particles){
    auto root = std::make_unique<QuadTreeNode>(bounding_square(particles));

    for(std::size_t index = 0; index < particles.size(); ++index){
        if(!root->insert(index, particles)){
            throw std::runtime_error("A particle fell outside the quadtree");
        }
    }

    return root;
}

Vec2 acceleration_from_node(
    std::size_t target_index,
    const std::vector<Particle>& particles,
    const QuadTreeNode& node,
    double theta){
    if(node.total_mass == 0.0){
        return Vec2{0.0, 0.0};
    }

    if(node.is_leaf()){
        if(node.particle_index == -1 ||
           node.particle_index == static_cast<int>(target_index)){
            return Vec2{0.0, 0.0};
        }

        return gravitational_acceleration(
            particles[target_index],
            particles[static_cast<std::size_t>(node.particle_index)]
        );
    }

    const Particle& target = particles[target_index];
    const double dx = node.centre_of_mass.x - target.position.x;
    const double dy = node.centre_of_mass.y - target.position.y;
    const double distance = std::sqrt(dx * dx + dy * dy);
    const double node_width = 2.0 * node.region.half_size;
    const bool contains_target = node.region.contains(target.position);

    if(!contains_target && distance > 0.0 && node_width / distance < theta){
        const Particle combined_particle{
            node.total_mass,
            node.centre_of_mass,
            Vec2{0.0, 0.0}
        };
        return gravitational_acceleration(target, combined_particle);
    }

    Vec2 acceleration{0.0, 0.0};
    for(const auto& child : node.children){
        const Vec2 contribution =
            acceleration_from_node(target_index, particles, *child, theta);
        acceleration.x += contribution.x;
        acceleration.y += contribution.y;
    }

    return acceleration;
}

std::vector<Vec2> barnes_hut_accelerations(
    const std::vector<Particle>& particles,
    double theta){
    if(theta <= 0.0){
        throw std::invalid_argument("Barnes-Hut theta must be positive");
    }

    if(particles.empty()){
        return {};
    }

    const std::unique_ptr<QuadTreeNode> root = build_quadtree(particles);
    std::vector<Vec2> accelerations(
        particles.size(),
        Vec2{0.0, 0.0}
    );

    for(std::size_t index = 0; index < particles.size(); ++index){
        accelerations[index] =
            acceleration_from_node(index, particles, *root, theta);
    }

    return accelerations;
}
