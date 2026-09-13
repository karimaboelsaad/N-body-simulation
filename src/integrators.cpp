#include "integrators.hpp"

void euler_step(
    std::vector<Particle>& particles,
    double dt,
    ForceMethod force_method){
    const std::vector<Vec2> accelerations =
        calculate_accelerations(particles, force_method);

    for(std::size_t index = 0; index < particles.size(); ++index){
        particles[index].position.x += particles[index].velocity.x * dt;
        particles[index].position.y += particles[index].velocity.y * dt;
        particles[index].velocity.x += accelerations[index].x * dt;
        particles[index].velocity.y += accelerations[index].y * dt;
    }
}

void velocity_verlet_step(
    std::vector<Particle>& particles,
    double dt,
    ForceMethod force_method){
    const std::vector<Vec2> old_accelerations =
        calculate_accelerations(particles, force_method);
    const double half_dt_squared = 0.5 * dt * dt;

    for(std::size_t index = 0; index < particles.size(); ++index){
        particles[index].position.x +=
            particles[index].velocity.x * dt +
            old_accelerations[index].x * half_dt_squared;
        particles[index].position.y +=
            particles[index].velocity.y * dt +
            old_accelerations[index].y * half_dt_squared;
    }

    const std::vector<Vec2> new_accelerations =
        calculate_accelerations(particles, force_method);

    for(std::size_t index = 0; index < particles.size(); ++index){
        particles[index].velocity.x +=
            0.5 * (old_accelerations[index].x +
                   new_accelerations[index].x) * dt;
        particles[index].velocity.y +=
            0.5 * (old_accelerations[index].y +
                   new_accelerations[index].y) * dt;
    }
}
