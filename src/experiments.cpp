#include "experiments.hpp"

#include "gravity.hpp"
#include "integrators.hpp"
#include "quadtree.hpp"
#include "scenarios.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double median(std::vector<double> values){
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

double acceleration_checksum(const std::vector<Vec2>& accelerations){
    double checksum = 0.0;
    for(const Vec2& acceleration : accelerations){
        checksum += acceleration.x + acceleration.y;
    }
    return checksum;
}

double measure_direct_runtime(
    const std::vector<Particle>& particles,
    std::size_t repetitions){
    std::vector<double> durations;
    durations.reserve(repetitions);
    double checksum = 0.0;

    for(std::size_t repetition = 0; repetition < repetitions; ++repetition){
        const auto start = Clock::now();
        const std::vector<Vec2> accelerations =
            direct_accelerations(particles);
        const auto finish = Clock::now();

        checksum += acceleration_checksum(accelerations);
        durations.push_back(
            std::chrono::duration<double, std::milli>(finish - start).count()
        );
    }

    if(!std::isfinite(checksum)){
        throw std::runtime_error("Invalid direct-force benchmark result");
    }

    return median(durations);
}

double measure_barnes_hut_runtime(
    const std::vector<Particle>& particles,
    double theta,
    std::size_t repetitions){
    std::vector<double> durations;
    durations.reserve(repetitions);
    double checksum = 0.0;

    for(std::size_t repetition = 0; repetition < repetitions; ++repetition){
        const auto start = Clock::now();
        const std::vector<Vec2> accelerations =
            barnes_hut_accelerations(particles, theta);
        const auto finish = Clock::now();

        checksum += acceleration_checksum(accelerations);
        durations.push_back(
            std::chrono::duration<double, std::milli>(finish - start).count()
        );
    }

    if(!std::isfinite(checksum)){
        throw std::runtime_error("Invalid Barnes-Hut benchmark result");
    }

    return median(durations);
}

double particle_separation(
    const Particle& first,
    const Particle& second){
    const double dx = second.position.x - first.position.x;
    const double dy = second.position.y - first.position.y;
    return std::sqrt(dx * dx + dy * dy);
}

double relative_difference(double value, double reference){
    if(std::abs(reference) < 1.0e-15){
        return 0.0;
    }

    return std::abs(value - reference) / std::abs(reference);
}

void write_integrator_row(
    std::ostream& output,
    double time,
    const std::vector<Particle>& euler_particles,
    const std::vector<Particle>& verlet_particles,
    double initial_energy,
    double initial_separation){
    const double euler_energy = total_energy(euler_particles);
    const double verlet_energy = total_energy(verlet_particles);
    const double euler_separation =
        particle_separation(euler_particles[0], euler_particles[1]);
    const double verlet_separation =
        particle_separation(verlet_particles[0], verlet_particles[1]);

    output << time << ','
           << euler_energy << ','
           << relative_energy_drift(euler_energy, initial_energy) << ','
           << relative_difference(euler_separation, initial_separation) << ','
           << verlet_energy << ','
           << relative_energy_drift(verlet_energy, initial_energy) << ','
           << relative_difference(verlet_separation, initial_separation) << '\n';
}

void run_integrator_accuracy_experiment(const std::string& output_path){
    constexpr double dt = 0.001;
    constexpr std::size_t number_of_steps = 20000;
    constexpr std::size_t output_interval = 100;

    const std::vector<Particle> initial_particles = make_two_body_orbit();
    std::vector<Particle> euler_particles = initial_particles;
    std::vector<Particle> verlet_particles = initial_particles;
    const double initial_energy = total_energy(initial_particles);
    const double initial_separation =
        particle_separation(initial_particles[0], initial_particles[1]);

    std::ofstream output{output_path};
    if(!output){
        throw std::runtime_error("Could not open " + output_path);
    }

    output << std::setprecision(12);
    output << "time,euler_energy,euler_energy_drift,euler_separation_drift,"
              "verlet_energy,verlet_energy_drift,verlet_separation_drift\n";
    write_integrator_row(
        output,
        0.0,
        euler_particles,
        verlet_particles,
        initial_energy,
        initial_separation
    );

    for(std::size_t step = 0; step < number_of_steps; ++step){
        euler_step(euler_particles, dt, ForceMethod::direct);
        velocity_verlet_step(verlet_particles, dt, ForceMethod::direct);

        if((step + 1) % output_interval == 0){
            const double time = static_cast<double>(step + 1) * dt;
            write_integrator_row(
                output,
                time,
                euler_particles,
                verlet_particles,
                initial_energy,
                initial_separation
            );
        }
    }

    const double final_euler_energy = total_energy(euler_particles);
    const double final_verlet_energy = total_energy(verlet_particles);
    const double final_euler_separation =
        particle_separation(euler_particles[0], euler_particles[1]);
    const double final_verlet_separation =
        particle_separation(verlet_particles[0], verlet_particles[1]);
    std::cerr << "Integrator accuracy results written to " << output_path << '\n';
    std::cerr << "Final Euler energy drift: "
              << relative_energy_drift(final_euler_energy, initial_energy) << '\n';
    std::cerr << "Final Verlet energy drift: "
              << relative_energy_drift(final_verlet_energy, initial_energy) << '\n';
    std::cerr << "Final Euler separation drift: "
              << relative_difference(final_euler_separation, initial_separation)
              << '\n';
    std::cerr << "Final Verlet separation drift: "
              << relative_difference(final_verlet_separation, initial_separation)
              << '\n';
}

void run_theta_accuracy_experiment(const std::string& output_path){
    constexpr std::size_t particle_count = 5000;
    constexpr std::size_t repetitions = 5;
    const std::vector<double> theta_values{0.2, 0.35, 0.5, 0.75, 1.0};
    const std::vector<Particle> particles =
        make_uniform_cloud(particle_count, 42);

    const std::vector<Vec2> exact = direct_accelerations(particles);
    const double direct_ms =
        measure_direct_runtime(particles, repetitions);

    std::ofstream output{output_path};
    if(!output){
        throw std::runtime_error("Could not open " + output_path);
    }

    output << std::setprecision(12);
    output << "theta,relative_acceleration_error,barnes_hut_ms,direct_ms,speedup\n";

    for(const double theta : theta_values){
        const std::vector<Vec2> approximate =
            barnes_hut_accelerations(particles, theta);
        const double barnes_hut_ms =
            measure_barnes_hut_runtime(particles, theta, repetitions);
        const double error =
            relative_acceleration_error(exact, approximate);

        output << theta << ','
               << error << ','
               << barnes_hut_ms << ','
               << direct_ms << ','
               << direct_ms / barnes_hut_ms << '\n';

        std::cerr << "theta=" << theta
                  << ": error=" << error
                  << ", Barnes-Hut=" << barnes_hut_ms << " ms\n";
    }

    std::cerr << "Theta accuracy results written to " << output_path << '\n';
}

} // namespace

void run_benchmark_experiment(const std::string& output_path){
    constexpr std::size_t repetitions = 7;
    const std::vector<std::size_t> particle_counts{
        100,
        250,
        500,
        1000,
        2000,
        5000,
        10000
    };

    std::ofstream output{output_path};
    if(!output){
        throw std::runtime_error("Could not open " + output_path);
    }

    output << std::setprecision(12);
    output << "particle_count,direct_ms,barnes_hut_ms,speedup,relative_error\n";

    for(const std::size_t particle_count : particle_counts){
        const std::vector<Particle> particles =
            make_uniform_cloud(particle_count, 42);

        const std::vector<Vec2> exact = direct_accelerations(particles);
        const std::vector<Vec2> approximate =
            barnes_hut_accelerations(particles);
        const double direct_ms =
            measure_direct_runtime(particles, repetitions);
        const double barnes_hut_ms =
            measure_barnes_hut_runtime(
                particles,
                barnes_hut_theta,
                repetitions
            );
        const double speedup = direct_ms / barnes_hut_ms;
        const double error =
            relative_acceleration_error(exact, approximate);

        output << particle_count << ','
               << direct_ms << ','
               << barnes_hut_ms << ','
               << speedup << ','
               << error << '\n';

        std::cerr << particle_count << " particles: direct="
                  << direct_ms << " ms, Barnes-Hut="
                  << barnes_hut_ms << " ms, speedup="
                  << speedup << "x, error=" << error << '\n';
    }

    std::cerr << "Benchmark results written to " << output_path << '\n';
}

void run_accuracy_experiments(
    const std::string& integrator_output_path,
    const std::string& theta_output_path){
    run_integrator_accuracy_experiment(integrator_output_path);
    run_theta_accuracy_experiment(theta_output_path);
}
