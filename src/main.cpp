#include "experiments.hpp"
#include "gravity.hpp"
#include "integrators.hpp"
#include "io.hpp"
#include "quadtree.hpp"
#include "scenarios.hpp"

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char* argv[]){
    const std::string scenario = argc > 1 ? argv[1] : "random";
    std::size_t particle_count = 0;
    const unsigned int seed = std::random_device{}();

    if(scenario == "benchmark"){
        if(argc > 2){
            std::cerr << "Usage: " << argv[0]
                      << " [random | figure-eight | benchmark | accuracy]\n";
            return 1;
        }

        run_benchmark_experiment("benchmark.csv");
        return 0;
    }

    if(scenario == "accuracy"){
        if(argc > 2){
            std::cerr << "Usage: " << argv[0]
                      << " [random | figure-eight | benchmark | accuracy]\n";
            return 1;
        }

        run_accuracy_experiments(
            "integrator_accuracy.csv",
            "theta_accuracy.csv"
        );
        return 0;
    }

    if(scenario == "random"){
        if(argc > 2){
            std::cerr << "Usage: " << argv[0]
                      << " [random | figure-eight | benchmark | accuracy]\n";
            return 1;
        }

        std::cerr << "Choose the number of bodies (1-100): " << std::flush;
        if(!(std::cin >> particle_count) ||
           particle_count < 1 || particle_count > 100){
            std::cerr << "\nPlease enter a whole number between 1 and 100.\n";
            return 1;
        }
    } else if(scenario != "figure-eight"){
        std::cerr << "Usage: " << argv[0]
                  << " [random | figure-eight | benchmark | accuracy]\n";
        return 1;
    } else if(argc > 2){
        std::cerr << "Usage: " << argv[0]
                  << " [random | figure-eight | benchmark | accuracy]\n";
        return 1;
    }

    const double dt = 0.001;
    const std::size_t number_of_steps =
        scenario == "random" ? 10000 : 20000;

    const std::vector<Particle> initial_particles =
        scenario == "random"
            ? make_random_system(particle_count, seed)
            : make_figure_eight();

    if(scenario == "random"){
        std::cerr << "Random system: " << particle_count
                  << " bodies (seed " << seed << ")\n";
    }

    const ForceMethod force_method =
        scenario == "random" ? ForceMethod::barnes_hut : ForceMethod::direct;

    if(force_method == ForceMethod::barnes_hut){
        const std::vector<Vec2> exact = direct_accelerations(initial_particles);
        const std::vector<Vec2> approximate =
            barnes_hut_accelerations(initial_particles);
        std::cerr << "Force method: Barnes-Hut (theta = "
                  << barnes_hut_theta << ")\n";
        std::cerr << "Initial relative acceleration error: "
                  << relative_acceleration_error(exact, approximate) << '\n';
    } else {
        std::cerr << "Force method: direct O(n^2)\n";
    }

    std::vector<Particle> euler_particles = initial_particles;
    std::vector<Particle> verlet_particles = initial_particles;
    const double initial_energy = total_energy(initial_particles);

    std::ofstream trajectory_file{"trajectory.csv"};
    if(!trajectory_file){
        std::cerr << "Could not open trajectory.csv\n";
        return 1;
    }

    trajectory_file << "time,particle_id,x,y\n";
    write_state(trajectory_file, 0.0, verlet_particles);

    std::cout << "time,euler_energy,euler_drift,verlet_energy,verlet_drift\n";
    std::cout << 0.0 << ','
              << initial_energy << ','
              << 0.0 << ','
              << initial_energy << ','
              << 0.0 << '\n';

    for(std::size_t step = 0; step < number_of_steps; ++step){
        euler_step(euler_particles, dt, force_method);
        velocity_verlet_step(verlet_particles, dt, force_method);

        if((step + 1) % 10 == 0){
            const double trajectory_time =
                static_cast<double>(step + 1) * dt;
            write_state(trajectory_file, trajectory_time, verlet_particles);
        }

        if((step + 1) % 100 == 0){
            const double time = static_cast<double>(step + 1) * dt;
            const double euler_energy = total_energy(euler_particles);
            const double verlet_energy = total_energy(verlet_particles);
            const double euler_drift =
                relative_energy_drift(euler_energy, initial_energy);
            const double verlet_drift =
                relative_energy_drift(verlet_energy, initial_energy);

            std::cout << time << ','
                      << euler_energy << ','
                      << euler_drift << ','
                      << verlet_energy << ','
                      << verlet_drift << '\n';
        }
    }

    return 0;
}
