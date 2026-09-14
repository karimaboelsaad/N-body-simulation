#pragma once

#include <string>

void run_benchmark_experiment(const std::string& output_path);

void run_accuracy_experiments(
    const std::string& integrator_output_path,
    const std::string& theta_output_path
);
