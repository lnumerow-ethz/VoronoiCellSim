#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"
#include <fstream>

#include <random>
#include <chrono>

class ExperimentRigidBody2D : public Experiment {
    std::mt19937 rng = std::mt19937(time(NULL));

    bool pre_optimized_0 = false;

    int frame = 0;
    int num_frame_iters = 0;
    std::chrono::high_resolution_clock::time_point iter_start_time;

   public:
    [[nodiscard]] std::string getName() const override { return "Rigid Body"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
