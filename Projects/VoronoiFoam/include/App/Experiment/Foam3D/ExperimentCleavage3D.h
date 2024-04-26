#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"

#include <random>
#include <chrono>

class ExperimentCleavage3D : public Experiment {
    std::mt19937 rng = std::mt19937(time(NULL));

    int frame = 0;
    int num_frame_iters = 0;
    std::chrono::high_resolution_clock::time_point iter_start_time;

    std::map<int, F> growing_cells;

   public:
    [[nodiscard]] std::string getName() const override { return "Cleavage"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
