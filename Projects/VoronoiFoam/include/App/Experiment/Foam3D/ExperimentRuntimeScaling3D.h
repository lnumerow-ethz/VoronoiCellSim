#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"

#include <chrono>

class ExperimentRuntimeScaling3D : public Experiment {
    int frame = 0;

    std::chrono::high_resolution_clock::time_point start_time;
    F num_sites_float = 10.0;

   public:
    [[nodiscard]] std::string getName() const override { return "Runtime Scaling"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
