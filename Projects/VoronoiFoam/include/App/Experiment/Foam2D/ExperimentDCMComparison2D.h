#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"

#include <chrono>

class ExperimentDCMComparison2D : public Experiment {
    int stage = 0;

    int frame = 0;
    int num_frame_iters = 0;
    std::chrono::high_resolution_clock::time_point iter_start_time;

   public:
    [[nodiscard]] std::string getName() const override { return "DCM Comparison"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
