#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"

#include <chrono>

class ExperimentCoarsening3D : public Experiment {
    int num_cells = 2000;

    int frame = 0;
    int num_frame_iters = 0;
    std::chrono::high_resolution_clock::time_point iter_start_time;

   private:
    void saveFrame(FoamSubApp* foam_sub_app);

   public:
    [[nodiscard]] std::string getName() const override { return "Coarsening"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
