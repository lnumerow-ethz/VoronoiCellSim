#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"

class ExperimentCoarsening2D : public Experiment {
    int num_cells = 100;

   public:
    [[nodiscard]] std::string getName() const override { return "Coarsening"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
