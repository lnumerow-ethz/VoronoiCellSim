#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"
#include "Projects/VoronoiFoam/include/Model/Energy/DynamicsObjective.h"
#include <fstream>

struct ConvergenceTest2DParams {
    int tessellation = 0;  /// Voronoi (0) or Power (1)
    DynamicsDiscretization discretization;

    F area_weight = 0;
    F perimeter_weight = 0;
    F centroid_weight = 0;
    F second_moment_weight = 0;
    F time_step = 0;
    int num_time_steps = 0;

    ConvergenceTest2DParams(int tessellation, DynamicsDiscretization discretization, F area_weight, F perimeter_weight,
                            F centroid_weight, F second_moment_weight, F time_step, int num_time_steps)
        : tessellation(tessellation),
          discretization(discretization),
          area_weight(area_weight),
          perimeter_weight(perimeter_weight),
          centroid_weight(centroid_weight),
          second_moment_weight(second_moment_weight),
          time_step(time_step),
          num_time_steps(num_time_steps) {}
};

class ExperimentConvergenceTest2D : public Experiment {
   private:
    int run_index = 0;
    std::vector<ConvergenceTest2DParams> run_params;

    int current_time_step = 0;
    std::stringstream output_string;

   public:
    [[nodiscard]] std::string getName() const override { return "Convergence Test"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;

   public:
    ExperimentConvergenceTest2D();
};
