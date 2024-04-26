#pragma once

#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"
#include <fstream>

#include <random>
#include <chrono>

class ExperimentImageMatch2D : public Experiment {
   private:
    void computeImageMatchObjective(FoamSubApp* foam_sub_app, F& objective);
    void computeImageMatchGradient(FoamSubApp* foam_sub_app, VectorXF& pLpc, VectorXF& pLpu);
    void computeMixedHessian(FoamSubApp* foam_sub_app, SparseMatrixF& d2Edy2, SparseMatrixF& d2Edydu);

   private:
    F dx, dy;
    std::vector<TessellationNode> target_nodes;

    MatrixXF bfgs_inverse_hessian;
    VectorXF g_prev;

    bool line_search_in_progress = false;
    F L_prev = 1e10;
    VectorXF y_prev;
    VectorXF u_prev;
    VectorXF search_direction;
    F alpha;

    int frame = 0;
    int num_frame_iters = 0;
    std::chrono::high_resolution_clock::time_point iter_start_time;

    int num_pre_steps = 0;

   private:
    void loadImageNetwork(FoamSubApp* foam_sub_app, std::string image_file_name, std::string network_file_name);

   public:
    [[nodiscard]] std::string getName() const override { return "Image Match"; };
    void setup(FoamSubApp* foam_sub_app) override;
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override;
};
