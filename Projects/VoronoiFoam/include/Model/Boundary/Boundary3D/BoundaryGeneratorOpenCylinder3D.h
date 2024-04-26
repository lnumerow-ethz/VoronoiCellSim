#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGenerator3D.h"
#include <set>

class BoundaryGeneratorOpenCylinder3D : public BoundaryGenerator3D {
   private:
    MatrixXI tri;
    MatrixXF original_points;
    std::set<int> top_edge_points;
    std::set<int> top_interior_points;
    VectorXI idx;  /// Vertex (flattened) index to param index.
    F z_base;
    F spring_constant;
    F barrier_weight = 1e-8;

   private:
    void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                         int order) const override;

    [[nodiscard]] bool checkValidNumBoundaryParams(int n_params) const override;

   public:
    [[nodiscard]] F computeEnergy(const Model &model) const override;

    void computeEnergyGradient(const Model &model, VectorXF &gradient) const override;

    void computeEnergyHessian(const Model &model, HessianF &hessian) const override;

   public:
    BoundaryGeneratorOpenCylinder3D(F spring_constant);

    void getDefaultParams(VectorXF &params);
};
