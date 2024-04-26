#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGenerator2D.h"

class BoundaryGeneratorRigidBody2D : public BoundaryGenerator2D {
   protected:
    MatrixXF vertex_boundary;
    MatrixXI edge_boundary;
    MatrixXF vertex_rigid_body;
    MatrixXI edge_rigid_body;

    F rigid_body_force_magnitude;

   private:
    void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                         int order) const override;

   private:
    [[nodiscard]] bool checkValidNumBoundaryParams(int n_params) const override;

   public:
    [[nodiscard]] F computeEnergy(const Model &model) const override;

    void computeEnergyGradient(const Model &model, VectorXF &gradient) const override;

    void computeEnergyHessian(const Model &model, HessianF &hessian) const override;

   public:
    BoundaryGeneratorRigidBody2D(MatrixXF vertex_boundary, MatrixXI edge_boundary, MatrixXF vertex_rigid_body,
                                 MatrixXI edge_rigid_body, F rigid_body_force_magnitude);
};
