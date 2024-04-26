#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorMesh3D.h"

class BoundaryGeneratorMeshLaplacian3D : public BoundaryGeneratorMesh3D {
   private:
    F spring_constant;

   public:
    [[nodiscard]] F computeEnergy(const Model &model) const override;

    void computeEnergyGradient(const Model &model, VectorXF &gradient) const override;

    void computeEnergyHessian(const Model &model, HessianF &hessian) const override;

   public:
    BoundaryGeneratorMeshLaplacian3D(const MatrixXI &tri, F spring_constant);
};
