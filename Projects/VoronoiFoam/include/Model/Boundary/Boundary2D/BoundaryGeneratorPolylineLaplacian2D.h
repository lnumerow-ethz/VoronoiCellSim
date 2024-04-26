#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorPolyline2D.h"

class BoundaryGeneratorPolylineLaplacian2D : public BoundaryGeneratorPolyline2D {
   private:
    F spring_constant;

   public:
    [[nodiscard]] F computeEnergy(const Model &model) const override;

    void computeEnergyGradient(const Model &model, VectorXF &gradient) const override;

    void computeEnergyHessian(const Model &model, HessianF &hessian) const override;

   public:
    BoundaryGeneratorPolylineLaplacian2D(const MatrixXI &edge, F spring_constant);
};
