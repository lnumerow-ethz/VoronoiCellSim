#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGenerator3D.h"

class BoundaryGeneratorMesh3D : public BoundaryGenerator3D {
   protected:
    MatrixXI tri;

   private:
    void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                         int order) const override;

    [[nodiscard]] bool checkValidNumBoundaryParams(int n_params) const override;

   public:
    explicit BoundaryGeneratorMesh3D(MatrixXI tri);
};
