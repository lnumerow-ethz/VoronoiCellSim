#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGenerator2D.h"

class BoundaryGeneratorPolyline2D : public BoundaryGenerator2D {
   protected:
    MatrixXI edge;

   private:
    void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                         int order) const override;

    [[nodiscard]] bool checkValidNumBoundaryParams(int n_params) const override;

   public:
    explicit BoundaryGeneratorPolyline2D(MatrixXI edge);
};
