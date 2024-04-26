#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGenerator2D.h"

class BoundaryGeneratorBox2D : public BoundaryGenerator2D {
   private:
    void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                         int order) const override;

    [[nodiscard]] bool checkValidNumBoundaryParams(int n_params) const override;
};
