#pragma once

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/PerCellFunctionFromSimplex.h"

/// Second moment of area about perpendicular axis passing through ***ORIGIN***.
/// Use parallel axis theorem to get moment about centroid.
class PerCellSecondMoment2D : public PerCellFunctionFromSimplex {
   private:
    void getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const override;
};
