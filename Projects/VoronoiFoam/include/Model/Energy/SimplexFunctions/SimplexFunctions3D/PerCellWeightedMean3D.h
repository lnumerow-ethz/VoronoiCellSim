#pragma once

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/PerCellFunctionFromSimplex.h"

class PerCellWeightedMeanX3D : public PerCellFunctionFromSimplex {
   private:
    void getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const override;
};

class PerCellWeightedMeanY3D : public PerCellFunctionFromSimplex {
   private:
    void getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const override;
};

class PerCellWeightedMeanZ3D : public PerCellFunctionFromSimplex {
   private:
    void getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const override;

    void getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const override;
};
