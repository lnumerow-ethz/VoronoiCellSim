#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellWeightedMean2D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellWeightedMeanY2D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];


    F unknown[1];

    unknown[0] = 0.1666666667e0 * (x0 * y1 - 0.1e1 * x1 * y0) * (y0 + y1);

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellWeightedMeanY2D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = y0 + y1;
    F t7 = 0.1666666667e0 * x0 * y1;
    F t9 = 0.1666666667e0 * x1 * y0;

    F unknown[4];

    unknown[0] = 0.1666666667e0 * t1 * y1;
    unknown[1] = -0.1666666667e0 * t1 * x1 + t7 - t9;
    unknown[2] = -0.1666666667e0 * t1 * y0;
    unknown[3] = 0.1666666667e0 * t1 * x0 + t7 - t9;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 4, 1);
    // clang-format on
}

void PerCellWeightedMeanY2D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = 0.1666666667e0 * y1;
    F t2 = 0.1666666667e0 * y0;
    F t4 = t2 + 0.3333333334e0 * y1;
    F t7 = -0.3333333334e0 * y0 - t1;
    F t10 = -0.1666666667e0 * x1 + 0.1666666667e0 * x0;

    F unknown[4][4];

    unknown[0][0] = 0.0e0;
    unknown[0][1] = t1;
    unknown[0][2] = 0.0e0;
    unknown[0][3] = t4;
    unknown[1][0] = t1;
    unknown[1][1] = -0.3333333334e0 * x1;
    unknown[1][2] = t7;
    unknown[1][3] = t10;
    unknown[2][0] = 0.0e0;
    unknown[2][1] = t7;
    unknown[2][2] = 0.0e0;
    unknown[2][3] = -t2;
    unknown[3][0] = t4;
    unknown[3][1] = t10;
    unknown[3][2] = -t2;
    unknown[3][3] = 0.3333333334e0 * x0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 4, 4);
    // clang-format on
}
