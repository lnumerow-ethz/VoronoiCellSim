#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellArea2D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellArea2D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];


    F unknown[1];

    unknown[0] = 0.5e0 * x0 * y1 - 0.5e0 * x1 * y0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellArea2D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];


    F unknown[4];

    unknown[0] = 0.5e0 * y1;
    unknown[1] = -0.5e0 * x1;
    unknown[2] = -0.5e0 * y0;
    unknown[3] = 0.5e0 * x0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 4, 1);
    // clang-format on
}

void PerCellArea2D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];


    F unknown[4][4];

    unknown[0][0] = 0.0e0;
    unknown[0][1] = 0.0e0;
    unknown[0][2] = 0.0e0;
    unknown[0][3] = 0.5e0;
    unknown[1][0] = 0.0e0;
    unknown[1][1] = 0.0e0;
    unknown[1][2] = -0.5e0;
    unknown[1][3] = 0.0e0;
    unknown[2][0] = 0.0e0;
    unknown[2][1] = -0.5e0;
    unknown[2][2] = 0.0e0;
    unknown[2][3] = 0.0e0;
    unknown[3][0] = 0.5e0;
    unknown[3][1] = 0.0e0;
    unknown[3][2] = 0.0e0;
    unknown[3][3] = 0.0e0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 4, 4);
    // clang-format on
}
