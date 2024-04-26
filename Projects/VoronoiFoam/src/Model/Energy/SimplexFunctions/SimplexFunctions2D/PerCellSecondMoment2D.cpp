#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellSecondMoment2D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellSecondMoment2D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t4 = x0 * x0;
    F t6 = x1 * x1;
    F t7 = y0 * y0;
    F t9 = y1 * y1;

    F unknown[1];

    unknown[0] = (x0 * x1 + y0 * y1 + t4 + t6 + t7 + t9) * (x0 * y1 - x1 * y0) / 0.12e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellSecondMoment2D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = x0 * x0;
    F t3 = x1 * x1;
    F t4 = y0 * y0;
    F t6 = y1 * y1;
    F t7 = x0 * x1 + y0 * y1 + t1 + t3 + t4 + t6;
    F t11 = x0 * y1 - x1 * y0;

    F unknown[4];

    unknown[0] = t7 * y1 / 0.12e2 + (0.2e1 * x0 + x1) * t11 / 0.12e2;
    unknown[1] = -t7 * x1 / 0.12e2 + (0.2e1 * y0 + y1) * t11 / 0.12e2;
    unknown[2] = -t7 * y0 / 0.12e2 + (x0 + 0.2e1 * x1) * t11 / 0.12e2;
    unknown[3] = t7 * x0 / 0.12e2 + (y0 + 0.2e1 * y1) * t11 / 0.12e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 4, 1);
    // clang-format on
}

void PerCellSecondMoment2D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t2 = 0.2e1 * x0 + x1;
    F t4 = x0 * y1;
    F t5 = x1 * y0;
    F t8 = 0.2e1 * y0 + y1;
    F t11 = -t2 * x1 + t8 * y1;
    F t13 = x0 + 0.2e1 * x1;
    F t16 = t13 * y1 - t2 * y0 + t4 - t5;
    F t17 = x0 * x0;
    F t18 = x0 * x1;
    F t19 = x1 * x1;
    F t20 = y0 * y0;
    F t21 = y0 * y1;
    F t22 = y1 * y1;
    F t24 = y0 + 0.2e1 * y1;
    F t27 = t2 * x0 + t24 * y1 + t17 + t18 + t19 + t20 + t21 + t22;
    F t32 = -t13 * x1 - t8 * y0 - t17 - t18 - t19 - t20 - t21 - t22;
    F t35 = -t24 * x1 + t8 * x0 + t4 - t5;
    F t40 = t13 * x0 - t24 * y0;

    F unknown[4][4];

    unknown[0][0] = t2 * y1 / 0.6e1 + t4 / 0.6e1 - t5 / 0.6e1;
    unknown[0][1] = t11 / 0.12e2;
    unknown[0][2] = t16 / 0.12e2;
    unknown[0][3] = t27 / 0.12e2;
    unknown[1][0] = t11 / 0.12e2;
    unknown[1][1] = -t8 * x1 / 0.6e1 + t4 / 0.6e1 - t5 / 0.6e1;
    unknown[1][2] = t32 / 0.12e2;
    unknown[1][3] = t35 / 0.12e2;
    unknown[2][0] = t16 / 0.12e2;
    unknown[2][1] = t32 / 0.12e2;
    unknown[2][2] = -t13 * y0 / 0.6e1 + t4 / 0.6e1 - t5 / 0.6e1;
    unknown[2][3] = t40 / 0.12e2;
    unknown[3][0] = t27 / 0.12e2;
    unknown[3][1] = t35 / 0.12e2;
    unknown[3][2] = t40 / 0.12e2;
    unknown[3][3] = t24 * x0 / 0.6e1 + t4 / 0.6e1 - t5 / 0.6e1;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 4, 4);
    // clang-format on
}
