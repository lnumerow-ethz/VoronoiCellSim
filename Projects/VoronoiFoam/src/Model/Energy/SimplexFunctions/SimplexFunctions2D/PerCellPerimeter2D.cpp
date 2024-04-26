#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellPerimeter2D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellPerimeter2D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = x0 * x0;
    F t4 = x1 * x1;
    F t5 = y0 * y0;
    F t8 = y1 * y1;
    F t10 = sqrt(-0.2e1 * x1 * x0 - 0.2e1 * y1 * y0 + t1 + t4 + t5 + t8);

    F unknown[1];

    unknown[0] = t10;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellPerimeter2D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = x0 * x0;
    F t4 = x1 * x1;
    F t5 = y0 * y0;
    F t8 = y1 * y1;
    F t10 = sqrt(-0.2e1 * x1 * x0 - 0.2e1 * y1 * y0 + t1 + t4 + t5 + t8);
    F t11 = 0.1e1 / t10;
    F t12 = x0 - x1;
    F t15 = y0 - y1;

    F unknown[4];

    unknown[0] = t12 * t11;
    unknown[1] = t15 * t11;
    unknown[2] = -t12 * t11;
    unknown[3] = -t15 * t11;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 4, 1);
    // clang-format on
}

void PerCellPerimeter2D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];

    F t1 = x0 * x0;
    F t4 = x1 * x1;
    F t5 = y0 * y0;
    F t8 = y1 * y1;
    F t9 = -0.2e1 * x1 * x0 - 0.2e1 * y1 * y0 + t1 + t4 + t5 + t8;
    F t10 = sqrt(t9);
    F t12 = 0.1e1 / t10 / t9;
    F t13 = x0 - x1;
    F t17 = 0.1e1 / t10;
    F t18 = -t13 * t13 * t12 + t17;
    F t19 = 0.2e1 * t13 * t12;
    F t20 = y0 - y1;
    F t22 = t20 * t19 / 0.2e1;
    F t25 = t13 * t19 / 0.2e1 - t17;
    F t27 = -t20 * t19 / 0.2e1;
    F t31 = -t20 * t20 * t12 + t17;
    F t32 = 0.2e1 * t20 * t12;
    F t34 = -t13 * t32 / 0.2e1;
    F t37 = t20 * t32 / 0.2e1 - t17;
    F t40 = t20 * t13 * t12;

    F unknown[4][4];

    unknown[0][0] = t18;
    unknown[0][1] = -t22;
    unknown[0][2] = t25;
    unknown[0][3] = -t27;
    unknown[1][0] = -t22;
    unknown[1][1] = t31;
    unknown[1][2] = -t34;
    unknown[1][3] = t37;
    unknown[2][0] = t25;
    unknown[2][1] = -t34;
    unknown[2][2] = t18;
    unknown[2][3] = -t40;
    unknown[3][0] = -t27;
    unknown[3][1] = t37;
    unknown[3][2] = -t40;
    unknown[3][3] = t31;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 4, 4);
    // clang-format on
}
