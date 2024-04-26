#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellWeightedMean3D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellWeightedMeanZ3D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F z0 = inputs[2];
    F x1 = inputs[3];
    F y1 = inputs[4];
    F z1 = inputs[5];
    F x2 = inputs[6];
    F y2 = inputs[7];
    F z2 = inputs[8];


    F unknown[1];

    unknown[0] = ((x1 * y2 - y1 * x2) * z0 + (-x0 * y2 + y0 * x2) * z1 + z2 * (y1 * x0 - x1 * y0)) * (z0 + z1 + z2) / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellWeightedMeanZ3D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F z0 = inputs[2];
    F x1 = inputs[3];
    F y1 = inputs[4];
    F z1 = inputs[5];
    F x2 = inputs[6];
    F y2 = inputs[7];
    F z2 = inputs[8];

    F t4 = z0 + z1 + z2;
    F t14 = x1 * y2 - y1 * x2;
    F t16 = z0 * t14;
    F t19 = -x0 * y2 + y0 * x2;
    F t20 = z1 * t19;
    F t23 = y1 * x0 - x1 * y0;
    F t24 = t23 * z2;

    F unknown[9];

    unknown[0] = t4 * (y1 * z2 - z1 * y2) / 0.24e2;
    unknown[1] = t4 * (-x1 * z2 + z1 * x2) / 0.24e2;
    unknown[2] = t4 * t14 / 0.24e2 + t16 / 0.24e2 + t20 / 0.24e2 + t24 / 0.24e2;
    unknown[3] = t4 * (-y0 * z2 + z0 * y2) / 0.24e2;
    unknown[4] = t4 * (x0 * z2 - z0 * x2) / 0.24e2;
    unknown[5] = t4 * t19 / 0.24e2 + t16 / 0.24e2 + t20 / 0.24e2 + t24 / 0.24e2;
    unknown[6] = t4 * (z1 * y0 - y1 * z0) / 0.24e2;
    unknown[7] = t4 * (-z1 * x0 + x1 * z0) / 0.24e2;
    unknown[8] = t4 * t23 / 0.24e2 + t16 / 0.24e2 + t20 / 0.24e2 + t24 / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 9, 1);
    // clang-format on
}

void PerCellWeightedMeanZ3D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F z0 = inputs[2];
    F x1 = inputs[3];
    F y1 = inputs[4];
    F z1 = inputs[5];
    F x2 = inputs[6];
    F y2 = inputs[7];
    F z2 = inputs[8];

    F t1 = y1 * z2;
    F t2 = z1 * y2;
    F t3 = t1 - t2;
    F t4 = z0 + z1 + z2;
    F t6 = t4 * z2 / 0.24e2;
    F t7 = t4 * y2;
    F t8 = -t7 + t1 - t2;
    F t10 = t4 * z1 / 0.24e2;
    F t11 = t4 * y1;
    F t12 = t11 + t1 - t2;
    F t13 = x1 * z2;
    F t14 = z1 * x2;
    F t15 = -t13 + t14;
    F t16 = t4 * x2;
    F t17 = t16 - t13 + t14;
    F t18 = t4 * x1;
    F t19 = -t18 - t13 + t14;
    F t20 = x1 * y2;
    F t21 = y1 * x2;
    F t23 = z0 * y2;
    F t24 = y0 * z2;
    F t25 = t7 + t23 - t24;
    F t26 = z0 * x2;
    F t27 = x0 * z2;
    F t28 = -t16 - t26 + t27;
    F t29 = x0 * y2;
    F t30 = y0 * x2;
    F t31 = -t29 + t20 + t30 - t21;
    F t32 = y1 * z0;
    F t33 = z1 * y0;
    F t34 = -t11 - t32 + t33;
    F t35 = x1 * z0;
    F t36 = z1 * x0;
    F t37 = t18 + t35 - t36;
    F t38 = y1 * x0;
    F t39 = x1 * y0;
    F t40 = t38 - t39 + t20 - t21;
    F t41 = -t24 + t23;
    F t43 = t4 * z0 / 0.24e2;
    F t44 = t4 * y0;
    F t45 = -t44 - t24 + t23;
    F t46 = -t26 + t27;
    F t47 = t4 * x0;
    F t48 = t47 + t27 - t26;
    F t50 = t44 - t32 + t33;
    F t51 = -t47 + t35 - t36;
    F t52 = t38 - t29 - t39 + t30;
    F t53 = t33 - t32;
    F t54 = -t36 + t35;

    F unknown[9][9];

    unknown[0][0] = 0.0e0;
    unknown[0][1] = 0.0e0;
    unknown[0][2] = t3 / 0.24e2;
    unknown[0][3] = 0.0e0;
    unknown[0][4] = t6;
    unknown[0][5] = t8 / 0.24e2;
    unknown[0][6] = 0.0e0;
    unknown[0][7] = -t10;
    unknown[0][8] = t12 / 0.24e2;
    unknown[1][0] = 0.0e0;
    unknown[1][1] = 0.0e0;
    unknown[1][2] = t15 / 0.24e2;
    unknown[1][3] = -t6;
    unknown[1][4] = 0.0e0;
    unknown[1][5] = t17 / 0.24e2;
    unknown[1][6] = t10;
    unknown[1][7] = 0.0e0;
    unknown[1][8] = t19 / 0.24e2;
    unknown[2][0] = t3 / 0.24e2;
    unknown[2][1] = t15 / 0.24e2;
    unknown[2][2] = t20 / 0.12e2 - t21 / 0.12e2;
    unknown[2][3] = t25 / 0.24e2;
    unknown[2][4] = t28 / 0.24e2;
    unknown[2][5] = t31 / 0.24e2;
    unknown[2][6] = t34 / 0.24e2;
    unknown[2][7] = t37 / 0.24e2;
    unknown[2][8] = t40 / 0.24e2;
    unknown[3][0] = 0.0e0;
    unknown[3][1] = -t6;
    unknown[3][2] = t25 / 0.24e2;
    unknown[3][3] = 0.0e0;
    unknown[3][4] = 0.0e0;
    unknown[3][5] = t41 / 0.24e2;
    unknown[3][6] = 0.0e0;
    unknown[3][7] = t43;
    unknown[3][8] = t45 / 0.24e2;
    unknown[4][0] = t6;
    unknown[4][1] = 0.0e0;
    unknown[4][2] = t28 / 0.24e2;
    unknown[4][3] = 0.0e0;
    unknown[4][4] = 0.0e0;
    unknown[4][5] = t46 / 0.24e2;
    unknown[4][6] = -t43;
    unknown[4][7] = 0.0e0;
    unknown[4][8] = t48 / 0.24e2;
    unknown[5][0] = t8 / 0.24e2;
    unknown[5][1] = t17 / 0.24e2;
    unknown[5][2] = t31 / 0.24e2;
    unknown[5][3] = t41 / 0.24e2;
    unknown[5][4] = t46 / 0.24e2;
    unknown[5][5] = -t29 / 0.12e2 + t30 / 0.12e2;
    unknown[5][6] = t50 / 0.24e2;
    unknown[5][7] = t51 / 0.24e2;
    unknown[5][8] = t52 / 0.24e2;
    unknown[6][0] = 0.0e0;
    unknown[6][1] = t10;
    unknown[6][2] = t34 / 0.24e2;
    unknown[6][3] = 0.0e0;
    unknown[6][4] = -t43;
    unknown[6][5] = t50 / 0.24e2;
    unknown[6][6] = 0.0e0;
    unknown[6][7] = 0.0e0;
    unknown[6][8] = t53 / 0.24e2;
    unknown[7][0] = -t10;
    unknown[7][1] = 0.0e0;
    unknown[7][2] = t37 / 0.24e2;
    unknown[7][3] = t43;
    unknown[7][4] = 0.0e0;
    unknown[7][5] = t51 / 0.24e2;
    unknown[7][6] = 0.0e0;
    unknown[7][7] = 0.0e0;
    unknown[7][8] = t54 / 0.24e2;
    unknown[8][0] = t12 / 0.24e2;
    unknown[8][1] = t19 / 0.24e2;
    unknown[8][2] = t40 / 0.24e2;
    unknown[8][3] = t45 / 0.24e2;
    unknown[8][4] = t48 / 0.24e2;
    unknown[8][5] = t52 / 0.24e2;
    unknown[8][6] = t53 / 0.24e2;
    unknown[8][7] = t54 / 0.24e2;
    unknown[8][8] = t38 / 0.12e2 - t39 / 0.12e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 9, 9);
    // clang-format on
}
