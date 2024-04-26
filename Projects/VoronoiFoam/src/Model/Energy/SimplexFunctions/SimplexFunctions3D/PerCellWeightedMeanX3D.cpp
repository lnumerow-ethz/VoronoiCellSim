#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellWeightedMean3D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellWeightedMeanX3D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
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

    unknown[0] = ((y1 * z2 - z1 * y2) * x0 + (-y0 * z2 + z0 * y2) * x1 + x2 * (z1 * y0 - y1 * z0)) * (x0 + x1 + x2) / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellWeightedMeanX3D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t3 = y1 * z2 - z1 * y2;
    F t4 = x0 + x1 + x2;
    F t6 = x0 * t3;
    F t9 = -y0 * z2 + z0 * y2;
    F t10 = x1 * t9;
    F t13 = z1 * y0 - y1 * z0;
    F t14 = t13 * x2;

    F unknown[9];

    unknown[0] = t4 * t3 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 + t6 / 0.24e2;
    unknown[1] = t4 * (-x1 * z2 + z1 * x2) / 0.24e2;
    unknown[2] = t4 * (x1 * y2 - y1 * x2) / 0.24e2;
    unknown[3] = t4 * t9 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 + t6 / 0.24e2;
    unknown[4] = t4 * (x0 * z2 - z0 * x2) / 0.24e2;
    unknown[5] = t4 * (-x0 * y2 + y0 * x2) / 0.24e2;
    unknown[6] = t4 * t13 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 + t6 / 0.24e2;
    unknown[7] = t4 * (-z1 * x0 + x1 * z0) / 0.24e2;
    unknown[8] = t4 * (y1 * x0 - x1 * y0) / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 9, 1);
    // clang-format on
}

void PerCellWeightedMeanX3D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
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
    F t4 = x1 * z2;
    F t5 = z1 * x2;
    F t6 = -t4 + t5;
    F t7 = x1 * y2;
    F t8 = y1 * x2;
    F t9 = t7 - t8;
    F t10 = y0 * z2;
    F t11 = z0 * y2;
    F t12 = -t10 + t1 + t11 - t2;
    F t13 = x0 + x1 + x2;
    F t14 = t13 * z2;
    F t15 = x0 * z2;
    F t16 = z0 * x2;
    F t17 = t14 + t15 - t16;
    F t18 = t13 * y2;
    F t19 = x0 * y2;
    F t20 = y0 * x2;
    F t21 = -t18 - t19 + t20;
    F t22 = z1 * y0;
    F t23 = y1 * z0;
    F t24 = t22 - t23 + t1 - t2;
    F t25 = t13 * z1;
    F t26 = z1 * x0;
    F t27 = x1 * z0;
    F t28 = -t25 - t26 + t27;
    F t29 = t13 * y1;
    F t30 = y1 * x0;
    F t31 = x1 * y0;
    F t32 = t29 + t30 - t31;
    F t33 = -t14 - t4 + t5;
    F t35 = t13 * x2 / 0.24e2;
    F t36 = t25 - t4 + t5;
    F t38 = t13 * x1 / 0.24e2;
    F t39 = t18 + t7 - t8;
    F t40 = -t29 + t7 - t8;
    F t42 = t15 - t16;
    F t43 = -t19 + t20;
    F t44 = t22 - t10 - t23 + t11;
    F t45 = t13 * z0;
    F t46 = t45 - t26 + t27;
    F t47 = t13 * y0;
    F t48 = -t47 + t30 - t31;
    F t49 = -t45 + t15 - t16;
    F t51 = t13 * x0 / 0.24e2;
    F t52 = t47 - t19 + t20;
    F t54 = -t26 + t27;
    F t55 = t30 - t31;

    F unknown[9][9];

    unknown[0][0] = t1 / 0.12e2 - t2 / 0.12e2;
    unknown[0][1] = t6 / 0.24e2;
    unknown[0][2] = t9 / 0.24e2;
    unknown[0][3] = t12 / 0.24e2;
    unknown[0][4] = t17 / 0.24e2;
    unknown[0][5] = t21 / 0.24e2;
    unknown[0][6] = t24 / 0.24e2;
    unknown[0][7] = t28 / 0.24e2;
    unknown[0][8] = t32 / 0.24e2;
    unknown[1][0] = t6 / 0.24e2;
    unknown[1][1] = 0.0e0;
    unknown[1][2] = 0.0e0;
    unknown[1][3] = t33 / 0.24e2;
    unknown[1][4] = 0.0e0;
    unknown[1][5] = t35;
    unknown[1][6] = t36 / 0.24e2;
    unknown[1][7] = 0.0e0;
    unknown[1][8] = -t38;
    unknown[2][0] = t9 / 0.24e2;
    unknown[2][1] = 0.0e0;
    unknown[2][2] = 0.0e0;
    unknown[2][3] = t39 / 0.24e2;
    unknown[2][4] = -t35;
    unknown[2][5] = 0.0e0;
    unknown[2][6] = t40 / 0.24e2;
    unknown[2][7] = t38;
    unknown[2][8] = 0.0e0;
    unknown[3][0] = t12 / 0.24e2;
    unknown[3][1] = t33 / 0.24e2;
    unknown[3][2] = t39 / 0.24e2;
    unknown[3][3] = -t10 / 0.12e2 + t11 / 0.12e2;
    unknown[3][4] = t42 / 0.24e2;
    unknown[3][5] = t43 / 0.24e2;
    unknown[3][6] = t44 / 0.24e2;
    unknown[3][7] = t46 / 0.24e2;
    unknown[3][8] = t48 / 0.24e2;
    unknown[4][0] = t17 / 0.24e2;
    unknown[4][1] = 0.0e0;
    unknown[4][2] = -t35;
    unknown[4][3] = t42 / 0.24e2;
    unknown[4][4] = 0.0e0;
    unknown[4][5] = 0.0e0;
    unknown[4][6] = t49 / 0.24e2;
    unknown[4][7] = 0.0e0;
    unknown[4][8] = t51;
    unknown[5][0] = t21 / 0.24e2;
    unknown[5][1] = t35;
    unknown[5][2] = 0.0e0;
    unknown[5][3] = t43 / 0.24e2;
    unknown[5][4] = 0.0e0;
    unknown[5][5] = 0.0e0;
    unknown[5][6] = t52 / 0.24e2;
    unknown[5][7] = -t51;
    unknown[5][8] = 0.0e0;
    unknown[6][0] = t24 / 0.24e2;
    unknown[6][1] = t36 / 0.24e2;
    unknown[6][2] = t40 / 0.24e2;
    unknown[6][3] = t44 / 0.24e2;
    unknown[6][4] = t49 / 0.24e2;
    unknown[6][5] = t52 / 0.24e2;
    unknown[6][6] = t22 / 0.12e2 - t23 / 0.12e2;
    unknown[6][7] = t54 / 0.24e2;
    unknown[6][8] = t55 / 0.24e2;
    unknown[7][0] = t28 / 0.24e2;
    unknown[7][1] = 0.0e0;
    unknown[7][2] = t38;
    unknown[7][3] = t46 / 0.24e2;
    unknown[7][4] = 0.0e0;
    unknown[7][5] = -t51;
    unknown[7][6] = t54 / 0.24e2;
    unknown[7][7] = 0.0e0;
    unknown[7][8] = 0.0e0;
    unknown[8][0] = t32 / 0.24e2;
    unknown[8][1] = -t38;
    unknown[8][2] = 0.0e0;
    unknown[8][3] = t48 / 0.24e2;
    unknown[8][4] = t51;
    unknown[8][5] = 0.0e0;
    unknown[8][6] = t55 / 0.24e2;
    unknown[8][7] = 0.0e0;
    unknown[8][8] = 0.0e0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 9, 9);
    // clang-format on
}
