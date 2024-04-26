#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellWeightedMean3D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellWeightedMeanY3D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
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

    unknown[0] = (y0 + y1 + y2) * ((-z1 * x0 + x1 * z0) * y2 + (x0 * z2 - z0 * x2) * y1 - y0 * (x1 * z2 - z1 * x2)) / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellWeightedMeanY3D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t1 = y0 + y1 + y2;
    F t9 = -z1 * x0 + x1 * z0;
    F t10 = y2 * t9;
    F t13 = x0 * z2 - z0 * x2;
    F t14 = y1 * t13;
    F t17 = x1 * z2 - z1 * x2;
    F t18 = t17 * y0;

    F unknown[9];

    unknown[0] = (y1 * z2 - z1 * y2) * t1 / 0.24e2;
    unknown[1] = -t17 * t1 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 - t18 / 0.24e2;
    unknown[2] = (x1 * y2 - y1 * x2) * t1 / 0.24e2;
    unknown[3] = (-y0 * z2 + z0 * y2) * t1 / 0.24e2;
    unknown[4] = t13 * t1 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 - t18 / 0.24e2;
    unknown[5] = (-x0 * y2 + y0 * x2) * t1 / 0.24e2;
    unknown[6] = (z1 * y0 - y1 * z0) * t1 / 0.24e2;
    unknown[7] = t9 * t1 / 0.24e2 + t10 / 0.24e2 + t14 / 0.24e2 - t18 / 0.24e2;
    unknown[8] = (y1 * x0 - x1 * y0) * t1 / 0.24e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 9, 1);
    // clang-format on
}

void PerCellWeightedMeanY3D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
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
    F t4 = y0 + y1 + y2;
    F t5 = t4 * z2;
    F t6 = t1 - t2 + t5;
    F t8 = t4 * y2 / 0.24e2;
    F t9 = t4 * z1;
    F t10 = t1 - t2 - t9;
    F t12 = t4 * y1 / 0.24e2;
    F t13 = x1 * z2;
    F t14 = z1 * x2;
    F t16 = x1 * y2;
    F t17 = y1 * x2;
    F t18 = t16 - t17;
    F t19 = z0 * y2;
    F t20 = y0 * z2;
    F t21 = t19 - t20 - t5;
    F t22 = x0 * z2;
    F t23 = z0 * x2;
    F t24 = t22 - t13 - t23 + t14;
    F t25 = x0 * y2;
    F t26 = y0 * x2;
    F t27 = t4 * x2;
    F t28 = -t25 + t26 + t27;
    F t29 = y1 * z0;
    F t30 = z1 * y0;
    F t31 = -t29 + t30 + t9;
    F t32 = z1 * x0;
    F t33 = x1 * z0;
    F t34 = -t32 + t33 - t13 + t14;
    F t35 = y1 * x0;
    F t36 = x1 * y0;
    F t37 = t4 * x1;
    F t38 = t35 - t36 - t37;
    F t39 = t16 - t17 - t27;
    F t40 = t16 - t17 + t37;
    F t41 = -t20 + t19;
    F t42 = t4 * z0;
    F t43 = -t20 + t19 + t42;
    F t45 = t4 * y0 / 0.24e2;
    F t47 = -t25 + t26;
    F t48 = -t29 + t30 - t42;
    F t49 = -t32 + t22 + t33 - t23;
    F t50 = t4 * x0;
    F t51 = t35 - t36 + t50;
    F t52 = -t25 + t26 - t50;
    F t53 = -t29 + t30;
    F t55 = -t36 + t35;

    F unknown[9][9];

    unknown[0][0] = 0.0e0;
    unknown[0][1] = t3 / 0.24e2;
    unknown[0][2] = 0.0e0;
    unknown[0][3] = 0.0e0;
    unknown[0][4] = t6 / 0.24e2;
    unknown[0][5] = -t8;
    unknown[0][6] = 0.0e0;
    unknown[0][7] = t10 / 0.24e2;
    unknown[0][8] = t12;
    unknown[1][0] = t3 / 0.24e2;
    unknown[1][1] = -t13 / 0.12e2 + t14 / 0.12e2;
    unknown[1][2] = t18 / 0.24e2;
    unknown[1][3] = t21 / 0.24e2;
    unknown[1][4] = t24 / 0.24e2;
    unknown[1][5] = t28 / 0.24e2;
    unknown[1][6] = t31 / 0.24e2;
    unknown[1][7] = t34 / 0.24e2;
    unknown[1][8] = t38 / 0.24e2;
    unknown[2][0] = 0.0e0;
    unknown[2][1] = t18 / 0.24e2;
    unknown[2][2] = 0.0e0;
    unknown[2][3] = t8;
    unknown[2][4] = t39 / 0.24e2;
    unknown[2][5] = 0.0e0;
    unknown[2][6] = -t12;
    unknown[2][7] = t40 / 0.24e2;
    unknown[2][8] = 0.0e0;
    unknown[3][0] = 0.0e0;
    unknown[3][1] = t21 / 0.24e2;
    unknown[3][2] = t8;
    unknown[3][3] = 0.0e0;
    unknown[3][4] = t41 / 0.24e2;
    unknown[3][5] = 0.0e0;
    unknown[3][6] = 0.0e0;
    unknown[3][7] = t43 / 0.24e2;
    unknown[3][8] = -t45;
    unknown[4][0] = t6 / 0.24e2;
    unknown[4][1] = t24 / 0.24e2;
    unknown[4][2] = t39 / 0.24e2;
    unknown[4][3] = t41 / 0.24e2;
    unknown[4][4] = t22 / 0.12e2 - t23 / 0.12e2;
    unknown[4][5] = t47 / 0.24e2;
    unknown[4][6] = t48 / 0.24e2;
    unknown[4][7] = t49 / 0.24e2;
    unknown[4][8] = t51 / 0.24e2;
    unknown[5][0] = -t8;
    unknown[5][1] = t28 / 0.24e2;
    unknown[5][2] = 0.0e0;
    unknown[5][3] = 0.0e0;
    unknown[5][4] = t47 / 0.24e2;
    unknown[5][5] = 0.0e0;
    unknown[5][6] = t45;
    unknown[5][7] = t52 / 0.24e2;
    unknown[5][8] = 0.0e0;
    unknown[6][0] = 0.0e0;
    unknown[6][1] = t31 / 0.24e2;
    unknown[6][2] = -t12;
    unknown[6][3] = 0.0e0;
    unknown[6][4] = t48 / 0.24e2;
    unknown[6][5] = t45;
    unknown[6][6] = 0.0e0;
    unknown[6][7] = t53 / 0.24e2;
    unknown[6][8] = 0.0e0;
    unknown[7][0] = t10 / 0.24e2;
    unknown[7][1] = t34 / 0.24e2;
    unknown[7][2] = t40 / 0.24e2;
    unknown[7][3] = t43 / 0.24e2;
    unknown[7][4] = t49 / 0.24e2;
    unknown[7][5] = t52 / 0.24e2;
    unknown[7][6] = t53 / 0.24e2;
    unknown[7][7] = t33 / 0.12e2 - t32 / 0.12e2;
    unknown[7][8] = t55 / 0.24e2;
    unknown[8][0] = t12;
    unknown[8][1] = t38 / 0.24e2;
    unknown[8][2] = 0.0e0;
    unknown[8][3] = -t45;
    unknown[8][4] = t51 / 0.24e2;
    unknown[8][5] = 0.0e0;
    unknown[8][6] = 0.0e0;
    unknown[8][7] = t55 / 0.24e2;
    unknown[8][8] = 0.0e0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 9, 9);
    // clang-format on
}
