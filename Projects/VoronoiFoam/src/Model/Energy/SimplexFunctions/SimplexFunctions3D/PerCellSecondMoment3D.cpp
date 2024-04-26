#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellSecondMoment3D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellSecondMoment3D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t14 = x0 * x0;
    F t17 = x1 * x1;
    F t19 = x2 * x2;
    F t20 = y0 * y0;
    F t23 = y1 * y1;
    F t25 = y2 * y2;
    F t26 = z0 * z0;
    F t29 = z1 * z1;
    F t31 = z2 * z2;
    F t32 = t14 + (x1 + x2) * x0 + t17 + x1 * x2 + t19 + t20 + (y1 + y2) * y0 + t23 + y1 * y2 + t25 + t26 + (z1 + z2) * z0 + t29 + z1 * z2 + t31;

    F unknown[1];

    unknown[0] = t32 * ((y1 * z2 - y2 * z1) * x0 + (-y0 * z2 + y2 * z0) * x1 + x2 * (y0 * z1 - y1 * z0)) / 0.60e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellSecondMoment3D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t3 = y1 * z2 - y2 * z1;
    F t4 = x0 * x0;
    F t7 = x1 * x1;
    F t9 = x2 * x2;
    F t10 = y0 * y0;
    F t13 = y1 * y1;
    F t15 = y2 * y2;
    F t16 = z0 * z0;
    F t19 = z1 * z1;
    F t21 = z2 * z2;
    F t22 = t4 + (x1 + x2) * x0 + t7 + x1 * x2 + t9 + t10 + (y1 + y2) * y0 + t13 + y1 * y2 + t15 + t16 + (z1 + z2) * z0 + t19 + z1 * z2 + t21;
    F t27 = -y0 * z2 + y2 * z0;
    F t31 = y0 * z1 - y1 * z0;
    F t33 = x1 * t27 + x0 * t3 + t31 * x2;

    F unknown[9];

    unknown[0] = t22 * t3 / 0.60e2 + (0.2e1 * x0 + x1 + x2) * t33 / 0.60e2;
    unknown[1] = t22 * (-x1 * z2 + x2 * z1) / 0.60e2 + (0.2e1 * y0 + y1 + y2) * t33 / 0.60e2;
    unknown[2] = t22 * (x1 * y2 - x2 * y1) / 0.60e2 + (0.2e1 * z0 + z1 + z2) * t33 / 0.60e2;
    unknown[3] = t22 * t27 / 0.60e2 + (x0 + 0.2e1 * x1 + x2) * t33 / 0.60e2;
    unknown[4] = t22 * (x0 * z2 - x2 * z0) / 0.60e2 + (y0 + 0.2e1 * y1 + y2) * t33 / 0.60e2;
    unknown[5] = t22 * (-x0 * y2 + x2 * y0) / 0.60e2 + (z0 + 0.2e1 * z1 + z2) * t33 / 0.60e2;
    unknown[6] = t22 * t31 / 0.60e2 + (x0 + x1 + 0.2e1 * x2) * t33 / 0.60e2;
    unknown[7] = t22 * (-x0 * z1 + x1 * z0) / 0.60e2 + (y0 + y1 + 0.2e1 * y2) * t33 / 0.60e2;
    unknown[8] = t22 * (x0 * y1 - x1 * y0) / 0.60e2 + (z0 + z1 + 0.2e1 * z2) * t33 / 0.60e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 9, 1);
    // clang-format on
}

void PerCellSecondMoment3D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t3 = y1 * z2 - y2 * z1;
    F t5 = 0.2e1 * x0 + x1 + x2;
    F t7 = x0 * t3;
    F t10 = -y0 * z2 + y2 * z0;
    F t11 = x1 * t10;
    F t14 = y0 * z1 - y1 * z0;
    F t15 = t14 * x2;
    F t18 = 0.2e1 * y0 + y1 + y2;
    F t22 = -x1 * z2 + x2 * z1;
    F t24 = t18 * t3 + t5 * t22;
    F t26 = 0.2e1 * z0 + z1 + z2;
    F t30 = x1 * y2 - x2 * y1;
    F t32 = t26 * t3 + t5 * t30;
    F t34 = x0 + 0.2e1 * x1 + x2;
    F t37 = t5 * t10 + t34 * t3 + t11 + t15 + t7;
    F t38 = x0 * x0;
    F t41 = x1 * x1;
    F t43 = x2 * x2;
    F t44 = y0 * y0;
    F t47 = y1 * y1;
    F t49 = y2 * y2;
    F t50 = z0 * z0;
    F t53 = z1 * z1;
    F t55 = z2 * z2;
    F t56 = t38 + (x1 + x2) * x0 + t41 + x1 * x2 + t43 + t44 + (y1 + y2) * y0 + t47 + y1 * y2 + t49 + t50 + (z1 + z2) * z0 + t53 + z1 * z2 + t55;
    F t57 = t56 * z2;
    F t59 = y0 + 0.2e1 * y1 + y2;
    F t63 = x0 * z2 - x2 * z0;
    F t65 = t59 * t3 + t5 * t63 + t57;
    F t66 = t56 * y2;
    F t68 = z0 + 0.2e1 * z1 + z2;
    F t72 = -x0 * y2 + x2 * y0;
    F t74 = t68 * t3 + t5 * t72 - t66;
    F t76 = x0 + x1 + 0.2e1 * x2;
    F t79 = t5 * t14 + t76 * t3 + t11 + t15 + t7;
    F t80 = t56 * z1;
    F t82 = y0 + y1 + 0.2e1 * y2;
    F t86 = -x0 * z1 + x1 * z0;
    F t88 = t82 * t3 + t5 * t86 - t80;
    F t89 = t56 * y1;
    F t91 = z0 + z1 + 0.2e1 * z2;
    F t95 = x0 * y1 - x1 * y0;
    F t97 = t91 * t3 + t5 * t95 + t89;
    F t102 = t18 * t30 + t26 * t22;
    F t105 = t18 * t10 + t34 * t22 - t57;
    F t108 = t18 * t63 + t59 * t22 + t11 + t15 + t7;
    F t109 = t56 * x2;
    F t112 = t18 * t72 + t68 * t22 + t109;
    F t115 = t18 * t14 + t76 * t22 + t80;
    F t118 = t18 * t86 + t22 * t82 + t11 + t15 + t7;
    F t119 = t56 * x1;
    F t122 = t18 * t95 + t91 * t22 - t119;
    F t127 = t26 * t10 + t34 * t30 + t66;
    F t130 = t26 * t63 + t59 * t30 - t109;
    F t133 = t26 * t72 + t68 * t30 + t11 + t15 + t7;
    F t136 = t26 * t14 + t76 * t30 - t89;
    F t139 = t26 * t86 + t82 * t30 + t119;
    F t142 = t26 * t95 + t91 * t30 + t11 + t15 + t7;
    F t147 = t59 * t10 + t34 * t63;
    F t150 = t68 * t10 + t34 * t72;
    F t153 = t76 * t10 + t34 * t14 + t11 + t15 + t7;
    F t154 = t56 * z0;
    F t157 = t82 * t10 + t34 * t86 + t154;
    F t158 = t56 * y0;
    F t161 = t91 * t10 + t34 * t95 - t158;
    F t166 = t59 * t72 + t68 * t63;
    F t169 = t59 * t14 + t76 * t63 - t154;
    F t172 = t59 * t86 + t82 * t63 + t11 + t15 + t7;
    F t173 = t56 * x0;
    F t176 = t59 * t95 + t91 * t63 + t173;
    F t181 = t68 * t14 + t76 * t72 + t158;
    F t184 = t68 * t86 + t82 * t72 - t173;
    F t187 = t68 * t95 + t91 * t72 + t11 + t15 + t7;
    F t192 = t82 * t14 + t76 * t86;
    F t195 = t91 * t14 + t76 * t95;
    F t200 = t82 * t95 + t91 * t86;

    F unknown[9][9];

    unknown[0][0] = t5 * t3 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[0][1] = t24 / 0.60e2;
    unknown[0][2] = t32 / 0.60e2;
    unknown[0][3] = t37 / 0.60e2;
    unknown[0][4] = t65 / 0.60e2;
    unknown[0][5] = t74 / 0.60e2;
    unknown[0][6] = t79 / 0.60e2;
    unknown[0][7] = t88 / 0.60e2;
    unknown[0][8] = t97 / 0.60e2;
    unknown[1][0] = t24 / 0.60e2;
    unknown[1][1] = t18 * t22 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[1][2] = t102 / 0.60e2;
    unknown[1][3] = t105 / 0.60e2;
    unknown[1][4] = t108 / 0.60e2;
    unknown[1][5] = t112 / 0.60e2;
    unknown[1][6] = t115 / 0.60e2;
    unknown[1][7] = t118 / 0.60e2;
    unknown[1][8] = t122 / 0.60e2;
    unknown[2][0] = t32 / 0.60e2;
    unknown[2][1] = t102 / 0.60e2;
    unknown[2][2] = t26 * t30 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[2][3] = t127 / 0.60e2;
    unknown[2][4] = t130 / 0.60e2;
    unknown[2][5] = t133 / 0.60e2;
    unknown[2][6] = t136 / 0.60e2;
    unknown[2][7] = t139 / 0.60e2;
    unknown[2][8] = t142 / 0.60e2;
    unknown[3][0] = t37 / 0.60e2;
    unknown[3][1] = t105 / 0.60e2;
    unknown[3][2] = t127 / 0.60e2;
    unknown[3][3] = t34 * t10 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[3][4] = t147 / 0.60e2;
    unknown[3][5] = t150 / 0.60e2;
    unknown[3][6] = t153 / 0.60e2;
    unknown[3][7] = t157 / 0.60e2;
    unknown[3][8] = t161 / 0.60e2;
    unknown[4][0] = t65 / 0.60e2;
    unknown[4][1] = t108 / 0.60e2;
    unknown[4][2] = t130 / 0.60e2;
    unknown[4][3] = t147 / 0.60e2;
    unknown[4][4] = t59 * t63 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[4][5] = t166 / 0.60e2;
    unknown[4][6] = t169 / 0.60e2;
    unknown[4][7] = t172 / 0.60e2;
    unknown[4][8] = t176 / 0.60e2;
    unknown[5][0] = t74 / 0.60e2;
    unknown[5][1] = t112 / 0.60e2;
    unknown[5][2] = t133 / 0.60e2;
    unknown[5][3] = t150 / 0.60e2;
    unknown[5][4] = t166 / 0.60e2;
    unknown[5][5] = t68 * t72 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[5][6] = t181 / 0.60e2;
    unknown[5][7] = t184 / 0.60e2;
    unknown[5][8] = t187 / 0.60e2;
    unknown[6][0] = t79 / 0.60e2;
    unknown[6][1] = t115 / 0.60e2;
    unknown[6][2] = t136 / 0.60e2;
    unknown[6][3] = t153 / 0.60e2;
    unknown[6][4] = t169 / 0.60e2;
    unknown[6][5] = t181 / 0.60e2;
    unknown[6][6] = t76 * t14 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[6][7] = t192 / 0.60e2;
    unknown[6][8] = t195 / 0.60e2;
    unknown[7][0] = t88 / 0.60e2;
    unknown[7][1] = t118 / 0.60e2;
    unknown[7][2] = t139 / 0.60e2;
    unknown[7][3] = t157 / 0.60e2;
    unknown[7][4] = t172 / 0.60e2;
    unknown[7][5] = t184 / 0.60e2;
    unknown[7][6] = t192 / 0.60e2;
    unknown[7][7] = t82 * t86 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;
    unknown[7][8] = t200 / 0.60e2;
    unknown[8][0] = t97 / 0.60e2;
    unknown[8][1] = t122 / 0.60e2;
    unknown[8][2] = t142 / 0.60e2;
    unknown[8][3] = t161 / 0.60e2;
    unknown[8][4] = t176 / 0.60e2;
    unknown[8][5] = t187 / 0.60e2;
    unknown[8][6] = t195 / 0.60e2;
    unknown[8][7] = t200 / 0.60e2;
    unknown[8][8] = t91 * t95 / 0.30e2 + t11 / 0.30e2 + t15 / 0.30e2 + t7 / 0.30e2;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 9, 9);
    // clang-format on
}
