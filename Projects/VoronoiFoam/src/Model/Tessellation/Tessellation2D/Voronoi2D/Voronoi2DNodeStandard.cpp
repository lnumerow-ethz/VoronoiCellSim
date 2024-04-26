#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Voronoi2D/Voronoi2D.h"
#include "CRLHelper/MapleHelper.h"

void Voronoi2D::computeNodePosStandard(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];
    F x2 = inputs[4];
    F y2 = inputs[5];

    F t1 = y1 - y2;
    F t2 = y0 * y0;
    F t4 = x1 * x1;
    F t5 = x2 * x2;
    F t6 = y1 * y1;
    F t7 = y2 * y2;
    F t8 = -t4 + t5 - t6 + t7;
    F t11 = x0 * x0;
    F t17 = -x1 + x2;

    F unknown[2];

    unknown[0] = 0.1e1 / (0.2e1 * y0 * t17 + 0.2e1 * y1 * (x0 - x2) - 0.2e1 * y2 * (-x1 + x0)) * (t2 * t1 + y0 * t8 + y2 * t6 + y1 * (t11 - t5 - t7) + y2 * (-t11 + t4));
    unknown[1] = 0.1e1 / (0.2e1 * x0 * t1 + 0.2e1 * x1 * (y2 - y0) + 0.2e1 * x2 * (y0 - y1)) * (t11 * t17 - x0 * t8 - x2 * t4 + x1 * (t5 - t2 + t7) + (t2 - t6) * x2);

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.pos, 2, 1);
    // clang-format on
}

void Voronoi2D::computeNodeGradStandard(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];
    F x2 = inputs[4];
    F y2 = inputs[5];

    F t1 = x0 * y1;
    F t2 = x0 * y2;
    F t4 = -x1 + x2;
    F t6 = x0 - x2;
    F t8 = -x1 + x0;
    F t11 = 0.2e1 * y0 * t4 + 0.2e1 * y1 * t6 - 0.2e1 * t8 * y2;
    F t12 = 0.1e1 / t11;
    F t14 = y1 - y2;
    F t15 = y0 * y0;
    F t17 = x1 * x1;
    F t18 = x2 * x2;
    F t19 = y1 * y1;
    F t20 = y2 * y2;
    F t21 = -t17 + t18 - t19 + t20;
    F t24 = x0 * x0;
    F t30 = t11 * t11;
    F t32 = 0.1e1 / t30 * (t15 * t14 + y0 * t21 + y2 * t19 + y1 * (t24 - t18 - t20) + y2 * (-t24 + t17));
    F t41 = x1 * y0;
    F t42 = x1 * y2;
    F t45 = y2 - y0;
    F t51 = 0.2e1 * y1 * y2;
    F t56 = x2 * y0;
    F t57 = x2 * y1;
    F t60 = y0 - y1;
    F t76 = 0.2e1 * x0 * t14 + 0.2e1 * x1 * t45 + 0.2e1 * t60 * x2;
    F t77 = 0.1e1 / t76;
    F t87 = t76 * t76;
    F t89 = 0.1e1 / t87 * (t24 * t4 - x0 * t21 - x2 * t17 + x1 * (t18 - t15 + t20) + (t15 - t19) * x2);
    F t99 = 0.2e1 * x1 * x2;

    F unknown[2][6];

    unknown[0][0] = 0.2e1 * t12 * (t1 - t2) - 0.2e1 * t14 * t32;
    unknown[0][1] = t12 * (0.2e1 * y0 * t14 - t17 + t18 - t19 + t20) - 0.2e1 * t4 * t32;
    unknown[0][2] = 0.2e1 * t12 * (-t41 + t42) - 0.2e1 * t45 * t32;
    unknown[0][3] = t12 * (-0.2e1 * y1 * y0 + t15 - t18 - t20 + t24 + t51) - 0.2e1 * t6 * t32;
    unknown[0][4] = 0.2e1 * t12 * (t56 - t57) - 0.2e1 * t60 * t32;
    unknown[0][5] = t12 * (0.2e1 * y2 * y0 - t15 + t17 + t19 - t24 - t51) + 0.2e1 * t8 * t32;
    unknown[1][0] = t77 * (0.2e1 * x0 * t4 + t17 - t18 + t19 - t20) - 0.2e1 * t14 * t89;
    unknown[1][1] = 0.2e1 * t77 * (-t41 + t56) - 0.2e1 * t4 * t89;
    unknown[1][2] = t77 * (0.2e1 * x1 * x0 - t15 + t18 + t20 - t24 - t99) - 0.2e1 * t45 * t89;
    unknown[1][3] = 0.2e1 * t77 * (t1 - t57) - 0.2e1 * t6 * t89;
    unknown[1][4] = t77 * (-0.2e1 * x2 * x0 + t15 - t17 - t19 + t24 + t99) - 0.2e1 * t60 * t89;
    unknown[1][5] = 0.2e1 * t77 * (-t2 + t42) + 0.2e1 * t8 * t89;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.grad, 2, 6);
    // clang-format on
}

void Voronoi2D::computeNodeHessStandard(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F x1 = inputs[2];
    F y1 = inputs[3];
    F x2 = inputs[4];
    F y2 = inputs[5];

    F t1 = y1 - y2;
    F t2 = -x1 + x2;
    F t4 = x0 - x2;
    F t6 = -x1 + x0;
    F t9 = 0.2e1 * y0 * t2 + 0.2e1 * y1 * t4 - 0.2e1 * t6 * y2;
    F t10 = 0.1e1 / t9;
    F t11 = 0.2e1 * t10 * t1;
    F t12 = x0 * y1;
    F t13 = x0 * y2;
    F t15 = t9 * t9;
    F t16 = 0.1e1 / t15;
    F t17 = 0.2e1 * t16 * (t12 - t13);
    F t20 = y0 * y0;
    F t22 = x1 * x1;
    F t23 = x2 * x2;
    F t24 = y1 * y1;
    F t25 = y2 * y2;
    F t26 = -t22 + t23 - t24 + t25;
    F t29 = x0 * x0;
    F t34 = t20 * t1 + y0 * t26 + y2 * t24 + y1 * (t29 - t23 - t25) + y2 * (-t29 + t22);
    F t37 = 0.1e1 / t9 / t15 * t34;
    F t38 = 0.4e1 * t1 * t1;
    F t46 = t16 * (0.2e1 * y0 * t1 - t22 + t23 - t24 + t25);
    F t48 = 0.4e1 * t2 * t1;
    F t51 = -0.2e1 * t1 * t46 - 0.2e1 * t2 * t17 + 0.2e1 * t48 * t37;
    F t52 = y2 - y0;
    F t54 = x1 * y0;
    F t55 = x1 * y2;
    F t57 = 0.2e1 * t16 * (-t54 + t55);
    F t59 = 0.4e1 * t52 * t1;
    F t62 = -0.2e1 * t1 * t57 - 0.2e1 * t52 * t17 + 0.2e1 * t59 * t37;
    F t64 = 0.2e1 * t10 * x0;
    F t69 = 0.2e1 * y1 * y2;
    F t71 = t16 * (-0.2e1 * y1 * y0 + t20 - t23 - t25 + t29 + t69);
    F t73 = 0.4e1 * t4 * t1;
    F t77 = 0.2e1 * t16 * t34;
    F t78 = -0.2e1 * t1 * t71 - 0.2e1 * t4 * t17 + 0.2e1 * t73 * t37 + t64 - t77;
    F t79 = y0 - y1;
    F t81 = x2 * y0;
    F t82 = x2 * y1;
    F t84 = 0.2e1 * t16 * (t81 - t82);
    F t86 = 0.4e1 * t79 * t1;
    F t89 = -0.2e1 * t1 * t84 - 0.2e1 * t79 * t17 + 0.2e1 * t86 * t37;
    F t94 = t16 * (0.2e1 * y0 * y2 - t20 + t22 + t24 - t29 - t69);
    F t96 = -0.4e1 * t6 * t1;
    F t99 = -0.2e1 * t1 * t94 + 0.2e1 * t6 * t17 + 0.2e1 * t96 * t37 - t64 + t77;
    F t102 = 0.4e1 * t2 * t2;
    F t107 = 0.2e1 * t10 * x1;
    F t110 = 0.4e1 * t52 * t2;
    F t113 = 0.2e1 * t110 * t37 - 0.2e1 * t2 * t57 - 0.2e1 * t52 * t46 - t107 + t77;
    F t114 = 0.2e1 * t10 * t79;
    F t117 = 0.4e1 * t4 * t2;
    F t120 = 0.2e1 * t117 * t37 - 0.2e1 * t2 * t71 - 0.2e1 * t4 * t46 + t114;
    F t122 = 0.2e1 * t10 * x2;
    F t125 = 0.4e1 * t79 * t2;
    F t128 = 0.2e1 * t125 * t37 - 0.2e1 * t2 * t84 - 0.2e1 * t79 * t46 + t122 - t77;
    F t129 = 0.2e1 * t10 * t52;
    F t132 = -0.4e1 * t6 * t2;
    F t135 = 0.2e1 * t132 * t37 - 0.2e1 * t2 * t94 + 0.2e1 * t6 * t46 + t129;
    F t138 = 0.4e1 * t52 * t52;
    F t144 = 0.4e1 * t4 * t52;
    F t147 = 0.2e1 * t144 * t37 - 0.2e1 * t4 * t57 - 0.2e1 * t52 * t71;
    F t150 = 0.4e1 * t79 * t52;
    F t153 = 0.2e1 * t150 * t37 - 0.2e1 * t52 * t84 - 0.2e1 * t79 * t57;
    F t156 = -0.4e1 * t6 * t52;
    F t159 = 0.2e1 * t156 * t37 - 0.2e1 * t52 * t94 + 0.2e1 * t6 * t57 + t107 - t77;
    F t162 = 0.4e1 * t4 * t4;
    F t168 = 0.4e1 * t79 * t4;
    F t171 = 0.2e1 * t168 * t37 - 0.2e1 * t4 * t84 - 0.2e1 * t79 * t71 - t122 + t77;
    F t174 = -0.4e1 * t6 * t4;
    F t177 = 0.2e1 * t174 * t37 - 0.2e1 * t4 * t94 + 0.2e1 * t6 * t71 + t11;
    F t180 = 0.4e1 * t79 * t79;
    F t186 = -0.4e1 * t6 * t79;
    F t189 = 0.2e1 * t186 * t37 + 0.2e1 * t6 * t84 - 0.2e1 * t79 * t94;
    F t192 = 0.4e1 * t6 * t6;
    F t200 = 0.2e1 * x0 * t1 + 0.2e1 * x1 * t52 + 0.2e1 * t79 * x2;
    F t201 = 0.1e1 / t200;
    F t202 = 0.2e1 * t201 * t2;
    F t206 = t200 * t200;
    F t207 = 0.1e1 / t206;
    F t208 = t207 * (0.2e1 * x0 * t2 + t22 - t23 + t24 - t25);
    F t218 = t29 * t2 - x0 * t26 - x2 * t22 + x1 * (t23 - t20 + t25) + (t20 - t24) * x2;
    F t221 = 0.1e1 / t200 / t206 * t218;
    F t227 = 0.2e1 * t207 * (-t54 + t81);
    F t231 = -0.2e1 * t1 * t227 - 0.2e1 * t2 * t208 + 0.2e1 * t48 * t221;
    F t232 = -0.2e1 * t201 * t6;
    F t237 = 0.2e1 * x1 * x2;
    F t239 = t207 * (0.2e1 * x0 * x1 - t20 + t23 - t237 + t25 - t29);
    F t243 = -0.2e1 * t1 * t239 - 0.2e1 * t52 * t208 + 0.2e1 * t59 * t221 + t232;
    F t245 = 0.2e1 * t201 * y1;
    F t248 = 0.2e1 * t207 * (t12 - t82);
    F t253 = 0.2e1 * t207 * t218;
    F t254 = -0.2e1 * t1 * t248 - 0.2e1 * t4 * t208 + 0.2e1 * t73 * t221 + t245 - t253;
    F t255 = 0.2e1 * t201 * t4;
    F t260 = t207 * (-0.2e1 * x0 * x2 + t20 - t22 + t237 - t24 + t29);
    F t264 = -0.2e1 * t1 * t260 - 0.2e1 * t79 * t208 + 0.2e1 * t86 * t221 + t255;
    F t266 = 0.2e1 * t201 * y2;
    F t269 = 0.2e1 * t207 * (-t13 + t55);
    F t273 = -0.2e1 * t1 * t269 + 0.2e1 * t6 * t208 + 0.2e1 * t96 * t221 + t253 - t266;
    F t280 = 0.2e1 * t201 * y0;
    F t285 = 0.2e1 * t110 * t221 - 0.2e1 * t2 * t239 - 0.2e1 * t52 * t227 + t253 - t280;
    F t290 = 0.2e1 * t117 * t221 - 0.2e1 * t2 * t248 - 0.2e1 * t4 * t227;
    F t295 = 0.2e1 * t125 * t221 - 0.2e1 * t2 * t260 - 0.2e1 * t79 * t227 - t253 + t280;
    F t300 = 0.2e1 * t132 * t221 - 0.2e1 * t2 * t269 + 0.2e1 * t6 * t227;
    F t310 = 0.2e1 * t144 * t221 - 0.2e1 * t4 * t239 - 0.2e1 * t52 * t248;
    F t315 = 0.2e1 * t150 * t221 - 0.2e1 * t79 * t239 - 0.2e1 * t52 * t260 + t202;
    F t320 = 0.2e1 * t156 * t221 + 0.2e1 * t6 * t239 - 0.2e1 * t52 * t269 - t253 + t266;
    F t330 = 0.2e1 * t168 * t221 - 0.2e1 * t79 * t248 - 0.2e1 * t4 * t260 - t245 + t253;
    F t335 = 0.2e1 * t174 * t221 + 0.2e1 * t6 * t248 - 0.2e1 * t4 * t269;
    F t345 = 0.2e1 * t186 * t221 + 0.2e1 * t6 * t260 - 0.2e1 * t79 * t269;

    F unknown[12][6];

    unknown[0][0] = -0.4e1 * t1 * t17 + 0.2e1 * t38 * t37 + t11;
    unknown[0][1] = t51;
    unknown[0][2] = t62;
    unknown[0][3] = t78;
    unknown[0][4] = t89;
    unknown[0][5] = t99;
    unknown[1][0] = t51;
    unknown[1][1] = 0.2e1 * t102 * t37 - 0.4e1 * t2 * t46 + t11;
    unknown[1][2] = t113;
    unknown[1][3] = t120;
    unknown[1][4] = t128;
    unknown[1][5] = t135;
    unknown[2][0] = t62;
    unknown[2][1] = t113;
    unknown[2][2] = 0.2e1 * t138 * t37 - 0.4e1 * t52 * t57 + t129;
    unknown[2][3] = t147;
    unknown[2][4] = t153;
    unknown[2][5] = t159;
    unknown[3][0] = t78;
    unknown[3][1] = t120;
    unknown[3][2] = t147;
    unknown[3][3] = 0.2e1 * t162 * t37 - 0.4e1 * t4 * t71 + t129;
    unknown[3][4] = t171;
    unknown[3][5] = t177;
    unknown[4][0] = t89;
    unknown[4][1] = t128;
    unknown[4][2] = t153;
    unknown[4][3] = t171;
    unknown[4][4] = 0.2e1 * t180 * t37 - 0.4e1 * t79 * t84 + t114;
    unknown[4][5] = t189;
    unknown[5][0] = t99;
    unknown[5][1] = t135;
    unknown[5][2] = t159;
    unknown[5][3] = t177;
    unknown[5][4] = t189;
    unknown[5][5] = 0.2e1 * t192 * t37 + 0.4e1 * t6 * t94 + t114;
    unknown[6][0] = -0.4e1 * t1 * t208 + 0.2e1 * t38 * t221 + t202;
    unknown[6][1] = t231;
    unknown[6][2] = t243;
    unknown[6][3] = t254;
    unknown[6][4] = t264;
    unknown[6][5] = t273;
    unknown[7][0] = t231;
    unknown[7][1] = 0.2e1 * t102 * t221 - 0.4e1 * t2 * t227 + t202;
    unknown[7][2] = t285;
    unknown[7][3] = t290;
    unknown[7][4] = t295;
    unknown[7][5] = t300;
    unknown[8][0] = t243;
    unknown[8][1] = t285;
    unknown[8][2] = 0.2e1 * t138 * t221 - 0.4e1 * t52 * t239 + t255;
    unknown[8][3] = t310;
    unknown[8][4] = t315;
    unknown[8][5] = t320;
    unknown[9][0] = t254;
    unknown[9][1] = t290;
    unknown[9][2] = t310;
    unknown[9][3] = 0.2e1 * t162 * t221 - 0.4e1 * t4 * t248 + t255;
    unknown[9][4] = t330;
    unknown[9][5] = t335;
    unknown[10][0] = t264;
    unknown[10][1] = t295;
    unknown[10][2] = t315;
    unknown[10][3] = t330;
    unknown[10][4] = 0.2e1 * t180 * t221 - 0.4e1 * t79 * t260 + t232;
    unknown[10][5] = t345;
    unknown[11][0] = t273;
    unknown[11][1] = t300;
    unknown[11][2] = t320;
    unknown[11][3] = t335;
    unknown[11][4] = t345;
    unknown[11][5] = 0.2e1 * t192 * t221 + 0.4e1 * t6 * t269 + t232;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.hess, 12, 6);
    // clang-format on
}
