#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Power2D/Power2D.h"
#include "CRLHelper/MapleHelper.h"

void Power2D::computeNodePosBEdge(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F xb0 = inputs[0];
    F yb0 = inputs[1];
    F xb1 = inputs[2];
    F yb1 = inputs[3];
    F x0 = inputs[4];
    F y0 = inputs[5];
    F w0 = inputs[6];
    F x1 = inputs[7];
    F y1 = inputs[8];
    F w1 = inputs[9];

    F t1 = x0 * x0;
    F t2 = x1 * x1;
    F t3 = y0 * y0;
    F t6 = y1 * y1;
    F t18 = x1 - x0;
    F t21 = -y1 + y0;

    F unknown[2];

    unknown[0] = 0.1e1 / (0.2e1 * xb0 * t18 - 0.2e1 * xb1 * t18 + 0.2e1 * (yb1 - yb0) * t21) * (xb0 * (0.2e1 * y0 * yb1 - 0.2e1 * y1 * yb1 - t1 + t2 - t3 + t6 + w0 - w1) - (0.2e1 * y0 * yb0 - 0.2e1 * y1 * yb0 - t1 + t2 - t3 + t6 + w0 - w1) * xb1);
    unknown[1] = 0.1e1 / (-0.2e1 * yb0 * t21 + 0.2e1 * yb1 * t21 + 0.2e1 * t18 * (-xb1 + xb0)) * (yb0 * (0.2e1 * x0 * xb1 - 0.2e1 * x1 * xb1 - t1 + t2 - t3 + t6 + w0 - w1) - (0.2e1 * x0 * xb0 - 0.2e1 * x1 * xb0 - t1 + t2 - t3 + t6 + w0 - w1) * yb1);

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.pos, 2, 1);
    // clang-format on
}

void Power2D::computeNodeGradBEdge(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F xb0 = inputs[0];
    F yb0 = inputs[1];
    F xb1 = inputs[2];
    F yb1 = inputs[3];
    F x0 = inputs[4];
    F y0 = inputs[5];
    F w0 = inputs[6];
    F x1 = inputs[7];
    F y1 = inputs[8];
    F w1 = inputs[9];

    F t1 = x0 * x0;
    F t2 = x1 * x1;
    F t3 = y0 * y0;
    F t4 = y0 * yb1;
    F t6 = y1 * y1;
    F t7 = y1 * yb1;
    F t9 = -t1 + t2 - t3 + 0.2e1 * t4 + t6 - 0.2e1 * t7 + w0 - w1;
    F t10 = x1 - x0;
    F t13 = -y1 + y0;
    F t14 = yb1 - yb0;
    F t17 = 0.2e1 * xb0 * t10 - 0.2e1 * xb1 * t10 + 0.2e1 * t14 * t13;
    F t18 = 0.1e1 / t17;
    F t21 = y0 * yb0;
    F t23 = y1 * yb0;
    F t25 = -t1 + t2 - t3 + 0.2e1 * t21 + t6 - 0.2e1 * t23 + w0 - w1;
    F t28 = t17 * t17;
    F t30 = 0.1e1 / t28 * (-t25 * xb1 + xb0 * t9);
    F t44 = x0 * xb0;
    F t45 = x0 * xb1;
    F t48 = xb1 - xb0;
    F t60 = x1 * xb0;
    F t61 = x1 * xb1;
    F t80 = -0.2e1 * t10 * t48 - 0.2e1 * yb0 * t13 + 0.2e1 * yb1 * t13;
    F t81 = 0.1e1 / t80;
    F t85 = -t1 + 0.2e1 * t45 + t2 - 0.2e1 * t61 - t3 + t6 + w0 - w1;
    F t89 = -t1 + 0.2e1 * t44 + t2 - 0.2e1 * t60 - t3 + t6 + w0 - w1;
    F t92 = t80 * t80;
    F t94 = 0.1e1 / t92 * (yb0 * t85 - t89 * yb1);

    F unknown[2][10];

    unknown[0][0] = -0.2e1 * t10 * t30 + t18 * t9;
    unknown[0][1] = -0.2e1 * t18 * t13 * xb1 + 0.2e1 * t13 * t30;
    unknown[0][2] = 0.2e1 * t10 * t30 - t18 * t25;
    unknown[0][3] = 0.2e1 * t18 * xb0 * t13 - 0.2e1 * t13 * t30;
    unknown[0][4] = 0.2e1 * t18 * (-t44 + t45) - 0.2e1 * t48 * t30;
    unknown[0][5] = t18 * (0.2e1 * xb0 * (-y0 + yb1) - 0.2e1 * (-y0 + yb0) * xb1) - 0.2e1 * t14 * t30;
    unknown[0][6] = -t18 * t48;
    unknown[0][7] = 0.2e1 * t18 * (t60 - t61) + 0.2e1 * t48 * t30;
    unknown[0][8] = t18 * (0.2e1 * xb0 * (y1 - yb1) - 0.2e1 * (y1 - yb0) * xb1) + 0.2e1 * t14 * t30;
    unknown[0][9] = t18 * t48;
    unknown[1][0] = 0.2e1 * t81 * t10 * yb1 - 0.2e1 * t10 * t94;
    unknown[1][1] = 0.2e1 * t13 * t94 + t81 * t85;
    unknown[1][2] = -0.2e1 * t81 * yb0 * t10 + 0.2e1 * t10 * t94;
    unknown[1][3] = -0.2e1 * t13 * t94 - t81 * t89;
    unknown[1][4] = t81 * (0.2e1 * yb0 * (-x0 + xb1) - 0.2e1 * (-x0 + xb0) * yb1) - 0.2e1 * t48 * t94;
    unknown[1][5] = 0.2e1 * t81 * (-t21 + t4) - 0.2e1 * t14 * t94;
    unknown[1][6] = -t81 * t14;
    unknown[1][7] = t81 * (0.2e1 * yb0 * (x1 - xb1) - 0.2e1 * (x1 - xb0) * yb1) + 0.2e1 * t48 * t94;
    unknown[1][8] = 0.2e1 * t81 * (t23 - t7) + 0.2e1 * t14 * t94;
    unknown[1][9] = t81 * t14;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.grad, 2, 10);
    // clang-format on
}

void Power2D::computeNodeHessBEdge(const VectorXF &inputs, NodeData &node_data) const {
    // clang-format off
    F xb0 = inputs[0];
    F yb0 = inputs[1];
    F xb1 = inputs[2];
    F yb1 = inputs[3];
    F x0 = inputs[4];
    F y0 = inputs[5];
    F w0 = inputs[6];
    F x1 = inputs[7];
    F y1 = inputs[8];
    F w1 = inputs[9];

    F t1 = x0 * x0;
    F t2 = x1 * x1;
    F t3 = y0 * y0;
    F t4 = y0 * yb1;
    F t6 = y1 * y1;
    F t7 = y1 * yb1;
    F t9 = -t1 + t2 - t3 + 0.2e1 * t4 + t6 - 0.2e1 * t7 + w0 - w1;
    F t10 = x1 - x0;
    F t13 = -y1 + y0;
    F t14 = yb1 - yb0;
    F t17 = 0.2e1 * xb0 * t10 - 0.2e1 * xb1 * t10 + 0.2e1 * t14 * t13;
    F t18 = t17 * t17;
    F t19 = 0.1e1 / t18;
    F t20 = t19 * t9;
    F t23 = y0 * yb0;
    F t25 = y1 * yb0;
    F t27 = -t1 + t2 - t3 + 0.2e1 * t23 + t6 - 0.2e1 * t25 + w0 - w1;
    F t29 = -t27 * xb1 + xb0 * t9;
    F t32 = 0.1e1 / t17 / t18 * t29;
    F t33 = 0.4e1 * t10 * t10;
    F t34 = t33 * t32;
    F t37 = 0.2e1 * t13 * xb1;
    F t38 = 0.2e1 * t10 * t19;
    F t40 = -0.4e1 * t13 * t10;
    F t42 = 0.2e1 * t40 * t32;
    F t43 = 0.2e1 * t13 * t20 + t38 * t37 + t42;
    F t45 = -t19 * t27;
    F t47 = -0.4e1 * t10 * t10;
    F t50 = 0.2e1 * t10 * t20 - 0.2e1 * t10 * t45 + 0.2e1 * t47 * t32;
    F t51 = 0.1e1 / t17;
    F t52 = 0.2e1 * t51 * t13;
    F t54 = 0.2e1 * xb0 * t13;
    F t56 = 0.4e1 * t13 * t10;
    F t58 = 0.2e1 * t56 * t32;
    F t59 = -0.2e1 * t13 * t20 - t38 * t54 + t52 + t58;
    F t61 = 0.2e1 * t51 * x0;
    F t62 = xb1 - xb0;
    F t64 = x0 * xb0;
    F t65 = x0 * xb1;
    F t67 = 0.2e1 * t19 * (-t64 + t65);
    F t69 = 0.4e1 * t62 * t10;
    F t71 = 0.2e1 * t69 * t32;
    F t73 = 0.2e1 * t19 * t29;
    F t74 = -0.2e1 * t10 * t67 - 0.2e1 * t62 * t20 - t61 + t71 + t73;
    F t75 = -y0 + yb1;
    F t79 = -y0 + yb0;
    F t82 = t19 * (0.2e1 * xb0 * t75 - 0.2e1 * t79 * xb1);
    F t84 = 0.4e1 * t14 * t10;
    F t86 = 0.2e1 * t84 * t32;
    F t87 = -0.2e1 * t10 * t82 - 0.2e1 * t14 * t20 + 0.2e1 * t51 * t75 + t86;
    F t88 = -t19 * t62;
    F t90 = -0.2e1 * t10 * t88 + t51;
    F t92 = 0.2e1 * t51 * x1;
    F t94 = x1 * xb0;
    F t95 = x1 * xb1;
    F t97 = 0.2e1 * t19 * (t94 - t95);
    F t99 = -0.4e1 * t62 * t10;
    F t101 = 0.2e1 * t99 * t32;
    F t102 = -0.2e1 * t10 * t97 + 0.2e1 * t62 * t20 + t101 - t73 + t92;
    F t103 = y1 - yb1;
    F t107 = y1 - yb0;
    F t110 = t19 * (0.2e1 * xb0 * t103 - 0.2e1 * t107 * xb1);
    F t112 = -0.4e1 * t14 * t10;
    F t114 = 0.2e1 * t112 * t32;
    F t115 = -0.2e1 * t10 * t110 + 0.2e1 * t51 * t103 + 0.2e1 * t14 * t20 + t114;
    F t116 = t19 * t62;
    F t118 = -0.2e1 * t10 * t116 - t51;
    F t119 = -0.2e1 * t13 * t19;
    F t121 = 0.4e1 * t13 * t13;
    F t122 = t121 * t32;
    F t124 = -0.2e1 * t10 * t19;
    F t127 = t124 * t37 + 0.2e1 * t13 * t45 - t52 + t58;
    F t131 = -0.4e1 * t13 * t13;
    F t134 = t19 * t121 * xb1 - t119 * t54 + 0.2e1 * t131 * t32;
    F t135 = 0.2e1 * t19 * t62;
    F t138 = -0.4e1 * t62 * t13;
    F t140 = 0.2e1 * t138 * t32;
    F t141 = 0.2e1 * t13 * t67 + t135 * t37 + t140;
    F t143 = 0.2e1 * t51 * xb1;
    F t144 = 0.2e1 * t14 * t19;
    F t147 = -0.4e1 * t14 * t13;
    F t149 = 0.2e1 * t147 * t32;
    F t150 = 0.2e1 * t13 * t82 + t144 * t37 - t143 + t149 + t73;
    F t151 = -0.2e1 * t13 * t88;
    F t152 = -0.2e1 * t19 * t62;
    F t155 = 0.4e1 * t62 * t13;
    F t157 = 0.2e1 * t155 * t32;
    F t158 = 0.2e1 * t13 * t97 + t152 * t37 + t157;
    F t159 = -0.2e1 * t14 * t19;
    F t162 = 0.4e1 * t14 * t13;
    F t164 = 0.2e1 * t162 * t32;
    F t165 = 0.2e1 * t13 * t110 + t159 * t37 + t143 + t164 - t73;
    F t166 = -0.2e1 * t13 * t116;
    F t171 = -t124 * t54 - 0.2e1 * t13 * t45 + t42;
    F t174 = 0.2e1 * t10 * t67 - 0.2e1 * t62 * t45 + t101 + t61 - t73;
    F t178 = 0.2e1 * t10 * t82 - 0.2e1 * t14 * t45 - 0.2e1 * t51 * t79 + t114;
    F t180 = 0.2e1 * t10 * t88 - t51;
    F t183 = 0.2e1 * t10 * t97 + 0.2e1 * t62 * t45 + t71 + t73 - t92;
    F t187 = 0.2e1 * t10 * t110 - 0.2e1 * t51 * t107 + 0.2e1 * t14 * t45 + t86;
    F t189 = 0.2e1 * t10 * t116 + t51;
    F t195 = -0.2e1 * t13 * t67 - t135 * t54 + t157;
    F t197 = 0.2e1 * t51 * xb0;
    F t200 = -0.2e1 * t13 * t82 - t144 * t54 + t164 + t197 - t73;
    F t201 = 0.2e1 * t13 * t88;
    F t204 = -0.2e1 * t13 * t97 - t152 * t54 + t140;
    F t207 = -0.2e1 * t13 * t110 - t159 * t54 + t149 - t197 + t73;
    F t208 = 0.2e1 * t13 * t116;
    F t209 = 0.2e1 * t51 * t62;
    F t212 = 0.4e1 * t62 * t62;
    F t214 = 0.2e1 * t212 * t32;
    F t218 = 0.4e1 * t14 * t62;
    F t220 = 0.2e1 * t218 * t32;
    F t221 = -0.2e1 * t14 * t67 - 0.2e1 * t62 * t82 + t220;
    F t222 = 0.2e1 * t62 * t88;
    F t225 = -0.4e1 * t62 * t62;
    F t228 = 0.2e1 * t225 * t32 + 0.2e1 * t62 * t67 - 0.2e1 * t62 * t97;
    F t231 = -0.4e1 * t14 * t62;
    F t233 = 0.2e1 * t231 * t32;
    F t234 = -0.2e1 * t62 * t110 + 0.2e1 * t14 * t67 + t233;
    F t235 = 0.2e1 * t62 * t116;
    F t238 = 0.4e1 * t14 * t14;
    F t240 = 0.2e1 * t238 * t32;
    F t242 = 0.2e1 * t14 * t88;
    F t245 = -0.2e1 * t14 * t97 + 0.2e1 * t62 * t82 + t233;
    F t248 = -0.4e1 * t14 * t14;
    F t251 = -0.2e1 * t14 * t110 + 0.2e1 * t14 * t82 + 0.2e1 * t248 * t32;
    F t252 = 0.2e1 * t14 * t116;
    F t253 = -0.2e1 * t62 * t88;
    F t254 = -0.2e1 * t14 * t88;
    F t255 = -0.2e1 * t51 * t62;
    F t261 = 0.2e1 * t62 * t110 + 0.2e1 * t14 * t97 + t220;
    F t262 = -0.2e1 * t62 * t116;
    F t266 = -0.2e1 * t14 * t116;
    F t267 = -0.2e1 * t10 * yb1;
    F t272 = -0.2e1 * t62 * t10 - 0.2e1 * yb0 * t13 + 0.2e1 * yb1 * t13;
    F t273 = t272 * t272;
    F t274 = 0.1e1 / t273;
    F t275 = 0.2e1 * t10 * t274;
    F t279 = -t1 + 0.2e1 * t65 + t2 - 0.2e1 * t95 - t3 + t6 + w0 - w1;
    F t283 = -t1 + 0.2e1 * t64 + t2 - 0.2e1 * t94 - t3 + t6 + w0 - w1;
    F t285 = yb0 * t279 - t283 * yb1;
    F t288 = 0.1e1 / t272 / t273 * t285;
    F t289 = t33 * t288;
    F t291 = -0.2e1 * t13 * t274;
    F t293 = t274 * t279;
    F t296 = 0.2e1 * t40 * t288;
    F t297 = -0.2e1 * t10 * t293 + t291 * t267 + t296;
    F t300 = -0.2e1 * yb0 * t10;
    F t304 = t274 * t33 * yb1 - t275 * t300 + 0.2e1 * t47 * t288;
    F t305 = 0.1e1 / t272;
    F t306 = -0.2e1 * t305 * t10;
    F t307 = 0.2e1 * t13 * t274;
    F t309 = -t274 * t283;
    F t312 = 0.2e1 * t56 * t288;
    F t313 = -0.2e1 * t10 * t309 + t307 * t267 - t306 + t312;
    F t315 = 0.2e1 * t305 * yb1;
    F t316 = 0.2e1 * t62 * t274;
    F t318 = -x0 + xb1;
    F t320 = -x0 + xb0;
    F t323 = t274 * (0.2e1 * yb0 * t318 - 0.2e1 * t320 * yb1);
    F t326 = 0.2e1 * t69 * t288;
    F t328 = 0.2e1 * t274 * t285;
    F t329 = -0.2e1 * t10 * t323 + t316 * t267 - t315 + t326 + t328;
    F t330 = 0.2e1 * t14 * t274;
    F t333 = 0.2e1 * t274 * (-t23 + t4);
    F t336 = 0.2e1 * t84 * t288;
    F t337 = -0.2e1 * t10 * t333 + t330 * t267 + t336;
    F t338 = -t14 * t274;
    F t339 = 0.2e1 * t10 * t338;
    F t340 = -0.2e1 * t62 * t274;
    F t342 = x1 - xb1;
    F t344 = x1 - xb0;
    F t347 = t274 * (0.2e1 * yb0 * t342 - 0.2e1 * t344 * yb1);
    F t350 = 0.2e1 * t99 * t288;
    F t351 = -0.2e1 * t10 * t347 + t340 * t267 + t315 - t328 + t350;
    F t352 = -0.2e1 * t14 * t274;
    F t355 = 0.2e1 * t274 * (t25 - t7);
    F t358 = 0.2e1 * t112 * t288;
    F t359 = -0.2e1 * t10 * t355 + t352 * t267 + t358;
    F t360 = t14 * t274;
    F t361 = 0.2e1 * t10 * t360;
    F t363 = t121 * t288;
    F t367 = 0.2e1 * t10 * t293 - t291 * t300 + t306 + t312;
    F t372 = -0.2e1 * t13 * t293 + 0.2e1 * t13 * t309 + 0.2e1 * t131 * t288;
    F t377 = 0.2e1 * t138 * t288;
    F t378 = 0.2e1 * t13 * t323 - 0.2e1 * t62 * t293 + 0.2e1 * t305 * t318 + t377;
    F t380 = 0.2e1 * t305 * y0;
    F t384 = 0.2e1 * t147 * t288;
    F t385 = 0.2e1 * t13 * t333 - 0.2e1 * t14 * t293 + t328 - t380 + t384;
    F t387 = 0.2e1 * t13 * t338 + t305;
    F t392 = 0.2e1 * t155 * t288;
    F t393 = 0.2e1 * t13 * t347 + 0.2e1 * t62 * t293 + 0.2e1 * t305 * t342 + t392;
    F t395 = 0.2e1 * t305 * y1;
    F t399 = 0.2e1 * t162 * t288;
    F t400 = 0.2e1 * t13 * t355 + 0.2e1 * t14 * t293 - t328 + t395 + t399;
    F t402 = 0.2e1 * t13 * t360 - t305;
    F t408 = 0.2e1 * t10 * t309 - t307 * t300 + t296;
    F t410 = 0.2e1 * t305 * yb0;
    F t413 = 0.2e1 * t10 * t323 - t316 * t300 - t328 + t350 + t410;
    F t416 = 0.2e1 * t10 * t333 - t330 * t300 + t358;
    F t417 = -0.2e1 * t10 * t338;
    F t420 = 0.2e1 * t10 * t347 - t340 * t300 + t326 + t328 - t410;
    F t423 = 0.2e1 * t10 * t355 - t352 * t300 + t336;
    F t424 = -0.2e1 * t10 * t360;
    F t430 = -0.2e1 * t13 * t323 - 0.2e1 * t305 * t320 - 0.2e1 * t62 * t309 + t392;
    F t433 = -0.2e1 * t13 * t333 - 0.2e1 * t14 * t309 - t328 + t380 + t399;
    F t435 = -0.2e1 * t13 * t338 - t305;
    F t439 = -0.2e1 * t13 * t347 - 0.2e1 * t305 * t344 + 0.2e1 * t62 * t309 + t377;
    F t442 = -0.2e1 * t13 * t355 + 0.2e1 * t14 * t309 + t328 + t384 - t395;
    F t444 = -0.2e1 * t13 * t360 + t305;
    F t445 = 0.2e1 * t305 * t14;
    F t449 = 0.2e1 * t212 * t288;
    F t454 = 0.2e1 * t218 * t288;
    F t455 = -0.2e1 * t14 * t323 - 0.2e1 * t62 * t333 + t454;
    F t456 = 0.2e1 * t62 * t338;
    F t461 = 0.2e1 * t225 * t288 + 0.2e1 * t62 * t323 - 0.2e1 * t62 * t347;
    F t465 = 0.2e1 * t231 * t288;
    F t466 = 0.2e1 * t14 * t323 - 0.2e1 * t62 * t355 + t465;
    F t467 = 0.2e1 * t62 * t360;
    F t471 = 0.2e1 * t238 * t288;
    F t473 = 0.2e1 * t14 * t338;
    F t476 = -0.2e1 * t14 * t347 + 0.2e1 * t62 * t333 + t465;
    F t481 = 0.2e1 * t14 * t333 - 0.2e1 * t14 * t355 + 0.2e1 * t248 * t288;
    F t482 = 0.2e1 * t14 * t360;
    F t483 = -0.2e1 * t62 * t338;
    F t484 = -0.2e1 * t14 * t338;
    F t485 = -0.2e1 * t305 * t14;
    F t491 = 0.2e1 * t14 * t347 + 0.2e1 * t62 * t355 + t454;
    F t492 = -0.2e1 * t62 * t360;
    F t496 = -0.2e1 * t14 * t360;

    F unknown[20][10];

    unknown[0][0] = -0.4e1 * t10 * t20 + 0.2e1 * t34;
    unknown[0][1] = t43;
    unknown[0][2] = t50;
    unknown[0][3] = t59;
    unknown[0][4] = t74;
    unknown[0][5] = t87;
    unknown[0][6] = t90;
    unknown[0][7] = t102;
    unknown[0][8] = t115;
    unknown[0][9] = t118;
    unknown[1][0] = t43;
    unknown[1][1] = 0.2e1 * t119 * t37 + 0.2e1 * t122;
    unknown[1][2] = t127;
    unknown[1][3] = t134;
    unknown[1][4] = t141;
    unknown[1][5] = t150;
    unknown[1][6] = -t151;
    unknown[1][7] = t158;
    unknown[1][8] = t165;
    unknown[1][9] = -t166;
    unknown[2][0] = t50;
    unknown[2][1] = t127;
    unknown[2][2] = 0.4e1 * t10 * t45 + 0.2e1 * t34;
    unknown[2][3] = t171;
    unknown[2][4] = t174;
    unknown[2][5] = t178;
    unknown[2][6] = t180;
    unknown[2][7] = t183;
    unknown[2][8] = t187;
    unknown[2][9] = t189;
    unknown[3][0] = t59;
    unknown[3][1] = t134;
    unknown[3][2] = t171;
    unknown[3][3] = -0.2e1 * t19 * xb0 * t121 + 0.2e1 * t122;
    unknown[3][4] = t195;
    unknown[3][5] = t200;
    unknown[3][6] = -t201;
    unknown[3][7] = t204;
    unknown[3][8] = t207;
    unknown[3][9] = -t208;
    unknown[4][0] = t74;
    unknown[4][1] = t141;
    unknown[4][2] = t174;
    unknown[4][3] = t195;
    unknown[4][4] = -0.4e1 * t62 * t67 + t209 + t214;
    unknown[4][5] = t221;
    unknown[4][6] = -t222;
    unknown[4][7] = t228;
    unknown[4][8] = t234;
    unknown[4][9] = -t235;
    unknown[5][0] = t87;
    unknown[5][1] = t150;
    unknown[5][2] = t178;
    unknown[5][3] = t200;
    unknown[5][4] = t221;
    unknown[5][5] = -0.4e1 * t14 * t82 + t209 + t240;
    unknown[5][6] = -t242;
    unknown[5][7] = t245;
    unknown[5][8] = t251;
    unknown[5][9] = -t252;
    unknown[6][0] = t90;
    unknown[6][1] = -t151;
    unknown[6][2] = t180;
    unknown[6][3] = -t201;
    unknown[6][4] = -t222;
    unknown[6][5] = -t242;
    unknown[6][6] = 0.0e0;
    unknown[6][7] = -t253;
    unknown[6][8] = -t254;
    unknown[6][9] = 0.0e0;
    unknown[7][0] = t102;
    unknown[7][1] = t158;
    unknown[7][2] = t183;
    unknown[7][3] = t204;
    unknown[7][4] = t228;
    unknown[7][5] = t245;
    unknown[7][6] = -t253;
    unknown[7][7] = 0.4e1 * t62 * t97 + t214 + t255;
    unknown[7][8] = t261;
    unknown[7][9] = -t262;
    unknown[8][0] = t115;
    unknown[8][1] = t165;
    unknown[8][2] = t187;
    unknown[8][3] = t207;
    unknown[8][4] = t234;
    unknown[8][5] = t251;
    unknown[8][6] = -t254;
    unknown[8][7] = t261;
    unknown[8][8] = 0.4e1 * t14 * t110 + t240 + t255;
    unknown[8][9] = -t266;
    unknown[9][0] = t118;
    unknown[9][1] = -t166;
    unknown[9][2] = t189;
    unknown[9][3] = -t208;
    unknown[9][4] = -t235;
    unknown[9][5] = -t252;
    unknown[9][6] = 0.0e0;
    unknown[9][7] = -t262;
    unknown[9][8] = -t266;
    unknown[9][9] = 0.0e0;
    unknown[10][0] = 0.2e1 * t275 * t267 + 0.2e1 * t289;
    unknown[10][1] = t297;
    unknown[10][2] = t304;
    unknown[10][3] = t313;
    unknown[10][4] = t329;
    unknown[10][5] = t337;
    unknown[10][6] = -t339;
    unknown[10][7] = t351;
    unknown[10][8] = t359;
    unknown[10][9] = -t361;
    unknown[11][0] = t297;
    unknown[11][1] = 0.4e1 * t13 * t293 + 0.2e1 * t363;
    unknown[11][2] = t367;
    unknown[11][3] = t372;
    unknown[11][4] = t378;
    unknown[11][5] = t385;
    unknown[11][6] = t387;
    unknown[11][7] = t393;
    unknown[11][8] = t400;
    unknown[11][9] = t402;
    unknown[12][0] = t304;
    unknown[12][1] = t367;
    unknown[12][2] = -0.2e1 * t274 * yb0 * t33 + 0.2e1 * t289;
    unknown[12][3] = t408;
    unknown[12][4] = t413;
    unknown[12][5] = t416;
    unknown[12][6] = -t417;
    unknown[12][7] = t420;
    unknown[12][8] = t423;
    unknown[12][9] = -t424;
    unknown[13][0] = t313;
    unknown[13][1] = t372;
    unknown[13][2] = t408;
    unknown[13][3] = -0.4e1 * t13 * t309 + 0.2e1 * t363;
    unknown[13][4] = t430;
    unknown[13][5] = t433;
    unknown[13][6] = t435;
    unknown[13][7] = t439;
    unknown[13][8] = t442;
    unknown[13][9] = t444;
    unknown[14][0] = t329;
    unknown[14][1] = t378;
    unknown[14][2] = t413;
    unknown[14][3] = t430;
    unknown[14][4] = -0.4e1 * t62 * t323 + t445 + t449;
    unknown[14][5] = t455;
    unknown[14][6] = -t456;
    unknown[14][7] = t461;
    unknown[14][8] = t466;
    unknown[14][9] = -t467;
    unknown[15][0] = t337;
    unknown[15][1] = t385;
    unknown[15][2] = t416;
    unknown[15][3] = t433;
    unknown[15][4] = t455;
    unknown[15][5] = -0.4e1 * t14 * t333 + t445 + t471;
    unknown[15][6] = -t473;
    unknown[15][7] = t476;
    unknown[15][8] = t481;
    unknown[15][9] = -t482;
    unknown[16][0] = -t339;
    unknown[16][1] = t387;
    unknown[16][2] = -t417;
    unknown[16][3] = t435;
    unknown[16][4] = -t456;
    unknown[16][5] = -t473;
    unknown[16][6] = 0.0e0;
    unknown[16][7] = -t483;
    unknown[16][8] = -t484;
    unknown[16][9] = 0.0e0;
    unknown[17][0] = t351;
    unknown[17][1] = t393;
    unknown[17][2] = t420;
    unknown[17][3] = t439;
    unknown[17][4] = t461;
    unknown[17][5] = t476;
    unknown[17][6] = -t483;
    unknown[17][7] = 0.4e1 * t62 * t347 + t449 + t485;
    unknown[17][8] = t491;
    unknown[17][9] = -t492;
    unknown[18][0] = t359;
    unknown[18][1] = t400;
    unknown[18][2] = t423;
    unknown[18][3] = t442;
    unknown[18][4] = t466;
    unknown[18][5] = t481;
    unknown[18][6] = -t484;
    unknown[18][7] = t491;
    unknown[18][8] = 0.4e1 * t14 * t355 + t471 + t485;
    unknown[18][9] = -t496;
    unknown[19][0] = -t361;
    unknown[19][1] = t402;
    unknown[19][2] = -t424;
    unknown[19][3] = t444;
    unknown[19][4] = -t467;
    unknown[19][5] = -t482;
    unknown[19][6] = 0.0e0;
    unknown[19][7] = -t492;
    unknown[19][8] = -t496;
    unknown[19][9] = 0.0e0;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.hess, 20, 10);
    // clang-format on
}
