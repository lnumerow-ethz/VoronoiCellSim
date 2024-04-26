#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellVolume3D.h"
#include "CRLHelper/MapleHelper.h"

void PerCellVolume3D::getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const {
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

    unknown[0] = (y1 * z2 - y2 * z1) * x0 / 0.6e1 + (-y0 * z2 + y2 * z0) * x1 / 0.6e1 + x2 * (y0 * z1 - y1 * z0) / 0.6e1;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.value, 1, 1);
    // clang-format on
}

void PerCellVolume3D::getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const {
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


    F unknown[9];

    unknown[0] = y1 * z2 / 0.6e1 - y2 * z1 / 0.6e1;
    unknown[1] = -z2 * x1 / 0.6e1 + x2 * z1 / 0.6e1;
    unknown[2] = y2 * x1 / 0.6e1 - x2 * y1 / 0.6e1;
    unknown[3] = -y0 * z2 / 0.6e1 + y2 * z0 / 0.6e1;
    unknown[4] = z2 * x0 / 0.6e1 - x2 * z0 / 0.6e1;
    unknown[5] = -y2 * x0 / 0.6e1 + x2 * y0 / 0.6e1;
    unknown[6] = y0 * z1 / 0.6e1 - y1 * z0 / 0.6e1;
    unknown[7] = -z1 * x0 / 0.6e1 + z0 * x1 / 0.6e1;
    unknown[8] = y1 * x0 / 0.6e1 - y0 * x1 / 0.6e1;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.gradient, 9, 1);
    // clang-format on
}

void PerCellVolume3D::getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const {
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

    F t1 = z2 / 0.6e1;
    F t2 = y2 / 0.6e1;
    F t3 = z1 / 0.6e1;
    F t4 = y1 / 0.6e1;
    F t5 = x2 / 0.6e1;
    F t6 = x1 / 0.6e1;
    F t7 = z0 / 0.6e1;
    F t8 = y0 / 0.6e1;
    F t9 = x0 / 0.6e1;

    F unknown[9][9];

    unknown[0][0] = 0.0e0;
    unknown[0][1] = 0.0e0;
    unknown[0][2] = 0.0e0;
    unknown[0][3] = 0.0e0;
    unknown[0][4] = t1;
    unknown[0][5] = -t2;
    unknown[0][6] = 0.0e0;
    unknown[0][7] = -t3;
    unknown[0][8] = t4;
    unknown[1][0] = 0.0e0;
    unknown[1][1] = 0.0e0;
    unknown[1][2] = 0.0e0;
    unknown[1][3] = -t1;
    unknown[1][4] = 0.0e0;
    unknown[1][5] = t5;
    unknown[1][6] = t3;
    unknown[1][7] = 0.0e0;
    unknown[1][8] = -t6;
    unknown[2][0] = 0.0e0;
    unknown[2][1] = 0.0e0;
    unknown[2][2] = 0.0e0;
    unknown[2][3] = t2;
    unknown[2][4] = -t5;
    unknown[2][5] = 0.0e0;
    unknown[2][6] = -t4;
    unknown[2][7] = t6;
    unknown[2][8] = 0.0e0;
    unknown[3][0] = 0.0e0;
    unknown[3][1] = -t1;
    unknown[3][2] = t2;
    unknown[3][3] = 0.0e0;
    unknown[3][4] = 0.0e0;
    unknown[3][5] = 0.0e0;
    unknown[3][6] = 0.0e0;
    unknown[3][7] = t7;
    unknown[3][8] = -t8;
    unknown[4][0] = t1;
    unknown[4][1] = 0.0e0;
    unknown[4][2] = -t5;
    unknown[4][3] = 0.0e0;
    unknown[4][4] = 0.0e0;
    unknown[4][5] = 0.0e0;
    unknown[4][6] = -t7;
    unknown[4][7] = 0.0e0;
    unknown[4][8] = t9;
    unknown[5][0] = -t2;
    unknown[5][1] = t5;
    unknown[5][2] = 0.0e0;
    unknown[5][3] = 0.0e0;
    unknown[5][4] = 0.0e0;
    unknown[5][5] = 0.0e0;
    unknown[5][6] = t8;
    unknown[5][7] = -t9;
    unknown[5][8] = 0.0e0;
    unknown[6][0] = 0.0e0;
    unknown[6][1] = t3;
    unknown[6][2] = -t4;
    unknown[6][3] = 0.0e0;
    unknown[6][4] = -t7;
    unknown[6][5] = t8;
    unknown[6][6] = 0.0e0;
    unknown[6][7] = 0.0e0;
    unknown[6][8] = 0.0e0;
    unknown[7][0] = -t3;
    unknown[7][1] = 0.0e0;
    unknown[7][2] = t6;
    unknown[7][3] = t7;
    unknown[7][4] = 0.0e0;
    unknown[7][5] = -t9;
    unknown[7][6] = 0.0e0;
    unknown[7][7] = 0.0e0;
    unknown[7][8] = 0.0e0;
    unknown[8][0] = t4;
    unknown[8][1] = -t6;
    unknown[8][2] = 0.0e0;
    unknown[8][3] = -t8;
    unknown[8][4] = t9;
    unknown[8][5] = 0.0e0;
    unknown[8][6] = 0.0e0;
    unknown[8][7] = 0.0e0;
    unknown[8][8] = 0.0e0;

    processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 9, 9);
    // clang-format on
}
