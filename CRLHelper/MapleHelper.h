#pragma once

#include "CRLHelper/VecMatDef.h"

void processMapleOutput(const F *in, F &out, int rows, int cols);

void processMapleOutput(const F *in, VectorXF &out, int rows, int cols);

void processMapleOutput(const F *in, MatrixXF &out, int rows, int cols);

void processMapleOutput(const F *in, std::vector<MatrixXF> &out, int rows, int cols);
