#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "CRLHelper/VecMatDef.h"

struct HessianF {
    SparseMatrixF A;
    MatrixXF U;
    MatrixXF V;

    HessianF(int rows) { setZero(rows); }

    HessianF() : HessianF(0) {}

    void setZero(int rows);

    static HessianF outerProduct(const VectorXF &u, const VectorXF &v, F mul = 1.0);

    MatrixXF evalDense() const;

    HessianF operator+(const HessianF &other) const;

    HessianF operator-(const HessianF &other) const;

    HessianF operator*(const F &other) const;

    HessianF operator/(const F &other) const;
};
