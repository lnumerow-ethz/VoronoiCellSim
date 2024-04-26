#include "CRLHelper/Hessian.h"
#include "CRLHelper/EigenHelper.h"

void HessianF::setZero(int rows) {
    A.resize(rows, rows);
    A.setZero();
    U.resize(rows, 0);
    V.resize(rows, 0);
}

HessianF HessianF::outerProduct(const VectorXF &u, const VectorXF &v, F mul) {
    HessianF result(u.rows());
    result.U = u * mul;
    result.V = v;
    return result;
}

MatrixXF HessianF::evalDense() const {
    return A + U * V.transpose();
}

HessianF HessianF::operator+(const HessianF &other) const {
    HessianF result(A.rows());
    result.A = A + other.A;
    result.U = EigenHelper::horzCat(U, other.U);
    result.V = EigenHelper::horzCat(V, other.V);
    return result;
}

HessianF HessianF::operator-(const HessianF &other) const {
    HessianF result(A.rows());
    result.A = A - other.A;
    result.U = EigenHelper::horzCat(U, other.U);
    result.V = EigenHelper::horzCat(V, (MatrixXF)-other.V);
    return result;
}

HessianF HessianF::operator*(const F &other) const {
    HessianF result(A.rows());
    result.A = A * other;
    result.U = U * other;
    result.V = V;
    return result;
}

HessianF HessianF::operator/(const F &other) const {
    HessianF result(A.rows());
    result.A = A / other;
    result.U = U / other;
    result.V = V;
    return result;
}
