#pragma once

#include "CRLHelper/VecMatDef.h"

namespace EigenHelper {
template <typename T, int n, int m>
Matrix<T, n, m> horzCat(const Matrix<T, n, m> &A, const Matrix<T, n, m> &B) {
    assert(A.rows() == B.rows());
    Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(A.rows(), A.cols() + B.cols());
    result.leftCols(A.cols()) = A;
    result.rightCols(B.cols()) = B;
    return result;
}

template <typename T, int n, int m>
Matrix<T, n, m> vertCat(const Matrix<T, n, m> &A, const Matrix<T, n, m> &B) {
    assert(A.cols() == B.cols());
    Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(A.rows() + B.rows(), A.cols());
    result.topRows(A.rows()) = A;
    result.bottomRows(B.rows()) = B;
    return result;
}
}  // namespace EigenHelper
