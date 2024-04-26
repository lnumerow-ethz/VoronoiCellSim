#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"

/// Arithmetic operators with for scalar values with derivatives. Types must have member variables:
///     F value;
///     VectorXF gradient;
///     HessianF hessian; // (SparseMatrixF A, MatrixXF U, MatrixXF V)
/// Return value copies any other member variables from first operand. It is assumed that either they are the same for
/// both operands or they do not matter.
namespace GradArithmetic {
template <typename T>
T add(const T &a, const T &b, int order) {
    T result = a;
    result.value = a.value + b.value;
    if (order >= 1) result.gradient = a.gradient + b.gradient;
    if (order >= 2) result.hessian = a.hessian + b.hessian;
    return result;
}

template <typename T>
T subtract(const T &a, const T &b, int order) {
    T result = a;
    result.value = a.value - b.value;
    if (order >= 1) result.gradient = a.gradient - b.gradient;
    if (order >= 2) result.hessian = a.hessian - b.hessian;
    return result;
}

template <typename T>
T multiply(const T &a, const T &b, int order) {
    T result = a;
    result.value = a.value * b.value;
    if (order >= 1) result.gradient = a.gradient * b.value + a.value * b.gradient;
    if (order >= 2)
        result.hessian = a.hessian * b.value + HessianF::outerProduct(a.gradient, b.gradient) +
                         HessianF::outerProduct(b.gradient, a.gradient) + b.hessian * a.value;
    return result;
}

template <typename T>
T divide(const T &a, const T &b, int order) {
    T result = a;
    result.value = a.value / b.value;
    if (order >= 1) result.gradient = a.gradient / b.value - a.value * b.gradient / pow(b.value, 2);
    if (order >= 2)
        result.hessian = a.hessian / b.value - HessianF::outerProduct(a.gradient, b.gradient, 1.0 / pow(b.value, 2)) -
                         HessianF::outerProduct(b.gradient, a.gradient, 1.0 / pow(b.value, 2)) +
                         HessianF::outerProduct(b.gradient, b.gradient, 2 * a.value / pow(b.value, 3)) -
                         b.hessian * a.value / pow(b.value, 2);
    return result;
}

template <typename T>
T add(const T &a, const F &b) {
    T result = a;
    result.value = a.value + b;
    return result;
}

template <typename T>
T subtract(const T &a, const F &b) {
    T result = a;
    result.value = a.value - b;
    return result;
}

template <typename T>
T multiply(const T &a, const F &b, int order) {
    T result = a;
    result.value = a.value * b;
    if (order >= 1) result.gradient = a.gradient * b;
    if (order >= 2) result.hessian = a.hessian * b;
    return result;
}

template <typename T>
T divide(const T &a, const F &b, int order) {
    T result = a;
    result.value = a.value / b;
    if (order >= 1) result.gradient = a.gradient / b;
    if (order >= 2) result.hessian = a.hessian / b;
    return result;
}

template <typename T>
T divide(const F &a, const T &b, int order) {
    T result = b;
    result.value = a / b.value;
    if (order >= 1) result.gradient = -a * b.gradient / pow(b.value, 2);
    if (order >= 2)
        result.hessian =
            HessianF::outerProduct(b.gradient, b.gradient, 2 * a / pow(b.value, 3)) - b.hessian * a / pow(b.value, 2);
    return result;
}

template <typename T>
T powf(const T &a, const F &b, int order) {
    T result = a;
    result.value = pow(a.value, b);
    if (order >= 1) result.gradient = pow(a.value, b - 1.0) * b * a.gradient;
    if (order >= 2)
        result.hessian = HessianF::outerProduct(a.gradient, a.gradient, pow(a.value, b - 2.0) * b * (b - 1.0)) +
                         a.hessian * pow(a.value, b - 1.0) * b;
    return result;
}

template <typename T>
T square(const T &a, int order) {
    return powf(a, 2.0, order);
}
}  // namespace GradArithmetic
