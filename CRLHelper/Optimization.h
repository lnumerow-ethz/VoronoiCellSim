#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"

#include "ThirdParty/LBFGSpp/include/LBFGSpp/BFGSMat.h"

namespace Optimization {
enum OptimizationStatus {
    CONVERGED,  //
    SUCCESS,
    FAILURE
};

enum Optimizer {
    GRADIENT_DESCENT,  //
    NEWTON,
    NEWTON_WOODBURY,
    // BFGS,
    LBFGS,
    NUM_OPTIMIZERS
};
static std::vector<std::string> optimizer_names = {"Gradient Descent",      //
                                                   "Newton",                //
                                                   "Newton with Woodbury",  //
                                                                            //    "BFGS",                  //
                                                   "LBFGS"};

/// Solve for x in (A + U * V^T)x = b. Return true on success.
bool linearSolve(const HessianF &hessian, const VectorXF &b, VectorXF &x);

/// Solve for x in (A + U * V^T)x = b, using Woodbury formula. Return true on success.
bool linearSolveWoodbury(const HessianF &hessian, const VectorXF &b, VectorXF &x);

/// These optimization routines take lambda functions for computing the objective value and its gradients.
/// Gradient and hessian functions also compute objective value and lower-order derivatives.
OptimizationStatus stepGradientDescent(VectorXF &y, const std::function<bool(const VectorXF &, F &)> &objectiveFunction,
                                       const std::function<bool(const VectorXF &, F &, VectorXF &)> &gradientFunction,
                                       F tolerance);

OptimizationStatus stepBFGS(VectorXF &y, MatrixXF &bfgs_inverse_hessian,
                            const std::function<bool(const VectorXF &, F &)> &objectiveFunction,
                            const std::function<bool(const VectorXF &, F &, VectorXF &)> &gradientFunction,
                            F tolerance);

OptimizationStatus stepLBFGS(VectorXF &y, LBFGSpp::BFGSMat<F> &lbfgs_mat,
                             const std::function<bool(const VectorXF &, F &)> &objectiveFunction,
                             const std::function<bool(const VectorXF &, F &, VectorXF &)> &gradientFunction,
                             F tolerance);

OptimizationStatus stepNewton(VectorXF &y, const std::function<bool(const VectorXF &, F &)> &objectiveFunction,
                              const std::function<bool(const VectorXF &, F &, VectorXF &, HessianF &)> &hessianFunction,
                              bool use_woodbury, F tolerance);

/// Find step along search direction dy that decreases the objective value. Update y accordingly.
bool lineSearch(VectorXF &y, const VectorXF &dy, const std::function<bool(const VectorXF &, F &)> &objectiveFunction,
                F initial_objective_value);

bool lineSearch(VectorXF &y, const VectorXF &dy, const std::function<bool(const VectorXF &, F &)> &objectiveFunction);

/// Check gradient for scalar function of vector argument.
bool checkGradient(const VectorXF &y, VectorXF &error, const std::function<bool(const VectorXF &, F &)> &func,
                   const std::function<bool(const VectorXF &, VectorXF &)> &grad_func, F epsilon, int print_level = 0);

/// Check hessian for scalar function of vector argument.
bool checkHessian(const VectorXF &y, MatrixXF &error,
                  const std::function<bool(const VectorXF &, VectorXF &)> &grad_func,
                  const std::function<bool(const VectorXF &, MatrixXF &)> &hess_func, F epsilon, int print_level = 0);
}  // namespace Optimization
