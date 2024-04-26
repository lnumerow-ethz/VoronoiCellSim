#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"

bool BoundaryGenerator::generateBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                         int order) const {
    if (!checkValidNumBoundaryParams((int)degrees_of_freedom.boundary_param.rows())) return false;

    computeBoundary(degrees_of_freedom, boundary_data, order);

    return checkValidBoundary(degrees_of_freedom, boundary_data);
}

F BoundaryGenerator::computeEnergy(const Model &model) const {
    return 0;
}

void BoundaryGenerator::computeEnergyGradient(const Model &model, VectorXF &gradient) const {
    /// Default implementation constructs zero vector of appropriate size.
    gradient = VectorXF::Zero(model.dimensions_ind->np);
}

void BoundaryGenerator::computeEnergyHessian(const Model &model, HessianF &hessian) const {
    /// Default implementation constructs zero matrices of appropriate size.
    hessian.setZero(model.dimensions_ind->np);
}
