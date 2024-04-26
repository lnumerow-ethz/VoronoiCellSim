#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorPolylineLaplacian2D.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

BoundaryGeneratorPolylineLaplacian2D::BoundaryGeneratorPolylineLaplacian2D(const MatrixXI &edge, F spring_constant)
    : BoundaryGeneratorPolyline2D(edge), spring_constant(spring_constant) {}

F BoundaryGeneratorPolylineLaplacian2D::computeEnergy(const Model &model) const {
    F energy = 0;

    for (int i = 0; i < edge.rows(); i++) {
        Vector2F v0 = model.boundary_data.v[edge(i, 0)].pos;
        Vector2F v1 = model.boundary_data.v[edge(i, 1)].pos;
        energy += spring_constant * (v1 - v0).squaredNorm();
    }

    return energy;
}

void BoundaryGeneratorPolylineLaplacian2D::computeEnergyGradient(const Model &model, VectorXF &gradient) const {
    gradient = VectorXF::Zero(model.dimensions_ind->np);

    for (int i = 0; i < edge.rows(); i++) {
        Vector2F v0 = model.boundary_data.v[edge(i, 0)].pos;
        Vector2F v1 = model.boundary_data.v[edge(i, 1)].pos;
        Vector2F diff = 2 * spring_constant * (v1 - v0);

        int ip0 = 2 * edge(i, 0);
        int ip1 = 2 * edge(i, 1);
        gradient.segment(ip0, 2) -= diff;
        gradient.segment(ip1, 2) += diff;
    }
}

void BoundaryGeneratorPolylineLaplacian2D::computeEnergyHessian(const Model &model, HessianF &hessian) const {
    TripletListF hessian_triplets;
    for (int i = 0; i < edge.rows(); i++) {
        F value = 2 * spring_constant;
        int ip0 = 2 * edge(i, 0);
        int ip1 = 2 * edge(i, 1);
        for (int coord = 0; coord < 2; coord++) {
            hessian_triplets.emplace_back(ip0 + coord, ip0 + coord, value);
            hessian_triplets.emplace_back(ip0 + coord, ip1 + coord, -value);
            hessian_triplets.emplace_back(ip1 + coord, ip0 + coord, -value);
            hessian_triplets.emplace_back(ip1 + coord, ip1 + coord, value);
        }
    }

    hessian.setZero(model.dimensions_ind->np);
    hessian.A.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
}
