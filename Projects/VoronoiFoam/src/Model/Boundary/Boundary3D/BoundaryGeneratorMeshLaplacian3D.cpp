#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorMeshLaplacian3D.h"

BoundaryGeneratorMeshLaplacian3D::BoundaryGeneratorMeshLaplacian3D(const MatrixXI &tri, F spring_constant)
    : BoundaryGeneratorMesh3D(tri), spring_constant(spring_constant) {}

F BoundaryGeneratorMeshLaplacian3D::computeEnergy(const Model &model) const {
    F energy = 0;

    for (int i = 0; i < tri.rows(); i++) {
        Vector3F v0 = model.boundary_data.v[tri(i, 0)].pos;
        Vector3F v1 = model.boundary_data.v[tri(i, 1)].pos;
        Vector3F v2 = model.boundary_data.v[tri(i, 2)].pos;
        energy += spring_constant * (v1 - v0).squaredNorm();
        energy += spring_constant * (v2 - v1).squaredNorm();
        energy += spring_constant * (v0 - v2).squaredNorm();
    }

    return energy;
}

void BoundaryGeneratorMeshLaplacian3D::computeEnergyGradient(const Model &model, VectorXF &gradient) const {
    gradient = VectorXF::Zero(model.dimensions_ind->np);

    for (int i = 0; i < tri.rows(); i++) {
        Vector3F v0 = model.boundary_data.v[tri(i, 0)].pos;
        Vector3F v1 = model.boundary_data.v[tri(i, 1)].pos;
        Vector3F v2 = model.boundary_data.v[tri(i, 2)].pos;
        Vector3F diff01 = 2 * spring_constant * (v1 - v0);
        Vector3F diff12 = 2 * spring_constant * (v2 - v1);
        Vector3F diff20 = 2 * spring_constant * (v0 - v2);

        int ip0 = 3 * tri(i, 0);
        int ip1 = 3 * tri(i, 1);
        int ip2 = 3 * tri(i, 2);
        gradient.segment(ip0, 3) += diff20 - diff01;
        gradient.segment(ip1, 3) += diff01 - diff12;
        gradient.segment(ip2, 3) += diff12 - diff20;
    }
}

void BoundaryGeneratorMeshLaplacian3D::computeEnergyHessian(const Model &model, HessianF &hessian) const {
    TripletListF hessian_triplets;
    for (int i = 0; i < tri.rows(); i++) {
        F value = 2 * spring_constant;
        int ip0 = 3 * tri(i, 0);
        int ip1 = 3 * tri(i, 1);
        int ip2 = 3 * tri(i, 2);
        for (int coord = 0; coord < 3; coord++) {
            hessian_triplets.emplace_back(ip0 + coord, ip0 + coord, 2 * value);
            hessian_triplets.emplace_back(ip0 + coord, ip1 + coord, -value);
            hessian_triplets.emplace_back(ip0 + coord, ip2 + coord, -value);
            hessian_triplets.emplace_back(ip1 + coord, ip0 + coord, -value);
            hessian_triplets.emplace_back(ip1 + coord, ip1 + coord, 2 * value);
            hessian_triplets.emplace_back(ip1 + coord, ip2 + coord, -value);
            hessian_triplets.emplace_back(ip2 + coord, ip0 + coord, -value);
            hessian_triplets.emplace_back(ip2 + coord, ip1 + coord, -value);
            hessian_triplets.emplace_back(ip2 + coord, ip2 + coord, 2 * value);
        }
    }

    hessian.setZero(model.dimensions_ind->np);
    hessian.A.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
}
