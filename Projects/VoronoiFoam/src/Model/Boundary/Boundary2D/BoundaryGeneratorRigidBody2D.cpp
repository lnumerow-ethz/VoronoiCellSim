#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorRigidBody2D.h"

BoundaryGeneratorRigidBody2D::BoundaryGeneratorRigidBody2D(MatrixXF vertex_boundary, MatrixXI edge_boundary,
                                                           MatrixXF vertex_rigid_body, MatrixXI edge_rigid_body,
                                                           F rigid_body_force_magnitude)
    : vertex_boundary(std::move(vertex_boundary)),
      edge_boundary(std::move(edge_boundary)),
      vertex_rigid_body(std::move(vertex_rigid_body)),
      edge_rigid_body(std::move(edge_rigid_body)),
      rigid_body_force_magnitude(std::move(rigid_body_force_magnitude)) {}

void BoundaryGeneratorRigidBody2D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                                   BoundaryData &boundary_data, int order) const {
    int dims_space = getDims();
    int n_param = 3;

    int n_vertices_boundary = vertex_boundary.rows();
    int n_vertices_rigid_body = vertex_rigid_body.rows();
    int n_vertices_total = n_vertices_boundary + n_vertices_rigid_body;

    int n_faces_boundary = edge_boundary.rows();
    int n_faces_rigid_body = edge_rigid_body.rows();
    int n_faces_total = n_faces_boundary + n_faces_rigid_body;

    boundary_data.v.resize(n_vertices_total);
    for (int iv = 0; iv < n_vertices_boundary; iv++) {
        BoundaryVertex &vertex = boundary_data.v[iv];
        vertex.pos = vertex_boundary.row(iv);
        if (order >= 1) {
            // Gradients are zero for fixed boundary vertices.
            vertex.grad.resize(dims_space, n_param);
        }
        if (order >= 2) {
            // Hessians are zero.
            vertex.hess.resize(dims_space);
            for (int id = 0; id < dims_space; id++) {
                vertex.hess[id].resize(n_param, n_param);
            }
        }
    }

    F dx = degrees_of_freedom.boundary_param(0);
    F dy = degrees_of_freedom.boundary_param(1);
    F theta = degrees_of_freedom.boundary_param(2);

    for (int iv = 0; iv < n_vertices_rigid_body; iv++) {
        BoundaryVertex &vertex = boundary_data.v[n_vertices_total - 1 - iv];  /// Reverse order from CCW to CW.
        F x = vertex_rigid_body(iv, 0);
        F y = vertex_rigid_body(iv, 1);

        vertex.pos = Vector2F(dx + x * cos(theta) - y * sin(theta), dy + y * cos(theta) + x * sin(theta));
        if (order >= 1) {
            vertex.grad.resize(dims_space, n_param);
            vertex.grad.coeffRef(0, 0) = 1;
            vertex.grad.coeffRef(1, 1) = 1;
            vertex.grad.coeffRef(0, 2) = -x * sin(theta) - y * cos(theta);
            vertex.grad.coeffRef(1, 2) = x * cos(theta) - y * sin(theta);
        }
        if (order >= 2) {
            vertex.hess.resize(dims_space);
            for (int id = 0; id < dims_space; id++) {
                vertex.hess[id].resize(n_param, n_param);
            }
            vertex.hess[0].coeffRef(2, 2) = -x * cos(theta) + y * sin(theta);
            vertex.hess[1].coeffRef(2, 2) = -x * sin(theta) - y * cos(theta);
        }
    }

    boundary_data.f.resize(n_faces_total);
    for (int i = 0; i < n_faces_boundary; i++) {
        boundary_data.f[i].vert_indices = edge_boundary.row(i);
    }
    for (int i = 0; i < n_faces_rigid_body; i++) {
        boundary_data.f[i + n_faces_boundary].vert_indices =
            edge_rigid_body.row(i).transpose() + Vector2I::Constant(n_vertices_boundary);
    }
}

bool BoundaryGeneratorRigidBody2D::checkValidNumBoundaryParams(int n_params) const {
    /// Translation (dx, dy) and rotation (theta) of rigid body.
    return (n_params == 3);
}

F BoundaryGeneratorRigidBody2D::computeEnergy(const Model &model) const {
    F dx = model.degrees_of_freedom.boundary_param(0);
    F dy = model.degrees_of_freedom.boundary_param(1);

    return -rigid_body_force_magnitude * dx + 0 * pow(dy - 0.0, 2.0);
}

void BoundaryGeneratorRigidBody2D::computeEnergyGradient(const Model &model, VectorXF &gradient) const {
    int np = model.dimensions_ind->np;
    if (np == 0) return;

    gradient = VectorXF::Zero(np);

    F dx = model.degrees_of_freedom.boundary_param(0);
    F dy = model.degrees_of_freedom.boundary_param(1);

    gradient(0) = -rigid_body_force_magnitude;
    gradient(1) = 0 * 2 * (dy - 0.0);
}

void BoundaryGeneratorRigidBody2D::computeEnergyHessian(const Model &model, HessianF &hessian) const {
    int np = model.dimensions_ind->np;
    if (np == 0) return;

    hessian.setZero(model.dimensions_ind->np);

    hessian.A.coeffRef(0, 0) = 0;
    hessian.A.coeffRef(1, 1) = 0 * 2;
}
