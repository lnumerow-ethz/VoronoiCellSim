#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

void BoundaryHelper::computeSimpleBoundaryVertices(const DegreesOfFreedom &degrees_of_freedom,
                                                   BoundaryData &boundary_data, int dims_space, int order) {
    int n_param = (int)degrees_of_freedom.boundary_param.rows();
    int n_b_vertices = n_param / dims_space;

    boundary_data.v.resize(n_b_vertices);
    for (int iv = 0; iv < n_b_vertices; iv++) {
        boundary_data.v[iv].pos = degrees_of_freedom.boundary_param.segment(iv * dims_space, dims_space);
        if (order >= 1) {
            boundary_data.v[iv].grad.resize(dims_space, n_param);
            for (int id = 0; id < dims_space; id++) {
                boundary_data.v[iv].grad.coeffRef(id, iv * dims_space + id) = 1.0;
            }
        }
        if (order >= 2) {
            // Hessians are zero.
            boundary_data.v[iv].hess.resize(dims_space);
            for (int id = 0; id < dims_space; id++) {
                boundary_data.v[iv].hess[id].resize(n_param, n_param);
            }
        }
    }
}

void BoundaryHelper::computeDimensionScaledBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                                    BoundaryData &boundary_data, int order,
                                                    const MatrixXF &unscaled_vertices, const MatrixXI &unscaled_faces) {
    int n_params = (int)degrees_of_freedom.boundary_param.rows();
    int dims_space = (int)unscaled_vertices.cols();

    VectorXF scale_lengths;
    if (n_params == 1) {
        /// Scale whole box by a constant.
        scale_lengths = VectorXF::Constant(dims_space, degrees_of_freedom.boundary_param[0]);
    } else {
        /// Scale each dimension by its corresponding parameter.
        scale_lengths = degrees_of_freedom.boundary_param;
    }

    boundary_data.f.resize(unscaled_faces.rows());
    boundary_data.v.resize(unscaled_vertices.rows());

    for (int i = 0; i < unscaled_faces.rows(); i++) {
        boundary_data.f[i].vert_indices = unscaled_faces.row(i);
    }

    for (int i = 0; i < unscaled_vertices.rows(); i++) {
        boundary_data.v[i].pos = unscaled_vertices.row(i).transpose().array() * scale_lengths.array();

        if (order >= 1) {
            VectorXF unscaled_vertex = unscaled_vertices.row(i);
            if (n_params == 1) {
                boundary_data.v[i].grad = unscaled_vertex.sparseView();
            } else {
                boundary_data.v[i].grad = unscaled_vertex.asDiagonal();
            }
        }
        /// Hessians are zero.
        if (order >= 2) boundary_data.v[i].hess.resize(dims_space);
    }
}

void BoundaryHelper::getBoundaryBoundingBox(const BoundaryData &boundary_data, VectorXF &lower_bounds,
                                            VectorXF &upper_bounds) {
    int dims_space = (int)boundary_data.v[0].pos.rows();

    lower_bounds = VectorXF::Constant(dims_space, 1e10);
    upper_bounds = VectorXF::Constant(dims_space, -1e10);
    for (BoundaryVertex vertex : boundary_data.v) {
        lower_bounds = lower_bounds.cwiseMin(vertex.pos);
        upper_bounds = upper_bounds.cwiseMax(vertex.pos);
    }
}
