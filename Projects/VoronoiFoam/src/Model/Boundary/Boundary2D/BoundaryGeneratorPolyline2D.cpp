#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorPolyline2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

BoundaryGeneratorPolyline2D::BoundaryGeneratorPolyline2D(MatrixXI edge) : edge(std::move(edge)) {
    assert(this->edge.cols() == BoundaryGeneratorPolyline2D::getDims());
}

void BoundaryGeneratorPolyline2D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                                  BoundaryData &boundary_data, int order) const {
    BoundaryHelper::computeSimpleBoundaryVertices(degrees_of_freedom, boundary_data, getDims(), order);

    boundary_data.f.resize(edge.rows());
    for (int ii = 0; ii < edge.rows(); ii++) {
        boundary_data.f[ii].vert_indices = edge.row(ii);
    }
}

bool BoundaryGeneratorPolyline2D::checkValidNumBoundaryParams(int n_params) const {
    /// Max vertex index in 'edge' should be number of vertices minus one.
    return (n_params == getDims() * (edge.maxCoeff() + 1));
}
