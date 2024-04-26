#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorMesh3D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

BoundaryGeneratorMesh3D::BoundaryGeneratorMesh3D(MatrixXI tri) : tri(std::move(tri)) {
    assert(this->tri.cols() == BoundaryGeneratorMesh3D::getDims());
}

void BoundaryGeneratorMesh3D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                              int order) const {
    BoundaryHelper::computeSimpleBoundaryVertices(degrees_of_freedom, boundary_data, getDims(), order);

    boundary_data.f.resize(tri.rows());
    for (int ii = 0; ii < tri.rows(); ii++) {
        boundary_data.f[ii].vert_indices = tri.row(ii);
    }
}

bool BoundaryGeneratorMesh3D::checkValidNumBoundaryParams(int n_params) const {
    /// Max vertex index in 'edge' should be number of vertices minus one.
    return (n_params == getDims() * (tri.maxCoeff() + 1));
}
