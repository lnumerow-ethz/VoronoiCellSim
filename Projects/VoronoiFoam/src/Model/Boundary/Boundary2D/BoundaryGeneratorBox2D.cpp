#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

static void getSquareMesh(MatrixXF &vertices, MatrixXI &faces) {
    vertices.resize(4, 2);

    vertices << -1, -1,  //
        1, -1,           //
        1, 1,            //
        -1, 1;

    faces.resize(4, 2);

    faces << 0, 1,  //
        1, 2,       //
        2, 3,       //
        3, 0;
}

void BoundaryGeneratorBox2D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                             int order) const {
    MatrixXF square_vertices;
    MatrixXI square_faces;
    getSquareMesh(square_vertices, square_faces);

    BoundaryHelper::computeDimensionScaledBoundary(degrees_of_freedom, boundary_data, order, square_vertices,
                                                   square_faces);
}

bool BoundaryGeneratorBox2D::checkValidNumBoundaryParams(int n_params) const {
    /// Separate side lengths in each dimension, or one uniform side length.
    return (n_params == getDims()) || (n_params == 1);
}
