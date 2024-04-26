#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorBox3D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

static void getCubeMesh(MatrixXF &vertices, MatrixXI &faces) {
    vertices.resize(8, 3);

    vertices << -1, -1, -1,  //
        -1, -1, 1,           //
        -1, 1, -1,           //
        -1, 1, 1,            //
        1, -1, -1,           //
        1, -1, 1,            //
        1, 1, -1,            //
        1, 1, 1;             //

    faces.resize(12, 3);

    faces << 0, 1, 2,  //
        2, 1, 3,       //
        0, 4, 1,       //
        1, 4, 5,       //
        0, 2, 4,       //
        4, 2, 6,       //
        4, 6, 5,       //
        5, 6, 7,       //
        1, 5, 3,       //
        3, 5, 7,       //
        2, 3, 6,       //
        6, 3, 7;
}

void BoundaryGeneratorBox3D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                             int order) const {
    MatrixXF cube_vertices;
    MatrixXI cube_faces;
    getCubeMesh(cube_vertices, cube_faces);

    BoundaryHelper::computeDimensionScaledBoundary(degrees_of_freedom, boundary_data, order, cube_vertices, cube_faces);
}

bool BoundaryGeneratorBox3D::checkValidNumBoundaryParams(int n_params) const {
    /// Separate side lengths in each dimension, or one uniform side length.
    return (n_params == getDims()) || (n_params == 1);
}
