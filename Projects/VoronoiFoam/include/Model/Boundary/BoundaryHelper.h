#pragma once

#include "CRLHelper/VecMatDef.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"

/// Some helper functions for BoundaryGenerator derived classes, in order to share code between 2D and 3D.
namespace BoundaryHelper {
/// Computes vertices for boundaries whose parameters are simply the vertex positions (p = v).
void computeSimpleBoundaryVertices(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                   int dims_space, int order);

/// Computes vertices for boundaries whose parameters are simply the vertex positions (p = v).
void computeDimensionScaledBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data, int order,
                                    const MatrixXF &unscaled_vertices, const MatrixXI &unscaled_faces);

/// Computes a bounding box (lower and upper bounds for each coordinate) for a given boundary.
void getBoundaryBoundingBox(const BoundaryData &boundary_data, VectorXF &lower_bounds, VectorXF &upper_bounds);

/// Computes the enclosed volume (3D) or area (2D) of the boundary mesh.
F getBoundaryEnclosedVolume(const BoundaryData &boundary_data);
}  // namespace BoundaryHelper
