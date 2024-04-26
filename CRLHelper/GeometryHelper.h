#pragma once

#include "VecMatDef.h"

namespace GeometryHelper {
Vector3F rotateVectorAroundAxisAndCenter(const Vector3F &vec, const Vector3F &axis, F angle_rad,
                                         const Vector3F &center);

Vector2F rotateVector2D(const Vector2F &vec, F angle_rad);

bool intersectLinesPointDirection2D(const Vector2F &p0, const Vector2F &d0, const Vector2F &p1, const Vector2F &d1,
                                    Vector2F &t);

/// Computes the winding number of an edge with respect to a point. Sum for all edges of a polygon to determine the
/// winding number for the polygon. Winding number is 2 * pi if point is in polygon.
F windingNumber2D(const Vector2F &point, const Vector2F &edge0, const Vector2F &edge1);

/// Computes the winding number of a triangular face with respect to a point. Sum for all faces of a mesh to
/// determine the winding number for the mesh. Winding number is 4 * pi if point is inside mesh.
F windingNumber3D(const Vector3F &point, const Vector3F &face0, const Vector3F &face1, const Vector3F &face2);

/// Calculates the signed area of the CCW-oriented triangle O-p0-p1.
F triangleAreaWithOrigin2D(const Vector2F &p0, const Vector2F &p1);

/// Calculates the signed volume of the tetrahedron formed by connecting the origin to the triangle p0-p1-p2,
/// oriented CCW from the outside (normal pointing outward).
F tetVolumeWithOrigin3D(const Vector3F &p0, const Vector3F &p1, const Vector3F &p2);
}  // namespace GeometryHelper
