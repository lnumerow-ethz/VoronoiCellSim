#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGenerator2D.h"

#include "CRLHelper/GeometryHelper.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::Segment_2 Segment_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

bool BoundaryGenerator2D::checkValidBoundary2D(const BoundaryData &boundary_data) {
    /// Construct CGAL arrangement from boundary.
    Arrangement_2 arr;
    for (const BoundaryFace &bface : boundary_data.f) {
        Vector2F v0 = boundary_data.v[bface.vert_indices(0)].pos;
        Vector2F v1 = boundary_data.v[bface.vert_indices(1)].pos;
        CGAL::insert_curve(arr, Segment_2(Point_2(v0.x(), v0.y()), Point_2(v1.x(), v1.y())));
    }

    /// Check if non-manifold, non-closed or self-intersecting. Vertex of degree != 2 indicates one of these.
    bool non_manifold_or_non_closed = false;
    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        if (vit->degree() != 2) {
            non_manifold_or_non_closed = true;
            break;
        }
    }
    if (non_manifold_or_non_closed) {
        return false;
    }

    /// Check if inverted (encloses negative area). Compute a point offset inward from the boundary and check that this
    /// point is actually in the interior. Currently check all points to robustly handle boundaries with holes.
    for (const auto &f : boundary_data.f) {
        Vector2F v0 = boundary_data.v[f.vert_indices(0)].pos;
        Vector2F v1 = boundary_data.v[f.vert_indices(1)].pos;
        Vector2F offset_point = (v0 + v1) / 2.0 + 1e-8 * GeometryHelper::rotateVector2D(v1 - v0, M_PI_2);
        bool inverted = !checkPointInBounds2D(boundary_data, offset_point);
        if (inverted) {
            return false;
        }
    }

    /// Passed all tests. Return true for a valid boundary.
    return true;
}

bool BoundaryGenerator2D::checkValidBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                             const BoundaryData &boundary_data) const {
    return checkValidBoundary2D(boundary_data);
}

bool BoundaryGenerator2D::checkPointInBounds2D(const BoundaryData &boundary_data, const VectorXF &point) {
    F winding_number = 0;
    for (const BoundaryFace &bface : boundary_data.f) {
        winding_number += GeometryHelper::windingNumber2D(point, boundary_data.v[bface.vert_indices(0)].pos,
                                                          boundary_data.v[bface.vert_indices(1)].pos);
    }

    /// Winding number is 2pi in bounds, 0 otherwise. Split the difference.
    return winding_number > M_PI && winding_number < 3 * M_PI;
}

bool BoundaryGenerator2D::checkPointInBounds(const BoundaryData &boundary_data, const VectorXF &point) const {
    return checkPointInBounds2D(boundary_data, point);
}

F BoundaryGenerator2D::computeEnclosedVolume(const BoundaryData &boundary_data) const {
    F enclosed_area = 0;

    for (const BoundaryFace &bface : boundary_data.f) {
        enclosed_area += GeometryHelper::triangleAreaWithOrigin2D(boundary_data.v[bface.vert_indices(0)].pos,
                                                                  boundary_data.v[bface.vert_indices(1)].pos);
    }

    return enclosed_area;
}
