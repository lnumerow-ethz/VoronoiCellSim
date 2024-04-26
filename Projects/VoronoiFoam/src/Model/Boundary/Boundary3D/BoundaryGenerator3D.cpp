#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGenerator3D.h"

#include "CRLHelper/GeometryHelper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;

bool BoundaryGenerator3D::checkValidBoundary3D(const BoundaryData &boundary_data) {
    /// Construct CGAL mesh from boundary.
    Mesh mesh;
    std::vector<Mesh::Vertex_index> vertex_indices;
    for (const BoundaryVertex &v : boundary_data.v) {
        vertex_indices.push_back(mesh.add_vertex(Point_3(v.pos(0), v.pos(1), v.pos(2))));
    }
    for (const BoundaryFace &f : boundary_data.f) {
        mesh.add_face(vertex_indices[f.vert_indices(0)], vertex_indices[f.vert_indices(1)],
                      vertex_indices[f.vert_indices(2)]);
    }

    /// Check if self-intersecting.
    bool self_intersects =
        PMP::does_self_intersect(mesh, PMP::parameters::vertex_point_map(mesh.points()).geom_traits(Kernel()));
    if (self_intersects) {
        return false;
    }

    /// Check if non-manifold.
    bool non_manifold = false;
    for (auto v : mesh.vertices()) {
        if (PMP::is_non_manifold_vertex(v, mesh)) {
            non_manifold = true;
            break;
        }
    }
    if (non_manifold) {
        return false;
    }

    /// Check if non-closed.
    bool non_closed = !CGAL::is_closed(mesh);
    if (non_closed) {
        return false;
    }

    /// Check if inverted (encloses negative volume).
    bool inverted = !PMP::is_outward_oriented(mesh);
    if (inverted) {
        return false;
    }

    /// Passed all tests. Return true for a valid mesh.
    return true;
}

bool BoundaryGenerator3D::checkValidBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                             const BoundaryData &boundary_data) const {
    return checkValidBoundary3D(boundary_data);
}

bool BoundaryGenerator3D::checkPointInBounds3D(const BoundaryData &boundary_data, const VectorXF &point) {
    F winding_number = 0;
    for (const BoundaryFace &bface : boundary_data.f) {
        winding_number += GeometryHelper::windingNumber3D(point, boundary_data.v[bface.vert_indices(0)].pos,
                                                          boundary_data.v[bface.vert_indices(1)].pos,
                                                          boundary_data.v[bface.vert_indices(2)].pos);
    }

    /// Winding number is 4pi in bounds, 0 otherwise. Split the difference.
    return winding_number > 2 * M_PI;
}

bool BoundaryGenerator3D::checkPointInBounds(const BoundaryData &boundary_data, const VectorXF &point) const {
    return checkPointInBounds3D(boundary_data, point);
}

F BoundaryGenerator3D::computeEnclosedVolume(const BoundaryData &boundary_data) const {
    F enclosed_volume = 0;

    for (const BoundaryFace &bface : boundary_data.f) {
        enclosed_volume += GeometryHelper::tetVolumeWithOrigin3D(boundary_data.v[bface.vert_indices(0)].pos,
                                                                 boundary_data.v[bface.vert_indices(1)].pos,
                                                                 boundary_data.v[bface.vert_indices(2)].pos);
    }

    return enclosed_volume;
}
