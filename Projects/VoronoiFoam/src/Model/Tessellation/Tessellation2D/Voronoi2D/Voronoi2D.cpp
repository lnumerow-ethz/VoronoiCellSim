#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Voronoi2D/Voronoi2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

#include "CRLHelper/GeometryHelper.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular;
typedef Regular::Vertex_handle Vertex_handle;
typedef K::Weighted_point_2 Weighted_point;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

static void constructBoundaryPolygons(Model &model, std::vector<std::vector<int>> &boundary_polygons) {
    /// Map start vertex of segment to segment index.
    std::map<int, int> start_vertex_to_bface_index;
    for (int i = 0; i < (int)model.boundary_data.f.size(); i++) {
        const BoundaryFace &bface = model.boundary_data.f[i];
        start_vertex_to_bface_index[bface.vert_indices[0]] = i;
    }

    std::vector<bool> visited(model.boundary_data.f.size());
    std::fill(visited.begin(), visited.end(), false);

    for (int i = 0; i < (int)model.boundary_data.f.size(); i++) {
        if (visited[i]) continue;

        boundary_polygons.emplace_back();

        int curr_face_index = i;
        do {
            /// Add segment to current boundary polygon.
            boundary_polygons.back().push_back(curr_face_index);
            visited[curr_face_index] = true;

            /// Move to next segment (segment whose start is this one's end).
            curr_face_index = start_vertex_to_bface_index.at(model.boundary_data.f[curr_face_index].vert_indices[1]);
        } while (curr_face_index != i);
    }
}

struct BoundaryCellIntersectingSegment {
    int bface_index = -1;
    bool enters_cell = false;
    Vertex_handle cell_entered_from_vertex_handle;
    bool leaves_cell = false;
    Vertex_handle cell_leaving_to_vertex_handle;
};

struct BoundaryCellVoronoiEdge {
    Vertex_handle cell_vertex_handle;
    Vertex_handle opposite_vertex_handle;
    Vertex_handle start_endpoint_gen_vertex_handle;
    Vertex_handle end_endpoint_gen_vertex_handle;

    Vector2F point;
    Vector2F direction;

    /// First value in pair indicates how far along the Voronoi edge the intersection occurs, for sorting.
    /// Second value is boundary face index.
    std::vector<std::pair<F, int>> boundary_intersections;
};

struct BoundaryCellConstructor {
    Vertex_handle vertex_handle;
    std::vector<BoundaryCellVoronoiEdge> voronoi_edges;
    std::vector<BoundaryCellIntersectingSegment> intersecting_boundary_segments;
};

static void getBisectorPointAndDirection(const Regular &rt, const Vertex_handle &site_vertex_handle,
                                         const Vertex_handle &other_vertex_handle, Vector2F &point,
                                         Vector2F &direction) {
    Weighted_point weighted_point_0 = site_vertex_handle->point();
    Weighted_point weighted_point_1 = other_vertex_handle->point();

    Vector2F edge_normal_direction = {weighted_point_1.x() - weighted_point_0.x(),
                                      weighted_point_1.y() - weighted_point_0.y()};

    direction = GeometryHelper::rotateVector2D(edge_normal_direction, M_PI_2);
    Weighted_point weighted_point_dummy({weighted_point_0.x() + direction.x(), weighted_point_0.y() + direction.y()},
                                        0);
    Point circumcenter = rt.weighted_circumcenter(weighted_point_0, weighted_point_1, weighted_point_dummy);
    point = {circumcenter.x(), circumcenter.y()};
}

/// Initialize boundary cell struct for a given site, computing all the adjacent Voronoi edges.
/// Initializes with no intersections.
static void initializeBoundaryCell(const Regular &rt, const Vertex_handle &site_vertex_handle,
                                   std::map<Vertex_handle, BoundaryCellConstructor> &vh_to_boundary_cell_map) {
    /// Skip if already initialized.
    if (vh_to_boundary_cell_map.find(site_vertex_handle) != vh_to_boundary_cell_map.end()) return;
    /// Emplace boundary cell struct in map.
    BoundaryCellConstructor &boundary_cell_constructor = vh_to_boundary_cell_map[site_vertex_handle];

    Regular::Face_circulator fc = rt.incident_faces(site_vertex_handle);
    if (fc.is_empty()) {  /// Handle degenerate cases where input has only one or two sites (RT is empty).
        int num_total_vertices = (int)rt.number_of_vertices();
        assert(num_total_vertices == 1 || num_total_vertices == 2);

        /// Find point and direction for the bisector, create Voronoi edge.
        if (num_total_vertices == 2) {
            boundary_cell_constructor.voronoi_edges.emplace_back();
            BoundaryCellVoronoiEdge &voronoi_edge = boundary_cell_constructor.voronoi_edges.back();

            auto vertex_iterator = rt.finite_vertices_begin();
            Vertex_handle other_vertex_handle = (site_vertex_handle == Vertex_handle(vertex_iterator))
                                                    ? Vertex_handle(++vertex_iterator)
                                                    : Vertex_handle(vertex_iterator);

            voronoi_edge.cell_vertex_handle = site_vertex_handle;
            voronoi_edge.opposite_vertex_handle = other_vertex_handle;
            voronoi_edge.start_endpoint_gen_vertex_handle = rt.infinite_vertex();
            voronoi_edge.end_endpoint_gen_vertex_handle = rt.infinite_vertex();

            getBisectorPointAndDirection(rt, site_vertex_handle, other_vertex_handle, voronoi_edge.point,
                                         voronoi_edge.direction);
        }

        /// Do nothing if only one vertex, as there are no Voronoi edges.
        return;
    }
    auto lambda_neighbor_site_handle = [&](const Regular::Face_circulator &fc) {
        return fc->vertex(Regular::ccw(fc->index(site_vertex_handle)));
    };
    Vertex_handle neighbor_site_handles[3] = {lambda_neighbor_site_handle(fc++), lambda_neighbor_site_handle(fc++),
                                              lambda_neighbor_site_handle(fc++)};

    Regular::Face_circulator done = fc;
    do {
        if (!rt.is_infinite(neighbor_site_handles[1])) {
            boundary_cell_constructor.voronoi_edges.emplace_back();
            BoundaryCellVoronoiEdge &voronoi_edge = boundary_cell_constructor.voronoi_edges.back();

            voronoi_edge.cell_vertex_handle = site_vertex_handle;
            voronoi_edge.opposite_vertex_handle = neighbor_site_handles[1];
            voronoi_edge.start_endpoint_gen_vertex_handle = neighbor_site_handles[0];
            voronoi_edge.end_endpoint_gen_vertex_handle = neighbor_site_handles[2];

            getBisectorPointAndDirection(rt, site_vertex_handle, voronoi_edge.opposite_vertex_handle,
                                         voronoi_edge.point, voronoi_edge.direction);
        }

        neighbor_site_handles[0] = neighbor_site_handles[1];
        neighbor_site_handles[1] = neighbor_site_handles[2];
        neighbor_site_handles[2] = lambda_neighbor_site_handle(fc++);
    } while (fc != done);  /// Stop once we're back at the first edge.
}

/// Returns vertex handle of site whose cell this segment terminates in.
static Vertex_handle addIntersectingSegment(const Model &model, const Vertex_handle &site_vertex_handle,
                                            int bface_index,
                                            std::map<Vertex_handle, BoundaryCellConstructor> &vh_to_boundary_cell_map) {
    VectorXI bface_vert_indices = model.boundary_data.f[bface_index].vert_indices;
    Vector2F bface_start_point = model.boundary_data.v[bface_vert_indices(0)].pos;
    Vector2F bface_end_point = model.boundary_data.v[bface_vert_indices(1)].pos;

    /// Boundary cell struct should already be initialized.
    BoundaryCellConstructor &boundary_cell = vh_to_boundary_cell_map.at(site_vertex_handle);

    /// Check the intersection of each voronoi edge with the boundary segment.
    /// Tuple is (t along boundary segment, Voronoi edge index, t along Voronoi edge).
    std::vector<std::tuple<F, int, F>> voronoi_edges_entering;
    std::vector<std::tuple<F, int, F>> voronoi_edges_leaving;
    for (int i = 0; i < (int)boundary_cell.voronoi_edges.size(); i++) {
        /// Compute the intersection point.
        Vector2F t;
        if (!GeometryHelper::intersectLinesPointDirection2D(bface_start_point, bface_end_point - bface_start_point,
                                                            boundary_cell.voronoi_edges[i].point,
                                                            boundary_cell.voronoi_edges[i].direction, t))
            continue;

        /// Determine if the boundary segment crosses the voronoi edge into (entering) or out of (leaving) the cell.
        bool orientation_entering =
            (bface_end_point - bface_start_point)
                .dot(GeometryHelper::rotateVector2D(boundary_cell.voronoi_edges[i].direction, M_PI_2)) > 0;
        if (orientation_entering) {
            voronoi_edges_entering.emplace_back(t(0), i, t(1));
        } else {
            voronoi_edges_leaving.emplace_back(t(0), i, t(1));
        }
    }

    std::sort(voronoi_edges_entering.begin(), voronoi_edges_entering.end());
    std::sort(voronoi_edges_leaving.begin(), voronoi_edges_leaving.end());

    boundary_cell.intersecting_boundary_segments.emplace_back();
    BoundaryCellIntersectingSegment &intersecting_segment = boundary_cell.intersecting_boundary_segments.back();

    intersecting_segment.bface_index = bface_index;
    if (!voronoi_edges_entering.empty()) {
        auto last_edge_entering = voronoi_edges_entering.back();
        if (std::get<0>(last_edge_entering) > 0.0) {
            assert(std::get<0>(last_edge_entering) < 1.0);
            intersecting_segment.enters_cell = true;
            intersecting_segment.cell_entered_from_vertex_handle =
                boundary_cell.voronoi_edges[std::get<1>(last_edge_entering)].opposite_vertex_handle;
            boundary_cell.voronoi_edges[std::get<1>(last_edge_entering)].boundary_intersections.emplace_back(
                std::get<2>(last_edge_entering), bface_index);
        }
    }
    if (!voronoi_edges_leaving.empty()) {
        auto first_edge_leaving = voronoi_edges_leaving.front();
        if (std::get<0>(first_edge_leaving) < 1.0) {
            assert(std::get<0>(first_edge_leaving) > 0.0);
            intersecting_segment.leaves_cell = true;
            intersecting_segment.cell_leaving_to_vertex_handle =
                boundary_cell.voronoi_edges[std::get<1>(first_edge_leaving)].opposite_vertex_handle;
            boundary_cell.voronoi_edges[std::get<1>(first_edge_leaving)].boundary_intersections.emplace_back(
                std::get<2>(first_edge_leaving), bface_index);

            return intersecting_segment.cell_leaving_to_vertex_handle;
        }
    }

    return site_vertex_handle;
}

static void traverseBoundaryPolygon(const Model &model, const std::vector<int> &boundary_polygon, const Regular &rt,
                                    std::map<Vertex_handle, BoundaryCellConstructor> &vh_to_boundary_cell_map) {
    const BoundaryVertex &start_bvertex =
        model.boundary_data.v[model.boundary_data.f[boundary_polygon[0]].vert_indices[0]];
    Point start_point(start_bvertex.pos.x(), start_bvertex.pos.y());

    Vertex_handle start_cell_handle = rt.nearest_power_vertex(start_point);
    initializeBoundaryCell(rt, start_cell_handle, vh_to_boundary_cell_map);

    auto curr_cell_handle = start_cell_handle;
    for (int bface_index : boundary_polygon) {
        auto end_cell_handle = addIntersectingSegment(model, curr_cell_handle, bface_index, vh_to_boundary_cell_map);
        while (end_cell_handle != curr_cell_handle) {
            curr_cell_handle = end_cell_handle;
            initializeBoundaryCell(rt, curr_cell_handle, vh_to_boundary_cell_map);
            end_cell_handle = addIntersectingSegment(model, curr_cell_handle, bface_index, vh_to_boundary_cell_map);
        }
    }
}

int getNodeIndexAndAddToMap2(std::map<TessellationNode, int> &node_index_map,
                             TessellationNode &node) {  /// TODO: refactor
    auto found = node_index_map.find(node);
    if (found == node_index_map.end()) {
        int previous_num_nodes = node_index_map.size();
        node_index_map.emplace(node, previous_num_nodes);
        return previous_num_nodes;
    } else {
        return found->second;
    }
}

bool Voronoi2D::computeCellsClippedVoronoi2D(Model &model, bool use_weights) {
    int n_sites = model.dimensions_ind->n_sites;

    /// Maintain a map from Weighted_point to index.
    std::map<Weighted_point, int> point_to_site_index;

    /// Construct the unclipped Voronoi diagram from sites.
    std::vector<Weighted_point> weighted_points;
    for (int i = 0; i < n_sites; i++) {
        const Site &site = model.degrees_of_freedom.sites[i];
        if (site.is_removed) continue;
        weighted_points.emplace_back(Point(site.pos.x(), site.pos.y()),
                                     (use_weights ? site.param[SITE_PARAM_POWER_WEIGHT] : 0));
        point_to_site_index[weighted_points.back()] = i;
    }
    Regular rt(weighted_points.begin(), weighted_points.end());

    if (rt.number_of_hidden_vertices() > 0) {
        /// Fail if a non-removed site has no cell.
        return false;
    }

    /// Construct boundary.
    std::vector<std::vector<int>> boundary_polygons;
    constructBoundaryPolygons(model, boundary_polygons);

    /// Store all tagged boundary cells in a map with a struct to keep track of boundary intersections.
    std::map<Vertex_handle, BoundaryCellConstructor> vh_to_boundary_cell_map;

    /// Traverse each boundary polygon and add the boundary-cell intersections.
    for (const std::vector<int> &boundary_polygon : boundary_polygons) {
        traverseBoundaryPolygon(model, boundary_polygon, rt, vh_to_boundary_cell_map);
    }

    /// Construct a map of node indices for future storage in model.
    std::map<TessellationNode, int> node_index_map;

    /// TODO: Clean up from here!!!
    model.cells.resize(n_sites);
    for (Vertex_handle site_vertex_handle : rt.finite_vertex_handles()) {
        int cell_index = point_to_site_index[site_vertex_handle->point()];

        if (vh_to_boundary_cell_map.find(site_vertex_handle) == vh_to_boundary_cell_map.end()) {
            /// Cell has no boundary intersections. Compute as normal Voronoi cell and skip completely if out of bounds.
            Regular::Face_circulator fc = rt.incident_faces(site_vertex_handle);
            auto lambda_neighbor_site_handle = [&](const Regular::Face_circulator &fc) {
                return fc->vertex(Regular::ccw(fc->index(site_vertex_handle)));
            };
            Vertex_handle neighbor_site_handles[3] = {lambda_neighbor_site_handle(fc++),
                                                      lambda_neighbor_site_handle(fc++),
                                                      lambda_neighbor_site_handle(fc++)};
            if (rt.is_infinite(neighbor_site_handles[0]) || rt.is_infinite(neighbor_site_handles[1]) ||
                rt.is_infinite(neighbor_site_handles[2])) {
                continue;
            }
            Point test_point = rt.weighted_circumcenter(neighbor_site_handles[0]->point(),
                                                        neighbor_site_handles[1]->point(), site_vertex_handle->point());
            if (!model.model_definition.boundary_generator->checkPointInBounds(
                    model.boundary_data, Vector2F(test_point.x(), test_point.y()))) {
                continue;
            }

            /// Should be a valid Voronoi cell within bounds. Iterate over faces and add to tessellation.
            Regular::Face_circulator done = fc;
            do {
                TessellationFace face;
                face.is_boundary_face = false;
                face.opposite_gen_index = point_to_site_index[neighbor_site_handles[1]->point()];

                for (int i = 0; i < 2; i++) {
                    TessellationNode node;
                    node.gen.resize(
                        3);  /// In Voronoi and similar 2D tessellations, 3 generator indices are sufficient.
                    node.type = NodeType::STANDARD;
                    node.gen(0) = point_to_site_index[site_vertex_handle->point()];
                    node.gen(1) = point_to_site_index[neighbor_site_handles[0 + i]->point()];
                    node.gen(2) = point_to_site_index[neighbor_site_handles[1 + i]->point()];
                    std::sort(node.gen.data(), node.gen.data() + 3);

                    face.node_indices.push_back(getNodeIndexAndAddToMap2(node_index_map, node));
                }

                model.cells[cell_index].faces.push_back(face);

                neighbor_site_handles[0] = neighbor_site_handles[1];
                neighbor_site_handles[1] = neighbor_site_handles[2];
                neighbor_site_handles[2] = lambda_neighbor_site_handle(fc++);
            } while (fc != done);  /// Stop once we're back at the first edge.
        } else {
            BoundaryCellConstructor &boundary_cell = vh_to_boundary_cell_map.at(site_vertex_handle);
            for (const auto &boundary_segment : boundary_cell.intersecting_boundary_segments) {
                TessellationFace face;
                face.is_boundary_face = true;
                face.opposite_gen_index = boundary_segment.bface_index;

                TessellationNode nodes[2];

                for (int i = 0; i < 2; i++) {
                    TessellationNode &node = nodes[i];
                    node.gen.resize(3);
                    node.type = NodeType::B_VERTEX;
                    node.gen(0) = model.boundary_data.f[boundary_segment.bface_index].vert_indices(i);
                    node.gen(1) = -1;
                    node.gen(2) = -1;
                }

                if (boundary_segment.enters_cell) {
                    TessellationNode &node = nodes[0];
                    node.type = NodeType::B_EDGE;
                    node.gen(0) = boundary_segment.bface_index;
                    node.gen(1) = cell_index;
                    node.gen(2) = point_to_site_index[boundary_segment.cell_entered_from_vertex_handle->point()];
                    std::sort(node.gen.data() + 1, node.gen.data() + 3);
                }
                if (boundary_segment.leaves_cell) {
                    TessellationNode &node = nodes[1];
                    node.type = NodeType::B_EDGE;
                    node.gen(0) = boundary_segment.bface_index;
                    node.gen(1) = cell_index;
                    node.gen(2) = point_to_site_index[boundary_segment.cell_leaving_to_vertex_handle->point()];
                    std::sort(node.gen.data() + 1, node.gen.data() + 3);
                }

                for (int i = 0; i < 2; i++) {
                    face.node_indices.push_back(getNodeIndexAndAddToMap2(node_index_map, nodes[i]));
                }

                model.cells[cell_index].faces.push_back(face);
            }

            /// If there is only one site, there are no Voronoi edges.
            if (boundary_cell.voronoi_edges.empty()) {
                continue;
            }

            auto &first_voronoi_edge = boundary_cell.voronoi_edges[0];
            bool in_bounds;
            if (rt.is_infinite(first_voronoi_edge.start_endpoint_gen_vertex_handle)) {
                in_bounds = false;
            } else {
                Point test_point = rt.weighted_circumcenter(
                    first_voronoi_edge.cell_vertex_handle->point(), first_voronoi_edge.opposite_vertex_handle->point(),
                    first_voronoi_edge.start_endpoint_gen_vertex_handle->point());
                in_bounds = model.model_definition.boundary_generator->checkPointInBounds(
                    model.boundary_data, Vector2F(test_point.x(), test_point.y()));
            }
            for (auto &voronoi_edge : boundary_cell.voronoi_edges) {
                int opposite_cell_index = point_to_site_index[voronoi_edge.opposite_vertex_handle->point()];

                std::sort(voronoi_edge.boundary_intersections.begin(), voronoi_edge.boundary_intersections.end());
                int num_segments = 1 + (int)voronoi_edge.boundary_intersections.size();

                for (int i = 0; i < num_segments; i++) {
                    if (in_bounds) {
                        TessellationFace face;
                        face.is_boundary_face = false;
                        face.opposite_gen_index = opposite_cell_index;

                        TessellationNode nodes[2];
                        nodes[0].gen.resize(3);
                        nodes[1].gen.resize(3);

                        if (i == 0) {
                            TessellationNode &node = nodes[0];
                            node.type = NodeType::STANDARD;
                            node.gen(0) = cell_index;
                            node.gen(1) = opposite_cell_index;
                            node.gen(2) = point_to_site_index[voronoi_edge.start_endpoint_gen_vertex_handle->point()];
                            std::sort(node.gen.data(), node.gen.data() + 3);
                        } else {
                            TessellationNode &node = nodes[0];
                            node.type = NodeType::B_EDGE;
                            node.gen(0) = voronoi_edge.boundary_intersections[i - 1].second;
                            node.gen(1) = cell_index;
                            node.gen(2) = opposite_cell_index;
                            std::sort(node.gen.data() + 1, node.gen.data() + 3);
                        }

                        if (i == num_segments - 1) {
                            TessellationNode &node = nodes[1];
                            node.type = NodeType::STANDARD;
                            node.gen(0) = cell_index;
                            node.gen(1) = opposite_cell_index;
                            node.gen(2) = point_to_site_index[voronoi_edge.end_endpoint_gen_vertex_handle->point()];
                            std::sort(node.gen.data(), node.gen.data() + 3);
                        } else {
                            TessellationNode &node = nodes[1];
                            node.type = NodeType::B_EDGE;
                            node.gen(0) = voronoi_edge.boundary_intersections[i].second;
                            node.gen(1) = cell_index;
                            node.gen(2) = opposite_cell_index;
                            std::sort(node.gen.data() + 1, node.gen.data() + 3);
                        }

                        for (int i = 0; i < 2; i++) {
                            face.node_indices.push_back(getNodeIndexAndAddToMap2(node_index_map, nodes[i]));
                        }

                        model.cells[cell_index].faces.push_back(face);
                    }

                    /// Pass through a boundary edge, unless we've reached the end of this Voronoi edge.
                    if (i != num_segments - 1) in_bounds = !in_bounds;
                }
            }
        }
    }

    /// Fail if a non-removed site has no cell.
    for (int i = 0; i < n_sites; i++) {
        if (!model.degrees_of_freedom.sites[i].is_removed && model.cells[i].faces.empty()) {
            return false;
        }
    }

    /// Construct model nodes vector.
    model.nodes.resize(node_index_map.size());
    for (const auto &node_index_pair : node_index_map) {
        model.nodes[node_index_pair.second] = node_index_pair.first;
    }

    return true;
}

bool Voronoi2D::computeCells(Model &model) const {
    return computeCellsClippedVoronoi2D(model, false);
}
