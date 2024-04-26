#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/Voronoi3D/Voronoi3D.h"

#include "CRLHelper/EigenHelper.h"

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/voronoi/convex_cell.h>

/// Add some far-away sites to circumvent geogram bug when all vertices are coplanar.
static void geogramBugWorkAround(const Model &model, VectorXF &sites_vector, VectorXF &weights_vector,
                                 bool use_weights) {
    F min_coord = 1e10;
    F max_coord = 1e-10;

    for (Site site : model.degrees_of_freedom.sites) {
        VectorXF pos = site.pos;
        min_coord = std::min(min_coord, site.pos.minCoeff());
        max_coord = std::max(max_coord, site.pos.maxCoeff());
    }
    for (BoundaryVertex v : model.boundary_data.v) {
        VectorXF pos = v.pos;
        min_coord = std::min(min_coord, v.pos.minCoeff());
        max_coord = std::max(max_coord, v.pos.maxCoeff());
    }

    /// Place additional sites far enough away to ensure these cells are completely outside of boundary.
    F extra_sites_x = max_coord + 2 * (max_coord - min_coord);
    F extra_sites_weight = use_weights ? weights_vector.minCoeff() : 0.0;

    VectorXF extra_sites_vector(3 * 3);
    extra_sites_vector << extra_sites_x, 0, 0,  //
        extra_sites_x, 0, 1,                    //
        extra_sites_x, 1, 0;
    VectorXF extra_weights_vector = VectorXF::Constant(3, extra_sites_weight);

    /// Append extra weights and sites vectors to the original vectors.
    sites_vector = EigenHelper::vertCat(sites_vector, extra_sites_vector);
    weights_vector = EigenHelper::vertCat(weights_vector, extra_weights_vector);
}

int getNodeIndexAndAddToMap(std::map<TessellationNode, int> &node_index_map,
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

bool Voronoi3D::computeCellsClippedVoronoi3D(Model &model, bool use_weights) {
    int n_sites = model.dimensions_ind->n_sites;
    int dims_space = model.dimensions_ind->dims_space;

    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");

    /// Construct a map of node indices for future storage in model.
    std::map<TessellationNode, int> node_index_map;

    /// Map sorted vertex index triplet to boundary face. For identifying boundary faces later.
    std::map<std::tuple<int, int, int>, int> tri_to_bface_index;
    for (int i = 0; i < model.dimensions_boundary->n_b_faces; i++) {
        VectorXI vert_indices = model.boundary_data.f[i].vert_indices;
        std::sort(vert_indices.data(), vert_indices.data() + 3);
        tri_to_bface_index[{vert_indices(0), vert_indices(1), vert_indices(2)}] = i;
    }

    /// Assemble boundary mesh in geogram data structure.
    GEO::Mesh boundary_mesh;
    for (const BoundaryVertex &b_vertex : model.boundary_data.v) {
        boundary_mesh.vertices.create_vertex(b_vertex.pos.data());
    }
    for (const BoundaryFace &b_face : model.boundary_data.f) {
        VectorXI vert_indices = b_face.vert_indices;
        /// Reverse the cyclic order of vertices here, since our code defines the normal direction as outward.
        boundary_mesh.facets.create_triangle(vert_indices(0), vert_indices(2), vert_indices(1));
    }
    boundary_mesh.facets.connect();

    /// Identify indices of non-removed cells.
    std::vector<int> non_removed_cell_indices;
    for (int i = 0; i < n_sites; i++) {
        if (!model.degrees_of_freedom.sites[i].is_removed) non_removed_cell_indices.push_back(i);
    }
    int num_non_removed_cells = non_removed_cell_indices.size();

    /// Assemble site data vector for geogram Delaunay.
    VectorXF weights_vector(num_non_removed_cells);
    VectorXF sites_vector(num_non_removed_cells * dims_space);
    for (int i = 0; i < num_non_removed_cells; i++) {
        const Site &site = model.degrees_of_freedom.sites[non_removed_cell_indices[i]];
        sites_vector.segment(i * dims_space, dims_space) = site.pos;
        weights_vector(i) = site.param[SITE_PARAM_POWER_WEIGHT];
    }

    /// Add some far-away sites to circumvent geogram bug when all vertices are coplanar.
    geogramBugWorkAround(model, sites_vector, weights_vector, use_weights);

    GEO::PeriodicDelaunay3d dual(false);
    dual.set_vertices(sites_vector.rows() / 3, sites_vector.data());
    if (use_weights) dual.set_weights(weights_vector.data());
    dual.set_stores_cicl(true);
    dual.compute();
    if (dual.has_empty_cells()) {
        /// Fail if a non-removed site has no cell.
        return false;
    }

    /// Geogram requires to tetrahedralize the bounding domain before clipping.
    GEO::mesh_tetrahedralize(boundary_mesh, false, false);

    /// After tetrahedralizing the domain, the correspondence between our boundary face indices and faces in the
    /// tetrahedralization is lost. This block of code reconstructs this mapping.
    /// The geogram mesh stores "Cell facets" with unique indices.
    std::vector<std::pair<int, int>> cell_facet_to_facet_pairs;
    /// boundary_mesh.cells are tetrahedra from a tetrahedralization of the bounding domain.
    for (size_t cc = 0; cc < boundary_mesh.cells.nb(); cc++) {
        for (size_t ff = 0; ff < boundary_mesh.cells.nb_facets(cc); ff++) {
            Vector3I vert_indices = {(int)boundary_mesh.cells.facet_vertex(cc, ff, 0),
                                     (int)boundary_mesh.cells.facet_vertex(cc, ff, 1),
                                     (int)boundary_mesh.cells.facet_vertex(cc, ff, 2)};
            std::sort(vert_indices.data(), vert_indices.data() + 3);
            std::tuple<int, int, int> verts_tuple = {vert_indices(0), vert_indices(1), vert_indices(2)};
            if (tri_to_bface_index.find(verts_tuple) != tri_to_bface_index.end()) {
                cell_facet_to_facet_pairs.emplace_back(boundary_mesh.cells.facet(cc, ff),
                                                       tri_to_bface_index.at(verts_tuple));
                continue;
            }
        }
    }
    std::map<int, int> cell_facet_to_facet_map(cell_facet_to_facet_pairs.begin(), cell_facet_to_facet_pairs.end());

    /// Clip the Voronoi diagram with the boundary mesh.
    GEO::RestrictedVoronoiDiagram_var rvd = GEO::RestrictedVoronoiDiagram::create(&dual, &boundary_mesh);
    rvd->set_volumetric(true);
    GEO::Mesh clippedMesh;
    GEO::BuildRVDMesh rvdMeshBuilder(clippedMesh);
    rvdMeshBuilder.set_simplify_internal_tet_facets(true);
    rvdMeshBuilder.set_simplify_voronoi_facets(true);
    rvdMeshBuilder.set_generate_ids(true);
    rvd->for_each_polyhedron(rvdMeshBuilder, true, true, false);  /// Parallel=true causes crash. :(

    model.cells.resize(n_sites);
    GEO::Attribute<int> vertex_gen(clippedMesh.vertices.attributes(), "vertex_gen");
    GEO::Attribute<int> facet_seed_id(clippedMesh.facets.attributes(), "facet_seed_id");

    for (GEO::index_t f : clippedMesh.facets) {
        TessellationFace face;

        /// Find all boundary vertices involved in the generation of this face. Used to determine if this face is
        /// on a boundary facet, and if so, which one.
        std::set<int> b_vert_index_set;

        for (GEO::index_t lv = 0; lv < clippedMesh.facets.nb_vertices(f); ++lv) {
            GEO::index_t v = clippedMesh.facets.vertex(f, lv);

            /// Our fork of geogram provides information about which inputs a tessellation vertex was generated by.
            int gen[6];
            for (int i = 0; i < 6; i++) {
                gen[i] = vertex_gen[v * 6 + i];
            }

            TessellationNode node;
            node.gen.resize(4);  /// In Voronoi and similar 3D tessellations, 4 generator indices are sufficient.
            if (gen[1] < 0 && gen[2] < 0 && gen[3] < 0) {
                node.type = NodeType::B_VERTEX;
                node.gen(0) = gen[4];
                node.gen(1) = -1;
                node.gen(2) = -1;
                node.gen(3) = -1;
                b_vert_index_set.insert(gen[4]);
            } else if (gen[1] < 0 && gen[2] < 0) {
                node.type = NodeType::B_EDGE;
                node.gen(0) = gen[4];
                node.gen(1) = gen[5];
                node.gen(2) = gen[0];
                node.gen(3) = gen[3] - 1;
                std::sort(node.gen.data(), node.gen.data() + 2);
                std::sort(node.gen.data() + 2, node.gen.data() + 4);
                b_vert_index_set.insert(gen[4]);
                b_vert_index_set.insert(gen[5]);
            } else if (gen[1] < 0) {
                /// Abort if this fails. Cause of this error is difficult to understand for now.
                if (cell_facet_to_facet_map.find(-gen[1] - 1) == cell_facet_to_facet_map.end()) {
                    assert(0);
                    return false;
                }
                node.type = NodeType::B_FACE;
                node.gen(0) = cell_facet_to_facet_map.at(-gen[1] - 1);
                node.gen(1) = gen[0];
                node.gen(2) = gen[2] - 1;
                node.gen(3) = gen[3] - 1;
                std::sort(node.gen.data() + 1, node.gen.data() + 4);
                const BoundaryFace &b_face = model.boundary_data.f[node.gen[0]];
                b_vert_index_set.insert(b_face.vert_indices(0));
                b_vert_index_set.insert(b_face.vert_indices(1));
                b_vert_index_set.insert(b_face.vert_indices(2));
            } else {
                assert(gen[1] > 0 && gen[2] > 0 && gen[3] > 0);
                node.type = NodeType::STANDARD;
                node.gen(0) = gen[0];
                node.gen(1) = gen[1] - 1;
                node.gen(2) = gen[2] - 1;
                node.gen(3) = gen[3] - 1;
                std::sort(node.gen.data(), node.gen.data() + 4);
            }

            switch (node.type) {
                case B_VERTEX:
                    break;
                case STANDARD:
                    node.gen(0) = non_removed_cell_indices[node.gen(0)];
                case B_FACE:
                    node.gen(1) = non_removed_cell_indices[node.gen(1)];
                case B_EDGE:
                    node.gen(2) = non_removed_cell_indices[node.gen(2)];
                    node.gen(3) = non_removed_cell_indices[node.gen(3)];
            }

            face.node_indices.push_back(getNodeIndexAndAddToMap(node_index_map, node));
        }

        /// Determine the bisector or boundary face containing this face.
        if (facet_seed_id[f] >= 0) {
            face.is_boundary_face = false;
            face.opposite_gen_index = facet_seed_id[f];
        } else {
            assert(b_vert_index_set.size() == 3);
            face.is_boundary_face = true;
            std::vector<int> bv(b_vert_index_set.begin(), b_vert_index_set.end());
            face.opposite_gen_index = tri_to_bface_index.at({bv[0], bv[1], bv[2]});
        }

        int cell_idx = vertex_gen[clippedMesh.facets.vertex(f, 0) * 6 + 0];
        cell_idx = non_removed_cell_indices[cell_idx];
        model.cells[cell_idx].faces.push_back(face);
    }

    /// Fail if a site has no cell.
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

bool Voronoi3D::computeCells(Model &model) const {
    return computeCellsClippedVoronoi3D(model, false);
}
