#include <igl/triangle/triangulate.h>

#include "CRLHelper/GeometryHelper.h"

#include "Projects/VoronoiFoam/include/App/Rendering/Rendering.h"

#include <iostream>
#include <fstream>
#include <set>

void Rendering::writeFile2D(const Model& model, std::string file_name) {
    std::ofstream file(file_name);

    int n_nodes = model.dimensions_tess->n_nodes;
    int n_sites = model.dimensions_ind->n_sites;

    // Write nodes
    file << n_nodes << "\n";
    for (int i = 0; i < n_nodes; i++) {
        Vector2F node = model.nodes[i].data->pos;
        file << node(0) << " " << node(1) << " " << 0 << "\n";
    }

    MatrixXF igl_V;
    MatrixXI igl_F;
    VectorXI cell_idx;
    triangulateVoronoiCells2D(model, igl_V, igl_F, cell_idx);

    // Write cells
    file << n_sites << "\n";
    for (const TessellationCell& cell : model.cells) {
        file << model.degrees_of_freedom.sites[cell.site_index].pos(0) << " "
             << model.degrees_of_freedom.sites[cell.site_index].pos(1) << " ";

        std::vector<Vector2I> edges;
        std::vector<Vector3I> triangles;

        for (const TessellationFace& face : cell.faces) {
            edges.emplace_back(face.node_indices[0], face.node_indices[1]);
        }
        for (int i = 0; i < (int)edges.size(); i++) {
            for (int j = i + 2; j < (int)edges.size(); j++) {
                if (edges[j](0) == edges[i](1)) {
                    std::swap(edges[j], edges[i + 1]);
                    break;
                }
            }
        }
        for (int i = 0; i < igl_F.rows(); i++) {
            if (cell_idx(i) == cell.site_index) {
                triangles.emplace_back(igl_F.row(i));
            }
        }

        file << edges.size() << " ";
        for (const auto& edge : edges) {
            file << edge(0) << " " << edge(1) << " ";
        }

        file << triangles.size() << " ";
        for (const auto& tri : triangles) {
            file << tri(0) << " " << tri(1) << " " << tri(2) << " ";
        }

        file << "\n";
    }

    file.close();
}

void Rendering::writeFile3D(const Model& model, std::string file_name, bool include_internal_faces,
                            bool include_boundary_edges) {
    std::ofstream file(file_name);

    int n_nodes = model.dimensions_tess->n_nodes;
    int n_sites = model.dimensions_ind->n_sites;

    // Write nodes
    file << n_nodes << "\n";
    for (int i = 0; i < n_nodes; i++) {
        Vector3F node = model.nodes[i].data->pos;
        file << node(0) << " " << node(1) << " " << node(2) << "\n";
    }

    // Write cells
    file << n_sites << "\n";
    for (const TessellationCell& cell : model.cells) {
        file << model.degrees_of_freedom.sites[cell.site_index].pos(0) << " "
             << model.degrees_of_freedom.sites[cell.site_index].pos(1) << " "
             << model.degrees_of_freedom.sites[cell.site_index].pos(2) << " ";

        std::vector<Vector2I> edges;
        std::vector<Vector3I> triangles;

        for (const TessellationFace& face : cell.faces) {
            if (face.is_boundary_face || include_internal_faces) {
                for (size_t i = 1; i < face.node_indices.size() - 1; i++) {
                    triangles.emplace_back(face.node_indices[0], face.node_indices[i], face.node_indices[i + 1]);
                }
            }
            if (!face.is_boundary_face || include_boundary_edges) {
                for (size_t i = 0; i < face.node_indices.size(); i++) {
                    edges.emplace_back(face.node_indices[i], face.node_indices[(i + 1) % face.node_indices.size()]);
                }
            }
        }

        file << edges.size() << " ";
        for (const auto& edge : edges) {
            file << edge(0) << " " << edge(1) << " ";
        }

        file << triangles.size() << " ";
        for (const auto& tri : triangles) {
            file << tri(0) << " " << tri(1) << " " << tri(2) << " ";
        }

        file << "\n";
    }

    file.close();
}

void Rendering::writeFile(const Model& model, std::string file_name) {
    int dims_space = model.dimensions_ind->dims_space;
    switch (dims_space) {
        case 2:
            writeFile2D(model, file_name);
            break;
        case 3:
            writeFile3D(model, file_name, true, true);
            break;
        default:
            assert(0);
            break;
    }
}

/// Obtain a triangulation of the cells for visualization. Output contains cell index of each triangle.
void Rendering::triangulateVoronoiCells2D(const Model& model, MatrixXF& igl_V, MatrixXI& igl_F, VectorXI& cell_idx) {
    int n_nodes = model.dimensions_tess->n_nodes;

    std::map<int, std::set<int>> node_idx_to_adjacent_sites_map;

    std::vector<std::pair<int, int>> edge_node_index_pairs;
    for (const TessellationCell& cell : model.cells) {
        for (const TessellationFace& face : cell.faces) {
            int node0_idx = face.node_indices[0];
            int node1_idx = face.node_indices[1];

            /// Each cell should be counted exactly once per vertex here, so don't repeat for node1_idx.
            node_idx_to_adjacent_sites_map[node0_idx].insert(cell.site_index);

            /// To avoid counting Voronoi edges twice, only add them if opposite cell index is higher than this one.
            if (face.is_boundary_face || cell.site_index < face.opposite_gen_index) {
                edge_node_index_pairs.emplace_back(node0_idx, node1_idx);
            }
        }
    }

    MatrixXF triangulate_P(n_nodes, 2);
    for (int i = 0; i < n_nodes; i++) {
        triangulate_P.row(i) = model.nodes[i].data->pos;
    }
    igl_V = MatrixXF::Zero(n_nodes, 3);
    igl_V.leftCols(2) = triangulate_P;

    MatrixXI triangulate_E(edge_node_index_pairs.size(), 2);
    for (int i = 0; i < (int)edge_node_index_pairs.size(); i++) {
        const auto& edge = edge_node_index_pairs[i];
        triangulate_E.row(i) = Vector2I(edge.first, edge.second);
    }

    /// Compute offset points from each boundary face to determine out-of-bounds markers, or "hole" points.
    MatrixXF triangulate_H(model.boundary_data.f.size(), 2);
    for (int i = 0; i < (int)model.boundary_data.f.size(); i++) {
        Vector2F v0 = model.boundary_data.v[model.boundary_data.f[i].vert_indices(0)].pos;
        Vector2F v1 = model.boundary_data.v[model.boundary_data.f[i].vert_indices(1)].pos;
        Vector2F offset_point = (v0 + v1) / 2.0 + 1e-8 * GeometryHelper::rotateVector2D(v1 - v0, -M_PI_2);
        triangulate_H.row(i) = offset_point;
    }

    MatrixXF discard_output_V;
    MatrixXI triangulate_F;
    igl::triangle::triangulate(triangulate_P, triangulate_E, triangulate_H, "Q", discard_output_V, igl_F);

    std::vector<int> cell_idx_vector;
    for (int i = 0; i < igl_F.rows(); i++) {
        auto sites_indices_0 = node_idx_to_adjacent_sites_map.at(igl_F(i, 0));
        auto sites_indices_1 = node_idx_to_adjacent_sites_map.at(igl_F(i, 1));
        auto sites_indices_2 = node_idx_to_adjacent_sites_map.at(igl_F(i, 2));

        /// Common site should be found unless this triangle is part of a hole in the domain.
        int common_site = -1;
        for (int site_index : sites_indices_0) {
            if (sites_indices_1.find(site_index) != sites_indices_1.end() &&
                sites_indices_2.find(site_index) != sites_indices_2.end()) {
                common_site = site_index;
                break;
            }
        }

        assert(common_site != -1);
        cell_idx_vector.push_back(common_site);
    }

    int num_output_triangles = (int)cell_idx_vector.size();
    cell_idx = Eigen::Map<VectorXI>(cell_idx_vector.data(), num_output_triangles);
}
