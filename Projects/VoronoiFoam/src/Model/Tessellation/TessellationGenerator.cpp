#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>
#include "CRLHelper/CRLTimer.h"

bool TessellationGenerator::generateTessellation(Model &model, int order) const {
    bool success;

    model.dimensions_ind = std::make_shared<DimensionsIndependent>(getDimensionsIndependent(model));

    success = model.model_definition.boundary_generator->generateBoundary(model.degrees_of_freedom, model.boundary_data,
                                                                          order);
    if (!success) {
        return false;
    }

    model.dimensions_boundary = std::make_shared<DimensionsBoundary>(getDimensionsBoundary(model));

    success = computeCells(model);
    if (!success) {
        return false;
    }

    assignCellIndices(model);
    computeAllNodesData(model, order);

    model.dimensions_tess = std::make_shared<DimensionsTessellation>(getDimensionsTessellation(model));

    if (order >= 1) {
        computeTessellationMatrices(model, order);
    }

    return success;
}

DimensionsIndependent TessellationGenerator::getDimensionsIndependent(Model &model) {
    DimensionsIndependent dimensions;

    dimensions.dims_space = (int)model.degrees_of_freedom.sites[0].pos.rows();
    dimensions.dims_site_params_free = (int)model.model_definition.site_free_param_indices.rows();
    dimensions.dims_site_dof = dimensions.dims_space + dimensions.dims_site_params_free;

    dimensions.dims_site_dof_in_tessellation = model.model_definition.tessellation_generator->getNumSiteDOFs();

    dimensions.n_sites = (int)model.degrees_of_freedom.sites.size();

    dimensions.n_b_params_total = (int)model.degrees_of_freedom.boundary_param.rows();
    dimensions.n_b_params_free = (int)model.model_definition.boundary_free_param_indices.rows();

    dimensions.nc = dimensions.dims_site_dof * dimensions.n_sites;
    dimensions.np = dimensions.n_b_params_free;

    return dimensions;
}

DimensionsBoundary TessellationGenerator::getDimensionsBoundary(Model &model) {
    DimensionsBoundary dimensions;

    dimensions.n_b_vertices = (int)model.boundary_data.v.size();
    dimensions.n_b_faces = (int)model.boundary_data.f.size();

    dimensions.nv = dimensions.n_b_vertices * model.dimensions_ind->dims_space;

    return dimensions;
}

DimensionsTessellation TessellationGenerator::getDimensionsTessellation(Model &model) {
    DimensionsTessellation dimensions;

    dimensions.n_nodes = (int)model.nodes.size();

    dimensions.nx = dimensions.n_nodes * (int)model.dimensions_ind->dims_space;

    return dimensions;
}

void TessellationGenerator::assignCellIndices(Model &model) {
    int ic = 0;
    for (TessellationCell &cell : model.cells) {
        cell.site_index = ic;
        ic++;

        int i = 0;
        for (const TessellationFace &face : cell.faces) {
            for (int n : face.node_indices) {
                /// If node does not already in node_indices_in_cell map, add it and increment index.
                if (cell.node_indices_in_cell.find(n) == cell.node_indices_in_cell.end()) {
                    cell.node_indices_in_cell[n] = i;
                    i++;
                }
            }
        }
    }
}

void TessellationGenerator::computeAllNodesData(Model &model, int order) const {
    for (TessellationNode &node : model.nodes) {
        node.data = std::make_shared<NodeData>(computeNodeData(model, node, order));
    }
}

void TessellationGenerator::getBoundaryTriplets(const Model &model, TessellationDerivativeTriplets &triplets,
                                                int order) {
    /// Derivatives stored in model.boundary_v are wrt all boundary params (not just free params p).
    /// Use this map to convert overall param indices to free param indices.
    VectorXI p_index_map = ModelHelper::getBoundaryParamReverseIndexVector(model);
    int dims_v = model.dimensions_ind->dims_space;

    for (int iv = 0; iv < model.dimensions_boundary->n_b_vertices; iv++) {
        const BoundaryVertex &vertex = model.boundary_data.v[iv];

        for (int i = 0; i < vertex.grad.outerSize(); i++) {
            for (SparseMatrixF::InnerIterator it(vertex.grad, i); it; ++it) {
                int row = (int)it.row();
                int col = p_index_map(it.col());
                if (col >= 0) {
                    triplets.dvdp.emplace_back(iv * dims_v + row, col, it.value());
                }
            }
        }

        if (order >= 2) {
            for (int i = 0; i < dims_v; i++) {
                for (int j = 0; j < vertex.hess[i].outerSize(); j++) {
                    for (SparseMatrixF::InnerIterator it(vertex.hess[i], j); it; ++it) {
                        int row = p_index_map(it.row());
                        int col = p_index_map(it.col());
                        if (col < row) continue;  /// Symmetric matrix, only store upper triangular.
                        if (row >= 0 && col >= 0) triplets.d2vdp2[iv * dims_v + i].emplace_back(row, col, it.value());
                    }
                }
            }
        }
    }
}

void TessellationGenerator::getNodeTriplets(const Model &model, TessellationDerivativeTriplets &triplets,
                                            int order) const {
    int dims_c = model.dimensions_ind->dims_site_dof;
    int dims_c_tess = model.dimensions_ind->dims_site_dof_in_tessellation;
    int dims_v = model.dimensions_ind->dims_space;
    int dims_x = model.dimensions_ind->dims_space;

    VectorXI c_index_map = ModelHelper::getTessellationDOFIndexToSiteDOFIndexVector(model);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, model.nodes.size(), 2000), [&](const tbb::blocked_range<size_t> &range) {
            TripletListF local_triplets_dxdc;
            TripletListF local_triplets_dxdv;

            for (size_t node_idx = range.begin(); node_idx != range.end(); ++node_idx) {
                const TessellationNode &node = model.nodes[node_idx];

                auto &data = *node.data;

                std::vector<NodeDependency> node_dependencies;
                getNodeDependencies(model, node, node_dependencies);

                for (const NodeDependency &dependency0 : node_dependencies) {
                    switch (dependency0.type) {
                        case NODE_GEN_SITE:
                            for (int i = 0; i < dims_x; i++) {
                                for (int j = 0; j < dims_c_tess; j++) {
                                    if (c_index_map(j) == -1) continue;
                                    local_triplets_dxdc.emplace_back(node_idx * dims_x + i,
                                                                     dependency0.overall_idx * dims_c + c_index_map(j),
                                                                     data.grad(i, dependency0.node_data_idx + j));
                                }
                            }

                            if (order >= 2) {
                                for (const NodeDependency &dependency1 : node_dependencies) {
                                    switch (dependency1.type) {
                                        case NODE_GEN_SITE:  /// Two site dependencies
                                            for (int i = 0; i < dims_x; i++) {
                                                for (int j = 0; j < dims_c_tess; j++) {
                                                    if (c_index_map(j) == -1) continue;
                                                    int row = dependency0.overall_idx * dims_c + c_index_map(j);
                                                    for (int k = 0; k < dims_c_tess; k++) {
                                                        if (c_index_map(k) == -1) continue;
                                                        int col = dependency1.overall_idx * dims_c + c_index_map(k);
                                                        if (col < row) continue;  /// Only store upper triangular.
                                                        triplets.d2xdc2[node_idx * dims_x + i].emplace_back(
                                                            dependency0.overall_idx * dims_c + c_index_map(j),
                                                            dependency1.overall_idx * dims_c + c_index_map(k),
                                                            data.hess[i](dependency0.node_data_idx + j,
                                                                         dependency1.node_data_idx + k));
                                                    }
                                                }
                                            }
                                            break;
                                        case NODE_GEN_BVERTEX:  /// One site and one boundary vertex dependency
                                            for (int i = 0; i < dims_x; i++) {
                                                for (int j = 0; j < dims_c_tess; j++) {
                                                    for (int k = 0; k < dims_v; k++) {
                                                        if (c_index_map(j) == -1) continue;
                                                        triplets.d2xdcdv[node_idx * dims_x + i].emplace_back(
                                                            dependency0.overall_idx * dims_c + c_index_map(j),
                                                            dependency1.overall_idx * dims_v + k,
                                                            data.hess[i](dependency0.node_data_idx + j,
                                                                         dependency1.node_data_idx + k));
                                                    }
                                                }
                                            }
                                            break;
                                    }
                                }
                            }

                            break;
                        case NODE_GEN_BVERTEX:
                            for (int i = 0; i < dims_x; i++) {
                                for (int j = 0; j < dims_v; j++) {
                                    local_triplets_dxdv.emplace_back(node_idx * dims_x + i,
                                                                     dependency0.overall_idx * dims_v + j,
                                                                     data.grad(i, dependency0.node_data_idx + j));
                                }
                            }

                            if (order >= 2) {
                                for (const NodeDependency &dependency1 : node_dependencies) {
                                    switch (dependency1.type) {
                                        case NODE_GEN_SITE:  /// Redundant case (one site and one boundary vertex), do
                                                             /// nothing.
                                            break;
                                        case NODE_GEN_BVERTEX:  /// Two boundary vertex dependencies.
                                            for (int i = 0; i < dims_x; i++) {
                                                for (int j = 0; j < dims_v; j++) {
                                                    int row = dependency0.overall_idx * dims_v + j;
                                                    for (int k = 0; k < dims_v; k++) {
                                                        int col = dependency1.overall_idx * dims_v + k;
                                                        if (col < row) continue;  /// Only store upper triangular.
                                                        triplets.d2xdv2[node_idx * dims_x + i].emplace_back(
                                                            row, col,
                                                            data.hess[i](dependency0.node_data_idx + j,
                                                                         dependency1.node_data_idx + k));
                                                    }
                                                }
                                            }
                                            break;
                                    }
                                }
                            }

                            break;
                    }
                }
            }
            triplets.dxdc.grow_by(local_triplets_dxdc.begin(), local_triplets_dxdc.end());
            triplets.dxdv.grow_by(local_triplets_dxdv.begin(), local_triplets_dxdv.end());
        });
}

void TessellationGenerator::constructMatricesFromTriplets(Model &model, TessellationDerivativeTriplets &triplets,
                                                          int order) {
    int nc = model.dimensions_ind->nc;
    int nx = model.dimensions_tess->nx;
    int nv = model.dimensions_boundary->nv;
    int np = model.dimensions_ind->np;

    model.derivative_matrices.dxdc.resize(nx, nc);
    model.derivative_matrices.dxdv.resize(nx, nv);
    model.derivative_matrices.dvdp.resize(nv, np);

    model.derivative_matrices.dxdc.setFromTriplets(triplets.dxdc.begin(), triplets.dxdc.end());
    model.derivative_matrices.dxdv.setFromTriplets(triplets.dxdv.begin(), triplets.dxdv.end());
    model.derivative_matrices.dvdp.setFromTriplets(triplets.dvdp.begin(), triplets.dvdp.end());

    if (order >= 2) {
        /// TODO: needless copy
        model.derivative_matrices.d2xdc2 = triplets.d2xdc2;
        model.derivative_matrices.d2xdcdv = triplets.d2xdcdv;
        model.derivative_matrices.d2xdv2 = triplets.d2xdv2;
        model.derivative_matrices.d2vdp2 = triplets.d2vdp2;
    }
}

void TessellationGenerator::computeTessellationMatrices(Model &model, int order) const {
    TessellationDerivativeTriplets triplets;
    if (order >= 2) {
        int nx = model.dimensions_tess->nx;
        int nv = model.dimensions_boundary->nv;

        triplets.d2xdc2.resize(nx);
        triplets.d2xdcdv.resize(nx);
        triplets.d2xdv2.resize(nx);
        triplets.d2vdp2.resize(nv);
    }

    getNodeTriplets(model, triplets, order);
    getBoundaryTriplets(model, triplets, order);
    constructMatricesFromTriplets(model, triplets, order);
}

VectorXF TessellationGenerator::getSiteDOFsInTessellation(const Site &site) const {
    int dims_space = site.pos.rows();

    VectorXF site_dof(getNumSiteDOFs());

    site_dof.head(dims_space) = site.pos;
    for (int i = 0; i < NUM_SITE_PARAM_TESSELLATION; i++) {
        int site_param_dof_index;
        if (getSiteParamIndex((SiteParamTessellation)i, site_param_dof_index)) {
            site_dof(site_param_dof_index) = site.param(i);
        }
    }

    return site_dof;
}

bool TessellationGenerator::hasSiteParam(SiteParamTessellation param_id) const {
    int site_param_index;
    return getSiteParamIndex(param_id, site_param_index);
}
