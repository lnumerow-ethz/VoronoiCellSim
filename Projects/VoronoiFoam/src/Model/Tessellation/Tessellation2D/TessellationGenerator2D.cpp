#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/TessellationGenerator2D.h"

NodeData TessellationGenerator2D::computeNodeData(const Model &model, const TessellationNode &node, int order) const {
    int dims_c = getNumSiteDOFs();
    int dims_v = model.dimensions_ind->dims_space;
    int dims_x = model.dimensions_ind->dims_space;

    /// Get inputs vector.
    VectorXF inputs;
    switch (node.type) {
        case STANDARD:
            inputs.resize(3 * dims_c);
            for (int i = 0; i < 3; i++) {
                inputs.segment(i * dims_c, dims_c) =
                    getSiteDOFsInTessellation(model.degrees_of_freedom.sites[node.gen[i]]);
            }
            break;
        case B_EDGE:
            /// Input order is boundary vertices first, then sites.
            inputs.resize(2 * dims_v + 2 * dims_c);
            for (int i = 0; i < 2; i++) {
                inputs.segment(i * dims_v, dims_v) =
                    model.boundary_data.v[model.boundary_data.f[node.gen[0]].vert_indices(i)].pos;
            }
            for (int i = 0; i < 2; i++) {
                inputs.segment(2 * dims_v + i * dims_c, dims_c) =
                    getSiteDOFsInTessellation(model.degrees_of_freedom.sites[node.gen[i + 1]]);
            }
            break;
        case B_VERTEX:
            inputs = model.boundary_data.v[node.gen[0]].pos;
            break;
        default:
            assert(false);
    }

    /// Initialize node data matrices according to size of input vector.
    int n_inputs = (int)inputs.rows();
    NodeData node_data(dims_x, n_inputs, order);

    /// Finally calculate node data.
    switch (node.type) {
        case STANDARD:
            computeNodePosStandard(inputs, node_data);
            if (order >= 1) computeNodeGradStandard(inputs, node_data);
            if (order >= 2) computeNodeHessStandard(inputs, node_data);
            break;
        case B_EDGE:
            computeNodePosBEdge(inputs, node_data);
            if (order >= 1) computeNodeGradBEdge(inputs, node_data);
            if (order >= 2) computeNodeHessBEdge(inputs, node_data);
            break;
        case B_VERTEX:
            /// This case is trivial.
            node_data.pos = inputs;
            if (order >= 1) node_data.grad.setIdentity();
            if (order >= 2) {
                for (int i = 0; i < dims_x; i++) {
                    node_data.hess[i].setZero();
                }
            }
            break;
        default:
            assert(false);
    }

    return node_data;
}

void TessellationGenerator2D::getNodeDependencies(const Model &model, const TessellationNode &node,
                                                  std::vector<NodeDependency> &node_dependencies) const {
    int dims_c = getNumSiteDOFs();
    int dims_v = model.dimensions_ind->dims_space;

    /// Magic numbers in the switch statement below correspond to the indices of elements
    /// in node.gen and are explained in the comments.
    switch (node.type) {
        case STANDARD:  /// node.gen is 3 sites.
            for (int i = 0; i < 3; i++) {
                node_dependencies.emplace_back(NODE_GEN_SITE, node.gen[i], i * dims_c);
            }
            break;
        case B_EDGE:  /// node.gen is [b_edge, site, site].
            for (int i = 0; i < 2; i++) {
                node_dependencies.emplace_back(NODE_GEN_BVERTEX, model.boundary_data.f[node.gen[0]].vert_indices(i),
                                               i * dims_v);
            }
            for (int i = 0; i < 2; i++) {
                node_dependencies.emplace_back(NODE_GEN_SITE, node.gen[i + 1], 2 * dims_v + i * dims_c);
            }
            break;
        case B_VERTEX:  /// node.gen is one b_vert.
            node_dependencies.emplace_back(NODE_GEN_BVERTEX, node.gen[0], 0);
            break;
        default:
            assert(false);
    }
}
