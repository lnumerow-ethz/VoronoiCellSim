#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/PerCellFunctionFromSimplex.h"

#include <iostream>

PerSimplexValue PerCellFunctionFromSimplex::getSimplexValue(const Model &model, const std::vector<int> &simplex_nodes,
                                                            int order) const {
    int dims_space = model.dimensions_ind->dims_space;

    /// A simplex consists of dims_space nodes (plus the origin, in a sense),
    /// and each node has dims_space coordinates.
    VectorXF inputs(dims_space * dims_space);
    for (int i = 0; i < dims_space; i++) {
        inputs.segment(i * dims_space, dims_space) = model.nodes[simplex_nodes[i]].data->pos;
    }

    PerSimplexValue simplex_value(dims_space);
    getSimplexValue(inputs, simplex_value);
    if (order >= 1) getSimplexGradient(inputs, simplex_value);
    if (order >= 2) getSimplexHessian(inputs, simplex_value);

    return simplex_value;
}

PerFaceValue PerCellFunctionFromSimplex::getFaceValue2D(const Model &model, const TessellationFace &face,
                                                        int order) const {
    const std::vector<int> simplex_nodes = {face.node_indices[0], face.node_indices[1]};
    PerSimplexValue simplex_value = getSimplexValue(model, simplex_nodes, order);

    PerFaceValue face_value;
    face_value.face_simplex_values.emplace_back(Vector2I(0, 1), simplex_value);

    return face_value;
}

PerFaceValue PerCellFunctionFromSimplex::getFaceValue3D(const Model &model, const TessellationFace &face,
                                                        int order) const {
    PerFaceValue face_value;

    for (int i = 1; i < (int)face.node_indices.size() - 1; i++) {
        Vector3I node_indices_in_face = {0, i, i + 1};
        const std::vector<int> simplex_nodes = {face.node_indices[node_indices_in_face[0]],
                                                face.node_indices[node_indices_in_face[1]],
                                                face.node_indices[node_indices_in_face[2]]};
        PerSimplexValue simplex_value = getSimplexValue(model, simplex_nodes, order);

        face_value.face_simplex_values.emplace_back(node_indices_in_face, simplex_value);
    }

    return face_value;
}

PerFaceValue PerCellFunctionFromSimplex::getFaceValue(const Model &model, const TessellationFace &face,
                                                      int order) const {
    int dims_space = model.dimensions_ind->dims_space;

    switch (dims_space) {
        case 2:
            return getFaceValue2D(model, face, order);
        case 3:
            return getFaceValue3D(model, face, order);
        default:
            return {};
    }
}

void PerCellFunctionFromSimplex::getValue(const Model &model, PerCellValue &cell_value) const {
    int dims_space = model.dimensions_ind->dims_space;
    int order = cell_value.order;

    const TessellationCell &cell = cell_value.cell;

    TripletListF hessian_triplets;
    for (TessellationFace face : cell.faces) {
        PerFaceValue face_value = getFaceValue(model, face, order);
        postProcessFaceValue(model, cell, face, order);

        for (const auto &face_simplex_value : face_value.face_simplex_values) {
            const VectorXI &node_indices_in_face = face_simplex_value.first;
            const PerSimplexValue &simplex_value = face_simplex_value.second;
            if (!simplex_value.allFinite()) continue;

            cell_value.value += simplex_value.value;

            if (order >= 1) {
                for (int i = 0; i < dims_space; i++) {
                    int node_index_in_cell_i = cell.node_indices_in_cell.at(face.node_indices[node_indices_in_face[i]]);
                    cell_value.gradient.segment(node_index_in_cell_i * dims_space, dims_space) +=
                        simplex_value.gradient.segment(i * dims_space, dims_space);
                }
            }

            if (order >= 2) {
                for (int i = 0; i < dims_space; i++) {
                    int node_index_in_cell_i = cell.node_indices_in_cell.at(face.node_indices[node_indices_in_face[i]]);
                    for (int j = 0; j < dims_space; j++) {
                        int node_index_in_cell_j =
                            cell.node_indices_in_cell.at(face.node_indices[node_indices_in_face[j]]);

                        // cell_value.hessian.block(node_index_in_cell_i * dims_space, node_index_in_cell_j *
                        // dims_space,
                        //                          dims_space, dims_space) +=
                        //     simplex_value.hessian.block(i * dims_space, j * dims_space, dims_space, dims_space);

                        for (int ii = 0; ii < dims_space; ii++) {
                            for (int jj = 0; jj < dims_space; jj++) {
                                hessian_triplets.emplace_back(
                                    node_index_in_cell_i * dims_space + ii, node_index_in_cell_j * dims_space + jj,
                                    simplex_value.hessian(i * dims_space + ii, j * dims_space + jj));
                            }
                        }
                    }
                }
            }
        }
    }
    cell_value.hessian.A.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
}
