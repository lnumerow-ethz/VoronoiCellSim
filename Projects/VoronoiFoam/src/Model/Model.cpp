#include <utility>

#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

Model::Model(ModelDefinition model_definition, DegreesOfFreedom degrees_of_freedom, int order, bool &success)
    : model_definition(std::move(model_definition)), degrees_of_freedom(std::move(degrees_of_freedom)) {
    /// The model_definition from the argument has been cleared by std::move, need to use this->model_definition.
    success = this->model_definition.tessellation_generator->generateTessellation(*this, order);
}

bool operator<(const TessellationNode &a, const TessellationNode &b) {
    if (a.type != b.type) return a.type < b.type;
    for (int i = 0; i < a.gen.rows(); i++) {
        if (a.gen[i] != b.gen[i]) return a.gen[i] < b.gen[i];
    }
    return false;
}

NodeData::NodeData(int dims_x, int n_inputs, int order) {
    pos.resize(dims_x);
    if (order >= 1) grad.resize(dims_x, n_inputs);
    if (order >= 2) {
        hess.resize(dims_x);
        for (int i = 0; i < dims_x; i++) {
            hess[i].resize(n_inputs, n_inputs);
        }
    }
}
