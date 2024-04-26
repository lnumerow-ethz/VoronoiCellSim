#pragma once

#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

class TessellationGenerator2D : public TessellationGenerator {
   private:
    [[nodiscard]] virtual NodeData computeNodeData(const Model &model, const TessellationNode &node, int order) const;

    virtual void getNodeDependencies(const Model &model, const TessellationNode &node,
                                     std::vector<NodeDependency> &node_dependencies) const;

   private:
    virtual void computeNodePosStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeGradStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeHessStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodePosBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeGradBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeHessBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;
};
