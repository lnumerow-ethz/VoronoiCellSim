#pragma once

#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

class TessellationGenerator3D : public TessellationGenerator {
   private:
    [[nodiscard]] NodeData computeNodeData(const Model &model, const TessellationNode &node, int order) const override;

    void getNodeDependencies(const Model &model, const TessellationNode &node,
                             std::vector<NodeDependency> &node_dependencies) const override;

   private:
    virtual void computeNodePosStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeGradStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeHessStandard(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodePosBFace(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeGradBFace(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeHessBFace(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodePosBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeGradBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;

    virtual void computeNodeHessBEdge(const VectorXF &inputs, NodeData &node_data) const = 0;
};
