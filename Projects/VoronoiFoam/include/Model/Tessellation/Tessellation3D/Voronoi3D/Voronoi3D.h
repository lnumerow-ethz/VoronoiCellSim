#pragma once

#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/TessellationGenerator3D.h"

class Voronoi3D : public TessellationGenerator3D {
   private:
    [[nodiscard]] bool computeCells(Model &model) const override;

   protected:
    [[nodiscard]] static bool computeCellsClippedVoronoi3D(Model &model, bool use_weights);

   private:
    void computeNodePosStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeGradStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeHessStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodePosBFace(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeGradBFace(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeHessBFace(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodePosBEdge(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeGradBEdge(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeHessBEdge(const VectorXF &inputs, NodeData &node_data) const override;

   protected:
    [[nodiscard]] int getNumSiteDOFs() const override { return 3; }

   public:
    [[nodiscard]] std::string getName() const override { return "Voronoi"; };
};
