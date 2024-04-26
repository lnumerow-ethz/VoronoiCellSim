#pragma once

#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/TessellationGenerator2D.h"

class Voronoi2D : public TessellationGenerator2D {
   private:
    [[nodiscard]] bool computeCells(Model &model) const override;

   protected:
    [[nodiscard]] static bool computeCellsClippedVoronoi2D(Model &model, bool use_weights);

   private:
    void computeNodePosStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeGradStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeHessStandard(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodePosBEdge(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeGradBEdge(const VectorXF &inputs, NodeData &node_data) const override;

    void computeNodeHessBEdge(const VectorXF &inputs, NodeData &node_data) const override;

   protected:
    [[nodiscard]] int getNumSiteDOFs() const override { return 2; }

   public:
    [[nodiscard]] std::string getName() const override { return "Voronoi"; };
};
