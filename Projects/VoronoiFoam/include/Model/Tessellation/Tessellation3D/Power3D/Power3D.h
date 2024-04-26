#pragma once

#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/Voronoi3D/Voronoi3D.h"

class Power3D : public Voronoi3D {
   private:
    [[nodiscard]] bool computeCells(Model &model) const override;

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
    [[nodiscard]] int getNumSiteDOFs() const override { return 4; }

   public:
    bool getSiteParamIndex(SiteParamTessellation param_id, int &index) const override;

    [[nodiscard]] std::string getName() const override { return "Power"; };
};
