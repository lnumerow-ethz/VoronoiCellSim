#pragma once

#include "CRLHelper/VecMatDef.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"

/// Triplets to construct sparse tessellation derivative matrices.
struct TessellationDerivativeTriplets {
    tbb::concurrent_vector<TripletF> dxdc;
    tbb::concurrent_vector<TripletF> dxdv;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdc2;  /// Only upper triangular coefficients.
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdcdv;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdv2;  /// Only upper triangular coefficients.
    tbb::concurrent_vector<TripletF> dvdp;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2vdp2;  /// Only upper triangular coefficients.
};

class TessellationGenerator {
   private:
    /// Returns coordinates and gradients to specified order for a single node.
    [[nodiscard]] virtual NodeData computeNodeData(const Model &model, const TessellationNode &node,
                                                   int order) const = 0;

    /// Assembles node_data_map from nodes in cells.
    void computeAllNodesData(Model &model, int order) const;

    /// Assembles gradient matrices for tessellation nodes with respect to system DOFs.
    void computeTessellationMatrices(Model &model, int order) const;

   protected:
    enum NodeGeneratorType {
        NODE_GEN_SITE,  //
        NODE_GEN_BVERTEX
    };

    struct NodeDependency {
        /// Is this a site or a boundary vertex.
        NodeGeneratorType type;
        /// The index of the element (site or boundary vertex) in the respective container in Model.
        int overall_idx;
        /// The starting row/column in NodeData gradient matrices corresponding to this dependency.
        int node_data_idx;

        NodeDependency(NodeGeneratorType type, int overall_idx, int node_data_idx)
            : type(type), overall_idx(overall_idx), node_data_idx(node_data_idx) {}
    };

   private:
    virtual void getNodeDependencies(const Model &model, const TessellationNode &node,
                                     std::vector<NodeDependency> &node_dependencies) const = 0;

    /// Add triplets for derivatives of tessellation nodes x wrt sites c and boundary vertices v.
    void getNodeTriplets(const Model &model, TessellationDerivativeTriplets &triplets, int order) const;

    /// Add triplets for derivatives of boundary vertices v wrt free boundary parameters p.
    static void getBoundaryTriplets(const Model &model, TessellationDerivativeTriplets &triplets, int order);

    static void constructMatricesFromTriplets(Model &model, TessellationDerivativeTriplets &triplets, int order);

   public:
    /// Should be the only public interface to this class, and is called once to fully populate Model given only the
    /// DOFs.
    [[nodiscard]] bool generateTessellation(Model &model, int order) const;

   private:
    /// Assign site_index and node_indices_in_cell for all cells.
    static void assignCellIndices(Model &model);

    [[nodiscard]] virtual bool computeCells(Model &model) const = 0;

   private:
    static DimensionsIndependent getDimensionsIndependent(Model &model);

    static DimensionsBoundary getDimensionsBoundary(Model &model);

    static DimensionsTessellation getDimensionsTessellation(Model &model);

   public:
    [[nodiscard]] virtual int getNumSiteDOFs() const = 0;

    [[nodiscard]] VectorXF getSiteDOFsInTessellation(const Site &site) const;

    virtual bool getSiteParamIndex(SiteParamTessellation param_id, int &index) const { return false; }

    bool hasSiteParam(SiteParamTessellation param_id) const;

    [[nodiscard]] virtual std::string getName() const = 0;
};
