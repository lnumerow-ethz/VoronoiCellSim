#pragma once

#include <utility>

#include "CRLHelper/VecMatDef.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>

/// Node type specifier.
enum NodeType {
    STANDARD,  //
    B_FACE,
    B_EDGE,
    B_VERTEX
};

/// Computed data for a node.
struct NodeData {
    VectorXF pos;
    MatrixXF grad;
    std::vector<MatrixXF> hess;

    NodeData(int dims_x, int n_inputs, int order);
};

/// This data defines a unique node for identification and sorting.
struct TessellationNode {
    NodeType type;
    VectorXI gen;
    std::shared_ptr<const NodeData> data;
};

bool operator<(const TessellationNode &a, const TessellationNode &b);

/// Face data structure. Faces are edges in 2D and facets in 3D. Face belongs to a cell.
/// Nodes oriented CCW (in 3D, normal direction faces out of the cell).
struct TessellationFace {
    bool is_boundary_face = false;
    int opposite_gen_index = -1;
    std::vector<int> node_indices;
};
/// Cell data structure. Faces unordered.
struct TessellationCell {
    int site_index;
    std::vector<TessellationFace> faces;

    /// Each node mapped to a unique index for cell energy gradient and hessian matrices.
    std::map<int, int> node_indices_in_cell;
};

/// Boundary vertex data.
struct BoundaryVertex {
    VectorXF pos;
    SparseMatrixF grad;
    std::vector<SparseMatrixF> hess;
};

/// Boundary face data structure. Vertices CCW (in 3D, normal direction faces out of domain).
struct BoundaryFace {
    VectorXI vert_indices;
};

/// Wrapper for boundary data.
struct BoundaryData {
    std::vector<BoundaryFace> f;
    std::vector<BoundaryVertex> v;
};

enum SiteParamTessellation {
    SITE_PARAM_POWER_WEIGHT,  //
    NUM_SITE_PARAM_TESSELLATION
};

enum SiteParamEnergy {
    SITE_PARAM_SIZE_TARGET = NUM_SITE_PARAM_TESSELLATION,  //
    SITE_PARAM_SURFACE_TARGET,
    NUM_SITE_PARAM_TOTAL
};

/// Site degrees of freedom.
struct Site {
    /// Allow removing sites without resizing everything. May need a cull function in the future.
    bool is_removed = false;

    VectorXF pos;
    VectorXF param = VectorXF::Zero(NUM_SITE_PARAM_TOTAL);

    Site() {
        param(SITE_PARAM_POWER_WEIGHT) = 0;
        param(SITE_PARAM_SIZE_TARGET) = 0;
        param(SITE_PARAM_SURFACE_TARGET) = 0;
    }
};

class TessellationGenerator;

class BoundaryGenerator;

class PerCellFunction;

/// Model hyperparameters.
struct ModelDefinition {
    /// Indices of free boundary parameters. Other parameters will be fixed during energy minimization etc.
    VectorXI boundary_free_param_indices;
    // Free site parameters. Other parameters will be fixed during energy minimization etc.
    VectorXI site_free_param_indices;

    std::shared_ptr<BoundaryGenerator> boundary_generator;
    std::shared_ptr<TessellationGenerator> tessellation_generator;
    std::shared_ptr<PerCellFunction> cell_energy_function;
};

/// System degrees of freedom.
struct DegreesOfFreedom {
    std::vector<Site> sites;
    VectorXF boundary_param;
};

/// System dimensions determined from DOF.
struct DimensionsIndependent {
    int dims_space = 0;
    int dims_site_params_free = 0;
    int dims_site_dof = 0;

    int dims_site_dof_in_tessellation = 0;

    int n_sites = 0;

    int n_b_params_total = 0;
    int n_b_params_free = 0;

    int nc = 0;
    /// The number of FREE boundary params.
    int np = 0;
};

/// System dimensions determined from boundary generation.
struct DimensionsBoundary {
    int n_b_faces = 0;
    int n_b_vertices = 0;

    int nv = 0;
};

/// System dimensions determined from tessellation generation.
struct DimensionsTessellation {
    int n_nodes = 0;

    int nx = 0;
};

/// Tessellation derivatives.
struct TessellationDerivativeMatrices {
    SparseMatrixF dxdc;
    SparseMatrixF dxdv;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdc2;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdcdv;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2xdv2;
    SparseMatrixF dvdp;
    tbb::concurrent_vector<tbb::concurrent_vector<TripletF>> d2vdp2;
};

/// Complete system model including independent variables (dofs) and dependent/computed data (tessellation).
struct Model {
    const ModelDefinition model_definition;
    const DegreesOfFreedom degrees_of_freedom;

    BoundaryData boundary_data;

    std::vector<TessellationCell> cells;
    std::vector<TessellationNode> nodes;

    /// These are pointers so that the dimensions objects will be nullptr until initialized.
    std::shared_ptr<DimensionsIndependent> dimensions_ind;
    std::shared_ptr<DimensionsBoundary> dimensions_boundary;
    std::shared_ptr<DimensionsTessellation> dimensions_tess;

    TessellationDerivativeMatrices derivative_matrices;

    /// 'success' will be assigned false if some part of model generation failed.
    Model(ModelDefinition model_definition, DegreesOfFreedom degrees_of_freedom, int order, bool &success);
};
