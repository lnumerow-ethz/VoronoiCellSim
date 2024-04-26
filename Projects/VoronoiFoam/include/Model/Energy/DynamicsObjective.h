#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

enum DynamicsDiscretization {
    FIRST_ORDER,  //
    SECOND_ORDER
};

/// Contains the additional data required for a model which evolves in time.
struct DynamicsData {
    DynamicsDiscretization discretization = FIRST_ORDER;
    F dt = DYNAMICS_DT;
    F tolerance = DYNAMICS_TOLERANCE;

    VectorXF mass_matrix_diagonal;
    VectorXF viscosity_matrix_diagonal;

    VectorXF prev_dof_1;
    VectorXF prev_dof_2;
    VectorXF prev_dof_3;
};

namespace DynamicsObjective {
bool computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                           const DynamicsData &dynamics_data, const VectorXF &y, F &energy);

bool computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                           const DynamicsData &dynamics_data, const VectorXF &y, F &energy,
                                           VectorXF &gradient);

bool computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                           const DynamicsData &dynamics_data, const VectorXF &y, F &energy,
                                           VectorXF &gradient, HessianF &hessian);

/// Wrapper functions to call either potential energy or dynamics depending on use_dynamics.
bool computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                               const DegreesOfFreedom &model_dof, const DynamicsData &dynamics_data,
                                               bool use_dynamics, const VectorXF &y, F &energy);

bool computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                               const DegreesOfFreedom &model_dof, const DynamicsData &dynamics_data,
                                               bool use_dynamics, const VectorXF &y, F &energy, VectorXF &gradient);

bool computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                               const DegreesOfFreedom &model_dof, const DynamicsData &dynamics_data,
                                               bool use_dynamics, const VectorXF &y, F &energy, VectorXF &gradient,
                                               HessianF &hessian);

void nextDynamicsStep(DynamicsData &dynamics_data, const VectorXF &y);
}  // namespace DynamicsObjective
