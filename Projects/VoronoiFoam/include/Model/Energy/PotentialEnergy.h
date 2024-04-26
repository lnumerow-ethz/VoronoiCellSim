#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"

namespace PotentialEnergy {
void computePotentialEnergy(const Model &model, F &energy, VectorXF &gradient, HessianF &hessian, int order);

void computePotentialEnergy(const Model &model, F &energy);

void computePotentialEnergyGradient(const Model &model, VectorXF &gradient);

void computePotentialEnergyHessian(const Model &model, HessianF &hessian);

bool computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                         const VectorXF &y, F &energy);

bool computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                         const VectorXF &y, F &energy, VectorXF &gradient);

bool computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                                         const VectorXF &y, F &energy, VectorXF &gradient, HessianF &hessian);
}  // namespace PotentialEnergy
