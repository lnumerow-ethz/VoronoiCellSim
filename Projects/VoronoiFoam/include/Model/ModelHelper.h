#pragma once

#include "Projects/VoronoiFoam/include/Model/Model.h"

namespace ModelHelper {
VectorXI getBoundaryParamReverseIndexVector(const Model &model);

void getDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof, VectorXF &dof_vector);

void setDOFFromVector(const ModelDefinition &model_definition, DegreesOfFreedom &model_dof, const VectorXF &dof_vector);

bool siteParamIdToSiteDOFIndex(const Model &model, int site_param_id, int &site_dof_index);

void setFreeSiteParam(ModelDefinition &model_definition, const bool site_param_free[]);

VectorXI getSiteParamReverseIndexVector(const Model &model);

VectorXI getTessellationDOFIndexToSiteDOFIndexVector(const Model &model);

/// Exclude sites from tessellation if they satisfy the condition for removal, currently size target below zero.
void markRemovedSites(const ModelDefinition &model_definition, DegreesOfFreedom &model_dof);
}  // namespace ModelHelper
