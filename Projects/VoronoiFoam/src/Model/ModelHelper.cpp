#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

VectorXI ModelHelper::getBoundaryParamReverseIndexVector(const Model &model) {
    VectorXI free_index = VectorXI::Constant(model.degrees_of_freedom.boundary_param.rows(), -1);
    for (int i = 0; i < model.model_definition.boundary_free_param_indices.rows(); i++) {
        free_index(model.model_definition.boundary_free_param_indices(i)) = i;
    }
    return free_index;
}

void ModelHelper::getDOFVector(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                               VectorXF &dof_vector) {
    int dims_space = model_definition.boundary_generator->getDims();
    int n_site_free_param = model_definition.site_free_param_indices.rows();
    int dims_c = n_site_free_param + dims_space;
    int n_sites = (int)model_dof.sites.size();
    int nc = n_sites * dims_c;
    int np = (int)model_definition.boundary_free_param_indices.rows();

    dof_vector.resize(nc + np);
    for (int i = 0; i < n_sites; i++) {
        dof_vector.segment(i * dims_c, dims_space) = model_dof.sites[i].pos;
        for (int j = 0; j < n_site_free_param; j++) {
            dof_vector(i * dims_c + dims_space + j) =
                model_dof.sites[i].param(model_definition.site_free_param_indices(j));
        }
    }
    for (int i = 0; i < np; i++) {
        dof_vector(nc + i) = model_dof.boundary_param(model_definition.boundary_free_param_indices(i));
    }
}

void ModelHelper::setDOFFromVector(const ModelDefinition &model_definition, DegreesOfFreedom &model_dof,
                                   const VectorXF &dof_vector) {
    int dims_space = model_definition.boundary_generator->getDims();
    int n_site_free_param = model_definition.site_free_param_indices.rows();
    int dims_c = n_site_free_param + dims_space;
    int n_sites = (int)model_dof.sites.size();
    int nc = n_sites * dims_c;
    int np = (int)model_definition.boundary_free_param_indices.rows();

    for (int i = 0; i < n_sites; i++) {
        model_dof.sites[i].pos = dof_vector.segment(i * dims_c, dims_space);
        for (int j = 0; j < n_site_free_param; j++) {
            model_dof.sites[i].param(model_definition.site_free_param_indices(j)) =
                dof_vector(i * dims_c + dims_space + j);
        }
    }
    for (int i = 0; i < np; i++) {
        model_dof.boundary_param(model_definition.boundary_free_param_indices(i)) = dof_vector(nc + i);
    }

    // markRemovedSites(model_definition, model_dof);
}

bool ModelHelper::siteParamIdToSiteDOFIndex(const Model &model, int site_param_id, int &site_dof_index) {
    int dims_space = model.dimensions_ind->dims_space;
    int n_site_free_param = model.dimensions_ind->dims_site_params_free;

    for (int i = 0; i < n_site_free_param; i++) {
        if (model.model_definition.site_free_param_indices(i) == site_param_id) {
            site_dof_index = dims_space + i;
            return true;
        }
    }
    return false;
}

void ModelHelper::setFreeSiteParam(ModelDefinition &model_definition, const bool site_param_free[]) {
    std::vector<int> free_param_indices;
    for (int i = 0; i < NUM_SITE_PARAM_TOTAL; i++) {
        if (i < NUM_SITE_PARAM_TESSELLATION &&
            !model_definition.tessellation_generator->hasSiteParam((SiteParamTessellation)i))
            continue;
        if (site_param_free[i]) free_param_indices.emplace_back(i);
    }
    model_definition.site_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), free_param_indices.size());
}

VectorXI ModelHelper::getSiteParamReverseIndexVector(const Model &model) {
    VectorXI free_index = VectorXI::Constant(NUM_SITE_PARAM_TOTAL, -1);
    for (int i = 0; i < model.model_definition.site_free_param_indices.rows(); i++) {
        free_index(model.model_definition.site_free_param_indices(i)) = i;
    }
    return free_index;
}

VectorXI ModelHelper::getTessellationDOFIndexToSiteDOFIndexVector(const Model &model) {
    int dims_space = model.dimensions_ind->dims_space;
    int dims_site_dof_in_tessellation = model.dimensions_ind->dims_site_dof_in_tessellation;

    VectorXI site_param_reverse_index_vector = getSiteParamReverseIndexVector(model);

    VectorXI c_index_map = VectorXI::Constant(model.dimensions_ind->dims_site_dof_in_tessellation, -1);
    for (int i = 0; i < dims_space; i++) {
        c_index_map(i) = i;
    }
    for (int i = 0; i < NUM_SITE_PARAM_TESSELLATION; i++) {
        if (site_param_reverse_index_vector(i) == -1) continue;
        int site_param_index;
        if (model.model_definition.tessellation_generator->getSiteParamIndex((SiteParamTessellation)i,
                                                                             site_param_index)) {
            c_index_map(site_param_index) = site_param_reverse_index_vector(i) + dims_space;
        }
    }

    return c_index_map;
}

void ModelHelper::markRemovedSites(const ModelDefinition &model_definition, DegreesOfFreedom &model_dof) {
    int n_sites = (int)model_dof.sites.size();
    for (int i = 0; i < n_sites; i++) {
        /// TODO: Sometimes already-removed sites reach this point with positive size target. Maybe because
        /// the Newton step perturbs null-space DOFs, but should investigate.
        if (model_dof.sites[i].param[SITE_PARAM_SIZE_TARGET] < 0) model_dof.sites[i].is_removed = true;
    }
}
