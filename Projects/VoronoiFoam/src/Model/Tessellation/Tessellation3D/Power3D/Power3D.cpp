#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/Power3D/Power3D.h"

bool Power3D::computeCells(Model &model) const {
    return Voronoi3D::computeCellsClippedVoronoi3D(model, true);
}

bool Power3D::getSiteParamIndex(SiteParamTessellation param_id, int &index) const {
    switch (param_id) {
        case SITE_PARAM_POWER_WEIGHT:
            index = 3;
            return true;
        default:
            return false;
    }
}
