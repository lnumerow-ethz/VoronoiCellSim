#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Power2D/Power2D.h"

bool Power2D::computeCells(Model &model) const {
    return Voronoi2D::computeCellsClippedVoronoi2D(model, true);
}

bool Power2D::getSiteParamIndex(SiteParamTessellation param_id, int &index) const {
    switch (param_id) {
        case SITE_PARAM_POWER_WEIGHT:
            index = 2;
            return true;
        default:
            return false;
    }
}
