#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergy3D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellVolume3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellSurfaceArea3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellWeightedMean3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellSecondMoment3D.h"

#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

void CellEnergy3D::getValue(const Model &model, PerCellValue &cell_value) const {
    const TessellationCell &cell = cell_value.cell;

    const Site &site = model.degrees_of_freedom.sites[cell.site_index];
    if (site.is_removed) {
        return;
    }

    int order = cell_value.order;

    PerCellValue volume = PerCellVolume3D().getValueWrapper(model, cell, order);
    PerCellValue surface_area = PerCellSurfaceArea3D().getValueWrapper(model, cell, order);
    PerCellValue weighted_mean_x = PerCellWeightedMeanX3D().getValueWrapper(model, cell, order);
    PerCellValue weighted_mean_y = PerCellWeightedMeanY3D().getValueWrapper(model, cell, order);
    PerCellValue weighted_mean_z = PerCellWeightedMeanZ3D().getValueWrapper(model, cell, order);
    PerCellValue second_moment = PerCellSecondMoment3D().getValueWrapper(model, cell, order);

    PerCellValue centroid_x = weighted_mean_x / volume;
    PerCellValue centroid_y = weighted_mean_y / volume;
    PerCellValue centroid_z = weighted_mean_z / volume;

    PerCellValue site_x = PerCellValue::sitePosCellValue(model, cell, order, 0);
    PerCellValue site_y = PerCellValue::sitePosCellValue(model, cell, order, 1);
    PerCellValue site_z = PerCellValue::sitePosCellValue(model, cell, order, 2);

    PerCellValue site_power_weight = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_POWER_WEIGHT);
    PerCellValue site_size_target = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_SIZE_TARGET);
    PerCellValue site_surface_target = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_SURFACE_TARGET);

    cell_value = (volume - site_size_target).square() * weights[VOLUME] +
                 (surface_area - site_surface_target).square() * weights[SURFACE_TARGET] +
                 surface_area * weights[SURFACE_MINIMIZATION] +
                 ((centroid_x - site_x).square() + (centroid_y - site_y).square() + (centroid_z - site_z).square()) *
                     weights[CENTROID] +
                 (second_moment - volume * (centroid_x.square() + centroid_y.square() + centroid_z.square())) *
                     weights[SECOND_MOMENT] +
                 (site_size_target / volume).square() * ENERGY_VOLUME_BARRIER_WEIGHT +
                 site_power_weight.square() * ENERGY_POWER_REGULARIZER_WEIGHT;
}

void CellEnergy3D::makeConfigMenu() {
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Volume", &weights[VOLUME], ENERGY_3D_STEP_VOLUME, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Surface Target", &weights[SURFACE_TARGET], ENERGY_3D_STEP_SURFACE_TARGET, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Surface Minimization", &weights[SURFACE_MINIMIZATION], ENERGY_3D_STEP_SURFACE_MINIMIZATION, 0,
                       "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Centroidal", &weights[CENTROID], ENERGY_3D_STEP_CENTROID, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Second Moment", &weights[SECOND_MOMENT], ENERGY_3D_STEP_SECOND_MOMENT, 0, "%.4f");
}
