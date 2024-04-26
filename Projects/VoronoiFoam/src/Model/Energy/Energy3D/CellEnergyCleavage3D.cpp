#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergyCleavage3D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellVolume3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellSurfaceArea3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellWeightedMean3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions3D/PerCellSecondMoment3D.h"

#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

void CellEnergyCleavage3D::getValue(const Model &model, PerCellValue &cell_value) const {
    const TessellationCell &cell = cell_value.cell;

    const Site &site = model.degrees_of_freedom.sites[cell.site_index];
    if (site.is_removed) {
        return;
    }

    int order = cell_value.order;

    PerCellValue volume = PerCellVolume3D().getValueWrapper(model, cell, order);
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

    PerCellValue site_size_target =
        PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_SIZE_TARGET).pow(1.0 / 3.0);

    cell_value = (volume / site_size_target.pow(3.0) - 1.0).square() * weights[VOLUME] +
                 ((centroid_x - site_x).square() + (centroid_y - site_y).square() + (centroid_z - site_z).square()) /
                     site_size_target.pow(2.0) * weights[CENTROID] +
                 (second_moment - volume * (centroid_x.square() + centroid_y.square() + centroid_z.square())) /
                     //  site_size_target.pow(5.0) * weights[SECOND_MOMENT] +
                     site_size_target.pow(4.0) * weights[SECOND_MOMENT] +
                 centroid_z / site_size_target.pow(1.0) * weights[GRAVITY] +
                 (site_size_target.pow(3.0) / volume).square() * ENERGY_VOLUME_BARRIER_WEIGHT;
}

void CellEnergyCleavage3D::makeConfigMenu() {
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Volume", &weights[VOLUME], ENERGY_3D_STEP_VOLUME, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Centroidal", &weights[CENTROID], ENERGY_3D_STEP_CENTROID, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Second Moment", &weights[SECOND_MOMENT], ENERGY_3D_STEP_SECOND_MOMENT, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Gravity", &weights[GRAVITY], ENERGY_3D_STEP_GRAVITY, 0, "%.4f");
}
