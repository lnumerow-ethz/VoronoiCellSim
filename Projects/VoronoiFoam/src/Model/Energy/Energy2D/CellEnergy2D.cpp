#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellArea2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellPerimeter2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellWeightedMean2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/SimplexFunctions/SimplexFunctions2D/PerCellSecondMoment2D.h"

#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"

void CellEnergy2D::getValue(const Model &model, PerCellValue &cell_value) const {
    const TessellationCell &cell = cell_value.cell;

    const Site &site = model.degrees_of_freedom.sites[cell.site_index];
    if (site.is_removed) {
        return;
    }

    int order = cell_value.order;

    PerCellValue area = PerCellArea2D().getValueWrapper(model, cell, order);
    PerCellValue perimeter = PerCellPerimeter2D().getValueWrapper(model, cell, order);
    PerCellValue weighted_mean_x = PerCellWeightedMeanX2D().getValueWrapper(model, cell, order);
    PerCellValue weighted_mean_y = PerCellWeightedMeanY2D().getValueWrapper(model, cell, order);
    PerCellValue second_moment = PerCellSecondMoment2D().getValueWrapper(model, cell, order);

    PerCellValue centroid_x = weighted_mean_x / area;
    PerCellValue centroid_y = weighted_mean_y / area;

    PerCellValue site_x = PerCellValue::sitePosCellValue(model, cell, order, 0);
    PerCellValue site_y = PerCellValue::sitePosCellValue(model, cell, order, 1);

    PerCellValue site_power_weight = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_POWER_WEIGHT);
    PerCellValue site_size_target = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_SIZE_TARGET);
    PerCellValue site_surface_target = PerCellValue::siteParamCellValue(model, cell, order, SITE_PARAM_SURFACE_TARGET);

    cell_value = (area - site_size_target).square() * weights[AREA] +
                 (perimeter - site_surface_target).square() * weights[PERIMETER_TARGET] +
                 perimeter * weights[PERIMETER_MINIMIZATION] +
                 ((centroid_x - site_x).square() + (centroid_y - site_y).square()) * weights[CENTROID] +
                 (second_moment - area * (centroid_x.square() + centroid_y.square())) * weights[SECOND_MOMENT] +
                 (site_size_target / area).square() * ENERGY_VOLUME_BARRIER_WEIGHT +
                 site_power_weight.square() * ENERGY_POWER_REGULARIZER_WEIGHT;
}

void CellEnergy2D::makeConfigMenu() {
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Area", &weights[AREA], ENERGY_2D_STEP_AREA, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Perimeter Target", &weights[PERIMETER_TARGET], ENERGY_2D_STEP_PERIMETER_TARGET, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Perimeter Minimization", &weights[PERIMETER_MINIMIZATION],
                       ENERGY_2D_STEP_PERIMETER_MINIMIZATION, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Centroidal", &weights[CENTROID], ENERGY_2D_STEP_CENTROID, 0, "%.4f");
    ImGui::SetNextItemWidth(ImGui::GetWindowWidth() * 0.5f);
    ImGui::InputDouble("Second Moment", &weights[SECOND_MOMENT], ENERGY_2D_STEP_SECOND_MOMENT, 0, "%.4f");
}
