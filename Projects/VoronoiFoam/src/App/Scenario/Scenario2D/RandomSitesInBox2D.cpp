#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorBox2D.h"

bool RandomSitesInBox2D::generateScenario(ModelDefinition &model_definition,
                                          DegreesOfFreedom &degrees_of_freedom) const {
    model_definition.boundary_generator = std::make_shared<BoundaryGeneratorBox2D>();

    int dims_space = model_definition.boundary_generator->getDims();
    int num_param = (dimensions_independent ? dims_space : 1);

    std::vector<int> free_param_indices;
    for (int i = 0; i < num_param; i++) {
        if (dimensions_free[i]) free_param_indices.emplace_back(i);
    }
    model_definition.boundary_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), (int)free_param_indices.size());
    degrees_of_freedom.boundary_param = Eigen::Map<const VectorXF>(dimension_defaults, num_param);

    return generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
}

void RandomSitesInBox2D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);

    ImGui::Checkbox("Independent Dimensions", &dimensions_independent);

    int dims_space = 2;
    int num_param = (dimensions_independent ? dims_space : 1);

    ImGui::Checkbox("##FreeDimension0", &dimensions_free[0]);
    for (int i = 1; i < num_param; i++) {
        if (i > 0) ImGui::SameLine();
        ImGui::Checkbox(("##FreeDimension" + std::to_string(i)).c_str(), &dimensions_free[i]);
    }
    ImGui::SameLine();
    ImGui::Text("Free Dimensions");
}
