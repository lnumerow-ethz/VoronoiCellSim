#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wnarrowing"
#pragma GCC diagnostic ignored "-Wreorder"

#include <igl/loop.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#pragma GCC diagnostic pop

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/OpenCylinder3D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorOpenCylinder3D.h"

bool OpenCylinder3D::generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const {
    std::shared_ptr<BoundaryGeneratorOpenCylinder3D> open_cylinder_boundary =
        std::make_shared<BoundaryGeneratorOpenCylinder3D>(spring_constant);
    open_cylinder_boundary->getDefaultParams(degrees_of_freedom.boundary_param);

    int num_param = degrees_of_freedom.boundary_param.rows();

    std::vector<int> free_param_indices;
    for (int i = 0; i < num_param; i++) {
        free_param_indices.emplace_back(i);
    }
    model_definition.boundary_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), (int)free_param_indices.size());

    model_definition.boundary_generator = open_cylinder_boundary;

    return generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
}

void OpenCylinder3D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);
    ImGui::InputDouble("Spring Constant", &spring_constant, 0.1, 0, "%.4f");
}
