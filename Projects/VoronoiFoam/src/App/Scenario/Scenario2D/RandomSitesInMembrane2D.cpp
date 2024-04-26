#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInMembrane2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorPolylineLaplacian2D.h"

bool RandomSitesInMembrane2D::generateScenario(ModelDefinition &model_definition,
                                               DegreesOfFreedom &degrees_of_freedom) const {
    MatrixXI edge(num_boundary_vertices, 2);
    VectorXF vertices(num_boundary_vertices * 2);
    for (int i = 0; i < num_boundary_vertices; i++) {
        edge(i, 0) = i;
        edge(i, 1) = (i + 1) % num_boundary_vertices;

        F angle = i * 2 * M_PI / num_boundary_vertices;
        vertices(2 * i + 0) = cos(angle);
        vertices(2 * i + 1) = sin(angle);
    }

    model_definition.boundary_generator = std::make_shared<BoundaryGeneratorPolylineLaplacian2D>(edge, spring_constant);

    std::vector<int> free_param_indices;
    for (int i = 0; i < vertices.rows(); i++) {
        free_param_indices.emplace_back(i);
    }
    model_definition.boundary_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), (int)free_param_indices.size());
    degrees_of_freedom.boundary_param = vertices;

    return generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
}

void RandomSitesInMembrane2D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);
    ImGui::InputInt("Boundary Vertices", &num_boundary_vertices, 1, 10);
    ImGui::InputDouble("Spring Constant", &spring_constant, 0.1, 0, "%.4f");
}
