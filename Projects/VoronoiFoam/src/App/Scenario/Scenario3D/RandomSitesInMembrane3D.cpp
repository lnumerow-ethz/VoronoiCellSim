#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wnarrowing"
#pragma GCC diagnostic ignored "-Wreorder"

#include <igl/loop.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#pragma GCC diagnostic pop

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/RandomSitesInMembrane3D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorMeshLaplacian3D.h"

bool RandomSitesInMembrane3D::generateScenario(ModelDefinition &model_definition,
                                               DegreesOfFreedom &degrees_of_freedom) const {
    MatrixXF V_ico(12, 3);
    MatrixXI F_ico(20, 3);
    F f = (1.0 + sqrt(5.0)) / 2;
    V_ico << -1, f, 0, 1, f, 0, -1, -f, 0, 1, -f, 0, 0, -1, f, 0, 1, f, 0, -1, -f, 0, 1, -f, f, 0, -1, f, 0, 1, -f, 0,
        -1, -f, 0, 1;
    F_ico << 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 11, 10, 2, 5, 11, 4, 1, 5, 9, 7, 1, 8, 10, 7, 6, 3, 9, 4,
        3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 9, 8, 1, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7;

    MatrixXF vertices_stacked;
    MatrixXI tri;
    igl::loop(V_ico, F_ico, vertices_stacked, tri, subdivision_level);
    VectorXF vertices(vertices_stacked.size());
    for (int i = 0; i < vertices_stacked.rows(); i++) {
        vertices.segment(i * 3, 3) = vertices_stacked.row(i);
    }

    model_definition.boundary_generator = std::make_shared<BoundaryGeneratorMeshLaplacian3D>(tri, spring_constant);

    std::vector<int> free_param_indices;
    for (int i = 0; i < vertices.rows(); i++) {
        free_param_indices.emplace_back(i);
    }
    model_definition.boundary_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), (int)free_param_indices.size());
    degrees_of_freedom.boundary_param = vertices;

    return generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
}

void RandomSitesInMembrane3D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);
    ImGui::InputInt("Subdivision Level", &subdivision_level, 1, 10);
    ImGui::InputDouble("Spring Constant", &spring_constant, 0.1, 0, "%.4f");
}
