#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RigidBody2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorRigidBody2D.h"

bool RigidBody2D::generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const {
    MatrixXF vertices_boundary(4, 2);
    MatrixXI edge_boundary(4, 2);
    vertices_boundary << -2, -1, 2, -1, 2, 1, -2, 1;
    edge_boundary << 0, 1, 1, 2, 2, 3, 3, 0;

    MatrixXF vertices_rigid_body(20, 2);
    MatrixXI edge_rigid_body(20, 2);
    vertices_rigid_body << 0, 0, -4, 4, -14, 10, -32, 17, -54, 21, -72, 21, -87, 35, -114, 25, -102, 14, -100, 7, -99,
        0, -100, -7, -102, -14, -114, -25, -87, -35, -72, -21, -54, -21, -32, -17, -14, -10, -4, -4;
    edge_rigid_body << 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14,
        15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 0;
    vertices_rigid_body.col(0) += VectorXF::Constant(20, 114 / 2.0);
    vertices_rigid_body = vertices_rigid_body * 1.0 / 200.0;

    // MatrixXF vertices_rigid_body(18, 2);
    // MatrixXI edge_rigid_body(18, 2);
    // vertices_rigid_body << 0, 0, -4, 4, -14, 10, -32, 17, -54, 21, -72, 21, -87, 35, -114, 25, -102, 14, -100, 7,
    // -99,
    //     0, -100, -7, -102, -14, -72, -21, -54, -21, -32, -17, -14, -10, -4, -4;
    // edge_rigid_body << 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
    // 14,
    //     15, 15, 16, 16, 17, 17, 0;
    // vertices_rigid_body.col(0) += VectorXF::Constant(18, 114 / 2.0);
    // vertices_rigid_body = vertices_rigid_body * 1.0 / 200.0;

    F rigid_body_force_magnitude = 0.035;

    model_definition.boundary_generator = std::make_shared<BoundaryGeneratorRigidBody2D>(
        vertices_boundary, edge_boundary, vertices_rigid_body, edge_rigid_body, rigid_body_force_magnitude);

    int num_param = 3;
    std::vector<int> free_param_indices;
    for (int i = 0; i < num_param; i++) {
        free_param_indices.emplace_back(i);
    }
    model_definition.boundary_free_param_indices =
        Eigen::Map<VectorXI>(free_param_indices.data(), (int)free_param_indices.size());
    degrees_of_freedom.boundary_param = Vector3F::Zero();
    degrees_of_freedom.boundary_param(0) = -1.25;

    return generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
}

void RigidBody2D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);
}
