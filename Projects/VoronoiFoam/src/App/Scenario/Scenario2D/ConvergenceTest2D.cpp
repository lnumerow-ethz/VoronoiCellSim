#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/ConvergenceTest2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorBox2D.h"

bool ConvergenceTest2D::generateScenario(ModelDefinition &model_definition,
                                         DegreesOfFreedom &degrees_of_freedom) const {
    model_definition.boundary_generator = std::make_shared<BoundaryGeneratorBox2D>();

    model_definition.boundary_free_param_indices.resize(0);
    degrees_of_freedom.boundary_param = Vector2F(1, 0.3);

    degrees_of_freedom.sites.resize(4);
    degrees_of_freedom.sites[0].pos = Vector2F(0.1, 0);
    degrees_of_freedom.sites[1].pos = Vector2F(0, 0.1001);
    degrees_of_freedom.sites[2].pos = Vector2F(-0.1, 0);
    degrees_of_freedom.sites[3].pos = Vector2F(0, -0.1001);
    for (int i = 0; i < 4; i++) {
        degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET) = 0.3;
    }

    return true;
}

void ConvergenceTest2D::makeConfigMenu() {}
