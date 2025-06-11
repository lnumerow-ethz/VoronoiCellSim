#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary2D/BoundaryGeneratorBox2D.h"

#include "CRLHelper/PeriodicHelper.h"

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
    degrees_of_freedom.boundary_param = VectorXF::Constant(num_param, 1.0);

    bool success = generateRandomSitesWithinBoundary(model_definition, degrees_of_freedom, num_sites);
    std::vector<Site> sites_tile = degrees_of_freedom.sites;
    MatrixXI tiles = PeriodicHelper::getPeriodicTiles(2, 2, true, true);
    degrees_of_freedom.sites.resize(tiles.rows() * num_sites);
    for (int i = 0; i < tiles.rows(); i++) {
        for (int j = 0; j < num_sites; j++) {
            degrees_of_freedom.sites[i * num_sites + j] = sites_tile[j];
            degrees_of_freedom.sites[i * num_sites + j].pos += Vector2F(2 * tiles(i, 0), 2 * tiles(i, 1));
        }
    }

    degrees_of_freedom.boundary_param = VectorXF::Constant(num_param, 7.0);
    return success;
}

void RandomSitesInBox2D::makeConfigMenu() {
    ImGui::InputInt("Number of Sites", &num_sites, 1, 10);
}
