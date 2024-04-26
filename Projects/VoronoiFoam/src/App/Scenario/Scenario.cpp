#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

#include <iostream>

bool Scenario::assignScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) {
    degrees_of_freedom.sites.clear();
    return generateScenario(model_definition, degrees_of_freedom);
}

bool Scenario::generateRandomSitesWithinBoundary(const ModelDefinition &model_definition,
                                                 DegreesOfFreedom &degrees_of_freedom, int num_sites) {
    BoundaryData boundary_data;
    bool valid_boundary = model_definition.boundary_generator->generateBoundary(degrees_of_freedom, boundary_data, 0);

    if (!valid_boundary) {
        return false;
    }

    int dims_space = model_definition.boundary_generator->getDims();
    VectorXF bbox_lower, bbox_upper;
    BoundaryHelper::getBoundaryBoundingBox(boundary_data, bbox_lower, bbox_upper);

    F cell_size_target = model_definition.boundary_generator->computeEnclosedVolume(boundary_data) / num_sites;

    srand(0);
    for (int i = 0; i < num_sites; i++) {
        degrees_of_freedom.sites.emplace_back();
        do {
            degrees_of_freedom.sites.back().pos = ((VectorXF::Random(dims_space) + VectorXF::Ones(dims_space)) / 2.0)
                                                      .cwiseProduct(bbox_upper - bbox_lower) +
                                                  bbox_lower;
            degrees_of_freedom.sites.back().param(SITE_PARAM_SIZE_TARGET) = cell_size_target;
        } while (!model_definition.boundary_generator->checkPointInBounds(boundary_data,
                                                                          degrees_of_freedom.sites.back().pos));
    }

    return true;
}
