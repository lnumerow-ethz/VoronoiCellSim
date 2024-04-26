#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class RandomSitesInMembrane2D : public Scenario {
   public:
    int num_sites = SCENARIO_SITES_IN_MEMBRANE_2D_NUM_SITES;
    int num_boundary_vertices = SCENARIO_SITES_IN_MEMBRANE_2D_NUM_VERTICES;
    F spring_constant = SCENARIO_SITES_IN_MEMBRANE_2D_SPRING_CONSTANT;

   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Random Sites in Membrane"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
