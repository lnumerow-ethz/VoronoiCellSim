#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class RandomSitesInBox2D : public Scenario {
   public:
    int num_sites = SCENARIO_SITES_IN_BOX_2D_NUM_SITES;

    bool dimensions_independent = false;
    bool dimensions_free[2] = {false, false};
    F dimension_defaults[2] = {1, 1};

   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Random Sites in Box"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
