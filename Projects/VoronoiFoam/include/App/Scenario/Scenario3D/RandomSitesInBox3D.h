#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class RandomSitesInBox3D : public Scenario {
   public:
    int num_sites = SCENARIO_SITES_IN_BOX_3D_NUM_SITES;

    bool dimensions_independent = false;
    bool dimensions_free[3] = {false, false, false};
    F dimension_defaults[3] = {1, 1, 1};

   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Random Sites in Box"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
