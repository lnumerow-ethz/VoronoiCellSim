#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class RigidBody2D : public Scenario {
   public:
    int num_sites = SCENARIO_SITES_IN_BOX_2D_NUM_SITES;

   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Rigid Body"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
