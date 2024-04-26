#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class OpenCylinder3D : public Scenario {
   public:
    int num_sites = SCENARIO_SITES_IN_MEMBRANE_3D_NUM_SITES;
    F spring_constant = SCENARIO_SITES_IN_MEMBRANE_3D_SPRING_CONSTANT;

   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Open-Top Cylinder"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
