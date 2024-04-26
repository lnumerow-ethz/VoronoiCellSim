#pragma once

#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class ConvergenceTest2D : public Scenario {
   public:
    void makeConfigMenu() override;

    [[nodiscard]] std::string getName() const override { return "Convergence Test"; };

   private:
    bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const override;
};
