#pragma once

#include "Projects/VoronoiFoam/include/Model/Model.h"

class Scenario {
   public:
    /// Wrapper for generateScenario to perform generic setup.
    bool assignScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom);

    virtual void makeConfigMenu() = 0;

    [[nodiscard]] virtual std::string getName() const = 0;

   private:
    virtual bool generateScenario(ModelDefinition &model_definition, DegreesOfFreedom &degrees_of_freedom) const = 0;

   protected:
    /// Returns true if valid boundary generated. Assigns equal cell size targets which sum to the total domain volume.
    static bool generateRandomSitesWithinBoundary(const ModelDefinition &model_definition,
                                                  DegreesOfFreedom &degrees_of_freedom, int num_sites);
};
