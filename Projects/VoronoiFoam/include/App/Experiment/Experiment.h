#pragma once

#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "CRLHelper/Optimization.h"

class FoamSubApp;

class Experiment {
   public:
    [[nodiscard]] virtual std::string getName() const = 0;
    virtual void setup(FoamSubApp* foam_sub_app) = 0;
    virtual void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) = 0;
};

/// Blank experiment class which does nothing, to be selected by default.
class NoExperimentSelected : public Experiment {
   public:
    [[nodiscard]] std::string getName() const override { return "--Select--"; };
    void setup(FoamSubApp* foam_sub_app) override{};
    void loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) override{};
};
