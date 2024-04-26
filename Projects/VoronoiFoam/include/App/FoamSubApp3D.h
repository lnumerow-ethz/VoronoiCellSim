#pragma once

#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "CRLHelper/Optimization.h"

class FoamSubApp3D : public FoamSubApp {
   public:
    void getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) override;

    void initializeSubApp() override;

    void mainLoop() override;

    void makeConfigWindow() override;

    void makeAnalysisWindow() override;

    [[nodiscard]] std::string getName() const override { return "Foam3D"; }

   public:
    bool callbackMouseMove(const CRLControlState &control_state, const Vector2F &mouse_pos,
                           const Vector2F &mouse_pos_delta) override;

    bool callbackMouseScroll(const CRLControlState &control_state, float t) override;
};
