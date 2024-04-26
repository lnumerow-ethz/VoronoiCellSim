#pragma once

#include <memory>
#include <vector>

#include "VecMatDef.h"

class CRLSubApp;

class CRLApp {
   public:
    /// SubApp list with index of selected subapp.
    std::pair<int, std::vector<std::shared_ptr<CRLSubApp>>> subapps_selector = {0, {}};

   public:
    void launch();

    std::shared_ptr<CRLSubApp> subApp() { return subapps_selector.second[subapps_selector.first]; }
};
