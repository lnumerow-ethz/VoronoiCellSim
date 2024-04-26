#include "CRLHelper/CRLApp.h"

#include <exception>
#include <iostream>

#include "Projects/VoronoiFoam/include/App/FoamSubApp2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp3D.h"

void customTerminate() {
    std::cerr << "Custom terminate handler" << std::endl;
    // Manually trigger a breakpoint, or log stack trace
    // For GDB or LLDB: __builtin_trap();
    abort();  // To ensure termination
}

int main() {
    std::set_terminate(customTerminate);

    CRLApp app;
    app.subapps_selector.second.push_back(std::make_shared<FoamSubApp2D>());
    app.subapps_selector.second.push_back(std::make_shared<FoamSubApp3D>());
    app.launch();

    // FoamSubApp3D sub_app;
    // sub_app.initializeSubApp();
    // sub_app.runExperiment(2);
    // while (true) {
    //     sub_app.mainLoop();
    // }

    // FoamSubApp2D sub_app;
    // sub_app.initializeSubApp();
    // sub_app.runExperiment(1);
    // while (true) {
    //     sub_app.mainLoop();
    // }

    return 0;
}
