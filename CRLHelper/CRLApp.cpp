#include "CRLApp.h"
#include "CRLSubApp.h"
#include "ImGuiHelpers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "polyscope/camera_view.h"

#include "ThirdParty/polyscope/deps/glfw/include/GLFW/glfw3.h"
#include "ThirdParty/polyscope/deps/imgui/imgui/imgui.h"

static CRLControlState getControlState() {
    CRLControlState control_state;

    GLFWwindow *window = glfwGetCurrentContext();

    control_state.modifiers[MOUSE_LEFT] = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    control_state.modifiers[MOUSE_RIGHT] = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    control_state.modifiers[MOUSE_MIDDLE] = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
    control_state.modifiers[KEY_ALT] =
        (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS);
    control_state.modifiers[KEY_SHIFT] = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                                          glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);
    control_state.modifiers[KEY_CTRL] = (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
                                         glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS);

    return control_state;
}

static void setCamera(const CRLCamera &camera) {
    Vector3F temp = camera.eye;
    glm::vec3 eye(temp(0), temp(1), temp(2));
    temp = camera.center;
    glm::vec3 center(temp(0), temp(1), temp(2));
    temp = camera.get_up_direction();
    glm::vec3 up(temp(0), temp(1), temp(2));

    /// Calculate the vertical FOV using the given window height
    F eye_center_distance = (camera.center - camera.eye).norm();
    F fovY_rad = 2.0 * std::atan(camera.height / eye_center_distance);
    F fovY_deg = fovY_rad * (180.0 / M_PI);

    polyscope::CameraParameters camera_params(polyscope::CameraIntrinsics::fromFoVDegVerticalAndAspect(fovY_deg, 1.0),
                                              polyscope::CameraExtrinsics::fromVectors(eye, center - eye, up));
    polyscope::view::setViewToCamera(camera_params);
}

static void initializePolyscope() {
    polyscope::options::autocenterStructures = false;
    polyscope::options::autoscaleStructures = false;
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = 1.;
    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{{-1., -1., -1.}, {1., 1., 1.}};
    polyscope::view::windowWidth = 1920;
    polyscope::view::windowHeight = 1080;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::groundPlaneHeightFactor = 0.;
    polyscope::options::shadowDarkness = 0.4;
    polyscope::options::buildGui = false;
    polyscope::options::ssaaFactor = 2;
    polyscope::init();
}

static void updateViewerData(const std::shared_ptr<CRLSubApp> &sub_app) {
    std::vector<CRLViewerData> viewer_data;
    CRLCamera camera;

    sub_app->getViewerData(viewer_data, camera);
    setCamera(camera);

    polyscope::removeAllStructures();
    for (int i = 0; i < (int)viewer_data.size(); i++) {
        CRLViewerData &curr_data = viewer_data[i];
        if (curr_data.mesh_f.rows() > 0) {
            auto mesh = polyscope::registerSurfaceMesh("Mesh" + std::to_string(i), curr_data.mesh_v, curr_data.mesh_f);
            auto color = mesh->addFaceColorQuantity("Color", curr_data.mesh_c);
            color->setEnabled(true);

            if (curr_data.show_lines) {
                mesh->setEdgeWidth(1.0);
            }
        }

        if (curr_data.points.rows() > 0) {
            auto points = polyscope::registerPointCloud("Points" + std::to_string(i), curr_data.points);
            auto color = points->addColorQuantity("Color", curr_data.points_c);
            color->setEnabled(true);

            points->setPointRadius(curr_data.point_size);
            if (curr_data.points_quad) {
                points->setPointRenderMode(polyscope::PointRenderMode::Quad);
            }
        }

        if (curr_data.lines_e.rows() > 0) {
            MatrixXF lines_v = curr_data.lines_v;
            MatrixXI lines_e = curr_data.lines_e;
            std::vector<bool> vertex_visited(lines_v.rows(), false);
            for (int ie = 0; ie < lines_e.rows(); ie++) {
                vertex_visited[lines_e(ie, 0)] = true;
                vertex_visited[lines_e(ie, 1)] = true;
            }
            for (int iv = 0; iv < lines_v.rows(); iv++) {
                if (!vertex_visited[iv]) {
                    lines_v.row(iv) = Vector3F::Constant(1e10);
                }
            }

            auto lines = polyscope::registerCurveNetwork("Lines" + std::to_string(i), lines_v, lines_e);
            auto color = lines->addEdgeColorQuantity("Color", curr_data.lines_c);
            color->setEnabled(true);

            lines->setRadius(curr_data.edge_width);
        }

        // TODO: Handle texture map data.
        // if (curr_data.mesh_f.rows() > 0) {
        //     viewer.data(i).set_mesh(curr_data.mesh_v, curr_data.mesh_f);
        //     if (curr_data.mesh_tex_uv.rows() == 0) {
        //         viewer.data(i).set_colors(curr_data.mesh_c);
        //     } else {
        //         viewer.data(i).show_texture = true;
        //         viewer.data(i).set_colors(Eigen::RowVector3d(1, 1, 1));
        //         viewer.data(i).set_texture(curr_data.mesh_tex_r.cast<unsigned char>(),
        //                                    curr_data.mesh_tex_g.cast<unsigned char>(),
        //                                    curr_data.mesh_tex_b.cast<unsigned char>());
        //         viewer.data(i).set_uv(curr_data.mesh_tex_uv);
        //     }
        // }

        if (curr_data.write_png) {
            polyscope::screenshot(curr_data.png_file_name);
        }
    }
}

void CRLApp::launch() {
    /// Initialize polyscope
    initializePolyscope();

    // /// Startup code for ImPlot.
    // if (!ImPlot::GetCurrentContext()) {
    //     auto ctx = ImPlot::CreateContext();
    //     ImPlot::SetCurrentContext(ctx);
    // }

    polyscope::state::userCallback = [&]() {
        /// Main menu window on screen left.
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x * 0.25, ImGui::GetIO().DisplaySize.y),
                                 ImGuiCond_Once);
        ImGui::Begin("Menu");

        /// Combo box to select SubApp.
        std::vector<std::string> subapp_names;
        for (const auto &subapp : subapps_selector.second) {
            subapp_names.push_back(subapp->getName());
        }
        if (ImGui::Combo("App", &subapps_selector.first, subapp_names)) {
            subApp()->initializeSubApp();
        }
        ImGui::Spacing();
        ImGui::Spacing();

        /// Remaining controls provided by SubApp.
        subApp()->makeConfigWindow();
        ImGui::End();

        /// Additional "Analysis" window on right, content provided by SubApp.
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.75, 0), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x * 0.25, ImGui::GetIO().DisplaySize.y),
                                 ImGuiCond_Once);
        ImGui::Begin("Analysis");
        subApp()->makeAnalysisWindow();
        ImGui::End();

        /// Handle mouse and keyboard controls.
        CRLControlState control_state = getControlState();
        ImVec2 im_mouse_pos = ImGui::GetIO().MousePos;
        ImVec2 im_mouse_delta = ImGui::GetIO().MouseDelta;
        Vector2F mouse_pos(im_mouse_pos.x, im_mouse_pos.y);
        Vector2F mouse_delta(im_mouse_delta.x, im_mouse_delta.y);
        F mouse_scroll = ImGui::GetIO().MouseWheel;

        if (!ImGui::GetIO().WantCaptureKeyboard) {
            for (int i = 0; i <= GLFW_KEY_LAST; i++) {
                if (ImGui::IsKeyPressed(i, false)) {
                    subApp()->callbackKeyPressed(control_state, i);
                }
            }
        }

        if (!ImGui::GetIO().WantCaptureMouse) {
            for (int i = 0; i < IM_ARRAYSIZE(ImGui::GetIO().MouseDown); i++) {
                if (ImGui::IsMouseClicked(i, false)) {
                    subApp()->callbackMouseDown(control_state, i, mouse_pos);
                }
                if (ImGui::IsMouseReleased(i)) {
                    subApp()->callbackMouseUp(control_state, i, mouse_pos);
                }
            }
            if (mouse_delta != Vector2F::Zero()) {
                subApp()->callbackMouseMove(control_state, mouse_pos, mouse_delta);
            }
            if (mouse_scroll != 0.0) {
                subApp()->callbackMouseScroll(control_state, mouse_scroll);
            }
        }

        /// Update display.
        updateViewerData(subApp());

        /// Main loop function provided by SubApp.
        subApp()->mainLoop();
    };

    /// Initialize default SubApp.
    subApp()->initializeSubApp();

    // /// Launch viewer.
    polyscope::show();
}
