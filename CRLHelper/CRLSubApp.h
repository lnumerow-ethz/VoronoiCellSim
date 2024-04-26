#pragma once

#include "VecMatDef.h"
#include <Eigen/Geometry>
#include "ThirdParty/polyscope/deps/glfw/include/GLFW/glfw3.h"  // For key

enum CRLModifier {
    MOUSE_LEFT,
    MOUSE_MIDDLE,
    MOUSE_RIGHT,
    KEY_SHIFT,
    KEY_CTRL,
    KEY_ALT,
    NUM_MODIFIERS  /// The integer value of this will be the number of modifiers.
};

struct CRLControlState {
    bool modifiers[CRLModifier::NUM_MODIFIERS];

    CRLControlState() : modifiers{} { std::fill(std::begin(modifiers), std::end(modifiers), false); }
};

struct CRLViewerData {
    MatrixXF points;
    MatrixXF points_c;
    MatrixXF lines_v;
    MatrixXI lines_e;
    MatrixXF lines_c;
    MatrixXF mesh_v;
    MatrixXI mesh_f;
    MatrixXF mesh_c;
    MatrixXF mesh_tex_uv;
    MatrixXF mesh_tex_r;
    MatrixXF mesh_tex_g;
    MatrixXF mesh_tex_b;

    bool show_lines = false;
    bool points_quad = false;
    F point_size = 0.01;
    F edge_width = 0.003;

    bool write_png = false;
    std::string png_file_name;
};

struct CRLCamera {
    Vector3F eye;
    Vector3F center;

    /// Half the window height, i.e. the distance (normal to camera direction) from center to top of window.
    F height;

   private:
    /// Private to require calling get_up_direction, which ensures a vector perpendicular to camera direction.
    Vector3F up;

   public:
    [[nodiscard]] Vector3F get_camera_direction() const { return (center - eye).normalized(); }

    [[nodiscard]] Vector3F get_right_direction() const { return get_camera_direction().cross(up).normalized(); }

    [[nodiscard]] Vector3F get_up_direction() const {
        return get_right_direction().cross(get_camera_direction()).normalized();
    }

    void set_up_direction(const Vector3F &up_direction) { up = up_direction; }
};

class CRLSubApp {
   public:
    /// Get data to display in 3D graphics window.
    virtual void getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &camera) = 0;

    /// Construct ImGui config menu on left of screen.
    virtual void makeConfigWindow() {}

    /// Construct ImGui analysis window on right of screen.
    virtual void makeAnalysisWindow() {}

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackKeyPressed(const CRLControlState &control_state, int key) { return true; }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseScroll(const CRLControlState &control_state, float t) { return true; }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseDown(const CRLControlState &control_state, int button, const Vector2F &mouse_pos) {
        return true;
    }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseUp(const CRLControlState &control_state, int button, const Vector2F &mouse_pos) {
        return true;
    }

    /// Control action callbacks return true if default UI behaviour should remain.
    virtual bool callbackMouseMove(const CRLControlState &control_state, const Vector2F &mouse_pos,
                                   const Vector2F &mouse_delta) {
        return true;
    }

    /// Main loop function.
    virtual void mainLoop() {}

    /// Called when switching from another SubApp or starting up the software.
    virtual void initializeSubApp() = 0;

    /// Get name of SubApp for selection in UI.
    [[nodiscard]] virtual std::string getName() const = 0;
};
