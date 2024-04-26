#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/FoamSubApp3D.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/Voronoi3D/Voronoi3D.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation3D/Power3D/Power3D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentCoarsening3D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentCleavage3D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentCleavageCylinder3D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam3D/ExperimentRuntimeScaling3D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/RandomSitesInBox3D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/RandomSitesInMembrane3D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario3D/OpenCylinder3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergy3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergyCleavage3D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergyCoarsening3D.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

#include "CRLHelper/CameraHelper.h"

#include <iostream>

void FoamSubApp3D::mainLoop() {
    FoamSubApp::mainLoop();
}

void FoamSubApp3D::makeConfigWindow() {
    FoamSubApp::makeConfigWindow();
}

void FoamSubApp3D::makeAnalysisWindow() {
    FoamSubApp::makeAnalysisWindow();
}

void FoamSubApp3D::getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) {
    bool model_generation_success;
    Model model(model_definition, degrees_of_freedom, 0, model_generation_success);
    if (!model_generation_success) {
        std::cout << "Broken model encountered in FoamSubApp3D::getViewerData." << std::endl;
        assert(0 && "Broken model encountered in FoamSubApp3D::getViewerData.");
    }

    int dims_space = model.dimensions_ind->dims_space;
    int n_nodes = model.dimensions_tess->n_nodes;
    int n_sites = model.dimensions_ind->n_sites;

    /// Initialize a single CRLViewerData struct.
    viewer_data.clear();
    viewer_data.emplace_back();
    CRLViewerData &viewer_data_tessellation = viewer_data[0];

    /// Set to generate image file if requested.
    if (write_png) {
        viewer_data_tessellation.write_png = true;
        viewer_data_tessellation.png_file_name = png_file_name;
        write_png = false;
    }

    /// Tessellation nodes are the vertices for both 'lines' and 'mesh'.
    viewer_data_tessellation.lines_v = MatrixXF::Zero(n_nodes, 3);
    for (int i = 0; i < n_nodes; i++) {
        /// The head(dims_space) allows this to work for 2D and 3D, leaving the z-coordinate as zero for 2D.
        viewer_data_tessellation.lines_v.row(i).head(dims_space) = model.nodes[i].data->pos;
    }
    viewer_data_tessellation.mesh_v = viewer_data_tessellation.lines_v;

    /// Count the numbers of lines and mesh triangles to draw.
    int n_lines = 0;
    int n_mesh_triangles = 0;
    for (int i = 0; i < n_sites; i++) {
        const TessellationCell &cell = model.cells[i];
        for (const TessellationFace &face : cell.faces) {
            if (face.is_boundary_face) {
                n_mesh_triangles += (int)face.node_indices.size() - 2;
            } else {
                n_lines += (int)face.node_indices.size();
            }
        }
    }

    /// Resize viewer data matrices accordingly.
    viewer_data_tessellation.lines_e.resize(n_lines, 2);
    viewer_data_tessellation.lines_c = MatrixXF::Zero(n_lines, 3);
    viewer_data_tessellation.mesh_f.resize(n_mesh_triangles, 3);
    viewer_data_tessellation.mesh_c.resize(n_mesh_triangles, 3);
    int line_index = 0;
    int triangle_index = 0;

    /// Generate random (but consistent between program instances) colors for the cells.
    srand(0);  // NOLINT (suppress Clang-Tidy warning about pseudo-randomness)
    MatrixXF colors(n_sites, 3);
    for (int i = 0; i < n_sites; i++) {
        colors.row(i) = 0.5 * (Vector3F::Random() + Vector3F::Ones());
    }

    /// Copy cell edges to 'lines' and faces (triangulated) to 'mesh'.
    for (const TessellationCell &cell : model.cells) {
        for (const TessellationFace &face : cell.faces) {
            if (face.is_boundary_face) {
                for (size_t i = 1; i < face.node_indices.size() - 1; i++) {
                    viewer_data_tessellation.mesh_f(triangle_index, 0) = face.node_indices[0];
                    viewer_data_tessellation.mesh_f(triangle_index, 1) = face.node_indices[i];
                    viewer_data_tessellation.mesh_f(triangle_index, 2) = face.node_indices[i + 1];
                    viewer_data_tessellation.mesh_c.row(triangle_index) = colors.row(cell.site_index);
                    triangle_index++;
                }
            } else {
                for (size_t i = 0; i < face.node_indices.size(); i++) {
                    viewer_data_tessellation.lines_e(line_index, 0) = face.node_indices[i];
                    viewer_data_tessellation.lines_e(line_index, 1) =
                        face.node_indices[(i + 1) % face.node_indices.size()];
                    line_index++;
                }
            }
        }
    }

    /// Copy camera state.
    viewer_camera = app_camera;
}

void FoamSubApp3D::initializeSubApp() {
    experiment_selector.second.clear();
    experiment_selector.second.emplace_back(std::make_shared<NoExperimentSelected>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentCoarsening3D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentCleavage3D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentCleavageCylinder3D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentRuntimeScaling3D>());

    scenario_selector.second.clear();
    scenario_selector.second.push_back(std::make_shared<RandomSitesInBox3D>());
    scenario_selector.second.push_back(std::make_shared<RandomSitesInMembrane3D>());
    scenario_selector.second.push_back(std::make_shared<OpenCylinder3D>());
    scenario_selector.second[scenario_selector.first]->assignScenario(model_definition, degrees_of_freedom);

    tessellation_selector.second.clear();
    tessellation_selector.second.push_back(std::make_shared<Voronoi3D>());
    tessellation_selector.second.push_back(std::make_shared<Power3D>());

    energy_selector.second.clear();
    energy_selector.second.emplace_back("Standard", std::make_shared<CellEnergy3D>());
    energy_selector.second.emplace_back("Cleavage", std::make_shared<CellEnergyCleavage3D>());
    energy_selector.second.emplace_back("Coarsening", std::make_shared<CellEnergyCoarsening3D>());

    model_definition.tessellation_generator = tessellation_selector.second[tessellation_selector.first];
    model_definition.cell_energy_function = energy_selector.second[energy_selector.first].second;

    /// Initialize camera.
    app_camera.eye = Vector3F(5, 0, 0);
    app_camera.center = Vector3F(0, 0, 0);
    app_camera.set_up_direction(Vector3F(0, 0, 1));
    app_camera.height = 1;

    /// Initialization from superclass.
    FoamSubApp::initializeSubApp();
}

bool FoamSubApp3D::callbackMouseMove(const CRLControlState &control_state, const Vector2F &mouse_pos,
                                     const Vector2F &mouse_delta) {
    if (control_state.modifiers[MOUSE_LEFT]) {
        if (control_state.modifiers[KEY_SHIFT]) {
            CameraHelper::panCameraFromDrag(app_camera, mouse_delta, CAMERA_PAN_SENSITIVITY);
        } else {
            CameraHelper::rotateCameraFromDrag(app_camera, mouse_delta, CAMERA_ROTATE_SENSITIVITY);
        }
    }

    return false;
}

bool FoamSubApp3D::callbackMouseScroll(const CRLControlState &control_state, float t) {
    Vector3F camera_dist = app_camera.eye - app_camera.center;
    app_camera.eye = app_camera.center + camera_dist * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);
    app_camera.height = app_camera.height * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);

    return false;
}
