#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/VoronoiFoam/include/App/FoamSubApp2D.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Voronoi2D/Voronoi2D.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/Tessellation2D/Power2D/Power2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentCoarsening2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentConvergenceTest2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentRigidBody2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentDCMComparison2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentImageMatch2D.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentRuntimeScaling2D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInMembrane2D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/ConvergenceTest2D.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RigidBody2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergyCleavage2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergyCoarsening2D.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"
#include "Projects/VoronoiFoam/include/App/Rendering/Rendering.h"

#include "CRLHelper/CameraHelper.h"

#include <iostream>

void FoamSubApp2D::mainLoop() {
    FoamSubApp::mainLoop();
}

void FoamSubApp2D::makeConfigWindow() {
    FoamSubApp::makeConfigWindow();
}

void FoamSubApp2D::makeAnalysisWindow() {
    FoamSubApp::makeAnalysisWindow();
}

void FoamSubApp2D::getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) {
    bool model_generation_success;
    Model model(model_definition, degrees_of_freedom, 0, model_generation_success);
    if (!model_generation_success) {
        std::cout << "Broken model encountered in FoamSubApp2D::getViewerData." << std::endl;
        assert(0 && "Broken model encountered in FoamSubApp2D::getViewerData.");
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

    /// Display sites as points.
    viewer_data_tessellation.points = MatrixXF::Zero(n_sites, 3);
    for (int i = 0; i < n_sites; i++) {
        viewer_data_tessellation.points.row(i).head(dims_space) = model.degrees_of_freedom.sites[i].pos;
        /// Move removed sites off screen.
        if (model.degrees_of_freedom.sites[i].is_removed)
            viewer_data_tessellation.points.row(i).head(dims_space).setConstant(1e10);
    }
    viewer_data_tessellation.points_c = MatrixXF::Zero(n_sites, 3);

    /// Tessellation nodes are the vertices for 'lines'.
    viewer_data_tessellation.lines_v = MatrixXF::Zero(n_nodes, 3);
    for (int i = 0; i < n_nodes; i++) {
        /// The head(dims_space) allows this to work for 2D and 3D, leaving the z-coordinate as zero for 2D.
        viewer_data_tessellation.lines_v.row(i).head(dims_space) = model.nodes[i].data->pos;
    }

    /// Count the numbers of lines to draw.
    int n_lines = 0;
    for (int i = 0; i < n_sites; i++) {
        const TessellationCell &cell = model.cells[i];
        n_lines += (int)cell.faces.size();
    }

    /// Resize viewer data matrices accordingly.
    viewer_data_tessellation.lines_e.resize(n_lines, 2);
    viewer_data_tessellation.lines_c = MatrixXF::Zero(n_lines, 3);
    int line_index = 0;

    /// Copy cell edges to 'lines'.
    for (const TessellationCell &cell : model.cells) {
        for (const TessellationFace &face : cell.faces) {
            viewer_data_tessellation.lines_e(line_index, 0) = face.node_indices[0];
            viewer_data_tessellation.lines_e(line_index, 1) = face.node_indices[1];
            line_index++;
        }
    }

    /// Generate random (but consistent between program instances) colors for the cells.
    srand(0);  // NOLINT (suppress Clang-Tidy warning about pseudo-randomness)
    MatrixXF colors(n_sites, 3);
    for (int i = 0; i < n_sites; i++) {
        colors.row(i) = 0.5 * (Vector3F::Random() + Vector3F::Ones());
    }

    /// Obtain a triangulation of the cells in 'mesh' and assign color to each triangle according to cell index.
    VectorXI triangulation_cell_indices;
    Rendering::triangulateVoronoiCells2D(model, viewer_data_tessellation.mesh_v, viewer_data_tessellation.mesh_f,
                                         triangulation_cell_indices);
    viewer_data_tessellation.mesh_c.resize(triangulation_cell_indices.rows(), 3);
    for (int i = 0; i < triangulation_cell_indices.rows(); i++) {
        viewer_data_tessellation.mesh_c.row(i) = colors.row(triangulation_cell_indices(i));
    }

    /// Copy camera state.
    viewer_camera = app_camera;
}

void FoamSubApp2D::initializeSubApp() {
    experiment_selector.second.clear();
    experiment_selector.second.emplace_back(std::make_shared<NoExperimentSelected>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentCoarsening2D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentConvergenceTest2D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentRigidBody2D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentDCMComparison2D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentImageMatch2D>());
    experiment_selector.second.emplace_back(std::make_shared<ExperimentRuntimeScaling2D>());

    scenario_selector.second.clear();
    scenario_selector.second.push_back(std::make_shared<RandomSitesInBox2D>());
    scenario_selector.second.push_back(std::make_shared<RandomSitesInMembrane2D>());
    scenario_selector.second.push_back(std::make_shared<ConvergenceTest2D>());
    scenario_selector.second.push_back(std::make_shared<RigidBody2D>());
    scenario_selector.second[scenario_selector.first]->assignScenario(model_definition, degrees_of_freedom);

    tessellation_selector.second.clear();
    tessellation_selector.second.push_back(std::make_shared<Voronoi2D>());
    tessellation_selector.second.push_back(std::make_shared<Power2D>());

    energy_selector.second.clear();
    energy_selector.second.emplace_back("Standard", std::make_shared<CellEnergy2D>());
    energy_selector.second.emplace_back("Cleavage", std::make_shared<CellEnergyCleavage2D>());
    energy_selector.second.emplace_back("Coarsening", std::make_shared<CellEnergyCoarsening2D>());

    model_definition.tessellation_generator = tessellation_selector.second[tessellation_selector.first];
    model_definition.cell_energy_function = energy_selector.second[energy_selector.first].second;

    /// Initialize camera.
    app_camera.eye = Vector3F(0, 0, 5);
    app_camera.center = Vector3F(0, 0, 0);
    app_camera.set_up_direction(Vector3F(0, 1, 0));
    app_camera.height = 1;

    /// Initialization from superclass.
    FoamSubApp::initializeSubApp();
}

bool FoamSubApp2D::callbackMouseMove(const CRLControlState &control_state, const Vector2F &mouse_pos,
                                     const Vector2F &mouse_delta) {
    if (control_state.modifiers[MOUSE_LEFT]) {
        CameraHelper::panCameraFromDrag(app_camera, mouse_delta, CAMERA_PAN_SENSITIVITY);
    }

    return false;
}

bool FoamSubApp2D::callbackMouseScroll(const CRLControlState &control_state, float t) {
    Vector3F camera_dist = app_camera.eye - app_camera.center;
    app_camera.eye = app_camera.center + camera_dist * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);
    app_camera.height = app_camera.height * (1.0 - CAMERA_ZOOM_SENSITIVITY * t);

    return false;
}
