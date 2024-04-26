#pragma once

#include "CRLHelper/CRLSubApp.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/Energy/DynamicsObjective.h"
#include "Projects/VoronoiFoam/include/App/Experiment/Experiment.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"
#include "CRLHelper/Optimization.h"

#include "ThirdParty/LBFGSpp/include/LBFGSpp/BFGSMat.h"

class FoamSubApp : public CRLSubApp {
   protected:
    /// Foam model
    ModelDefinition model_definition;
    DegreesOfFreedom degrees_of_freedom;
    std::pair<int, std::vector<std::shared_ptr<TessellationGenerator>>> tessellation_selector = {0, {}};
    std::pair<int, std::vector<std::pair<std::string, std::shared_ptr<PerCellFunction>>>> energy_selector = {0, {}};
    bool site_param_free[NUM_SITE_PARAM_TOTAL] = {true, false, false};

    /// Dynamics
    DynamicsData dynamics_data;
    bool use_dynamics = false;
    F dynamics_mass_pos = DYNAMICS_MASS;
    std::vector<F> dynamics_mass_param = std::vector<F>(NUM_SITE_PARAM_TOTAL, DYNAMICS_MASS);
    F dynamics_mass_boundary = DYNAMICS_MASS;
    F dynamics_viscosity_pos = DYNAMICS_VISCOSITY;
    std::vector<F> dynamics_viscosity_param = std::vector<F>(NUM_SITE_PARAM_TOTAL, DYNAMICS_VISCOSITY);
    F dynamics_viscosity_boundary = DYNAMICS_VISCOSITY;

    /// Scenario and preset experiments
    std::pair<int, std::vector<std::shared_ptr<Experiment>>> experiment_selector = {0, {}};
    std::pair<int, std::vector<std::shared_ptr<Scenario>>> scenario_selector = {0, {}};

    /// Camera state
    CRLCamera app_camera;

    /// Optimizer
    Optimization::Optimizer optimizer = Optimization::NEWTON;
    bool optimize = false;
    MatrixXF bfgs_inverse_hessian;
    LBFGSpp::BFGSMat<F> lbfgs_mat;
    int lbfgs_m = 6;

    /// Analysis
    bool check_gradient_print_all = false;
    bool check_hessian_print_all = false;
    int check_gradient_epsilon_exponent = -4;
    int check_hessian_epsilon_exponent = -4;

    /// Output
    bool write_png = false;
    std::string png_file_name;
    bool write_rendering_file = false;
    std::string render_file_name;

   protected:
    void clearState();
    void resetOptimization();
    void applyAllSettings();
    void setDynamicsMatrices();

   private:
    void optimizationConverged(const VectorXF &dof_vector);

    void checkEnergyGradients(int order, F epsilon, bool print_all);

   public:
    Optimization::OptimizationStatus energyMinimizationStep();

    void initializeSubApp() override;

    void mainLoop() override;

    void makeConfigWindow() override;

    void makeAnalysisWindow() override;

    bool callbackKeyPressed(const CRLControlState &control_state, int key) override;

   public:
    friend class ExperimentCoarsening2D;
    friend class ExperimentConvergenceTest2D;
    friend class ExperimentRigidBody2D;
    friend class ExperimentDCMComparison2D;
    friend class ExperimentImageMatch2D;
    friend class ExperimentRuntimeScaling2D;
    friend class ExperimentCoarsening3D;
    friend class ExperimentCleavage3D;
    friend class ExperimentCleavageCylinder3D;
    friend class ExperimentVideoDemo3D;
    friend class ExperimentRuntimeScaling3D;

   public:
    void runExperiment(int experiment);

    void generateRenderOutput();
};
