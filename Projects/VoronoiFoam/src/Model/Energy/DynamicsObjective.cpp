#include "Projects/VoronoiFoam/include/Model/Energy/DynamicsObjective.h"
#include "Projects/VoronoiFoam/include/Model/Energy/PotentialEnergy.h"

#include <iostream>

static VectorXF get_a_1(const VectorXF &y, const DynamicsData &dynamics_data) {
    return (y - 2 * dynamics_data.prev_dof_1 + dynamics_data.prev_dof_2) / pow(dynamics_data.dt, 2.0);
}
static F get_dady_1(const DynamicsData &dynamics_data) {
    return 1.0 / pow(dynamics_data.dt, 2.0);
}

static VectorXF get_a_2(const VectorXF &y, const DynamicsData &dynamics_data) {
    return 0.5 * (3 * y - 7 * dynamics_data.prev_dof_1 + 5 * dynamics_data.prev_dof_2 - dynamics_data.prev_dof_3) /
           pow(dynamics_data.dt, 2.0);
}
static F get_dady_2(const DynamicsData &dynamics_data) {
    return 1.5 / pow(dynamics_data.dt, 2.0);
}

static VectorXF get_v(const VectorXF &y, const DynamicsData &dynamics_data) {
    return (y - dynamics_data.prev_dof_1) / dynamics_data.dt;
}
static F get_dvdy(const DynamicsData &dynamics_data) {
    return 1.0 / dynamics_data.dt;
}

static VectorXF get_a(const VectorXF &y, const DynamicsData &dynamics_data) {
    VectorXF a;
    switch (dynamics_data.discretization) {
        case FIRST_ORDER:
            a = get_a_1(y, dynamics_data);
            break;
        case SECOND_ORDER:
            a = get_a_2(y, dynamics_data);
            break;
    }
    return a;
}
static F get_dady(const DynamicsData &dynamics_data) {
    F dady;
    switch (dynamics_data.discretization) {
        case FIRST_ORDER:
            dady = get_dady_1(dynamics_data);
            break;
        case SECOND_ORDER:
            dady = get_dady_2(dynamics_data);
            break;
    }
    return dady;
}

static void addMomentumAndViscosity(const VectorXF &y, const DynamicsData &dynamics_data, F &energy) {
    VectorXF a = get_a(y, dynamics_data);
    F dady = get_dady(dynamics_data);
    VectorXF v = get_v(y, dynamics_data);
    F dvdy = get_dvdy(dynamics_data);
    energy += 0.5 / dady * a.transpose() * dynamics_data.mass_matrix_diagonal.asDiagonal() * a;
    energy += 0.5 / dvdy * v.transpose() * dynamics_data.viscosity_matrix_diagonal.asDiagonal() * v;
}
static void addMomentumAndViscosity(const VectorXF &y, const DynamicsData &dynamics_data, F &energy,
                                    VectorXF &gradient) {
    addMomentumAndViscosity(y, dynamics_data, energy);

    VectorXF a = get_a(y, dynamics_data);
    VectorXF v = get_v(y, dynamics_data);
    gradient += dynamics_data.mass_matrix_diagonal.asDiagonal() * a;
    gradient += dynamics_data.viscosity_matrix_diagonal.asDiagonal() * v;
}
static void addMomentumAndViscosity(const VectorXF &y, const DynamicsData &dynamics_data, F &energy, VectorXF &gradient,
                                    HessianF &hessian) {
    addMomentumAndViscosity(y, dynamics_data, energy, gradient);

    F dady = get_dady(dynamics_data);
    F dvdy = get_dvdy(dynamics_data);
    hessian.A += dady * dynamics_data.mass_matrix_diagonal.asDiagonal().toDenseMatrix().sparseView();
    hessian.A += dvdy * dynamics_data.viscosity_matrix_diagonal.asDiagonal().toDenseMatrix().sparseView();
}

static void adamsMoulton(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                         const DynamicsData &dynamics_data, const VectorXF &y, F &energy, VectorXF &gradient,
                         HessianF &hessian) {
    F prev_energy;
    VectorXF prev_gradient;
    PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, dynamics_data.prev_dof_1,
                                                         prev_energy, prev_gradient);

    energy = 0.5 * energy + 0.5 * y.transpose() * prev_gradient;
    gradient *= 0.5;
    gradient += 0.5 * prev_gradient;
    hessian.A *= 0.5;
    hessian.U *= 0.5;
}

static void adamsMoulton(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                         const DynamicsData &dynamics_data, const VectorXF &y, F &energy, VectorXF &gradient) {
    HessianF temp_hessian;
    adamsMoulton(model_definition, model_dof, dynamics_data, y, energy, gradient, temp_hessian);
}

static void adamsMoulton(const ModelDefinition &model_definition, const DegreesOfFreedom &model_dof,
                         const DynamicsData &dynamics_data, const VectorXF &y, F &energy) {
    VectorXF temp_gradient = 0 * dynamics_data.prev_dof_1;
    HessianF temp_hessian;
    adamsMoulton(model_definition, model_dof, dynamics_data, y, energy, temp_gradient, temp_hessian);
}

bool DynamicsObjective::computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition,
                                                              const DegreesOfFreedom &model_dof,
                                                              const DynamicsData &dynamics_data, const VectorXF &y,
                                                              F &energy) {
    bool success = PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy);

    // if (dynamics_data.discretization == SECOND_ORDER)
    //     adamsMoulton(model_definition, model_dof, dynamics_data, y, energy);
    addMomentumAndViscosity(y, dynamics_data, energy);

    return success;
}

bool DynamicsObjective::computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition,
                                                              const DegreesOfFreedom &model_dof,
                                                              const DynamicsData &dynamics_data, const VectorXF &y,
                                                              F &energy, VectorXF &gradient) {
    bool success =
        PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy, gradient);

    // if (dynamics_data.discretization == SECOND_ORDER)
    //     adamsMoulton(model_definition, model_dof, dynamics_data, y, energy, gradient);
    addMomentumAndViscosity(y, dynamics_data, energy, gradient);

    return success;
}

bool DynamicsObjective::computeDynamicsObjectiveFromDOFVector(const ModelDefinition &model_definition,
                                                              const DegreesOfFreedom &model_dof,
                                                              const DynamicsData &dynamics_data, const VectorXF &y,
                                                              F &energy, VectorXF &gradient, HessianF &hessian) {
    bool success =
        PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy, gradient, hessian);

    // if (dynamics_data.discretization == SECOND_ORDER)
    //     adamsMoulton(model_definition, model_dof, dynamics_data, y, energy, gradient, A, U);
    addMomentumAndViscosity(y, dynamics_data, energy, gradient, hessian);

    return success;
}

bool DynamicsObjective::computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                                                  const DegreesOfFreedom &model_dof,
                                                                  const DynamicsData &dynamics_data, bool use_dynamics,
                                                                  const VectorXF &y, F &energy) {
    return use_dynamics ? computeDynamicsObjectiveFromDOFVector(model_definition, model_dof, dynamics_data, y, energy)
                        : PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy);
}

bool DynamicsObjective::computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                                                  const DegreesOfFreedom &model_dof,
                                                                  const DynamicsData &dynamics_data, bool use_dynamics,
                                                                  const VectorXF &y, F &energy, VectorXF &gradient) {
    return use_dynamics
               ? computeDynamicsObjectiveFromDOFVector(model_definition, model_dof, dynamics_data, y, energy, gradient)
               : PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy, gradient);
}

bool DynamicsObjective::computePotentialEnergyOrDynamicsObjective(const ModelDefinition &model_definition,
                                                                  const DegreesOfFreedom &model_dof,
                                                                  const DynamicsData &dynamics_data, bool use_dynamics,
                                                                  const VectorXF &y, F &energy, VectorXF &gradient,
                                                                  HessianF &hessian) {
    return use_dynamics ? computeDynamicsObjectiveFromDOFVector(model_definition, model_dof, dynamics_data, y, energy,
                                                                gradient, hessian)
                        : PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, model_dof, y, energy,
                                                                               gradient, hessian);
}

void DynamicsObjective::nextDynamicsStep(DynamicsData &dynamics_data, const VectorXF &y) {
    dynamics_data.prev_dof_3 = dynamics_data.prev_dof_2;
    dynamics_data.prev_dof_2 = dynamics_data.prev_dof_1;
    dynamics_data.prev_dof_1 = y;
}
