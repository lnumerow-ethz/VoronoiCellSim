#pragma once

#define CAMERA_ROTATE_SENSITIVITY (1.0 / 300.0)
#define CAMERA_PAN_SENSITIVITY (1.0 / 300.0)
#define CAMERA_ZOOM_SENSITIVITY (1.0 / 50.0)

#define SCENARIO_SITES_IN_BOX_2D_NUM_SITES (10)
#define SCENARIO_SITES_IN_MEMBRANE_2D_NUM_SITES (10)
#define SCENARIO_SITES_IN_MEMBRANE_2D_NUM_VERTICES (40)
#define SCENARIO_SITES_IN_MEMBRANE_2D_SPRING_CONSTANT (0.01)
#define SCENARIO_SITES_IN_BOX_3D_NUM_SITES (5)
#define SCENARIO_SITES_IN_MEMBRANE_3D_NUM_SITES (10)
#define SCENARIO_SITES_IN_MEMBRANE_3D_SUBDIV_LEVEL (0)
#define SCENARIO_SITES_IN_MEMBRANE_3D_SPRING_CONSTANT (0.01)

#define ENERGY_2D_WEIGHT_AREA (0.05)
#define ENERGY_2D_WEIGHT_PERIMETER_TARGET (0.001)
#define ENERGY_2D_WEIGHT_PERIMETER_MINIMIZATION (0.0)
#define ENERGY_2D_WEIGHT_CENTROID (0.1)
#define ENERGY_2D_WEIGHT_SECOND_MOMENT (0.0)
#define ENERGY_2D_STEP_AREA (0.01)
#define ENERGY_2D_STEP_PERIMETER_TARGET (0.001)
#define ENERGY_2D_STEP_PERIMETER_MINIMIZATION (0.001)
#define ENERGY_2D_STEP_CENTROID (0.1)
#define ENERGY_2D_STEP_SECOND_MOMENT (0.1)

#define ENERGY_3D_WEIGHT_VOLUME (0.2)
#define ENERGY_3D_WEIGHT_SURFACE_TARGET (0.01)
#define ENERGY_3D_WEIGHT_SURFACE_MINIMIZATION (0.0)
#define ENERGY_3D_WEIGHT_CENTROID (0.1)
#define ENERGY_3D_WEIGHT_SECOND_MOMENT (0.0)
#define ENERGY_3D_WEIGHT_GRAVITY (0.0)
#define ENERGY_3D_STEP_VOLUME (0.1)
#define ENERGY_3D_STEP_SURFACE_TARGET (0.01)
#define ENERGY_3D_STEP_SURFACE_MINIMIZATION (0.01)
#define ENERGY_3D_STEP_CENTROID (0.1)
#define ENERGY_3D_STEP_SECOND_MOMENT (0.1)
#define ENERGY_3D_STEP_GRAVITY (0.1)

#define ENERGY_VOLUME_BARRIER_WEIGHT (1e-5)
#define ENERGY_POWER_REGULARIZER_WEIGHT (0.0)

#define DYNAMICS_TOLERANCE (1e-10)
#define DYNAMICS_DT (0.01)
#define DYNAMICS_MASS (0.0)
#define DYNAMICS_VISCOSITY (0.01)
#define DYNAMICS_STEP_DT (0.01)
#define DYNAMICS_STEP_MASS (0.01)
#define DYNAMICS_STEP_VISCOSITY (0.01)

#define DEFAULT_CONVERGENCE_TOLERANCE (1e-10)