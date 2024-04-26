#pragma once

#include "Projects/VoronoiFoam/include/Model/Energy/Energy3D/CellEnergy3D.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class CellEnergyCleavage3D : public PerCellFunction {
   public:
    F weights[NUM_ENERGY_TERMS] = {ENERGY_3D_WEIGHT_VOLUME,                //
                                   ENERGY_3D_WEIGHT_SURFACE_TARGET,        //
                                   ENERGY_3D_WEIGHT_SURFACE_MINIMIZATION,  //
                                   ENERGY_3D_WEIGHT_CENTROID,              //
                                   ENERGY_3D_WEIGHT_SECOND_MOMENT,         //
                                   ENERGY_3D_WEIGHT_GRAVITY};

   public:
    void getValue(const Model &model, PerCellValue &cell_value) const override;

    void makeConfigMenu() override;
};
