#pragma once

#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

class CellEnergyCoarsening2D : public PerCellFunction {
   public:
    F weights[NUM_ENERGY_TERMS] = {ENERGY_2D_WEIGHT_AREA,                    //
                                   ENERGY_2D_WEIGHT_PERIMETER_TARGET,        //
                                   ENERGY_2D_WEIGHT_PERIMETER_MINIMIZATION,  //
                                   ENERGY_2D_WEIGHT_CENTROID,                //
                                   ENERGY_2D_WEIGHT_SECOND_MOMENT};

   public:
    void getValue(const Model &model, PerCellValue &cell_value) const override;

    void makeConfigMenu() override;
};
