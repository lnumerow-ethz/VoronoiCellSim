#pragma once

#include "Projects/VoronoiFoam/include/Model/Energy/PerCellFunction.h"
#include "Projects/VoronoiFoam/include/Config/Config.h"

enum EnergyTermWeight2D {
    AREA,  //
    PERIMETER_TARGET,
    PERIMETER_MINIMIZATION,
    CENTROID,
    SECOND_MOMENT,
    NUM_ENERGY_TERMS
};

class CellEnergy2D : public PerCellFunction {
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
