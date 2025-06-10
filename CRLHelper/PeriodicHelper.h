#pragma once

#include "CRLHelper/VecMatDef.h"

namespace PeriodicHelper {

MatrixXI getPeriodicTiles(int depth, int dims, bool include_zero = false, bool include_negative = false);
int getTileIndex(const VectorXI &tile_coords, int depth, bool include_zero = false, bool include_negative = false);

}  // namespace PeriodicHelper