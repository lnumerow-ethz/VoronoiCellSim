#pragma once

#include "Projects/VoronoiFoam/include/Model/Model.h"

namespace Rendering {
void writeFile2D(const Model &model, std::string file_name);
void writeFile3D(const Model &model, std::string file_name, bool include_internal_faces = false, bool include_boundary_edges = false);
void writeFile(const Model &model, std::string file_name);

void triangulateVoronoiCells2D(const Model &model, MatrixXF &igl_V, MatrixXI &igl_F, VectorXI &cell_idx);
}  // namespace Rendering
