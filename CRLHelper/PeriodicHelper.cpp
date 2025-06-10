#include "CRLHelper/PeriodicHelper.h"

#include <iostream>

MatrixXI PeriodicHelper::getPeriodicTiles(int depth, int dims, bool include_zero, bool include_negative) {
    int grid_size = 2 * depth + 1;
    int num_pairs_one_side = (pow(grid_size, dims) - 1) / 2;
    int num_pairs_total = num_pairs_one_side * (include_negative ? 2 : 1) + (include_zero ? 1 : 0);
    MatrixXI tiles = MatrixXI::Zero(num_pairs_total, dims);

    if (dims > 0) {
        MatrixXI tiles_one_side = MatrixXI::Zero(num_pairs_one_side, dims);

        int i = 0;
        VectorXI current_tile = VectorXI::Zero(dims);
        current_tile(dims - 1) = 1;

        while (i < num_pairs_one_side) {
            tiles_one_side.row(i) = current_tile;

            for (int j = dims - 1; j >= 0; j--) {
                if (current_tile(j) == depth) {
                    current_tile(j) = -depth;
                } else {
                    current_tile(j)++;
                    break;
                }
            }

            i++;
        }

        int curr_row = 0;
        if (include_zero) {
            tiles.row(curr_row) = VectorXI::Zero(dims);
            curr_row++;
        }
        tiles.middleRows(curr_row, num_pairs_one_side) = tiles_one_side;
        curr_row += num_pairs_one_side;
        if (include_negative) {
            tiles.bottomRows(num_pairs_one_side) = -tiles_one_side;
        }
    }

    return tiles;
}

int PeriodicHelper::getTileIndex(const VectorXI &tile_coords, int depth, bool include_zero, bool include_negative) {
    int dims = tile_coords.rows();

    MatrixXI tiles = getPeriodicTiles(depth, dims, include_zero, include_negative);
    for (int i = 0; i < tiles.rows(); i++) {
        VectorXI tile = tiles.row(i);
        if (tile == tile_coords) {
            return i;
        }
    }

    return -1;
}