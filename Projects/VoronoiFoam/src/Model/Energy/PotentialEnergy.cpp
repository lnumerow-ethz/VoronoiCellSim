#include "Projects/VoronoiFoam/include/Model/Energy/PotentialEnergy.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"
#include "Projects/VoronoiFoam/include/Model/Energy/PerCellFunction.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#define TBB_SUPPRESS_DEPRECATED_MESSAGES (true)

#include <tbb/tbb.h>

#include <iostream>
#include "CRLHelper/CRLTimer.h"

void addTripletsMulTriangular(TripletListF &v, const TripletListF &M, double mul) {
    TripletListF temp;
    for (auto triplet : M) {
        int row = triplet.row();
        int col = triplet.col();
        if (row <= col) temp.emplace_back(row, col, mul * triplet.value());
    }
    v.insert(v.end(), temp.begin(), temp.end());
}

template <typename TemplateTripletList1>
void addTripletsMul(TripletListF &v, const TemplateTripletList1 &M, double mul) {
    TripletListF temp;
    for (auto triplet : M) {
        temp.emplace_back(triplet.row(), triplet.col(), mul * triplet.value());
    }
    v.insert(v.end(), temp.begin(), temp.end());
}

void addTripletsBlockOffsetTriangular(TripletListF &v, const SparseMatrixF &M, int row_offset, int col_offset) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            int row = it.row() + row_offset;
            int col = it.col() + col_offset;
            if (row <= col) temp.emplace_back(row, col, it.value());
        }
    }
    v.insert(v.end(), temp.begin(), temp.end());
}

void addTripletsBlockOffset(TripletListF &v, const SparseMatrixF &M, int row_offset, int col_offset) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            temp.emplace_back(it.row() + row_offset, it.col() + col_offset, it.value());
        }
    }
    v.insert(v.end(), temp.begin(), temp.end());
}

void addTripletsBlockOffset(tbb::concurrent_vector<TripletF> &v, const SparseMatrixF &M, int row_offset,
                            int col_offset) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            temp.emplace_back(it.row() + row_offset, it.col() + col_offset, it.value());
        }
    }
    v.grow_by(temp.begin(), temp.end());
}

void addTriplets(tbb::concurrent_vector<TripletF> &v, const SparseMatrixF &M) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            temp.emplace_back(it.row(), it.col(), it.value());
        }
    }
    v.grow_by(temp.begin(), temp.end());
}

static void addTripletsParallel_mul_A_B_C(tbb::concurrent_vector<TripletF> &triplets, const SparseMatrixF &A,
                                          const SparseMatrixF &B, const SparseMatrixF &C) {
    int num_rows = A.rows();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_rows, 1000), [&](const tbb::blocked_range<size_t> &range) {
        Eigen::SparseMatrix<F, 1> temp = A.middleRows(range.begin(), range.end() - range.begin()) * B;
        Eigen::SparseMatrix<F, 1> rows = temp * C;
        // Eigen::SparseMatrix<F, 1> rows = A.middleRows(range.begin(), range.end() - range.begin()) * B * C;

        TripletListF local_triplets;
        addTripletsBlockOffset(local_triplets, rows, range.begin(), 0);

        triplets.grow_by(local_triplets.begin(), local_triplets.end());
    });
}

static void addTripletsParallel_mul_A_B_C_blockOffset(tbb::concurrent_vector<TripletF> &triplets,
                                                      const SparseMatrixF &A, const SparseMatrixF &B,
                                                      const SparseMatrixF &C, int row_offset, int col_offset) {
    int num_rows = A.rows();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_rows, 1000), [&](const tbb::blocked_range<size_t> &range) {
        Eigen::SparseMatrix<F, 1> temp = A.middleRows(range.begin(), range.end() - range.begin()) * B;
        Eigen::SparseMatrix<F, 1> rows = temp * C;
        // Eigen::SparseMatrix<F, 1> rows = A.middleRows(range.begin(), range.end() - range.begin()) * B * C;

        TripletListF local_triplets;
        addTripletsBlockOffset(local_triplets, rows, range.begin() + row_offset, col_offset);

        triplets.grow_by(local_triplets.begin(), local_triplets.end());
    });
}

static void addTripletsParallel_mul_A_B_C_triangular(tbb::concurrent_vector<TripletF> &triplets, const SparseMatrixF &A,
                                                     const SparseMatrixF &B, const SparseMatrixF &C) {
    int num_rows = A.rows();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_rows, 1000), [&](const tbb::blocked_range<size_t> &range) {
        Eigen::SparseMatrix<F, 1> temp = A.middleRows(range.begin(), range.end() - range.begin()) * B;
        Eigen::SparseMatrix<F, 1> rows = temp * C;
        // Eigen::SparseMatrix<F, 1> rows = A.middleRows(range.begin(), range.end() - range.begin()) * B * C;

        TripletListF local_triplets;
        addTripletsBlockOffsetTriangular(local_triplets, rows, range.begin(), 0);

        triplets.grow_by(local_triplets.begin(), local_triplets.end());
    });
}

static void addTripletsParallel_mul_A_B_C_blockOffsetTriangular(tbb::concurrent_vector<TripletF> &triplets,
                                                                const SparseMatrixF &A, const SparseMatrixF &B,
                                                                const SparseMatrixF &C, int row_offset,
                                                                int col_offset) {
    int num_rows = A.rows();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_rows, 1000), [&](const tbb::blocked_range<size_t> &range) {
        Eigen::SparseMatrix<F, 1> temp = A.middleRows(range.begin(), range.end() - range.begin()) * B;
        Eigen::SparseMatrix<F, 1> rows = temp * C;
        // Eigen::SparseMatrix<F, 1> rows = A.middleRows(range.begin(), range.end() - range.begin()) * B * C;

        TripletListF local_triplets;
        addTripletsBlockOffsetTriangular(local_triplets, rows, range.begin() + row_offset, col_offset);

        triplets.grow_by(local_triplets.begin(), local_triplets.end());
    });
}

void PotentialEnergy::computePotentialEnergy(const Model &model, F &energy, VectorXF &gradient, HessianF &hessian,
                                             int order) {
    int nc = model.dimensions_ind->nc;
    int np = model.dimensions_ind->np;
    int nx = model.dimensions_tess->nx;
    int nv = model.dimensions_boundary->nv;

    int dims_c = model.dimensions_ind->dims_site_dof;
    int dims_x = model.dimensions_ind->dims_space;

    std::vector<F> cell_energies(model.dimensions_ind->n_sites);
    tbb::concurrent_vector<TripletF> triplets_dFdx;
    tbb::concurrent_vector<TripletF> triplets_dFdc;
    tbb::concurrent_vector<TripletF> triplets_d2Fdc2;  /// Only upper triangular coefficients.
    tbb::concurrent_vector<TripletF> triplets_d2Fdcdx;
    tbb::concurrent_vector<TripletF> triplets_d2Fdx2;  /// Only upper triangular coefficients.

    CRLTimer timer;

    /// This usage of parallel_for (with blocked range and local vectors) doesn't make much any difference on my
    /// laptop, but it's way faster on the cluster.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, model.cells.size()), [&](const tbb::blocked_range<size_t> &range) {
        TripletListF local_triplets_dFdx;
        TripletListF local_triplets_dFdc;
        TripletListF local_triplets_d2Fdc2;
        TripletListF local_triplets_d2Fdcdx;
        TripletListF local_triplets_d2Fdx2;

        for (size_t i = range.begin(); i != range.end(); ++i) {
            const TessellationCell &cell = model.cells[i];

            PerCellValue cell_energy_value(model, cell, order);
            model.model_definition.cell_energy_function->getValue(model, cell_energy_value);

            cell_energies[cell.site_index] = cell_energy_value.value;

            /// Not necessary when order == 0, but putting here for clarity.
            int num_nodes_in_cell = cell.node_indices_in_cell.size();
            VectorXI nodes_in_cell(num_nodes_in_cell);
            for (const auto &n : cell.node_indices_in_cell) {
                nodes_in_cell(n.second) = n.first;
            }

            if (order >= 1) {
                for (int i = 0; i < num_nodes_in_cell; i++) {
                    for (int j = 0; j < dims_x; j++) {
                        local_triplets_dFdx.emplace_back(nodes_in_cell(i) * dims_x + j, 0,
                                                         cell_energy_value.gradient(i * dims_x + j));
                    }
                }
                for (int i = 0; i < dims_c; i++) {
                    local_triplets_dFdc.emplace_back(cell.site_index * dims_c + i, 0,
                                                     cell_energy_value.gradient(num_nodes_in_cell * dims_x + i));
                }
            }

            if (order >= 2) {
                MatrixXF cell_energy_hessian = cell_energy_value.hessian.evalDense();
                for (int ii = 0; ii < cell_energy_hessian.rows(); ii++) {
                    for (int jj = 0; jj < cell_energy_hessian.cols(); jj++) {
                        int node0_index_in_cell = std::min(ii / dims_x, num_nodes_in_cell);
                        int node1_index_in_cell = std::min(jj / dims_x, num_nodes_in_cell);
                        int i = ii - dims_x * node0_index_in_cell;
                        int j = jj - dims_x * node1_index_in_cell;

                        if (node0_index_in_cell < num_nodes_in_cell) {
                            int row = nodes_in_cell(node0_index_in_cell) * dims_x + i;
                            if (node1_index_in_cell < num_nodes_in_cell) {
                                int col = nodes_in_cell(node1_index_in_cell) * dims_x + j;
                                if (col < row) continue;

                                local_triplets_d2Fdx2.emplace_back(
                                    row, col,
                                    cell_energy_hessian.coeff(node0_index_in_cell * dims_x + i,
                                                              node1_index_in_cell * dims_x + j));
                            }
                        } else {
                            int row = cell.site_index * dims_c + i;
                            if (node1_index_in_cell < num_nodes_in_cell) {
                                int col = nodes_in_cell(node1_index_in_cell) * dims_x + j;
                                local_triplets_d2Fdcdx.emplace_back(
                                    row, col,
                                    cell_energy_hessian.coeff(num_nodes_in_cell * dims_x + i,
                                                              node1_index_in_cell * dims_x + j));
                            } else {
                                int col = cell.site_index * dims_c + j;
                                if (col < row) continue;
                                // d2Fdc2 added directly
                                local_triplets_d2Fdc2.emplace_back(
                                    row, col,
                                    cell_energy_hessian.coeff(num_nodes_in_cell * dims_x + i,
                                                              num_nodes_in_cell * dims_x + j));
                            }
                        }
                    }
                }
            }
        }

        // Safely append the local vector to the concurrent vector
        triplets_dFdx.grow_by(local_triplets_dFdx.begin(), local_triplets_dFdx.end());
        triplets_dFdc.grow_by(local_triplets_dFdc.begin(), local_triplets_dFdc.end());
        triplets_d2Fdc2.grow_by(local_triplets_d2Fdc2.begin(), local_triplets_d2Fdc2.end());
        triplets_d2Fdcdx.grow_by(local_triplets_d2Fdcdx.begin(), local_triplets_d2Fdcdx.end());
        triplets_d2Fdx2.grow_by(local_triplets_d2Fdx2.begin(), local_triplets_d2Fdx2.end());
    });

    energy = 0;
    for (F e : cell_energies) {
        energy += e;
    }
    energy += model.model_definition.boundary_generator->computeEnergy(model);

    if (order >= 1) {
        SparseMatrixF dFdx_sparse(nx, 1);
        dFdx_sparse.setFromTriplets(triplets_dFdx.begin(), triplets_dFdx.end());
        SparseMatrixF dFdc_sparse(nc, 1);
        dFdc_sparse.setFromTriplets(triplets_dFdc.begin(), triplets_dFdc.end());
        VectorXF dFdx = dFdx_sparse;
        VectorXF dFdc = dFdc_sparse;
        VectorXF dFdp;
        model.model_definition.boundary_generator->computeEnergyGradient(model, dFdp);

        gradient.resize(nc + np);
        gradient.head(nc) = dFdx.transpose() * model.derivative_matrices.dxdc + dFdc.transpose();
        gradient.tail(np) =
            dFdx.transpose() * model.derivative_matrices.dxdv * model.derivative_matrices.dvdp + dFdp.transpose();
    }

    if (order >= 2) {
        SparseMatrixF dFdx(nx, 1);
        SparseMatrixF d2Fdc2(nc, nc);
        Eigen::SparseMatrix<F, Eigen::Lower> d2Fdcdx(nc, nx);
        SparseMatrixF d2Fdx2(nx, nx);

        /// Not running these in tbb::task_group because d2Fdx2 takes virtually all the time.
        dFdx.setFromTriplets(triplets_dFdx.begin(), triplets_dFdx.end());
        d2Fdc2.setFromTriplets(triplets_d2Fdc2.begin(), triplets_d2Fdc2.end());
        d2Fdcdx.setFromTriplets(triplets_d2Fdcdx.begin(), triplets_d2Fdcdx.end());
        d2Fdx2.setFromTriplets(triplets_d2Fdx2.begin(), triplets_d2Fdx2.end());

        HessianF d2Fdp2;
        model.model_definition.boundary_generator->computeEnergyHessian(model, d2Fdp2);

        tbb::concurrent_vector<TripletF> triplets_sum_dFdx_d2xdc2;
        tbb::concurrent_vector<TripletF> triplets_sum_dFdx_d2xdcdv;
        tbb::concurrent_vector<TripletF> triplets_sum_dFdx_d2xdv2;
        tbb::concurrent_vector<TripletF> triplets_sum_dFdx_dxdv_d2vdp2;

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nx), [&](const tbb::blocked_range<size_t> &range) {
            TripletListF local_triplets_sum_dFdx_d2xdc2;
            TripletListF local_triplets_sum_dFdx_d2xdcdv;
            TripletListF local_triplets_sum_dFdx_d2xdv2;

            for (size_t i = range.begin(); i != range.end(); ++i) {
                addTripletsMul(local_triplets_sum_dFdx_d2xdc2, model.derivative_matrices.d2xdc2[i], dFdx.coeff(i, 0));
                addTripletsMul(local_triplets_sum_dFdx_d2xdcdv, model.derivative_matrices.d2xdcdv[i], dFdx.coeff(i, 0));
                addTripletsMul(local_triplets_sum_dFdx_d2xdv2, model.derivative_matrices.d2xdv2[i], dFdx.coeff(i, 0));
            }

            triplets_sum_dFdx_d2xdc2.grow_by(local_triplets_sum_dFdx_d2xdc2.begin(),
                                             local_triplets_sum_dFdx_d2xdc2.end());
            triplets_sum_dFdx_d2xdcdv.grow_by(local_triplets_sum_dFdx_d2xdcdv.begin(),
                                              local_triplets_sum_dFdx_d2xdcdv.end());
            triplets_sum_dFdx_d2xdv2.grow_by(local_triplets_sum_dFdx_d2xdv2.begin(),
                                             local_triplets_sum_dFdx_d2xdv2.end());
        });

        VectorXF temp_dFdx_dxdv = (dFdx.transpose() * model.derivative_matrices.dxdv).transpose();

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nv), [&](const tbb::blocked_range<size_t> &range) {
            TripletListF local_triplets_sum_dFdx_dxdv_d2vdp2;

            for (size_t i = range.begin(); i != range.end(); ++i) {
                addTripletsMul(local_triplets_sum_dFdx_dxdv_d2vdp2, model.derivative_matrices.d2vdp2[i],
                               temp_dFdx_dxdv(i));
            }

            triplets_sum_dFdx_dxdv_d2vdp2.grow_by(local_triplets_sum_dFdx_dxdv_d2vdp2.begin(),
                                                  local_triplets_sum_dFdx_dxdv_d2vdp2.end());
        });

        SparseMatrixF sum_dFdx_d2xdc2(nc, nc);
        Eigen::SparseMatrix<F, Eigen::Lower> sum_dFdx_d2xdcdv(nc, nv);
        SparseMatrixF sum_dFdx_d2xdv2(nv, nv);
        SparseMatrixF sum_dFdx_dxdv_d2vdp2(np, np);

        sum_dFdx_d2xdcdv.setFromTriplets(triplets_sum_dFdx_d2xdcdv.begin(), triplets_sum_dFdx_d2xdcdv.end());
        sum_dFdx_d2xdv2.setFromTriplets(triplets_sum_dFdx_d2xdv2.begin(), triplets_sum_dFdx_d2xdv2.end());
        sum_dFdx_dxdv_d2vdp2.setFromTriplets(triplets_sum_dFdx_dxdv_d2vdp2.begin(),
                                             triplets_sum_dFdx_dxdv_d2vdp2.end());

        SparseMatrixF dxdc = model.derivative_matrices.dxdc;
        Eigen::SparseMatrix<F, Eigen::Lower> dxdcT = dxdc.transpose();
        SparseMatrixF d2Fdxdc = d2Fdcdx.transpose();
        SparseMatrixF dvdp = model.derivative_matrices.dvdp;
        Eigen::SparseMatrix<F, Eigen::Lower> dvdpT = dvdp.transpose();
        SparseMatrixF dxdv_dvdp = model.derivative_matrices.dxdv * dvdp;
        Eigen::SparseMatrix<F, Eigen::Lower> dxdv_dvdpT = dxdv_dvdp.transpose();

        SparseMatrixF temp;
        temp = d2Fdx2;
        d2Fdx2 = temp.selfadjointView<Eigen::Upper>();
        temp = sum_dFdx_d2xdv2;
        sum_dFdx_d2xdv2 = temp.selfadjointView<Eigen::Upper>();

        tbb::concurrent_vector<TripletF> triplets_A;
        tbb::task_group group;
        group.run([&] {
            sum_dFdx_d2xdc2.setFromTriplets(triplets_sum_dFdx_d2xdc2.begin(), triplets_sum_dFdx_d2xdc2.end());
            addTriplets(triplets_A, sum_dFdx_d2xdc2.triangularView<Eigen::Upper>());
        });
        group.run([&] { addTripletsParallel_mul_A_B_C_triangular(triplets_A, dxdcT, d2Fdx2, dxdc); });
        group.run([&] { addTriplets(triplets_A, (dxdcT * d2Fdxdc).triangularView<Eigen::Upper>()); });
        group.run([&] { addTriplets(triplets_A, (d2Fdcdx * dxdc).triangularView<Eigen::Upper>()); });
        group.run([&] { addTriplets(triplets_A, d2Fdc2.triangularView<Eigen::Upper>()); });

        if (model.dimensions_ind->np > 0) {
            group.run([&] { addTripletsBlockOffset(triplets_A, sum_dFdx_d2xdcdv * dvdp, 0, nc); });
            group.run([&] { addTripletsParallel_mul_A_B_C_blockOffset(triplets_A, dxdcT, d2Fdx2, dxdv_dvdp, 0, nc); });
            group.run([&] { addTripletsBlockOffset(triplets_A, d2Fdcdx * dxdv_dvdp, 0, nc); });

            group.run([&] {
                addTripletsParallel_mul_A_B_C_blockOffsetTriangular(triplets_A, dxdv_dvdpT, d2Fdx2, dxdv_dvdp, nc, nc);
            });
            group.run([&] {
                addTripletsParallel_mul_A_B_C_blockOffsetTriangular(triplets_A, dvdpT, sum_dFdx_d2xdv2, dvdp, nc, nc);
            });
            group.run([&] {
                addTripletsBlockOffset(triplets_A, sum_dFdx_dxdv_d2vdp2.triangularView<Eigen::Upper>(), nc, nc);
            });
            group.run([&] { addTripletsBlockOffset(triplets_A, d2Fdp2.A.triangularView<Eigen::Upper>(), nc, nc); });
        }

        group.wait();

        hessian.setZero(nc + np);
        hessian.A.setFromTriplets(triplets_A.begin(), triplets_A.end());

        hessian.U = MatrixXF::Zero(nc + np, d2Fdp2.U.cols());
        hessian.U.bottomRows(np) = d2Fdp2.U;
        hessian.V = MatrixXF::Zero(nc + np, d2Fdp2.V.cols());
        hessian.V.bottomRows(np) = d2Fdp2.V;
    }
}

void PotentialEnergy::computePotentialEnergy(const Model &model, F &energy) {
    VectorXF gradient;
    HessianF hessian;
    computePotentialEnergy(model, energy, gradient, hessian, 0);
}

void PotentialEnergy::computePotentialEnergyGradient(const Model &model, VectorXF &gradient) {
    F energy;
    HessianF hessian;
    computePotentialEnergy(model, energy, gradient, hessian, 1);
}

void PotentialEnergy::computePotentialEnergyHessian(const Model &model, HessianF &hessian) {
    F energy;
    VectorXF gradient;
    computePotentialEnergy(model, energy, gradient, hessian, 2);
}

bool PotentialEnergy::computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition,
                                                          const DegreesOfFreedom &model_dof, const VectorXF &y,
                                                          F &energy) {
    int order = 0;

    DegreesOfFreedom temp_dof = model_dof;
    ModelHelper::setDOFFromVector(model_definition, temp_dof, y);

    bool model_generation_success;
    Model model(model_definition, temp_dof, order, model_generation_success);
    if (!model_generation_success) {
        return false;
    }

    VectorXF gradient;
    HessianF hessian;
    computePotentialEnergy(model, energy, gradient, hessian, order);

    return true;
}

bool PotentialEnergy::computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition,
                                                          const DegreesOfFreedom &model_dof, const VectorXF &y,
                                                          F &energy, VectorXF &gradient) {
    int order = 1;

    DegreesOfFreedom temp_dof = model_dof;
    ModelHelper::setDOFFromVector(model_definition, temp_dof, y);

    bool model_generation_success;
    Model model(model_definition, temp_dof, order, model_generation_success);
    if (!model_generation_success) {
        return false;
    }

    HessianF hessian;
    computePotentialEnergy(model, energy, gradient, hessian, order);

    return true;
}

bool PotentialEnergy::computePotentialEnergyFromDOFVector(const ModelDefinition &model_definition,
                                                          const DegreesOfFreedom &model_dof, const VectorXF &y,
                                                          F &energy, VectorXF &gradient, HessianF &hessian) {
    int order = 2;

    DegreesOfFreedom temp_dof = model_dof;
    ModelHelper::setDOFFromVector(model_definition, temp_dof, y);

    bool model_generation_success;
    Model model(model_definition, temp_dof, order, model_generation_success);
    if (!model_generation_success) {
        return false;
    }

    computePotentialEnergy(model, energy, gradient, hessian, order);

    return true;
}
