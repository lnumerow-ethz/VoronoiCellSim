#include <igl/readOBJ.h>

#include "Projects/VoronoiFoam/include/Model/Boundary/Boundary3D/BoundaryGeneratorOpenCylinder3D.h"
#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryHelper.h"

BoundaryGeneratorOpenCylinder3D::BoundaryGeneratorOpenCylinder3D(F spring_constant) : spring_constant(spring_constant) {
    MatrixXF TC, N, V2;
    MatrixXI FTC, FN, F2;
    igl::readOBJ("resource/cylinder.obj", original_points, TC, N, tri, FTC, FN);

    z_base = 1e10;
    for (int i = 0; i < original_points.rows(); i++) {
        z_base = std::min(z_base, original_points(i, 2));
    }

    for (int i = 0; i < original_points.rows(); i++) {
        if (original_points(i, 2) == z_base) continue;

        F r = original_points.row(i).segment<2>(0).norm();
        original_points(i, 2) = z_base + 0.001 + 0.35 / exp(10 * r * r);
    }

    for (int i = 0; i < tri.rows(); i++) {
        if (original_points(tri(i, 0), 2) == z_base || original_points(tri(i, 1), 2) == z_base ||
            original_points(tri(i, 2), 2) == z_base) {
            for (int j = 0; j < 3; j++) {
                if (original_points(tri(i, j), 2) != z_base) {
                    top_edge_points.insert(tri(i, j));
                }
            }
        }
    }
    for (int i = 0; i < original_points.rows(); i++) {
        if (original_points(i, 2) != z_base && top_edge_points.find(i) == top_edge_points.end()) {
            top_interior_points.insert(i);
        }
    }

    idx = VectorXI::Constant(3 * original_points.rows(), -1);

    int ip = 0;
    for (int i : top_edge_points) {
        idx(i * 3 + 2) = ip;
        ip++;
    }
    for (int i : top_interior_points) {
        idx(i * 3 + 0) = ip + 0;
        idx(i * 3 + 1) = ip + 1;
        idx(i * 3 + 2) = ip + 2;
        ip += 3;
    }
}

void BoundaryGeneratorOpenCylinder3D::getDefaultParams(VectorXF &params) {
    int n_params = top_edge_points.size() + top_interior_points.size() * 3;
    params.resize(n_params);
    for (int i = 0; i < original_points.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (idx(i * 3 + j) >= 0) {
                params(idx(i * 3 + j)) = original_points(i, j);
            }
        }
    }
}

void BoundaryGeneratorOpenCylinder3D::computeBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                                      BoundaryData &boundary_data, int order) const {
    int n_points = original_points.rows();
    int dims_space = getDims();
    int n_params = top_edge_points.size() + top_interior_points.size() * 3;

    boundary_data.v.resize(n_points);
    for (int iv = 0; iv < n_points; iv++) {
        boundary_data.v[iv].pos = original_points.row(iv);
        if (order >= 1) {
            boundary_data.v[iv].grad.resize(dims_space, n_params);
        }
        if (order >= 2) {
            // Hessians are zero.
            boundary_data.v[iv].hess.resize(dims_space);
            for (int id = 0; id < dims_space; id++) {
                boundary_data.v[iv].hess[id].resize(n_params, n_params);
            }
        }
    }

    for (int i = 0; i < original_points.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (idx(i * 3 + j) >= 0) {
                boundary_data.v[i].pos(j) = degrees_of_freedom.boundary_param(idx(i * 3 + j));
                if (order >= 1) boundary_data.v[i].grad.coeffRef(j, idx(i * 3 + j)) = 1.0;
            }
        }
    }

    boundary_data.f.resize(tri.rows());
    for (int ii = 0; ii < tri.rows(); ii++) {
        boundary_data.f[ii].vert_indices = tri.row(ii);
    }
}

bool BoundaryGeneratorOpenCylinder3D::checkValidNumBoundaryParams(int n_params) const {
    return n_params == (int)top_edge_points.size() + (int)top_interior_points.size() * 3;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;

F BoundaryGeneratorOpenCylinder3D::computeEnergy(const Model &model) const {
    F energy = 0;

    for (int i = 0; i < tri.rows(); i++) {
        if (top_edge_points.find(tri(i, 0)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 0)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 1)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 1)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 2)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 2)) == top_interior_points.end())
            continue;

        Vector3F v0 = model.boundary_data.v[tri(i, 0)].pos;
        Vector3F v1 = model.boundary_data.v[tri(i, 1)].pos;
        Vector3F v2 = model.boundary_data.v[tri(i, 2)].pos;
        energy += spring_constant * (v1 - v0).squaredNorm();
        energy += spring_constant * (v2 - v1).squaredNorm();
        energy += spring_constant * (v0 - v2).squaredNorm();
    }

    for (int i = 0; i < (int)model.boundary_data.v.size(); i++) {
        if (top_edge_points.find(i) == top_edge_points.end() &&
            top_interior_points.find(i) == top_interior_points.end())
            continue;

        int ip = idx(i * 3 + 2);
        if (ip >= 0) {
            energy += barrier_weight / pow((model.degrees_of_freedom.boundary_param(ip) - z_base), 2.0);
        }
    }

    return energy;
}

void BoundaryGeneratorOpenCylinder3D::computeEnergyGradient(const Model &model, VectorXF &gradient) const {
    gradient = VectorXF::Zero(model.dimensions_ind->np);

    for (int i = 0; i < tri.rows(); i++) {
        if (top_edge_points.find(tri(i, 0)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 0)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 1)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 1)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 2)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 2)) == top_interior_points.end())
            continue;

        Vector3F v0 = model.boundary_data.v[tri(i, 0)].pos;
        Vector3F v1 = model.boundary_data.v[tri(i, 1)].pos;
        Vector3F v2 = model.boundary_data.v[tri(i, 2)].pos;
        Vector3F diff01 = 2 * spring_constant * (v1 - v0);
        Vector3F diff12 = 2 * spring_constant * (v2 - v1);
        Vector3F diff20 = 2 * spring_constant * (v0 - v2);

        int ip0 = 3 * tri(i, 0);
        int ip1 = 3 * tri(i, 1);
        int ip2 = 3 * tri(i, 2);
        for (int coord = 0; coord < 3; coord++) {
            if (idx(ip0 + coord) >= 0) gradient(idx(ip0 + coord)) += diff20(coord) - diff01(coord);
            if (idx(ip1 + coord) >= 0) gradient(idx(ip1 + coord)) += diff01(coord) - diff12(coord);
            if (idx(ip2 + coord) >= 0) gradient(idx(ip2 + coord)) += diff12(coord) - diff20(coord);
        }
    }

    for (int i = 0; i < (int)model.boundary_data.v.size(); i++) {
        if (top_edge_points.find(i) == top_edge_points.end() &&
            top_interior_points.find(i) == top_interior_points.end())
            continue;

        int ip = idx(i * 3 + 2);
        if (ip >= 0) {
            gradient(ip) -= barrier_weight * 2.0 / pow((model.degrees_of_freedom.boundary_param(ip) - z_base), 3.0);
        }
    }
}

void BoundaryGeneratorOpenCylinder3D::computeEnergyHessian(const Model &model, HessianF &hessian) const {
    TripletListF hessian_triplets;
    for (int i = 0; i < tri.rows(); i++) {
        if (top_edge_points.find(tri(i, 0)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 0)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 1)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 1)) == top_interior_points.end())
            continue;
        if (top_edge_points.find(tri(i, 2)) == top_edge_points.end() &&
            top_interior_points.find(tri(i, 2)) == top_interior_points.end())
            continue;

        F value = 2 * spring_constant;
        int ip0 = 3 * tri(i, 0);
        int ip1 = 3 * tri(i, 1);
        int ip2 = 3 * tri(i, 2);
        for (int coord = 0; coord < 3; coord++) {
            if (idx(ip0 + coord) >= 0) {
                hessian_triplets.emplace_back(idx(ip0 + coord), idx(ip0 + coord), 2 * value);
                if (idx(ip1 + coord) >= 0) hessian_triplets.emplace_back(idx(ip0 + coord), idx(ip1 + coord), -value);
                if (idx(ip2 + coord) >= 0) hessian_triplets.emplace_back(idx(ip0 + coord), idx(ip2 + coord), -value);
            }
            if (idx(ip1 + coord) >= 0) {
                if (idx(ip0 + coord) >= 0) hessian_triplets.emplace_back(idx(ip1 + coord), idx(ip0 + coord), -value);
                hessian_triplets.emplace_back(idx(ip1 + coord), idx(ip1 + coord), 2 * value);
                if (idx(ip2 + coord) >= 0) hessian_triplets.emplace_back(idx(ip1 + coord), idx(ip2 + coord), -value);
            }
            if (idx(ip2 + coord) >= 0) {
                if (idx(ip0 + coord) >= 0) hessian_triplets.emplace_back(idx(ip2 + coord), idx(ip0 + coord), -value);
                if (idx(ip1 + coord) >= 0) hessian_triplets.emplace_back(idx(ip2 + coord), idx(ip1 + coord), -value);
                hessian_triplets.emplace_back(idx(ip2 + coord), idx(ip2 + coord), 2 * value);
            }
        }
    }

    for (int i = 0; i < (int)model.boundary_data.v.size(); i++) {
        if (top_edge_points.find(i) == top_edge_points.end() &&
            top_interior_points.find(i) == top_interior_points.end())
            continue;

        int ip = idx(i * 3 + 2);
        if (ip >= 0) {
            hessian_triplets.emplace_back(
                ip, ip, barrier_weight * 6.0 / pow((model.degrees_of_freedom.boundary_param(ip) - z_base), 4.0));
        }
    }

    hessian.setZero(model.dimensions_ind->np);
    hessian.A.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
}
