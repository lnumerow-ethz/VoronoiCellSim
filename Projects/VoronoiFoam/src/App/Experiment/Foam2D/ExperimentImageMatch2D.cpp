#include "ThirdParty/polyscope/deps/stb/stb_image.h"

#include "Projects/VoronoiFoam/include/App/Experiment/Foam2D/ExperimentImageMatch2D.h"
#include "Projects/VoronoiFoam/include/App/FoamSubApp.h"
#include "Projects/VoronoiFoam/include/App/Scenario/Scenario2D/RandomSitesInBox2D.h"
#include "Projects/VoronoiFoam/include/Model/Energy/Energy2D/CellEnergy2D.h"

#include "Projects/VoronoiFoam/include/Model/Energy/PotentialEnergy.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"

#include "CRLHelper/GeometryHelper.h"
#include "CRLHelper/MapleHelper.h"

#include <Eigen/SparseLU>
#include <Eigen/CholmodSupport>

#include <iostream>

typedef Eigen::CholmodSupernodalLLT<SparseMatrixF, Eigen::Lower> CholmodSolver;
typedef Eigen::SparseLU<SparseMatrixF> SimplicialSolver;

static void computeNodePosStandard(const VectorXF& inputs, NodeData& node_data) {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F w0 = inputs[2];
    F x1 = inputs[3];
    F y1 = inputs[4];
    F w1 = inputs[5];
    F x2 = inputs[6];
    F y2 = inputs[7];
    F w2 = inputs[8];

    F t1 = y1 - y2;
    F t2 = y0 * y0;
    F t4 = x1 * x1;
    F t5 = x2 * x2;
    F t6 = y1 * y1;
    F t7 = y2 * y2;
    F t8 = -t4 + t5 - t6 + t7 + w1 - w2;
    F t11 = x0 * x0;
    F t17 = -x1 + x2;

    F unknown[2];

    unknown[0] = 0.1e1 / (0.2e1 * y0 * t17 + 0.2e1 * y1 * (x0 - x2) - 0.2e1 * y2 * (x0 - x1)) * (t2 * t1 + y0 * t8 + y2 * t6 + y1 * (t11 - t5 - t7 - w0 + w2) + (-t11 + t4 + w0 - w1) * y2);
    unknown[1] = 0.1e1 / (0.2e1 * x0 * t1 + 0.2e1 * x1 * (-y0 + y2) + 0.2e1 * x2 * (-y1 + y0)) * (t11 * t17 - x0 * t8 - x2 * t4 + x1 * (t5 - t2 + t7 + w0 - w2) - (-t2 + t6 + w0 - w1) * x2);

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.pos, 2, 1);
    // clang-format on
}

static void computeNodeGradStandard(const VectorXF& inputs, NodeData& node_data) {
    // clang-format off
    F x0 = inputs[0];
    F y0 = inputs[1];
    F w0 = inputs[2];
    F x1 = inputs[3];
    F y1 = inputs[4];
    F w1 = inputs[5];
    F x2 = inputs[6];
    F y2 = inputs[7];
    F w2 = inputs[8];

    F t1 = x0 * y1;
    F t2 = y2 * x0;
    F t4 = -x1 + x2;
    F t6 = x0 - x2;
    F t8 = x0 - x1;
    F t11 = 0.2e1 * y0 * t4 + 0.2e1 * y1 * t6 - 0.2e1 * t8 * y2;
    F t12 = 0.1e1 / t11;
    F t14 = y1 - y2;
    F t15 = y0 * y0;
    F t17 = x1 * x1;
    F t18 = x2 * x2;
    F t19 = y1 * y1;
    F t20 = y2 * y2;
    F t21 = -t17 + t18 - t19 + t20 + w1 - w2;
    F t24 = x0 * x0;
    F t30 = t11 * t11;
    F t32 = 0.1e1 / t30 * (t15 * t14 + y0 * t21 + y2 * t19 + y1 * (t24 - t18 - t20 - w0 + w2) + (-t24 + t17 + w0 - w1) * y2);
    F t42 = x1 * y0;
    F t43 = y2 * x1;
    F t46 = -y0 + y2;
    F t52 = 0.2e1 * y1 * y2;
    F t58 = x2 * y0;
    F t59 = x2 * y1;
    F t62 = -y1 + y0;
    F t79 = 0.2e1 * x0 * t14 + 0.2e1 * x1 * t46 + 0.2e1 * t62 * x2;
    F t80 = 0.1e1 / t79;
    F t90 = t79 * t79;
    F t92 = 0.1e1 / t90 * (t24 * t4 - x0 * t21 - x2 * t17 + x1 * (t18 - t15 + t20 + w0 - w2) - (-t15 + t19 + w0 - w1) * x2);
    F t103 = 0.2e1 * x1 * x2;

    F unknown[2][9];

    unknown[0][0] = 0.2e1 * t12 * (t1 - t2) - 0.2e1 * t14 * t32;
    unknown[0][1] = t12 * (0.2e1 * y0 * t14 - t17 + t18 - t19 + t20 + w1 - w2) - 0.2e1 * t4 * t32;
    unknown[0][2] = -t12 * t14;
    unknown[0][3] = 0.2e1 * t12 * (-t42 + t43) - 0.2e1 * t46 * t32;
    unknown[0][4] = t12 * (-0.2e1 * y1 * y0 + t15 - t18 - t20 + t24 + t52 - w0 + w2) - 0.2e1 * t6 * t32;
    unknown[0][5] = -t12 * t46;
    unknown[0][6] = 0.2e1 * t12 * (t58 - t59) - 0.2e1 * t62 * t32;
    unknown[0][7] = t12 * (0.2e1 * y2 * y0 - t15 + t17 + t19 - t24 - t52 + w0 - w1) + 0.2e1 * t8 * t32;
    unknown[0][8] = -t12 * t62;
    unknown[1][0] = t80 * (0.2e1 * x0 * t4 + t17 - t18 + t19 - t20 - w1 + w2) - 0.2e1 * t14 * t92;
    unknown[1][1] = 0.2e1 * t80 * (-t42 + t58) - 0.2e1 * t4 * t92;
    unknown[1][2] = -t80 * t4;
    unknown[1][3] = t80 * (0.2e1 * x1 * x0 - t103 - t15 + t18 + t20 - t24 + w0 - w2) - 0.2e1 * t46 * t92;
    unknown[1][4] = 0.2e1 * t80 * (t1 - t59) - 0.2e1 * t6 * t92;
    unknown[1][5] = -t80 * t6;
    unknown[1][6] = t80 * (-0.2e1 * x2 * x0 + t103 + t15 - t17 - t19 + t24 - w0 + w1) - 0.2e1 * t62 * t92;
    unknown[1][7] = 0.2e1 * t80 * (-t2 + t43) + 0.2e1 * t8 * t92;
    unknown[1][8] = t80 * t8;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.grad, 2, 9);
    // clang-format on
}

static void computeNodePosBEdge(const VectorXF& inputs, NodeData& node_data) {
    // clang-format off
    F xb0 = inputs[0];
    F yb0 = inputs[1];
    F xb1 = inputs[2];
    F yb1 = inputs[3];
    F x0 = inputs[4];
    F y0 = inputs[5];
    F w0 = inputs[6];
    F x1 = inputs[7];
    F y1 = inputs[8];
    F w1 = inputs[9];

    F t1 = x0 * x0;
    F t2 = x1 * x1;
    F t3 = y0 * y0;
    F t6 = y1 * y1;
    F t18 = x1 - x0;
    F t21 = -y1 + y0;

    F unknown[2];

    unknown[0] = 0.1e1 / (0.2e1 * xb0 * t18 - 0.2e1 * xb1 * t18 + 0.2e1 * (yb1 - yb0) * t21) * (xb0 * (0.2e1 * y0 * yb1 - 0.2e1 * y1 * yb1 - t1 + t2 - t3 + t6 + w0 - w1) - (0.2e1 * y0 * yb0 - 0.2e1 * y1 * yb0 - t1 + t2 - t3 + t6 + w0 - w1) * xb1);
    unknown[1] = 0.1e1 / (-0.2e1 * yb0 * t21 + 0.2e1 * yb1 * t21 + 0.2e1 * t18 * (-xb1 + xb0)) * (yb0 * (0.2e1 * x0 * xb1 - 0.2e1 * x1 * xb1 - t1 + t2 - t3 + t6 + w0 - w1) - (0.2e1 * x0 * xb0 - 0.2e1 * x1 * xb0 - t1 + t2 - t3 + t6 + w0 - w1) * yb1);

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.pos, 2, 1);
    // clang-format on
}

static void computeNodeGradBEdge(const VectorXF& inputs, NodeData& node_data) {
    // clang-format off
    F xb0 = inputs[0];
    F yb0 = inputs[1];
    F xb1 = inputs[2];
    F yb1 = inputs[3];
    F x0 = inputs[4];
    F y0 = inputs[5];
    F w0 = inputs[6];
    F x1 = inputs[7];
    F y1 = inputs[8];
    F w1 = inputs[9];

    F t1 = x0 * x0;
    F t2 = x1 * x1;
    F t3 = y0 * y0;
    F t4 = y0 * yb1;
    F t6 = y1 * y1;
    F t7 = y1 * yb1;
    F t9 = -t1 + t2 - t3 + 0.2e1 * t4 + t6 - 0.2e1 * t7 + w0 - w1;
    F t10 = x1 - x0;
    F t13 = -y1 + y0;
    F t14 = yb1 - yb0;
    F t17 = 0.2e1 * xb0 * t10 - 0.2e1 * xb1 * t10 + 0.2e1 * t14 * t13;
    F t18 = 0.1e1 / t17;
    F t21 = y0 * yb0;
    F t23 = y1 * yb0;
    F t25 = -t1 + t2 - t3 + 0.2e1 * t21 + t6 - 0.2e1 * t23 + w0 - w1;
    F t28 = t17 * t17;
    F t30 = 0.1e1 / t28 * (-t25 * xb1 + xb0 * t9);
    F t44 = x0 * xb0;
    F t45 = x0 * xb1;
    F t48 = xb1 - xb0;
    F t60 = x1 * xb0;
    F t61 = x1 * xb1;
    F t80 = -0.2e1 * t10 * t48 - 0.2e1 * yb0 * t13 + 0.2e1 * yb1 * t13;
    F t81 = 0.1e1 / t80;
    F t85 = -t1 + 0.2e1 * t45 + t2 - 0.2e1 * t61 - t3 + t6 + w0 - w1;
    F t89 = -t1 + 0.2e1 * t44 + t2 - 0.2e1 * t60 - t3 + t6 + w0 - w1;
    F t92 = t80 * t80;
    F t94 = 0.1e1 / t92 * (yb0 * t85 - t89 * yb1);

    F unknown[2][10];

    unknown[0][0] = -0.2e1 * t10 * t30 + t18 * t9;
    unknown[0][1] = -0.2e1 * t18 * t13 * xb1 + 0.2e1 * t13 * t30;
    unknown[0][2] = 0.2e1 * t10 * t30 - t18 * t25;
    unknown[0][3] = 0.2e1 * t18 * xb0 * t13 - 0.2e1 * t13 * t30;
    unknown[0][4] = 0.2e1 * t18 * (-t44 + t45) - 0.2e1 * t48 * t30;
    unknown[0][5] = t18 * (0.2e1 * xb0 * (-y0 + yb1) - 0.2e1 * (-y0 + yb0) * xb1) - 0.2e1 * t14 * t30;
    unknown[0][6] = -t18 * t48;
    unknown[0][7] = 0.2e1 * t18 * (t60 - t61) + 0.2e1 * t48 * t30;
    unknown[0][8] = t18 * (0.2e1 * xb0 * (y1 - yb1) - 0.2e1 * (y1 - yb0) * xb1) + 0.2e1 * t14 * t30;
    unknown[0][9] = t18 * t48;
    unknown[1][0] = 0.2e1 * t81 * t10 * yb1 - 0.2e1 * t10 * t94;
    unknown[1][1] = 0.2e1 * t13 * t94 + t81 * t85;
    unknown[1][2] = -0.2e1 * t81 * yb0 * t10 + 0.2e1 * t10 * t94;
    unknown[1][3] = -0.2e1 * t13 * t94 - t81 * t89;
    unknown[1][4] = t81 * (0.2e1 * yb0 * (-x0 + xb1) - 0.2e1 * (-x0 + xb0) * yb1) - 0.2e1 * t48 * t94;
    unknown[1][5] = 0.2e1 * t81 * (-t21 + t4) - 0.2e1 * t14 * t94;
    unknown[1][6] = -t81 * t14;
    unknown[1][7] = t81 * (0.2e1 * yb0 * (x1 - xb1) - 0.2e1 * (x1 - xb0) * yb1) + 0.2e1 * t48 * t94;
    unknown[1][8] = 0.2e1 * t81 * (t23 - t7) + 0.2e1 * t14 * t94;
    unknown[1][9] = t81 * t14;

    processMapleOutput(reinterpret_cast<F *>(unknown), node_data.grad, 2, 10);
    // clang-format on
}

void ExperimentImageMatch2D::computeImageMatchObjective(FoamSubApp* foam_sub_app, F& objective) {
    objective = 0;

    VectorXF c;
    ModelHelper::getDOFVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, c);

    VectorXF u(foam_sub_app->degrees_of_freedom.sites.size());
    for (int i = 0; i < u.rows(); i++) {
        u(i) = foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET);
        objective += 1e-10 / pow(u(i), 2.0);
    }

    for (TessellationNode& node : target_nodes) {
        if (node.type == NodeType::STANDARD) {
            VectorXF inputs(9);
            for (int i = 0; i < 3; i++) {
                inputs.segment<3>(i * 3) = c.segment<3>(node.gen[i] * 3);
            }

            NodeData data(2, 9, 0);
            computeNodePosStandard(inputs, data);
            objective += (data.pos - node.data->pos).squaredNorm();
        } else if (node.type == NodeType::B_EDGE) {
            VectorXF inputs(10);
            Vector4F boundary_verts;
            if (node.gen[0] == 0) {
                boundary_verts << -dx, -dy, dx, -dy;
            }
            if (node.gen[0] == 1) {
                boundary_verts << dx, -dy, dx, dy;
            }
            if (node.gen[0] == 2) {
                boundary_verts << dx, dy, -dx, dy;
            }
            if (node.gen[0] == 3) {
                boundary_verts << -dx, dy, -dx, -dy;
            }
            inputs.segment<4>(0) = boundary_verts;
            for (int i = 0; i < 2; i++) {
                inputs.segment<3>(4 + i * 3) = c.segment<3>(node.gen[i + 1] * 3);
            }

            NodeData data(2, 10, 0);
            computeNodePosBEdge(inputs, data);
            objective += (data.pos - node.data->pos).squaredNorm();
        }
    }
}

void ExperimentImageMatch2D::computeImageMatchGradient(FoamSubApp* foam_sub_app, VectorXF& pLpc, VectorXF& pLpu) {
    VectorXF u(foam_sub_app->degrees_of_freedom.sites.size());
    pLpu = 0 * u;
    for (int i = 0; i < u.rows(); i++) {
        u(i) = foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET);
        pLpu(i) += -2.0 * 1e-10 / pow(u(i), 3.0);
    }

    VectorXF c;
    ModelHelper::getDOFVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, c);

    pLpc = 0 * c;

    for (TessellationNode& node : target_nodes) {
        if (node.type == NodeType::STANDARD) {
            VectorXF inputs(9);
            for (int i = 0; i < 3; i++) {
                inputs.segment<3>(i * 3) = c.segment<3>(node.gen[i] * 3);
            }

            NodeData data(2, 9, 1);
            computeNodePosStandard(inputs, data);
            computeNodeGradStandard(inputs, data);

            VectorXF this_grad = 2 * (data.pos - node.data->pos).transpose() * data.grad;
            for (int i = 0; i < 3; i++) {
                pLpc.segment<3>(node.gen[i] * 3) += this_grad.segment<3>(i * 3);
            }
        } else if (node.type == NodeType::B_EDGE) {
            VectorXF inputs(10);
            Vector4F boundary_verts;
            if (node.gen[0] == 0) {
                boundary_verts << -dx, -dy, dx, -dy;
            }
            if (node.gen[0] == 1) {
                boundary_verts << dx, -dy, dx, dy;
            }
            if (node.gen[0] == 2) {
                boundary_verts << dx, dy, -dx, dy;
            }
            if (node.gen[0] == 3) {
                boundary_verts << -dx, dy, -dx, -dy;
            }
            inputs.segment<4>(0) = boundary_verts;
            for (int i = 0; i < 2; i++) {
                inputs.segment<3>(4 + i * 3) = c.segment<3>(node.gen[i + 1] * 3);
            }

            NodeData data(2, 10, 1);
            computeNodePosBEdge(inputs, data);
            computeNodeGradBEdge(inputs, data);

            VectorXF this_grad = 2 * (data.pos - node.data->pos).transpose() * data.grad;
            for (int i = 0; i < 2; i++) {
                pLpc.segment<3>(node.gen[i + 1] * 3) += this_grad.segment<3>(4 + i * 3);
            }
        }
    }
}

void ExperimentImageMatch2D::computeMixedHessian(FoamSubApp* foam_sub_app, SparseMatrixF& d2Edy2,
                                                 SparseMatrixF& d2Edydu) {
    ModelDefinition model_definition = foam_sub_app->model_definition;
    DegreesOfFreedom degrees_of_freedom = foam_sub_app->degrees_of_freedom;
    model_definition.site_free_param_indices = Vector2I(SITE_PARAM_POWER_WEIGHT, SITE_PARAM_SIZE_TARGET);

    VectorXF y;
    ModelHelper::getDOFVector(model_definition, degrees_of_freedom, y);

    F energy;
    VectorXF grad;
    HessianF hessian;
    PotentialEnergy::computePotentialEnergyFromDOFVector(model_definition, degrees_of_freedom, y, energy, grad,
                                                         hessian);

    assert(hessian.U.size() == 0);

    TripletListF triplets_d2Edy2;
    TripletListF triplets_d2Edydu;
    for (int i = 0; i < hessian.A.outerSize(); i++) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(hessian.A, i); it; ++it) {
            int row = it.row();
            int col = it.col();

            int iv0 = row / 4;
            int iv1 = col / 4;
            int ic0 = row % 4;
            int ic1 = col % 4;

            if (ic0 == 3) continue;
            if (ic1 == 3) {
                triplets_d2Edydu.emplace_back(iv0 * 3 + ic0, iv1, it.value());
            } else {
                triplets_d2Edy2.emplace_back(iv0 * 3 + ic0, iv1 * 3 + ic1, it.value());
            }
        }
    }

    int n = foam_sub_app->degrees_of_freedom.sites.size();
    d2Edy2.resize(n * 3, n * 3);
    d2Edydu.resize(n * 3, n);
    d2Edy2.setFromTriplets(triplets_d2Edy2.begin(), triplets_d2Edy2.end());
    d2Edydu.setFromTriplets(triplets_d2Edydu.begin(), triplets_d2Edydu.end());
}

void ExperimentImageMatch2D::loadImageNetwork(FoamSubApp* foam_sub_app, std::string image_file_name,
                                              std::string network_file_name) {
    int width, height, channels;
    stbi_image_free(stbi_load(image_file_name.c_str(), &width, &height, &channels, 0));

    F scale = height * 0.5;  /// This number of pixels corresponds to 1 in system dimensions.
    dx = width * 0.5 / scale;
    dy = height * 0.5 / scale;
    foam_sub_app->degrees_of_freedom.boundary_param = Vector2F(dx, dy);

    std::ifstream network_file(network_file_name);
    std::vector<Vector2F> vertices;
    std::vector<std::vector<int>> cells;

    if (!network_file) {
        std::cerr << "Error opening file." << std::endl;
    }

    std::string line;
    bool readVertices = false;

    while (std::getline(network_file, line)) {
        if (line.empty()) {
            readVertices = false;
            continue;
        }

        if (line.find("Vertices:") != std::string::npos) {
            readVertices = true;
            continue;
        }

        if (line.find("Cells:") != std::string::npos) {
            readVertices = false;
            continue;
        }

        if (readVertices) {
            std::istringstream iss(line);
            char openParen, comma;
            int x, y;
            iss >> openParen >> x >> comma >> y >> std::ws;
            vertices.emplace_back((x - width * 0.5) / scale, -(y - height * 0.5) / scale);
        } else {
            std::istringstream iss(line);
            std::string cellLabel;
            char comma;
            iss >> cellLabel >> cellLabel >> std::ws;
            std::vector<int> cell;
            int vertexIndex;
            while (iss >> vertexIndex) {
                cell.push_back(vertexIndex);
                iss >> comma >> std::ws;
            }
            cells.push_back(cell);
        }
    }

    std::vector<std::vector<int>> vertex_gen(vertices.size());

    for (int i = 0; i < (int)cells.size(); i++) {
        for (int vert_index : cells[i]) {
            vertex_gen[vert_index].push_back(i);
        }
    }
    for (int i = 0; i < (int)vertices.size(); i++) {
        std::sort(vertex_gen[i].begin(), vertex_gen[i].end());

        if (vertex_gen[i].size() < 2) continue;

        bool is_bdry = false;
        if (vertex_gen[i].size() == 2) {
            is_bdry = true;

            F d_right = fabs(vertices[i].x() - dx);
            F d_left = fabs(vertices[i].x() + dx);
            F d_top = fabs(vertices[i].y() - dy);
            F d_bottom = fabs(vertices[i].y() + dy);
            F d_min = std::min({d_right, d_left, d_top, d_bottom});

            int b_edge_index = -1;
            if (d_bottom == d_min) b_edge_index = 0;
            if (d_right == d_min) b_edge_index = 1;
            if (d_top == d_min) b_edge_index = 2;
            if (d_left == d_min) b_edge_index = 3;

            assert(b_edge_index >= 0);
            vertex_gen[i].insert(vertex_gen[i].begin(), b_edge_index);
        } else {
            assert(vertex_gen[i].size() == 3);
        }

        target_nodes.emplace_back();
        TessellationNode& node = target_nodes.back();
        node.type = is_bdry ? NodeType::B_EDGE : NodeType::STANDARD;
        node.gen = Eigen::Map<Vector3I>(vertex_gen[i].data());

        NodeData data(2, 0, 0);
        data.pos = vertices[i];
        node.data = std::make_shared<NodeData>(data);
    }

    // Print cells
    for (const std::vector<int>& cell : cells) {
        foam_sub_app->degrees_of_freedom.sites.emplace_back();
        Site& site = foam_sub_app->degrees_of_freedom.sites.back();

        F area = 0;
        F wmx = 0;
        F wmy = 0;
        for (int i = 0; i < (int)cell.size(); i++) {
            Vector2F v0 = vertices[cell[i]];
            Vector2F v1 = vertices[cell[(i + 1) % cell.size()]];
            area += GeometryHelper::triangleAreaWithOrigin2D(v1, v0);  // Reverse order because files are saved CW.

            F x0 = v0(0);
            F y0 = v0(1);
            F x1 = v1(0);
            F y1 = v1(1);
            wmx -= 0.1666666667e0 * (x0 * y1 - 0.1e1 * x1 * y0) * (x0 + x1);
            wmy -= 0.1666666667e0 * (x0 * y1 - 0.1e1 * x1 * y0) * (y0 + y1);
        }

        site.param(SITE_PARAM_SIZE_TARGET) = std::max(area, 0.005);
        site.pos = Vector2F(wmx, wmy) / area;
    }

    bfgs_inverse_hessian = MatrixXF::Identity(cells.size(), cells.size());
}

void ExperimentImageMatch2D::setup(FoamSubApp* foam_sub_app) {
    foam_sub_app->tessellation_selector.first = 1;  /// Power diagram

    foam_sub_app->scenario_selector.first = 0;  /// Sites in Box
    auto scenario = std::dynamic_pointer_cast<RandomSitesInBox2D>(
        foam_sub_app->scenario_selector.second[foam_sub_app->scenario_selector.first]);
    scenario->num_sites = 0;  /// Don't generate random sites

    foam_sub_app->energy_selector.first = 0;
    auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
        foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
    energy->weights[AREA] = 1.0;
    energy->weights[PERIMETER_TARGET] = 0;
    energy->weights[PERIMETER_MINIMIZATION] = 0.0001;
    energy->weights[CENTROID] = 0.01;
    energy->weights[SECOND_MOMENT] = 0;

    std::fill(std::begin(foam_sub_app->site_param_free), std::end(foam_sub_app->site_param_free), false);
    foam_sub_app->site_param_free[SITE_PARAM_POWER_WEIGHT] = true;

    foam_sub_app->use_dynamics = true;
    foam_sub_app->applyAllSettings();

    loadImageNetwork(foam_sub_app, "resource/ronan.jpg", "resource/ronan.txt");
    foam_sub_app->clearState();
}

void ExperimentImageMatch2D::loop(FoamSubApp* foam_sub_app, Optimization::OptimizationStatus optimization_status) {
    if (optimization_status == Optimization::CONVERGED) {
        /// Outer optimization loop
        if (foam_sub_app->use_dynamics) {
            num_pre_steps++;

            if (num_pre_steps > 20) {
                auto energy = std::dynamic_pointer_cast<CellEnergy2D>(
                    foam_sub_app->energy_selector.second[foam_sub_app->energy_selector.first].second);
                energy->weights[PERIMETER_MINIMIZATION] = 0.005;
                energy->weights[CENTROID] = 0.0001;

                iter_start_time = std::chrono::high_resolution_clock::now();
                foam_sub_app->use_dynamics = false;
            }

            return;
        }

        if (frame == 0) {
            foam_sub_app->write_png = true;
            foam_sub_app->png_file_name = "image_match_frame_" + std::to_string(frame) + ".png";
            foam_sub_app->write_rendering_file = true;
            foam_sub_app->render_file_name = "render_image_match_frame_" + std::to_string(frame) + ".txt";

            auto iter_end_time = std::chrono::high_resolution_clock::now();
            auto iter_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);

            num_frame_iters++;

            std::string frame_info_file_name = "image_match_frame_info_" + std::to_string(frame) + ".txt";
            std::ofstream info_file(frame_info_file_name);

            info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << " " << L_prev << "\n";
            info_file.close();

            frame++;
        }

        if (line_search_in_progress) {
            F L;
            computeImageMatchObjective(foam_sub_app, L);
            std::cout << "Line search: " << L << " " << L_prev << std::endl;

            if (L < L_prev) {
                foam_sub_app->write_png = true;
                foam_sub_app->png_file_name = "image_match_frame_" + std::to_string(frame) + ".png";
                foam_sub_app->write_rendering_file = true;
                foam_sub_app->render_file_name = "render_image_match_frame_" + std::to_string(frame) + ".txt";

                auto iter_end_time = std::chrono::high_resolution_clock::now();
                auto iter_time_ms =
                    std::chrono::duration_cast<std::chrono::milliseconds>(iter_end_time - iter_start_time);

                num_frame_iters++;

                std::string frame_info_file_name = "image_match_frame_info_" + std::to_string(frame) + ".txt";
                std::ofstream info_file(frame_info_file_name);

                info_file << frame << " " << iter_time_ms.count() << " " << num_frame_iters << " " << L_prev << "\n";
                info_file.close();

                frame++;

                L_prev = L;
                line_search_in_progress = false;

                VectorXF u_new(foam_sub_app->degrees_of_freedom.sites.size());
                for (int i = 0; i < u_new.rows(); i++) {
                    u_new(i) = foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET);
                }

                SparseMatrixF d2Edy2, d2Edydu;
                computeMixedHessian(foam_sub_app, d2Edy2, d2Edydu);

                computeImageMatchObjective(foam_sub_app, L_prev);
                VectorXF dLdy, pLpu;
                computeImageMatchGradient(foam_sub_app, dLdy, pLpu);

                SimplicialSolver solver;
                solver.compute(d2Edy2);

                MatrixXF dydu = solver.solve(-d2Edydu.toDense());

                VectorXF dLdu = dLdy.transpose() * dydu + pLpu.transpose();
                VectorXF g_new = dLdu;
                std::cout << "Gradient norm: " << g_new.norm() << std::endl;

                VectorXF deltaY = u_new - u_prev;
                VectorXF deltaG = g_new - g_prev;

                if (deltaG.norm() > 1e-10 && deltaY.norm() > 1e-10) {
                    MatrixXF I = MatrixXF::Identity(u_new.rows(), u_new.rows());
                    MatrixXF A = (I - (deltaY * deltaG.transpose()) / (deltaG.transpose() * deltaY));
                    MatrixXF B = (I - (deltaG * deltaY.transpose()) / (deltaG.transpose() * deltaY));
                    MatrixXF C = (deltaY * deltaY.transpose()) / (deltaG.transpose() * deltaY);
                    bfgs_inverse_hessian = A * bfgs_inverse_hessian * B + C;
                }
            } else {
                alpha *= 0.5;
            }

            if (alpha < 1e-15) {
                line_search_in_progress = false;
                // bfgs_inverse_hessian = MatrixXF::Identity(u_prev.rows(), u_prev.rows());
            }
        }

        if (!line_search_in_progress) {
            ModelHelper::getDOFVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, y_prev);
            u_prev.resize(foam_sub_app->degrees_of_freedom.sites.size());
            for (int i = 0; i < u_prev.rows(); i++) {
                u_prev(i) = foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET);
            }

            SparseMatrixF d2Edy2, d2Edydu;
            computeMixedHessian(foam_sub_app, d2Edy2, d2Edydu);

            F L_curr;
            computeImageMatchObjective(foam_sub_app, L_curr);
            L_prev = std::min(L_prev, L_curr);
            VectorXF dLdy, pLpu;
            computeImageMatchGradient(foam_sub_app, dLdy, pLpu);

            SimplicialSolver solver;
            solver.compute(d2Edy2);

            MatrixXF dydu = solver.solve(-d2Edydu.toDense());

            VectorXF dLdu = dLdy.transpose() * dydu + pLpu.transpose();
            g_prev = dLdu;

            search_direction = -bfgs_inverse_hessian * dLdu;
            if (search_direction.dot(dLdu) >= 0) {
                std::cout << "Invalid search direction in inverse problem." << std::endl;
            }
            search_direction = search_direction.normalized() * 1e-2;

            alpha = 1.0;
            while ((u_prev + alpha * search_direction).minCoeff() < 0) {
                alpha *= 0.5;
            }

            line_search_in_progress = true;
        }

        ModelHelper::setDOFFromVector(foam_sub_app->model_definition, foam_sub_app->degrees_of_freedom, y_prev);
        VectorXF u_iter = u_prev + alpha * search_direction;
        for (int i = 0; i < u_prev.rows(); i++) {
            foam_sub_app->degrees_of_freedom.sites[i].param(SITE_PARAM_SIZE_TARGET) = u_iter(i);
        }
        foam_sub_app->optimize = true;  // Turn optimization back on.
    }
}
