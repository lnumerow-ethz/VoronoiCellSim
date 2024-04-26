#include "CRLHelper/MapleHelper.h"

typedef Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXFRowMajor;

void processMapleOutput(const F *in, F &out, int dim1, int dim2) {
    /// Suppress unused parameter warning.
    (void)(dim1);
    (void)(dim2);
    /// Just copy the numerical value here.
    out = *in;
}

void processMapleOutput(const F *in, VectorXF &out, int dim1, int dim2) {
    /// Suppress unused parameter warning.
    (void)(dim2);
    /// Construct map then copy data on assignment to 'out'.
    out = Eigen::Map<const VectorXF>(in, dim1);
}

void processMapleOutput(const F *in, MatrixXF &out, int dim1, int dim2) {
    /// Converts row major Maple output to column major Eigen matrix.
    out = Eigen::Map<const MatrixXFRowMajor>(in, dim1, dim2);
}

void processMapleOutput(const F *in, std::vector<MatrixXF> &out, int dim1, int dim2) {
    /// Input is assumed to be Hessians (square matrices) concatenated vertically.
    int num_hessians = dim1 / dim2;
    auto stacked_hessians = Eigen::Map<const MatrixXFRowMajor>(in, dim1, dim2);

    out.resize(num_hessians);
    for (int i = 0; i < num_hessians; i++) {
        out[i] = stacked_hessians.middleRows(i * dim2, dim2);
    }
}
