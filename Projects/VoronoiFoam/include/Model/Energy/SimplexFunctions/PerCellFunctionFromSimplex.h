#pragma once

#include "CRLHelper/VecMatDef.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"
#include "Projects/VoronoiFoam/include/Model/Energy/PerCellFunction.h"

struct PerSimplexValue {
    F value = 0;
    VectorXF gradient;
    MatrixXF hessian;

    PerSimplexValue(int dims) {
        int n_inputs = dims * dims;
        gradient = VectorXF::Zero(n_inputs);
        hessian = MatrixXF::Zero(n_inputs, n_inputs);
    }

    [[nodiscard]] bool allFinite() const { return std::isfinite(value) && gradient.allFinite() && hessian.allFinite(); }
};

struct PerFaceValue {
    /// VectorXI is node indices in face for the corresponding simplex.
    std::vector<std::pair<VectorXI, PerSimplexValue>> face_simplex_values;
};

/// Many cell functions (e.g. Volume) can be computed as a sum of per-simplex values.
class PerCellFunctionFromSimplex : public PerCellFunction {
   private:
    virtual void getSimplexValue(const VectorXF &inputs, PerSimplexValue &value) const = 0;

    virtual void getSimplexGradient(const VectorXF &inputs, PerSimplexValue &value) const = 0;

    virtual void getSimplexHessian(const VectorXF &inputs, PerSimplexValue &value) const = 0;

   public:
    void getValue(const Model &model, PerCellValue &cell_value) const override;

   private:
    /// Does nothing by default. To be overridden for effects e.g. adhesion which depend on per-face information.
    virtual void postProcessFaceValue(const Model &model, const TessellationCell &cell, const TessellationFace &face,
                                      int order) const {}

    [[nodiscard]] PerFaceValue getFaceValue(const Model &model, const TessellationFace &face, int order) const;

    [[nodiscard]] PerFaceValue getFaceValue2D(const Model &model, const TessellationFace &face, int order) const;

    [[nodiscard]] PerFaceValue getFaceValue3D(const Model &model, const TessellationFace &face, int order) const;

    [[nodiscard]] PerSimplexValue getSimplexValue(const Model &model, const std::vector<int> &simplex_nodes,
                                                  int order) const;
};
