#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"

struct PerCellValue {
    TessellationCell cell;

    F value = 0;
    VectorXF gradient;
    HessianF hessian;

    int order;

    PerCellValue(const Model &model, const TessellationCell &cell, int order) : cell(cell), order(order) {
        int dims_x = model.dimensions_ind->dims_space;
        int dims_c = model.dimensions_ind->dims_site_dof;

        int n_vars = (int)cell.node_indices_in_cell.size() * dims_x + dims_c;

        if (order >= 1) gradient = VectorXF::Zero(n_vars);
        if (order >= 2) {
            hessian.setZero(n_vars);
        }
    }

    PerCellValue operator+(const PerCellValue &other) const;

    PerCellValue operator-(const PerCellValue &other) const;

    PerCellValue operator*(const PerCellValue &other) const;

    PerCellValue operator/(const PerCellValue &other) const;

    PerCellValue operator+(const F &other) const;

    PerCellValue operator-(const F &other) const;

    PerCellValue operator*(const F &other) const;

    PerCellValue operator/(const F &other) const;

    friend PerCellValue operator/(const F &a, const PerCellValue &b);

    [[nodiscard]] PerCellValue pow(const F &other) const;

    [[nodiscard]] PerCellValue square() const;

    static PerCellValue sitePosCellValue(const Model &model, const TessellationCell &cell, int order,
                                         int site_coord_index);

    static PerCellValue siteParamCellValue(const Model &model, const TessellationCell &cell, int order,
                                           int site_param_id);
};

class PerCellFunction {
   public:
    /// Computes derivatives up to and including the specified order.
    virtual void getValue(const Model &model, PerCellValue &cell_value) const = 0;

    /// Wrapper for virtual getValue function which allows calling it without first declaring cell_value variable.
    [[nodiscard]] PerCellValue getValueWrapper(const Model &model, const TessellationCell &cell, int order) const;

    /// Allow derived classes to provide a config menu. Seems like the cleanest way to do this despite needing ImGui
    /// outside of the App code directory.
    virtual void makeConfigMenu() {}
};
