#include "Projects/VoronoiFoam/include/Model/Energy/PerCellFunction.h"
#include "Projects/VoronoiFoam/include/Model/Tessellation/TessellationGenerator.h"
#include "Projects/VoronoiFoam/include/Model/ModelHelper.h"
#include "CRLHelper/GradArithmetic.h"

PerCellValue PerCellFunction::getValueWrapper(const Model &model, const TessellationCell &cell, int order) const {
    PerCellValue value(model, cell, order);
    getValue(model, value);
    return value;
}

PerCellValue PerCellValue::sitePosCellValue(const Model &model, const TessellationCell &cell, int order,
                                            int site_coord_index) {
    PerCellValue value(model, cell, order);
    int dims_site_dof = model.dimensions_ind->dims_site_dof;

    value.value = model.degrees_of_freedom.sites[cell.site_index].pos(site_coord_index);

    if (order >= 1) value.gradient.tail(dims_site_dof)(site_coord_index) = 1;
    return value;
}

PerCellValue PerCellValue::siteParamCellValue(const Model &model, const TessellationCell &cell, int order,
                                              int site_param_id) {
    PerCellValue value(model, cell, order);
    int dims_site_dof = model.dimensions_ind->dims_site_dof;

    value.value = model.degrees_of_freedom.sites[cell.site_index].param(site_param_id);

    int site_dof_index;
    if (order >= 1 && ModelHelper::siteParamIdToSiteDOFIndex(model, site_param_id, site_dof_index)) {
        value.gradient.tail(dims_site_dof)(site_dof_index) = 1;
    }

    return value;
}

PerCellValue PerCellValue::operator+(const PerCellValue &other) const {
    return GradArithmetic::add(*this, other, order);
}

PerCellValue PerCellValue::operator-(const PerCellValue &other) const {
    return GradArithmetic::subtract(*this, other, order);
}

PerCellValue PerCellValue::operator*(const PerCellValue &other) const {
    return GradArithmetic::multiply(*this, other, order);
}

PerCellValue PerCellValue::operator/(const PerCellValue &other) const {
    return GradArithmetic::divide(*this, other, order);
}

PerCellValue PerCellValue::operator+(const F &other) const {
    return GradArithmetic::add(*this, other);
}

PerCellValue PerCellValue::operator-(const F &other) const {
    return GradArithmetic::subtract(*this, other);
}

PerCellValue PerCellValue::operator*(const F &other) const {
    return GradArithmetic::multiply(*this, other, order);
}

PerCellValue PerCellValue::operator/(const F &other) const {
    return GradArithmetic::divide(*this, other, order);
}

PerCellValue PerCellValue::pow(const F &other) const {
    return GradArithmetic::powf(*this, other, order);
}

PerCellValue PerCellValue::square() const {
    return GradArithmetic::square(*this, order);
}

PerCellValue operator/(const F &a, const PerCellValue &b) {
    return GradArithmetic::divide(a, b, b.order);
}
