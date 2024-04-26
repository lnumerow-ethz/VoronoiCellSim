#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"
#include "Projects/VoronoiFoam/include/Model/Model.h"

class BoundaryGenerator {
   public:
    [[nodiscard]] virtual F computeEnergy(const Model &model) const;

    virtual void computeEnergyGradient(const Model &model, VectorXF &gradient) const;

    virtual void computeEnergyHessian(const Model &model, HessianF &hessian) const;

    /// Wrapper for virtual function computeBoundary which also performs validity checks.
    /// Returns false if boundary is invalid. Computes boundary mesh and derivatives to specified order (max 2).
    [[nodiscard]] bool generateBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                        int order) const;

   private:
    /// Computes boundary mesh and vertex derivatives with respect to all params (including fixed params).
    virtual void computeBoundary(const DegreesOfFreedom &degrees_of_freedom, BoundaryData &boundary_data,
                                 int order) const = 0;

    /// Returns false if boundary is invalid.
    [[nodiscard]] virtual bool checkValidBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                                  const BoundaryData &boundary_data) const = 0;

    /// Checks if number of boundary params is valid.
    [[nodiscard]] virtual bool checkValidNumBoundaryParams(int n_params) const = 0;

   public:
    [[nodiscard]] virtual int getDims() const = 0;

    [[nodiscard]] virtual bool checkPointInBounds(const BoundaryData &boundary_data, const VectorXF &point) const = 0;

    [[nodiscard]] virtual F computeEnclosedVolume(const BoundaryData &boundary_data) const = 0;
};
