#pragma once

#include "Projects/VoronoiFoam/include/Model/Boundary/BoundaryGenerator.h"

class BoundaryGenerator3D : public BoundaryGenerator {
   protected:
    /// Returns false if boundary is self-intersecting, non-closed, non-manifold, or encloses negative volume.
    static bool checkValidBoundary3D(const BoundaryData &boundary_data);

    static bool checkPointInBounds3D(const BoundaryData &boundary_data, const VectorXF &point);

   private:
    [[nodiscard]] bool checkValidBoundary(const DegreesOfFreedom &degrees_of_freedom,
                                          const BoundaryData &boundary_data) const override;

   public:
    [[nodiscard]] int getDims() const override { return 3; }

    [[nodiscard]] bool checkPointInBounds(const BoundaryData &boundary_data, const VectorXF &point) const override;

    [[nodiscard]] F computeEnclosedVolume(const BoundaryData &boundary_data) const override;
};
