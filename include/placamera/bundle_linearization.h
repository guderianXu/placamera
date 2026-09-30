#pragma once

#include "placamera/frame_camera.h"
#include "placamera/linescan_numeric_state.h"

#include <array>
#include <vector>

namespace placamera
{

    /**
     * Projection and first-order geometry exposed to a bundle solver.
     * Jacobians are row-major with two rows (sample, line). pointJacobian is
     * 2x3 and poseJacobian is 2x6 in rotation-then-translation order.
     * parameterJacobian is 2 x layout.parameterCount(), with columns indexed
     * by the offsets in layout.
     */
    struct BundleProjectionLinearization
    {
        Projection projection;
        std::array<double, 6> pointJacobian{};
        std::array<double, 12> poseJacobian{};
        OptimizationLayout layout;
        std::vector<double> parameterJacobian;
    };

    EvaluationResult<BundleProjectionLinearization>
    linearizeForBundle(const FramePinholeModel& model,
                       const GroundCoordinate& ground,
                       const FrameOptimizationSelection& selection = FrameOptimizationSelection::all());

    EvaluationResult<BundleProjectionLinearization>
    linearizeForBundle(const LineScanModel& model,
                       const GroundCoordinate& ground,
                       double observationLine,
                       const LineScanOptimizationSelection& selection = {});

} // namespace placamera
