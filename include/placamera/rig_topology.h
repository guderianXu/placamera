#pragma once

#include "placamera/camera_topology.h"
#include "placamera/frame_numeric_state.h"

#include <cstddef>
#include <vector>

namespace placamera
{

    /** Validate one complete rig topology against the number of bound frame cameras. */
    Result<void> validateRigTopology(const RigTopology& topology, std::size_t cameraCount);

    /** Apply a local six-degree-of-freedom pose increment to a rigid pose. */
    bool applyPoseDelta(RotationMatrix* rotation, Vector3* center, const std::array<double, 6>& delta) noexcept;

    /** Compose bound frame-camera poses from rig captures and shared sensor extrinsics. */
    Result<std::vector<FramePinholeNumericState>>
    composeRigCameras(const std::vector<FramePinholeNumericState>& cameraModels, const RigTopology& topology);

} // namespace placamera
