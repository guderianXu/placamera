#include "placamera/types.h"

#include "../internal/plamatrix_rotation.h"

#include <cmath>
#include <utility>

namespace placamera
{

    Pose Pose::create(FrameId frame, Vector3 center, RotationMatrix cameraToWorldRotation)
    {
        for (const double value : center)
        {
            if (!std::isfinite(value))
            {
                throw CameraValidationError(CameraErrorCode::InvalidPose, "pose center must contain finite values");
            }
        }
        for (const double value : cameraToWorldRotation)
        {
            if (!std::isfinite(value))
            {
                throw CameraValidationError(CameraErrorCode::InvalidPose, "pose rotation must contain finite values");
            }
        }

        if (!internal::validRotation(cameraToWorldRotation))
        {
            throw CameraValidationError(CameraErrorCode::InvalidPose,
                                        "pose rotation must be a proper orthonormal matrix");
        }

        return Pose{std::move(frame), center, cameraToWorldRotation};
    }

} // namespace placamera
