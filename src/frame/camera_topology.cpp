#include "placamera/camera_topology.h"

#include "../internal/plamatrix_rotation.h"

#include <cmath>

namespace placamera
{
    namespace
    {

        bool finiteVector(const Vector3& value) noexcept
        {
            return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
        }

    } // namespace

    bool RollingShutterMotion::isIdentity() const noexcept
    {
        return translation == Vector3{} && rotationVector == Vector3{};
    }

    Result<void> validateSensorMount(const SensorMountState& mount, const CameraDefinitionId* definitionId)
    {
        if (!finiteVector(mount.translation) || !internal::validRotation(mount.rotation))
        {
            return Result<void>::failure(CameraErrorCode::InvalidPose,
                                         "sensor mount requires a finite translation and proper rotation");
        }
        if (definitionId && mount.masterSensorId && *definitionId == *mount.masterSensorId)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "sensor mount cannot name its own definition as master");
        }
        return Result<void>::success();
    }

    Result<void> validateCameraAcquisition(const CameraAcquisitionState& acquisition,
                                           const CameraInstanceId* instanceId)
    {
        if (acquisition.role != CameraRole::Regular && acquisition.role != CameraRole::Keyframe)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument, "camera acquisition role is invalid");
        }
        if (acquisition.rollingShutterMode != RollingShutterMode::Disabled &&
            acquisition.rollingShutterMode != RollingShutterMode::Regularized &&
            acquisition.rollingShutterMode != RollingShutterMode::Full)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument, "rolling-shutter mode is invalid");
        }
        if (!finiteVector(acquisition.rollingShutter.translation) ||
            !finiteVector(acquisition.rollingShutter.rotationVector))
        {
            return Result<void>::failure(CameraErrorCode::InvalidPose,
                                         "rolling-shutter motion must contain finite values");
        }
        if (instanceId && acquisition.masterCameraId && *instanceId == *acquisition.masterCameraId)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "camera acquisition cannot name its own instance as master");
        }
        if (acquisition.rollingShutterMode == RollingShutterMode::Disabled && !acquisition.rollingShutter.isIdentity())
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "disabled rolling shutter cannot carry non-zero motion");
        }
        if (acquisition.rollingShutterMode == RollingShutterMode::Regularized &&
            (acquisition.rollingShutter.translation[2] != 0.0 ||
             acquisition.rollingShutter.rotationVector != Vector3{}))
        {
            return Result<void>::failure(
                CameraErrorCode::InvalidArgument,
                "regularized rolling shutter only permits x/y translation; use Full for z or rotation");
        }
        return Result<void>::success();
    }

} // namespace placamera
