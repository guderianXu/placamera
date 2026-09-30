#pragma once

#include "placamera/result.h"
#include "placamera/types.h"

#include <cstddef>
#include <optional>
#include <array>
#include <string>
#include <vector>

namespace placamera
{

    enum class CameraRole
    {
        Regular,
        Keyframe,
    };

    enum class RollingShutterMode
    {
        Disabled,
        Regularized,
        Full,
    };

    struct RollingShutterMotion
    {
        Vector3 translation{};
        Vector3 rotationVector{};

        bool isIdentity() const noexcept;
    };

    /** Sensor-to-master mount state, independent of the sensor projection family. */
    struct SensorMountState
    {
        std::optional<CameraDefinitionId> masterSensorId;
        Vector3 translation{};
        RotationMatrix rotation{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
        bool fixedTranslation = true;
        bool fixedRotation = true;
    };

    /** Per-image acquisition topology, independent of the sensor projection family. */
    struct CameraAcquisitionState
    {
        CameraRole role = CameraRole::Regular;
        std::optional<CaptureGroupId> captureGroupId;
        std::optional<CameraInstanceId> masterCameraId;
        std::size_t layerIndex = 0;
        RollingShutterMode rollingShutterMode = RollingShutterMode::Disabled;
        RollingShutterMotion rollingShutter;
        bool rollingShutterInitialized = false;
    };

    /** One rig capture pose shared by the sensors acquired at the same time. */
    struct RigCapture
    {
        int rigId = -1;
        int captureId = -1;
        std::array<double, 9> rigToWorldRotation{{1.0, 0.0, 0.0,
                                                    0.0, 1.0, 0.0,
                                                    0.0, 0.0, 1.0}};
        std::array<double, 3> rigCenterInWorld{{0.0, 0.0, 0.0}};
        bool fixedPose = false;
    };

    /** One sensor-to-rig extrinsic state shared by all captures of a rig. */
    struct RigSensor
    {
        int rigId = -1;
        int sensorId = -1;
        std::array<double, 9> cameraToRigRotation{{1.0, 0.0, 0.0,
                                                     0.0, 1.0, 0.0,
                                                     0.0, 0.0, 1.0}};
        std::array<double, 3> cameraCenterInRig{{0.0, 0.0, 0.0}};
        bool fixedExtrinsic = true;
    };

    /** Binds one raster instance to a capture and a sensor in a rig. */
    struct RigCameraBinding
    {
        int cameraIndex = -1;
        int rigId = -1;
        int captureId = -1;
        int sensorId = -1;
    };

    /** Complete capture/sensor topology consumed by bundle adjustment. */
    struct RigTopology
    {
        std::vector<RigCapture> captures;
        std::vector<RigSensor> sensors;
        std::vector<RigCameraBinding> cameraBindings;

        bool empty() const noexcept
        {
            return captures.empty() && sensors.empty() && cameraBindings.empty();
        }
    };

    Result<void> validateSensorMount(const SensorMountState& mount, const CameraDefinitionId* definitionId = nullptr);
    Result<void> validateCameraAcquisition(const CameraAcquisitionState& acquisition,
                                           const CameraInstanceId* instanceId = nullptr);

} // namespace placamera
