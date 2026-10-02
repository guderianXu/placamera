#include "placamera/rig_topology.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace placamera
{
    namespace
    {
        using Key = std::pair<int, int>;

        bool finite(const Vector3& value) noexcept
        {
            return std::all_of(value.begin(), value.end(), [](double component) { return std::isfinite(component); });
        }

        bool validRotation(const RotationMatrix& rotation) noexcept
        {
            if (!std::all_of(rotation.begin(), rotation.end(), [](double value) { return std::isfinite(value); }))
            {
                return false;
            }
            const auto dotRow = [&rotation](int first, int second)
            {
                return rotation[static_cast<std::size_t>(first)] * rotation[static_cast<std::size_t>(second)] +
                       rotation[static_cast<std::size_t>(first + 1)] * rotation[static_cast<std::size_t>(second + 1)] +
                       rotation[static_cast<std::size_t>(first + 2)] * rotation[static_cast<std::size_t>(second + 2)];
            };
            const double determinant = rotation[0] * (rotation[4] * rotation[8] - rotation[5] * rotation[7]) -
                                       rotation[1] * (rotation[3] * rotation[8] - rotation[5] * rotation[6]) +
                                       rotation[2] * (rotation[3] * rotation[7] - rotation[4] * rotation[6]);
            constexpr double tolerance = 1.0e-8;
            return std::abs(dotRow(0, 0) - 1.0) <= tolerance && std::abs(dotRow(3, 3) - 1.0) <= tolerance &&
                   std::abs(dotRow(6, 6) - 1.0) <= tolerance && std::abs(dotRow(0, 3)) <= tolerance &&
                   std::abs(dotRow(0, 6)) <= tolerance && std::abs(dotRow(3, 6)) <= tolerance &&
                   std::abs(determinant - 1.0) <= tolerance;
        }

        RotationMatrix multiply(const RotationMatrix& left, const RotationMatrix& right) noexcept
        {
            RotationMatrix product{};
            for (int row = 0; row < 3; ++row)
            {
                for (int column = 0; column < 3; ++column)
                {
                    for (int inner = 0; inner < 3; ++inner)
                    {
                        product[static_cast<std::size_t>(row * 3 + column)] +=
                            left[static_cast<std::size_t>(row * 3 + inner)] *
                            right[static_cast<std::size_t>(inner * 3 + column)];
                    }
                }
            }
            return product;
        }

        Vector3 transform(const RotationMatrix& rotation, const Vector3& point) noexcept
        {
            Vector3 result{};
            for (int row = 0; row < 3; ++row)
            {
                for (int column = 0; column < 3; ++column)
                {
                    result[static_cast<std::size_t>(row)] +=
                        rotation[static_cast<std::size_t>(row * 3 + column)] * point[static_cast<std::size_t>(column)];
                }
            }
            return result;
        }

        RotationMatrix rotationFromDelta(const std::array<double, 6>& delta) noexcept
        {
            const double wx = delta[0];
            const double wy = delta[1];
            const double wz = delta[2];
            const double thetaSquared = wx * wx + wy * wy + wz * wz;
            if (thetaSquared < 1.0e-20)
            {
                return {1.0, -wz, wy, wz, 1.0, -wx, -wy, wx, 1.0};
            }
            const double theta = std::sqrt(thetaSquared);
            const double sineOverTheta = std::sin(theta) / theta;
            const double oneMinusCosineOverThetaSquared = (1.0 - std::cos(theta)) / thetaSquared;
            return {1.0 - oneMinusCosineOverThetaSquared * (wy * wy + wz * wz),
                    oneMinusCosineOverThetaSquared * wx * wy - sineOverTheta * wz,
                    oneMinusCosineOverThetaSquared * wx * wz + sineOverTheta * wy,
                    oneMinusCosineOverThetaSquared * wx * wy + sineOverTheta * wz,
                    1.0 - oneMinusCosineOverThetaSquared * (wx * wx + wz * wz),
                    oneMinusCosineOverThetaSquared * wy * wz - sineOverTheta * wx,
                    oneMinusCosineOverThetaSquared * wx * wz - sineOverTheta * wy,
                    oneMinusCosineOverThetaSquared * wy * wz + sineOverTheta * wx,
                    1.0 - oneMinusCosineOverThetaSquared * (wx * wx + wy * wy)};
        }

    } // namespace

    bool applyPoseDelta(RotationMatrix* rotation, Vector3* center, const std::array<double, 6>& delta) noexcept
    {
        if (!rotation || !center || !validRotation(*rotation) || !finite(*center) ||
            !std::all_of(delta.begin(), delta.end(), [](double value) { return std::isfinite(value); }))
        {
            return false;
        }
        const RotationMatrix increment = rotationFromDelta(delta);
        *rotation = multiply(increment, *rotation);
        for (int axis = 0; axis < 3; ++axis)
        {
            (*center)[static_cast<std::size_t>(axis)] += delta[static_cast<std::size_t>(axis + 3)];
        }
        return validRotation(*rotation) && finite(*center);
    }

    Result<void> validateRigTopology(const RigTopology& topology, std::size_t cameraCount)
    {
        if (topology.empty())
        {
            return Result<void>::success();
        }
        if (topology.captures.empty() || topology.sensors.empty() || topology.cameraBindings.size() != cameraCount)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "rig topology requires captures, sensors, and exactly one binding per camera");
        }
        std::set<Key> captureKeys;
        for (const RigCapture& capture : topology.captures)
        {
            if (capture.rigId < 0 || capture.captureId < 0 || !validRotation(capture.rigToWorldRotation) ||
                !finite(capture.rigCenterInWorld) || !captureKeys.emplace(capture.rigId, capture.captureId).second)
            {
                return Result<void>::failure(CameraErrorCode::InvalidPose,
                                             "rig captures must have unique non-negative ids and finite rigid poses");
            }
        }
        std::set<Key> sensorKeys;
        std::set<int> rigsWithFixedSensor;
        for (const RigSensor& sensor : topology.sensors)
        {
            if (sensor.rigId < 0 || sensor.sensorId < 0 || !validRotation(sensor.cameraToRigRotation) ||
                !finite(sensor.cameraCenterInRig) || !sensorKeys.emplace(sensor.rigId, sensor.sensorId).second)
            {
                return Result<void>::failure(CameraErrorCode::InvalidPose,
                                             "rig sensors must have unique non-negative ids and finite extrinsics");
            }
            if (sensor.fixedExtrinsic)
            {
                rigsWithFixedSensor.insert(sensor.rigId);
            }
        }
        for (const RigCapture& capture : topology.captures)
        {
            if (!rigsWithFixedSensor.count(capture.rigId))
            {
                return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                             "each rig requires at least one fixed sensor extrinsic");
            }
        }
        std::set<int> cameraIndices;
        for (const RigCameraBinding& binding : topology.cameraBindings)
        {
            if (binding.cameraIndex < 0 || static_cast<std::size_t>(binding.cameraIndex) >= cameraCount ||
                binding.rigId < 0 || binding.captureId < 0 || binding.sensorId < 0 ||
                !cameraIndices.insert(binding.cameraIndex).second ||
                !captureKeys.count({binding.rigId, binding.captureId}) ||
                !sensorKeys.count({binding.rigId, binding.sensorId}))
            {
                return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                             "rig camera bindings must be unique and reference existing capture and sensor");
            }
        }
        if (cameraIndices.size() != cameraCount)
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "rig camera bindings must cover every camera exactly once");
        }
        return Result<void>::success();
    }

    Result<std::vector<FramePinholeNumericState>>
    composeRigCameras(const std::vector<FramePinholeNumericState>& cameraModels, const RigTopology& topology)
    {
        const auto valid = validateRigTopology(topology, cameraModels.size());
        if (!valid)
        {
            return Result<std::vector<FramePinholeNumericState>>::failure(valid.error());
        }
        if (topology.empty())
        {
            return Result<std::vector<FramePinholeNumericState>>::success(cameraModels);
        }
        std::map<Key, const RigCapture*> captures;
        std::map<Key, const RigSensor*> sensors;
        for (const RigCapture& capture : topology.captures)
        {
            captures[{capture.rigId, capture.captureId}] = &capture;
        }
        for (const RigSensor& sensor : topology.sensors)
        {
            sensors[{sensor.rigId, sensor.sensorId}] = &sensor;
        }
        std::vector<FramePinholeNumericState> result = cameraModels;
        for (const RigCameraBinding& binding : topology.cameraBindings)
        {
            const RigCapture& capture = *captures.at({binding.rigId, binding.captureId});
            const RigSensor& sensor = *sensors.at({binding.rigId, binding.sensorId});
            auto& camera = result[static_cast<std::size_t>(binding.cameraIndex)];
            const RotationMatrix rotation = multiply(capture.rigToWorldRotation, sensor.cameraToRigRotation);
            const Vector3 offset = transform(capture.rigToWorldRotation, sensor.cameraCenterInRig);
            Vector3 center = capture.rigCenterInWorld;
            for (int axis = 0; axis < 3; ++axis)
            {
                center[static_cast<std::size_t>(axis)] += offset[static_cast<std::size_t>(axis)];
            }
            camera.setPose(Pose::create(camera.groundFrame(), center, rotation));
        }
        return Result<std::vector<FramePinholeNumericState>>::success(std::move(result));
    }

} // namespace placamera
