#include "placamera/linescan_camera.h"

#include "../internal/plamatrix_rotation.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace placamera
{

    namespace
    {

        template <typename Sample> std::size_t lowerInterval(const std::vector<Sample>& samples, double seconds)
        {
            const auto upper = std::upper_bound(samples.begin(),
                                                samples.end(),
                                                seconds,
                                                [](double candidate, const Sample& sample)
                                                { return candidate < sample.time.seconds; });
            if (upper == samples.begin())
            {
                return 0;
            }
            return std::min<std::size_t>(samples.size() - 2,
                                         static_cast<std::size_t>(std::distance(samples.begin(), upper) - 1));
        }

        bool supports(const std::vector<QuaternionTrajectorySample>& samples, double seconds) noexcept
        {
            return samples.size() >= 2 && seconds >= samples.front().time.seconds - 1.0e-9 &&
                   seconds <= samples.back().time.seconds + 1.0e-9;
        }

        RotationMatrix rotationAt(const FrameRotationTrajectory& trajectory, double seconds)
        {
            const std::size_t index = lowerInterval(trajectory.samples, seconds);
            const QuaternionTrajectorySample& first = trajectory.samples[index];
            const QuaternionTrajectorySample& second = trajectory.samples[index + 1];
            const double fraction =
                std::clamp((seconds - first.time.seconds) / (second.time.seconds - first.time.seconds), 0.0, 1.0);
            const auto first_quaternion = internal::quaternionFromScalarFirst(first.scalarFirst);
            const auto second_quaternion = internal::quaternionFromScalarFirst(second.scalarFirst);
            return internal::multiply(
                trajectory.constantRotation,
                internal::rotationFromQuaternion(first_quaternion.slerp(fraction, second_quaternion)));
        }

        void validateQuaternionTrajectory(const FrameRotationTrajectory& trajectory, TimeScale scale, const char* label)
        {
            if (trajectory.samples.size() < 2)
            {
                throw CameraValidationError(CameraErrorCode::InvalidTime,
                                            std::string(label) + " requires at least two samples");
            }
            (void)Pose::create(FrameId("trajectory-validation"), {0.0, 0.0, 0.0}, trajectory.constantRotation);
            for (std::size_t index = 0; index < trajectory.samples.size(); ++index)
            {
                const QuaternionTrajectorySample& sample = trajectory.samples[index];
                if (sample.time.scale != scale || !std::isfinite(sample.time.seconds) ||
                    (index > 0 && sample.time.seconds <= trajectory.samples[index - 1].time.seconds))
                {
                    throw CameraValidationError(CameraErrorCode::InvalidTime,
                                                std::string(label) + " times are invalid");
                }
                if (!internal::validQuaternion(sample.scalarFirst))
                {
                    throw CameraValidationError(CameraErrorCode::InvalidPose,
                                                std::string(label) + " contains an invalid quaternion");
                }
            }
        }

        EvaluationResult<Pose> makePose(FrameId frame, Vector3 center, RotationMatrix rotation)
        {
            try
            {
                return EvaluationResult<Pose>::success(Pose::create(std::move(frame), center, rotation));
            }
            catch (const CameraValidationError& error)
            {
                return EvaluationResult<Pose>::failure(error.code(), error.what());
            }
        }

    } // namespace

    TrajectorySample::TrajectorySample(TimeReference sampleTime,
                                       Vector3 sampleCenter,
                                       RotationMatrix sampleCameraToWorldRotation,
                                       TrajectoryKnotConstraints sampleConstraints)
        : time(sampleTime), center(sampleCenter), cameraToWorldRotation(sampleCameraToWorldRotation),
          constraints(std::move(sampleConstraints))
    {
    }

    LineScanTrajectory LineScanTrajectory::create(std::vector<TrajectorySample> samples)
    {
        if (samples.size() < 2)
        {
            throw CameraValidationError(CameraErrorCode::InvalidTime,
                                        "line-scan trajectory requires at least two samples");
        }
        const TimeScale scale = samples.front().time.scale;
        for (std::size_t index = 0; index < samples.size(); ++index)
        {
            const TrajectorySample& sample = samples[index];
            if (sample.time.scale != scale || !std::isfinite(sample.time.seconds) ||
                (index > 0 && sample.time.seconds <= samples[index - 1].time.seconds))
            {
                throw CameraValidationError(CameraErrorCode::InvalidTime,
                                            "line-scan trajectory times must be finite, same-scale, and increasing");
            }
            (void)Pose::create(FrameId("trajectory-validation"), sample.center, sample.cameraToWorldRotation);
            const auto valid_sigma = [](const std::optional<Vector3>& sigma)
            {
                return !sigma || std::all_of(sigma->begin(),
                                             sigma->end(),
                                             [](double value) { return std::isfinite(value) && value > 0.0; });
            };
            if (!valid_sigma(sample.constraints.positionSigmaMeters) ||
                !valid_sigma(sample.constraints.rotationSigmaRadians))
            {
                throw CameraValidationError(CameraErrorCode::InvalidModelState,
                                            "line-scan knot prior sigmas must be finite and positive");
            }
        }
        return LineScanTrajectory(std::move(samples));
    }

    LineScanTrajectory LineScanTrajectory::createFrameComposed(FrameComposedTrajectory trajectory)
    {
        if (trajectory.inertialStates.size() < 2)
        {
            throw CameraValidationError(CameraErrorCode::InvalidTime,
                                        "frame-composed trajectory requires at least two states");
        }
        const TimeScale scale = trajectory.inertialStates.front().time.scale;
        for (std::size_t index = 0; index < trajectory.inertialStates.size(); ++index)
        {
            const TranslationalStateSample& state = trajectory.inertialStates[index];
            const bool finite_state = std::all_of(state.positionMeters.begin(),
                                                  state.positionMeters.end(),
                                                  [](double value) { return std::isfinite(value); }) &&
                                      std::all_of(state.velocityMetersPerSecond.begin(),
                                                  state.velocityMetersPerSecond.end(),
                                                  [](double value) { return std::isfinite(value); });
            if (state.time.scale != scale || !std::isfinite(state.time.seconds) || !finite_state ||
                (index > 0 && state.time.seconds <= trajectory.inertialStates[index - 1].time.seconds))
            {
                throw CameraValidationError(CameraErrorCode::InvalidTime,
                                            "frame-composed translational states are invalid");
            }
        }
        validateQuaternionTrajectory(trajectory.inertialToWorld, scale, "inertial-to-world trajectory");
        validateQuaternionTrajectory(trajectory.inertialToSensor, scale, "inertial-to-sensor trajectory");
        return LineScanTrajectory(std::move(trajectory));
    }

    const std::vector<TrajectorySample>& LineScanTrajectory::samples() const noexcept
    {
        return _samples;
    }

    const FrameComposedTrajectory* LineScanTrajectory::frameComposed() const noexcept
    {
        return _frameComposed ? &*_frameComposed : nullptr;
    }

    TimeScale LineScanTrajectory::timeScale() const noexcept
    {
        return _timeScale;
    }

    EvaluationResult<Pose> LineScanTrajectory::poseAt(TimeReference time, const FrameId& frame) const
    {
        if (time.scale != _timeScale || !std::isfinite(time.seconds))
        {
            return EvaluationResult<Pose>::failure(CameraErrorCode::InvalidTime,
                                                   "line-scan pose time uses the wrong time scale");
        }
        if (_frameComposed)
        {
            const FrameComposedTrajectory& trajectory = *_frameComposed;
            if (time.seconds < trajectory.inertialStates.front().time.seconds ||
                time.seconds > trajectory.inertialStates.back().time.seconds ||
                !supports(trajectory.inertialToWorld.samples, time.seconds) ||
                !supports(trajectory.inertialToSensor.samples, time.seconds))
            {
                return EvaluationResult<Pose>::failure(CameraErrorCode::OutsideModelDomain,
                                                       "line-scan pose time is outside composed trajectory support");
            }

            const std::size_t state_index = lowerInterval(trajectory.inertialStates, time.seconds);
            const TranslationalStateSample& first = trajectory.inertialStates[state_index];
            const TranslationalStateSample& second = trajectory.inertialStates[state_index + 1];
            const double duration = second.time.seconds - first.time.seconds;
            const double fraction = std::clamp((time.seconds - first.time.seconds) / duration, 0.0, 1.0);
            const double h00 = 2.0 * fraction * fraction * fraction - 3.0 * fraction * fraction + 1.0;
            const double h10 = fraction * fraction * fraction - 2.0 * fraction * fraction + fraction;
            const double h01 = -2.0 * fraction * fraction * fraction + 3.0 * fraction * fraction;
            const double h11 = fraction * fraction * fraction - fraction * fraction;
            Vector3 inertial_center{};
            for (int axis = 0; axis < 3; ++axis)
            {
                inertial_center[static_cast<std::size_t>(axis)] =
                    h00 * first.positionMeters[static_cast<std::size_t>(axis)] +
                    h10 * duration * first.velocityMetersPerSecond[static_cast<std::size_t>(axis)] +
                    h01 * second.positionMeters[static_cast<std::size_t>(axis)] +
                    h11 * duration * second.velocityMetersPerSecond[static_cast<std::size_t>(axis)];
            }
            const RotationMatrix inertial_to_world = rotationAt(trajectory.inertialToWorld, time.seconds);
            const RotationMatrix inertial_to_sensor = rotationAt(trajectory.inertialToSensor, time.seconds);
            return makePose(frame,
                            internal::multiply(inertial_to_world, inertial_center),
                            internal::multiply(inertial_to_world, internal::transposeArray(inertial_to_sensor)));
        }

        if (time.seconds < _samples.front().time.seconds || time.seconds > _samples.back().time.seconds)
        {
            return EvaluationResult<Pose>::failure(CameraErrorCode::OutsideModelDomain,
                                                   "line-scan pose time is outside trajectory support");
        }
        const auto upper = std::upper_bound(_samples.begin(),
                                            _samples.end(),
                                            time.seconds,
                                            [](double seconds, const TrajectorySample& sample)
                                            { return seconds < sample.time.seconds; });
        const std::size_t upper_index =
            upper == _samples.begin() ? 0 : static_cast<std::size_t>(std::distance(_samples.begin(), upper) - 1);
        if (upper_index >= _samples.size() - 1)
        {
            const TrajectorySample& last = _samples.back();
            return makePose(frame, last.center, last.cameraToWorldRotation);
        }

        const TrajectorySample& first = _samples[upper_index];
        const TrajectorySample& second = _samples[upper_index + 1];
        const double fraction = (time.seconds - first.time.seconds) / (second.time.seconds - first.time.seconds);
        const auto orientation = internal::quaternionFromRotation(first.cameraToWorldRotation)
                                     .slerp(std::clamp(fraction, 0.0, 1.0),
                                            internal::quaternionFromRotation(second.cameraToWorldRotation));
        const Vector3 center{first.center[0] + fraction * (second.center[0] - first.center[0]),
                             first.center[1] + fraction * (second.center[1] - first.center[1]),
                             first.center[2] + fraction * (second.center[2] - first.center[2])};
        return makePose(frame, center, internal::rotationFromQuaternion(orientation));
    }

    LineScanTrajectory::LineScanTrajectory(std::vector<TrajectorySample> samples)
        : _samples(std::move(samples)), _timeScale(_samples.front().time.scale)
    {
    }

    LineScanTrajectory::LineScanTrajectory(FrameComposedTrajectory trajectory)
        : _frameComposed(std::move(trajectory)), _timeScale(_frameComposed->inertialStates.front().time.scale)
    {
    }

} // namespace placamera
