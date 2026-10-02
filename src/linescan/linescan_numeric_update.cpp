#include "placamera/linescan_numeric_state.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <utility>

namespace placamera
{
    namespace
    {
        RotationMatrix multiplyRotation(const RotationMatrix& first, const RotationMatrix& second) noexcept
        {
            RotationMatrix result{};
            for (int row = 0; row < 3; ++row)
            {
                for (int column = 0; column < 3; ++column)
                {
                    for (int inner = 0; inner < 3; ++inner)
                    {
                        result[static_cast<std::size_t>(row * 3 + column)] +=
                            first[static_cast<std::size_t>(row * 3 + inner)] *
                            second[static_cast<std::size_t>(inner * 3 + column)];
                    }
                }
            }
            return result;
        }

        RotationMatrix angleAxisRotation(const Vector3& rotation) noexcept
        {
            const double angle_squared =
                rotation[0] * rotation[0] + rotation[1] * rotation[1] + rotation[2] * rotation[2];
            if (angle_squared < 1.0e-20)
            {
                return {1.0, -rotation[2], rotation[1], rotation[2], 1.0, -rotation[0], -rotation[1], rotation[0], 1.0};
            }
            const double angle = std::sqrt(angle_squared);
            const double sine_over_angle = std::sin(angle) / angle;
            const double one_minus_cosine_over_angle_squared = (1.0 - std::cos(angle)) / angle_squared;
            const double x = rotation[0];
            const double y = rotation[1];
            const double z = rotation[2];
            return {1.0 - one_minus_cosine_over_angle_squared * (y * y + z * z),
                    one_minus_cosine_over_angle_squared * x * y - sine_over_angle * z,
                    one_minus_cosine_over_angle_squared * x * z + sine_over_angle * y,
                    one_minus_cosine_over_angle_squared * x * y + sine_over_angle * z,
                    1.0 - one_minus_cosine_over_angle_squared * (x * x + z * z),
                    one_minus_cosine_over_angle_squared * y * z - sine_over_angle * x,
                    one_minus_cosine_over_angle_squared * x * z - sine_over_angle * y,
                    one_minus_cosine_over_angle_squared * y * z + sine_over_angle * x,
                    1.0 - one_minus_cosine_over_angle_squared * (x * x + y * y)};
        }

        bool validOptics(const LineScanOptics& optics) noexcept
        {
            if (!std::isfinite(optics.focalLengthMillimeters) || !(optics.focalLengthMillimeters > 0.0) ||
                !std::isfinite(optics.distortionK1))
            {
                return false;
            }
            if (!std::isfinite(optics.samplePitchMillimeters) || !std::isfinite(optics.principalSample))
            {
                return false;
            }
            if (optics.completeCalibration)
            {
                const MetashapeCalibration& calibration = *optics.completeCalibration;
                const double values[] = {calibration.f,
                                         calibration.cx,
                                         calibration.cy,
                                         calibration.b1,
                                         calibration.b2,
                                         calibration.k1,
                                         calibration.k2,
                                         calibration.k3,
                                         calibration.k4,
                                         calibration.p1,
                                         calibration.p2,
                                         calibration.p3,
                                         calibration.p4};
                if (!std::all_of(
                        std::begin(values), std::end(values), [](double value) { return std::isfinite(value); }) ||
                    !(calibration.f > 0.0) || !(calibration.f + calibration.b1 > 0.0))
                {
                    return false;
                }
                if (calibration.principalPointDecomposition)
                {
                    const PrincipalPointDecomposition& principal = *calibration.principalPointDecomposition;
                    const double principal_values[] = {
                        principal.imageCenterX, principal.imageCenterY, principal.cxOffset, principal.cyOffset};
                    if (!std::all_of(std::begin(principal_values),
                                     std::end(principal_values),
                                     [](double value) { return std::isfinite(value); }))
                    {
                        return false;
                    }
                }
            }
            if (!optics.detectorGeometry)
            {
                return optics.completeCalibration || optics.samplePitchMillimeters > 0.0;
            }
            if (optics.completeCalibration && !(optics.samplePitchMillimeters > 0.0))
            {
                return false;
            }
            const auto& detector = *optics.detectorGeometry;
            const double determinant = detector.focalToPixelLines[1] * detector.focalToPixelSamples[2] -
                                       detector.focalToPixelLines[2] * detector.focalToPixelSamples[1];
            const double values[] = {detector.detectorSampleSumming,
                                     detector.detectorLineSumming,
                                     detector.detectorSampleOrigin,
                                     detector.detectorLineOrigin,
                                     detector.startingDetectorSample,
                                     detector.startingDetectorLine,
                                     detector.focalToPixelSamples[0],
                                     detector.focalToPixelSamples[1],
                                     detector.focalToPixelSamples[2],
                                     detector.focalToPixelLines[0],
                                     detector.focalToPixelLines[1],
                                     detector.focalToPixelLines[2]};
            return std::all_of(
                       std::begin(values), std::end(values), [](double value) { return std::isfinite(value); }) &&
                   detector.detectorSampleSumming > 0.0 && detector.detectorLineSumming > 0.0 &&
                   std::abs(determinant) >= 1.0e-15;
        }

    } // namespace
    Result<void> LineScanNumericState::applyOptimizationDelta(std::span<const double> delta)
    {
        if (delta.size() != _layout.parameterCount())
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "line-scan numeric update size does not match its cached layout");
        }
        if (!std::all_of(delta.begin(), delta.end(), [](double value) { return std::isfinite(value); }))
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "line-scan numeric update must contain only finite values");
        }

        LineScanOptics next_optics = _optics;
        LineScanTrajectoryBias next_bias = _bias;
        std::size_t cursor = 0;
        for (const auto& sample : _samples)
        {
            if (_selection.knotPositions && !sample.constraints.positionFixed)
            {
                cursor += 3;
            }
            if (_selection.knotRotations && !sample.constraints.rotationFixed)
            {
                cursor += 3;
            }
        }
        if (_selection.globalTranslation)
        {
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
                next_bias.translationMeters[axis] += delta[cursor++];
            }
        }
        if (_selection.globalRotation)
        {
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
                next_bias.rotationVectorRadians[axis] += delta[cursor++];
            }
        }
        if (_selection.timeOffset)
        {
            next_bias.timeOffsetSeconds += delta[cursor++];
        }

        bool definition_changed = false;
        const auto scalar = [&](bool selected, double* value)
        {
            if (selected)
            {
                definition_changed = definition_changed || delta[cursor] != 0.0;
                *value += delta[cursor++];
            }
        };
        if (_selection.calibration.any())
        {
            MetashapeCalibration& calibration = *next_optics.completeCalibration;
            scalar(_selection.calibration.f, &calibration.f);
            const auto principal_scalar = [&](bool selected, bool x_axis)
            {
                if (!selected)
                {
                    return;
                }
                const double value = delta[cursor++];
                definition_changed = definition_changed || value != 0.0;
                if (value == 0.0)
                {
                    return;
                }
                double& absolute = x_axis ? calibration.cx : calibration.cy;
                if (!calibration.principalPointDecomposition)
                {
                    absolute += value;
                    return;
                }
                PrincipalPointDecomposition& principal = *calibration.principalPointDecomposition;
                double& offset = x_axis ? principal.cxOffset : principal.cyOffset;
                offset += value;
                absolute = (x_axis ? principal.imageCenterX : principal.imageCenterY) + offset;
            };
            principal_scalar(_selection.calibration.cx, true);
            principal_scalar(_selection.calibration.cy, false);
            scalar(_selection.calibration.b1, &calibration.b1);
            scalar(_selection.calibration.b2, &calibration.b2);
            scalar(_selection.calibration.k1, &calibration.k1);
            scalar(_selection.calibration.k2, &calibration.k2);
            scalar(_selection.calibration.k3, &calibration.k3);
            scalar(_selection.calibration.k4, &calibration.k4);
            scalar(_selection.calibration.p1, &calibration.p1);
            scalar(_selection.calibration.p2, &calibration.p2);
            scalar(_selection.calibration.p3, &calibration.p3);
            scalar(_selection.calibration.p4, &calibration.p4);
        }
        scalar(_selection.detector.focalLength, &next_optics.focalLengthMillimeters);
        scalar(_selection.detector.distortionK1, &next_optics.distortionK1);
        scalar(_selection.detector.samplePitch, &next_optics.samplePitchMillimeters);
        scalar(_selection.detector.principalSample, &next_optics.principalSample);
        if (_selection.detector.anyDetectorGeometryParameter())
        {
            auto& detector = *next_optics.detectorGeometry;
            scalar(_selection.detector.detectorSampleSumming, &detector.detectorSampleSumming);
            scalar(_selection.detector.detectorLineSumming, &detector.detectorLineSumming);
            scalar(_selection.detector.detectorSampleOrigin, &detector.detectorSampleOrigin);
            scalar(_selection.detector.detectorLineOrigin, &detector.detectorLineOrigin);
            scalar(_selection.detector.startingDetectorSample, &detector.startingDetectorSample);
            scalar(_selection.detector.startingDetectorLine, &detector.startingDetectorLine);
            if (_selection.detector.focalToPixelSamples)
            {
                for (double& value : detector.focalToPixelSamples)
                {
                    scalar(true, &value);
                }
            }
            if (_selection.detector.focalToPixelLines)
            {
                for (double& value : detector.focalToPixelLines)
                {
                    scalar(true, &value);
                }
            }
        }
        if (!validOptics(next_optics))
        {
            return Result<void>::failure(CameraErrorCode::InvalidIntrinsics,
                                         "line-scan numeric update produced invalid detector optics");
        }

        cursor = 0;
        for (std::size_t index = 0; index < _samples.size(); ++index)
        {
            const auto& constraints = _samples[index].constraints;
            if (_selection.knotPositions && !constraints.positionFixed)
            {
                for (std::size_t axis = 0; axis < 3; ++axis)
                {
                    _positionDeltas[index][axis] += delta[cursor++];
                    _samples[index].center[axis] = _nominalSamples[index].center[axis] + _positionDeltas[index][axis];
                }
            }
            if (_selection.knotRotations && !constraints.rotationFixed)
            {
                for (double& value : _rotationDeltas[index])
                {
                    value += delta[cursor++];
                }
                _samples[index].cameraToWorldRotation = multiplyRotation(angleAxisRotation(_rotationDeltas[index]),
                                                                         _nominalSamples[index].cameraToWorldRotation);
            }
        }
        _bias = next_bias;
        _optics = next_optics;
        _definitionDirty = _definitionDirty || definition_changed;
        return Result<void>::success();
    }

} // namespace placamera
