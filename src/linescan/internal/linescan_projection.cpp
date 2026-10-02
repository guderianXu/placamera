#include "linescan_projection.h"

#include "../../internal/plamatrix_rotation.h"

#include <algorithm>
#include <cmath>

namespace placamera::internal
{

    namespace
    {

        bool normalize(Vector3* vector) noexcept
        {
            const double norm = std::hypot((*vector)[0], std::hypot((*vector)[1], (*vector)[2]));
            if (!std::isfinite(norm) || !(norm > 0.0))
            {
                return false;
            }
            for (double& value : *vector)
            {
                value /= norm;
            }
            return true;
        }

        bool validOptions(const EvaluationOptions& options) noexcept
        {
            return std::isfinite(options.desiredPrecisionPixels) && options.desiredPrecisionPixels > 0.0 &&
                   options.maximumIterations > 0;
        }

        double minimumCoordinate(LineScanPixelConvention convention) noexcept
        {
            return convention == LineScanPixelConvention::PixelCenter ? 0.5 : 0.0;
        }

        bool insideImage(const LineScanModel& model, const ImageCoordinate& image) noexcept
        {
            const double minimum = minimumCoordinate(model.lineScanDefinition().pixelConvention());
            return image.sample >= minimum && image.line >= minimum &&
                   image.sample <= minimum + static_cast<double>(model.imageSize().samples - 1) &&
                   image.line <= minimum + static_cast<double>(model.imageSize().lines - 1);
        }

        EvaluationResult<Pose> poseForLine(const LineScanModel& model, double line, const LineScanTrajectoryBias& bias)
        {
            const auto nominal_time = model.timeForLine(line);
            if (!nominal_time)
            {
                return EvaluationResult<Pose>::failure(nominal_time.errorCode(), nominal_time.message());
            }

            const TimeReference biased_time{nominal_time.value().scale,
                                            nominal_time.value().seconds + bias.timeOffsetSeconds};
            const auto nominal_pose = model.trajectory().poseAt(biased_time, model.groundFrame());
            if (!nominal_pose)
            {
                return nominal_pose;
            }

            Vector3 center = nominal_pose.value().center;
            for (int axis = 0; axis < 3; ++axis)
            {
                center[static_cast<std::size_t>(axis)] += bias.translationMeters[static_cast<std::size_t>(axis)];
            }
            try
            {
                return EvaluationResult<Pose>::success(
                    Pose::create(model.groundFrame(),
                                 center,
                                 multiply(angleAxisRotation(bias.rotationVectorRadians),
                                          nominal_pose.value().cameraToWorldRotation)));
            }
            catch (const CameraValidationError& error)
            {
                return EvaluationResult<Pose>::failure(error.code(), error.what());
            }
        }

    } // namespace

    EvaluationResult<LineScanProjectionDetails> projectLineScanAtLine(const LineScanModel& model,
                                                                      const GroundCoordinate& ground,
                                                                      double line,
                                                                      const EvaluationOptions& options)
    {
        return projectLineScanAtLine(model, ground, line, model.trajectoryBias(), options);
    }

    EvaluationResult<LineScanProjectionDetails> projectLineScanAtLine(const LineScanModel& model,
                                                                      const GroundCoordinate& ground,
                                                                      double line,
                                                                      const LineScanTrajectoryBias& bias,
                                                                      const EvaluationOptions& options)
    {
        if (!validOptions(options) || !std::isfinite(line) || !std::isfinite(bias.timeOffsetSeconds) ||
            !std::all_of(bias.translationMeters.begin(),
                         bias.translationMeters.end(),
                         [](double value) { return std::isfinite(value); }) ||
            !std::all_of(bias.rotationVectorRadians.begin(),
                         bias.rotationVectorRadians.end(),
                         [](double value) { return std::isfinite(value); }) ||
            !std::all_of(
                ground.position.begin(), ground.position.end(), [](double value) { return std::isfinite(value); }))
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::InvalidArgument,
                "line-scan projection inputs and options must be finite and positive");
        }
        if (ground.frame != model.groundFrame())
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::FrameMismatch, "ground coordinate frame does not match the line-scan model");
        }

        const auto time = model.timeForLine(line);
        const auto pose = poseForLine(model, line, bias);
        if (!time || !pose)
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::OutsideModelDomain, "line is outside line-scan timing or trajectory support");
        }

        const Vector3 offset{ground.position[0] - pose.value().center[0],
                             ground.position[1] - pose.value().center[1],
                             ground.position[2] - pose.value().center[2]};
        const Vector3 sensor = transposeMultiply(pose.value().cameraToWorldRotation, offset);
        if (!std::isfinite(sensor[2]) || sensor[2] <= 1.0e-9)
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::OutsideModelDomain, "ground coordinate is not in front of the line-scan sensor");
        }

        const double focal_length = model.lineScanDefinition().optics().focalLengthMillimeters;
        const FocalPlaneCoordinate focal{focal_length * sensor[0] / sensor[2], focal_length * sensor[1] / sensor[2]};
        const auto detector = lineScanUndistortedFocalToPixel(model.lineScanDefinition(), focal, options);
        if (!detector)
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(detector.errorCode(), detector.message());
        }

        const TimeReference acquisition_time{time.value().scale, time.value().seconds + bias.timeOffsetSeconds};
        const Projection projection{ImageCoordinate{detector.value().sample, line}, sensor[2], acquisition_time};
        const LineScanProjectionDetails details{projection, detector.value().lineResidualPixels, focal};
        return EvaluationResult<LineScanProjectionDetails>::success(details,
                                                                    std::abs(detector.value().lineResidualPixels));
    }

    EvaluationResult<LineScanProjectionDetails>
    projectLineScan(const LineScanModel& model, const GroundCoordinate& ground, const EvaluationOptions& options)
    {
        if (!validOptions(options))
        {
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::InvalidArgument, "line-scan evaluation options must be finite and positive");
        }

        const double minimum_line = minimumCoordinate(model.lineScanDefinition().pixelConvention());
        const double maximum_line = minimum_line + static_cast<double>(model.imageSize().lines - 1);
        double line0 = 0.5 * (minimum_line + maximum_line);
        double line1 = std::min(maximum_line, line0 + 0.1);
        auto projection0 = projectLineScanAtLine(model, ground, line0, options);
        if (!projection0)
        {
            return projection0;
        }
        if (line1 == line0)
        {
            if (std::abs(projection0.value().lineResidualPixels) <= options.desiredPrecisionPixels)
            {
                return projection0;
            }
            return EvaluationResult<LineScanProjectionDetails>::failure(
                CameraErrorCode::NonConvergence, "single-line sensor residual exceeds the requested precision");
        }
        auto projection1 = projectLineScanAtLine(model, ground, line1, options);
        if (!projection1)
        {
            return projection1;
        }

        for (int iteration = 0; iteration < options.maximumIterations; ++iteration)
        {
            if (std::abs(projection1.value().lineResidualPixels) <= options.desiredPrecisionPixels)
            {
                if (options.requireInsideImage && !insideImage(model, projection1.value().projection.image))
                {
                    return EvaluationResult<LineScanProjectionDetails>::failure(
                        CameraErrorCode::OutsideModelDomain, "line-scan projection lies outside the image");
                }
                return projection1;
            }

            const double denominator = projection1.value().lineResidualPixels - projection0.value().lineResidualPixels;
            if (std::abs(denominator) < 1.0e-15)
            {
                break;
            }
            const double next_line =
                std::clamp(line1 - projection1.value().lineResidualPixels * (line1 - line0) / denominator,
                           minimum_line,
                           maximum_line);
            line0 = line1;
            projection0 = projection1;
            line1 = next_line;
            projection1 = projectLineScanAtLine(model, ground, line1, options);
            if (!projection1)
            {
                return projection1;
            }
        }

        return EvaluationResult<LineScanProjectionDetails>::failure(
            CameraErrorCode::NonConvergence, "line-scan projection did not reach the requested line precision");
    }

    EvaluationResult<ImagingLocus>
    lineScanImagingLocus(const LineScanModel& model, const ImageCoordinate& image, const EvaluationOptions& options)
    {
        if (!validOptions(options) || !std::isfinite(image.sample) || !std::isfinite(image.line))
        {
            return EvaluationResult<ImagingLocus>::failure(
                CameraErrorCode::InvalidArgument, "line-scan image coordinate and options must be finite and positive");
        }
        if (options.requireInsideImage && !insideImage(model, image))
        {
            return EvaluationResult<ImagingLocus>::failure(CameraErrorCode::OutsideModelDomain,
                                                           "image coordinate lies outside the line-scan image");
        }

        const auto nominal_time = model.timeForLine(image.line);
        const auto pose = poseForLine(model, image.line, model.trajectoryBias());
        const auto focal = lineScanPixelToUndistortedFocal(model.lineScanDefinition(), image.sample, options);
        if (!nominal_time || !pose || !focal)
        {
            return EvaluationResult<ImagingLocus>::failure(
                CameraErrorCode::OutsideModelDomain, "line-scan ray is outside timing, trajectory, or optical support");
        }

        Vector3 sensor_direction{focal.value().xMillimeters,
                                 focal.value().yMillimeters,
                                 model.lineScanDefinition().optics().focalLengthMillimeters};
        if (!normalize(&sensor_direction))
        {
            return EvaluationResult<ImagingLocus>::failure(CameraErrorCode::OutsideModelDomain,
                                                           "line-scan sensor ray is invalid");
        }
        Vector3 direction = multiply(pose.value().cameraToWorldRotation, sensor_direction);
        if (!normalize(&direction))
        {
            return EvaluationResult<ImagingLocus>::failure(CameraErrorCode::OutsideModelDomain,
                                                           "line-scan world ray is invalid");
        }

        const TimeReference acquisition_time{nominal_time.value().scale,
                                             nominal_time.value().seconds + model.trajectoryBias().timeOffsetSeconds};
        const ImagingLocus locus{
            GroundCoordinate{model.groundFrame(), pose.value().center}, direction, acquisition_time};
        return EvaluationResult<ImagingLocus>::success(locus);
    }

} // namespace placamera::internal
