#include "rpc_math.h"

#include <algorithm>
#include <cmath>

namespace placamera
{

    namespace
    {

        constexpr double kDenominatorEpsilon = 1.0e-14;
        bool finite(double value) noexcept
        {
            return std::isfinite(value);
        }

        RpcCoefficients polynomialTerms(double longitude, double latitude, double height) noexcept
        {
            const double longitude2 = longitude * longitude;
            const double latitude2 = latitude * latitude;
            const double height2 = height * height;
            return {{1.0,
                     longitude,
                     latitude,
                     height,
                     longitude * latitude,
                     longitude * height,
                     latitude * height,
                     longitude2,
                     latitude2,
                     height2,
                     longitude * latitude * height,
                     longitude2 * longitude,
                     longitude * latitude2,
                     longitude * height2,
                     longitude2 * latitude,
                     latitude2 * latitude,
                     latitude * height2,
                     longitude2 * height,
                     latitude2 * height,
                     height2 * height}};
        }

        std::array<RpcCoefficients, 3>
        polynomialTermDerivatives(double longitude, double latitude, double height) noexcept
        {
            const double longitude2 = longitude * longitude;
            const double latitude2 = latitude * latitude;
            const double height2 = height * height;
            std::array<RpcCoefficients, 3> derivatives{};
            derivatives[0] = {{0.0, 1.0, 0.0, 0.0, latitude, height, 0.0, 2.0 * longitude, 0.0, 0.0,
                               latitude * height, 3.0 * longitude2, latitude2, height2,
                               2.0 * longitude * latitude, 0.0, 0.0, 2.0 * longitude * height, 0.0, 0.0}};
            derivatives[1] = {{0.0, 0.0, 1.0, 0.0, longitude, 0.0, height, 0.0, 2.0 * latitude, 0.0,
                               longitude * height, 0.0, 2.0 * longitude * latitude, 0.0, longitude2,
                               3.0 * latitude2, height2, 0.0, 2.0 * latitude * height, 0.0}};
            derivatives[2] = {{0.0, 0.0, 0.0, 1.0, 0.0, longitude, latitude, 0.0, 0.0, 2.0 * height,
                               longitude * latitude, 0.0, 0.0, 2.0 * longitude * height, 0.0, 0.0,
                               2.0 * latitude * height, longitude2, latitude2, 3.0 * height2}};
            return derivatives;
        }

        double dot(const RpcCoefficients& first, const RpcCoefficients& second) noexcept
        {
            double result = 0.0;
            for (std::size_t index = 0; index < first.size(); ++index)
            {
                result += first[index] * second[index];
            }
            return result;
        }

        double derivativeDot(const RpcCoefficients& coefficients,
                             const std::array<RpcCoefficients, 3>& derivatives,
                             int axis) noexcept
        {
            return dot(coefficients, derivatives[static_cast<std::size_t>(axis)]);
        }

        double longitudeDifference(double longitude, double reference) noexcept
        {
            double difference = longitude - reference;
            if (difference < -270.0)
            {
                difference += 360.0;
            }
            else if (difference > 270.0)
            {
                difference -= 360.0;
            }
            return difference;
        }

        ImageCoordinate applyImageCorrection(const ImageCoordinate& uncorrected,
                                             const RpcParameters& parameters,
                                             const RpcImageCorrection& correction) noexcept
        {
            const double normalized_sample = (uncorrected.sample - parameters.sampleOffset) / parameters.sampleScale;
            const double normalized_line = (uncorrected.line - parameters.lineOffset) / parameters.lineScale;
            return {uncorrected.sample + correction.sampleOffsetPixels +
                        correction.sampleSamplePixels * normalized_sample +
                        correction.sampleLinePixels * normalized_line,
                    uncorrected.line + correction.lineOffsetPixels + correction.lineSamplePixels * normalized_sample +
                        correction.lineLinePixels * normalized_line};
        }

        ImageCoordinate applyGroundCorrection(const ImageCoordinate& uncorrected,
                                              const RpcParameters& parameters,
                                              const GeodeticCoordinate& ground,
                                              const RpcGroundCorrection& correction) noexcept
        {
            const double longitude_delta = ground.longitudeDegrees - parameters.longitudeOffset;
            const double latitude_delta = ground.latitudeDegrees - parameters.latitudeOffset;
            const double height_delta = ground.heightMeters - parameters.heightOffset;
            return {uncorrected.sample + correction.sampleOffsetPixels +
                        correction.sampleLongitudePixelsPerDegree * longitude_delta +
                        correction.sampleLatitudePixelsPerDegree * latitude_delta +
                        correction.sampleHeightPixelsPerMeter * height_delta,
                    uncorrected.line + correction.lineOffsetPixels +
                        correction.lineLongitudePixelsPerDegree * longitude_delta +
                        correction.lineLatitudePixelsPerDegree * latitude_delta +
                        correction.lineHeightPixelsPerMeter * height_delta};
        }

        bool evaluate(const RpcDefinition& definition,
                      const RpcCorrection& correction,
                      double normalizedLongitude,
                      double normalizedLatitude,
                      double normalizedHeight,
                      bool applyCorrection,
                      ImageCoordinate* image) noexcept
        {
            if (!image || !finite(normalizedLongitude) || !finite(normalizedLatitude) || !finite(normalizedHeight))
            {
                return false;
            }

            const RpcParameters& parameters = definition.parameters();
            const RpcCoefficients terms = polynomialTerms(normalizedLongitude, normalizedLatitude, normalizedHeight);
            const double line_denominator = dot(parameters.lineDenominator, terms);
            const double sample_denominator = dot(parameters.sampleDenominator, terms);
            if (std::abs(line_denominator) < kDenominatorEpsilon || std::abs(sample_denominator) < kDenominatorEpsilon)
            {
                return false;
            }

            const ImageCoordinate uncorrected{
                parameters.sampleOffset +
                    parameters.sampleScale * dot(parameters.sampleNumerator, terms) / sample_denominator,
                parameters.lineOffset + parameters.lineScale * dot(parameters.lineNumerator, terms) / line_denominator};
            if (!finite(uncorrected.sample) || !finite(uncorrected.line))
            {
                return false;
            }
            if (!applyCorrection)
            {
                *image = uncorrected;
            }
            else if (const RpcImageCorrection* image_correction = correction.normalizedImageValue())
            {
                *image = applyImageCorrection(uncorrected, parameters, *image_correction);
            }
            else
            {
                const GeodeticCoordinate ground{
                    parameters.longitudeOffset + normalizedLongitude * parameters.longitudeScale,
                    parameters.latitudeOffset + normalizedLatitude * parameters.latitudeScale,
                    parameters.heightOffset + normalizedHeight * parameters.heightScale};
                *image = applyGroundCorrection(uncorrected, parameters, ground, *correction.groundCoordinateValue());
            }
            return finite(image->sample) && finite(image->line);
        }

    } // namespace

    namespace internal
    {

        EvaluationResult<ImageCoordinate> projectRpc(const RpcDefinition& definition,
                                                     const RpcCorrection& correction,
                                                     const GeodeticCoordinate& ground,
                                                     bool applyCorrection)
        {
            if (!finite(ground.longitudeDegrees) || !finite(ground.latitudeDegrees) || !finite(ground.heightMeters) ||
                ground.latitudeDegrees < -90.0 || ground.latitudeDegrees > 90.0)
            {
                return EvaluationResult<ImageCoordinate>::failure(CameraErrorCode::InvalidArgument,
                                                                  "RPC geodetic coordinate is invalid");
            }

            const RpcParameters& parameters = definition.parameters();
            const double normalized_longitude =
                longitudeDifference(ground.longitudeDegrees, parameters.longitudeOffset) / parameters.longitudeScale;
            const double normalized_latitude =
                (ground.latitudeDegrees - parameters.latitudeOffset) / parameters.latitudeScale;
            const double normalized_height = (ground.heightMeters - parameters.heightOffset) / parameters.heightScale;
            ImageCoordinate image;
            if (!evaluate(definition,
                          correction,
                          normalized_longitude,
                          normalized_latitude,
                          normalized_height,
                          applyCorrection,
                          &image))
            {
                return EvaluationResult<ImageCoordinate>::failure(
                    CameraErrorCode::OutsideModelDomain,
                    "RPC polynomial is singular or non-finite at the requested coordinate");
            }
            return EvaluationResult<ImageCoordinate>::success(image);
        }

        EvaluationResult<RpcProjectionJacobian> jacobianRpc(const RpcDefinition& definition,
                                                            const RpcCorrection& correction,
                                                            const GeodeticCoordinate& ground,
                                                            bool applyCorrection)
        {
            if (!finite(ground.longitudeDegrees) || !finite(ground.latitudeDegrees) ||
                !finite(ground.heightMeters) || ground.latitudeDegrees < -90.0 || ground.latitudeDegrees > 90.0)
            {
                return EvaluationResult<RpcProjectionJacobian>::failure(
                    CameraErrorCode::InvalidArgument, "RPC geodetic coordinate is invalid");
            }

            const RpcParameters& parameters = definition.parameters();
            const double normalized_longitude =
                longitudeDifference(ground.longitudeDegrees, parameters.longitudeOffset) / parameters.longitudeScale;
            const double normalized_latitude =
                (ground.latitudeDegrees - parameters.latitudeOffset) / parameters.latitudeScale;
            const double normalized_height =
                (ground.heightMeters - parameters.heightOffset) / parameters.heightScale;
            if (!finite(normalized_longitude) || !finite(normalized_latitude) || !finite(normalized_height))
            {
                return EvaluationResult<RpcProjectionJacobian>::failure(
                    CameraErrorCode::OutsideModelDomain, "RPC normalization produced a non-finite coordinate");
            }

            const RpcCoefficients terms =
                polynomialTerms(normalized_longitude, normalized_latitude, normalized_height);
            const auto derivatives =
                polynomialTermDerivatives(normalized_longitude, normalized_latitude, normalized_height);
            const double line_numerator = dot(parameters.lineNumerator, terms);
            const double line_denominator = dot(parameters.lineDenominator, terms);
            const double sample_numerator = dot(parameters.sampleNumerator, terms);
            const double sample_denominator = dot(parameters.sampleDenominator, terms);
            if (!finite(line_numerator) || !finite(line_denominator) || !finite(sample_numerator) ||
                !finite(sample_denominator) || std::abs(line_denominator) < kDenominatorEpsilon ||
                std::abs(sample_denominator) < kDenominatorEpsilon)
            {
                return EvaluationResult<RpcProjectionJacobian>::failure(
                    CameraErrorCode::OutsideModelDomain, "RPC polynomial is singular or non-finite");
            }

            std::array<double, 3> sample{};
            std::array<double, 3> line{};
            std::array<double, 3> normalized_sample{};
            std::array<double, 3> normalized_line{};
            const std::array<double, 3> normalizationScales{
                parameters.longitudeScale, parameters.latitudeScale, parameters.heightScale};
            for (int axis = 0; axis < 3; ++axis)
            {
                const double sample_derivative =
                    (derivativeDot(parameters.sampleNumerator, derivatives, axis) * sample_denominator -
                     sample_numerator * derivativeDot(parameters.sampleDenominator, derivatives, axis)) /
                    (sample_denominator * sample_denominator);
                const double line_derivative =
                    (derivativeDot(parameters.lineNumerator, derivatives, axis) * line_denominator -
                     line_numerator * derivativeDot(parameters.lineDenominator, derivatives, axis)) /
                    (line_denominator * line_denominator);
                normalized_sample[axis] = sample_derivative / normalizationScales[static_cast<std::size_t>(axis)];
                normalized_line[axis] = line_derivative / normalizationScales[static_cast<std::size_t>(axis)];
                sample[axis] = parameters.sampleScale * normalized_sample[axis];
                line[axis] = parameters.lineScale * normalized_line[axis];
            }

            if (applyCorrection)
            {
                if (const RpcImageCorrection* image_correction = correction.normalizedImageValue())
                {
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        sample[axis] += image_correction->sampleSamplePixels * normalized_sample[axis] +
                                        image_correction->sampleLinePixels * normalized_line[axis];
                        line[axis] += image_correction->lineSamplePixels * normalized_sample[axis] +
                                      image_correction->lineLinePixels * normalized_line[axis];
                    }
                }
                else if (const RpcGroundCorrection* ground_correction = correction.groundCoordinateValue())
                {
                    sample[0] += ground_correction->sampleLongitudePixelsPerDegree;
                    sample[1] += ground_correction->sampleLatitudePixelsPerDegree;
                    sample[2] += ground_correction->sampleHeightPixelsPerMeter;
                    line[0] += ground_correction->lineLongitudePixelsPerDegree;
                    line[1] += ground_correction->lineLatitudePixelsPerDegree;
                    line[2] += ground_correction->lineHeightPixelsPerMeter;
                }
            }

            for (int axis = 0; axis < 3; ++axis)
            {
                if (!finite(sample[axis]) || !finite(line[axis]))
                {
                    return EvaluationResult<RpcProjectionJacobian>::failure(
                        CameraErrorCode::OutsideModelDomain, "RPC projection Jacobian is non-finite");
                }
            }
            return EvaluationResult<RpcProjectionJacobian>::success(RpcProjectionJacobian{sample, line});
        }

        EvaluationResult<GeodeticCoordinate> invertRpcAtHeight(const RpcDefinition& definition,
                                                               const RpcCorrection& correction,
                                                               const ImageCoordinate& image,
                                                               double ellipsoidalHeightMeters,
                                                               const EvaluationOptions& options)
        {
            if (!finite(image.sample) || !finite(image.line) || !finite(ellipsoidalHeightMeters) ||
                !finite(options.desiredPrecisionPixels) || options.desiredPrecisionPixels <= 0.0 ||
                options.maximumIterations <= 0)
            {
                return EvaluationResult<GeodeticCoordinate>::failure(
                    CameraErrorCode::InvalidArgument,
                    "RPC inverse inputs and evaluation options must be finite and positive");
            }

            const RpcParameters& parameters = definition.parameters();
            const double normalized_height =
                (ellipsoidalHeightMeters - parameters.heightOffset) / parameters.heightScale;
            double normalized_longitude = 0.0;
            double normalized_latitude = 0.0;
            for (int iteration = 0; iteration < options.maximumIterations; ++iteration)
            {
                ImageCoordinate current;
                if (!evaluate(definition,
                              correction,
                              normalized_longitude,
                              normalized_latitude,
                              normalized_height,
                              true,
                              &current))
                {
                    break;
                }
                const double sample_residual = image.sample - current.sample;
                const double line_residual = image.line - current.line;
                const double current_error = std::hypot(sample_residual, line_residual);
                if (current_error <= options.desiredPrecisionPixels)
                {
                    const GeodeticCoordinate ground{
                        parameters.longitudeOffset + normalized_longitude * parameters.longitudeScale,
                        parameters.latitudeOffset + normalized_latitude * parameters.latitudeScale,
                        ellipsoidalHeightMeters};
                    if (ground.latitudeDegrees < -90.0 || ground.latitudeDegrees > 90.0)
                    {
                        break;
                    }
                    return EvaluationResult<GeodeticCoordinate>::success(ground, current_error);
                }

                const GeodeticCoordinate current_ground{
                    parameters.longitudeOffset + normalized_longitude * parameters.longitudeScale,
                    parameters.latitudeOffset + normalized_latitude * parameters.latitudeScale,
                    ellipsoidalHeightMeters};
                const auto jacobian = jacobianRpc(definition, correction, current_ground, true);
                if (!jacobian)
                {
                    break;
                }

                const double d_sample_d_longitude =
                    jacobian.value().sample[0] * parameters.longitudeScale;
                const double d_line_d_longitude = jacobian.value().line[0] * parameters.longitudeScale;
                const double d_sample_d_latitude = jacobian.value().sample[1] * parameters.latitudeScale;
                const double d_line_d_latitude = jacobian.value().line[1] * parameters.latitudeScale;
                const double determinant =
                    d_sample_d_longitude * d_line_d_latitude - d_sample_d_latitude * d_line_d_longitude;
                if (!finite(determinant) || std::abs(determinant) < 1.0e-12)
                {
                    break;
                }

                const double longitude_update =
                    (sample_residual * d_line_d_latitude - line_residual * d_sample_d_latitude) / determinant;
                const double latitude_update =
                    (d_sample_d_longitude * line_residual - d_line_d_longitude * sample_residual) / determinant;
                bool accepted = false;
                double step_scale = 1.0;
                for (int line_search = 0; line_search < 10; ++line_search)
                {
                    const double trial_longitude = normalized_longitude + step_scale * longitude_update;
                    const double trial_latitude = normalized_latitude + step_scale * latitude_update;
                    ImageCoordinate trial;
                    if (evaluate(
                            definition, correction, trial_longitude, trial_latitude, normalized_height, true, &trial) &&
                        std::hypot(image.sample - trial.sample, image.line - trial.line) < current_error)
                    {
                        normalized_longitude = trial_longitude;
                        normalized_latitude = trial_latitude;
                        accepted = true;
                        break;
                    }
                    step_scale *= 0.5;
                }
                if (!accepted || !finite(normalized_longitude) || !finite(normalized_latitude))
                {
                    break;
                }
            }

            return EvaluationResult<GeodeticCoordinate>::failure(
                CameraErrorCode::NonConvergence, "RPC inverse did not reach the requested pixel precision");
        }

    } // namespace internal

} // namespace placamera
