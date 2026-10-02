#include "placamera/bundle_linearization.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace placamera
{
    namespace
    {

        constexpr double kPointStepMeters = 1.0e-6;
        constexpr double kRotationStepRadians = 1.0e-7;
        constexpr double kTranslationStepMeters = 1.0e-6;
        constexpr double kIntrinsicStep = 1.0e-5;
        constexpr double kDistortionStep = 1.0e-7;
        constexpr double kTimeStepSeconds = 1.0e-7;

        double stepFor(const OptimizationParameterBlock& block) noexcept
        {
            switch (block.kind)
            {
            case OptimizationParameterKind::RotationVector:
                return kRotationStepRadians;
            case OptimizationParameterKind::Translation:
                return kTranslationStepMeters;
            case OptimizationParameterKind::Intrinsics:
                return kIntrinsicStep;
            case OptimizationParameterKind::Distortion:
                return kDistortionStep;
            case OptimizationParameterKind::TimeOffset:
                return kTimeStepSeconds;
            case OptimizationParameterKind::ImageCorrection:
                return kIntrinsicStep;
            }
            return kDistortionStep;
        }

        EvaluationResult<Projection>
        project(const LineScanNumericState& state, const GroundCoordinate& ground, double line)
        {
            const auto result = state.projectAtLine(ground, line);
            if (!result)
            {
                return EvaluationResult<Projection>::failure(result.errorCode(), result.message());
            }
            return EvaluationResult<Projection>::success(result.value().projection, result.achievedPrecisionPixels());
        }

        EvaluationResult<Projection> projectWithDelta(const LineScanNumericState& state,
                                                      const GroundCoordinate& ground,
                                                      double line,
                                                      std::size_t column,
                                                      double delta)
        {
            LineScanNumericState perturbed = state;
            std::vector<double> update(state.optimizationLayout().parameterCount(), 0.0);
            update[column] = delta;
            const auto applied = perturbed.applyOptimizationDelta(update);
            if (!applied)
            {
                return EvaluationResult<Projection>::failure(applied.errorCode(), applied.message());
            }
            return project(perturbed, ground, line);
        }

        EvaluationResult<BundleProjectionLinearization> linearizeState(const LineScanNumericState& state,
                                                                       const LineScanModel& model,
                                                                       const GroundCoordinate& ground,
                                                                       double line)
        {
            if (ground.frame != state.groundFrame())
            {
                return EvaluationResult<BundleProjectionLinearization>::failure(
                    CameraErrorCode::FrameMismatch, "bundle point frame does not match the line-scan camera");
            }
            const auto base = project(state, ground, line);
            if (!base)
            {
                return EvaluationResult<BundleProjectionLinearization>::failure(base.errorCode(), base.message());
            }

            BundleProjectionLinearization output;
            output.projection = base.value();
            output.layout = state.optimizationLayout();
            output.parameterJacobian.assign(2U * output.layout.parameterCount(), 0.0);

            for (std::size_t axis = 0; axis < ground.position.size(); ++axis)
            {
                GroundCoordinate plus_point = ground;
                GroundCoordinate minus_point = ground;
                const double step = kPointStepMeters * std::max(1.0, std::abs(ground.position[axis]));
                plus_point.position[axis] += step;
                minus_point.position[axis] -= step;
                const auto plus = project(state, plus_point, line);
                const auto minus = project(state, minus_point, line);
                if (!plus || !minus)
                {
                    const auto& failed = !plus ? plus : minus;
                    return EvaluationResult<BundleProjectionLinearization>::failure(failed.errorCode(),
                                                                                    failed.message());
                }
                output.pointJacobian[axis] = (plus.value().image.sample - minus.value().image.sample) / (2.0 * step);
                output.pointJacobian[3U + axis] = (plus.value().image.line - minus.value().image.line) / (2.0 * step);
            }

            const std::size_t parameter_count = output.layout.parameterCount();
            for (const auto& block : output.layout.blocks)
            {
                const double step = stepFor(block);
                for (std::size_t component = 0; component < block.size; ++component)
                {
                    const std::size_t column = block.offset + component;
                    const auto plus = projectWithDelta(state, ground, line, column, step);
                    const auto minus = projectWithDelta(state, ground, line, column, -step);
                    if (!plus || !minus)
                    {
                        const auto& failed = !plus ? plus : minus;
                        return EvaluationResult<BundleProjectionLinearization>::failure(failed.errorCode(),
                                                                                        failed.message());
                    }
                    output.parameterJacobian[column] =
                        (plus.value().image.sample - minus.value().image.sample) / (2.0 * step);
                    output.parameterJacobian[parameter_count + column] =
                        (plus.value().image.line - minus.value().image.line) / (2.0 * step);
                }
            }

            LineScanOptimizationSelection pose_selection;
            pose_selection.knotPositions = false;
            pose_selection.knotRotations = false;
            pose_selection.globalTranslation = true;
            pose_selection.globalRotation = true;
            pose_selection.timeOffset = false;
            const auto pose_state_result = LineScanNumericState::fromModel(model, pose_selection);
            if (!pose_state_result)
            {
                return EvaluationResult<BundleProjectionLinearization>::failure(pose_state_result.errorCode(),
                                                                                pose_state_result.message());
            }
            {
                const auto& pose_state = pose_state_result.value();
                for (const auto& block : pose_state.optimizationLayout().blocks)
                {
                    for (std::size_t component = 0; component < block.size; ++component)
                    {
                        const std::size_t column = block.offset + component;
                        const double step = stepFor(block);
                        const auto plus = projectWithDelta(pose_state, ground, line, column, step);
                        const auto minus = projectWithDelta(pose_state, ground, line, column, -step);
                        if (!plus || !minus)
                        {
                            const auto& failed = !plus ? plus : minus;
                            return EvaluationResult<BundleProjectionLinearization>::failure(failed.errorCode(),
                                                                                            failed.message());
                        }
                        const double sample_derivative =
                            (plus.value().image.sample - minus.value().image.sample) / (2.0 * step);
                        const double line_derivative =
                            (plus.value().image.line - minus.value().image.line) / (2.0 * step);
                        const std::size_t pose_offset = block.name == "trajectory.rotation" ? 0U : 3U;
                        output.poseJacobian[pose_offset + component] = sample_derivative;
                        output.poseJacobian[6U + pose_offset + component] = line_derivative;
                    }
                }
            }
            return EvaluationResult<BundleProjectionLinearization>::success(std::move(output));
        }

    } // namespace

    EvaluationResult<BundleProjectionLinearization> linearizeForBundle(const LineScanModel& model,
                                                                       const GroundCoordinate& ground,
                                                                       double observationLine,
                                                                       const LineScanOptimizationSelection& selection)
    {
        if (!std::isfinite(observationLine))
        {
            return EvaluationResult<BundleProjectionLinearization>::failure(
                CameraErrorCode::InvalidArgument, "line-scan observation line must be finite");
        }
        const auto state = LineScanNumericState::fromModel(model, selection);
        if (!state)
        {
            return EvaluationResult<BundleProjectionLinearization>::failure(state.errorCode(), state.message());
        }
        return linearizeState(state.value(), model, ground, observationLine);
    }

} // namespace placamera
