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
        constexpr double kIntrinsicStepPixels = 1.0e-5;
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
                return kIntrinsicStepPixels;
            case OptimizationParameterKind::Distortion:
                return kDistortionStep;
            case OptimizationParameterKind::TimeOffset:
                return kTimeStepSeconds;
            case OptimizationParameterKind::ImageCorrection:
                return kIntrinsicStepPixels;
            }
            return kDistortionStep;
        }

        EvaluationResult<Projection> project(const FramePinholeModel& model, const GroundCoordinate& ground)
        {
            return model.groundToImage(ground);
        }

        EvaluationResult<Projection> perturbedProjection(const FramePinholeModel& model,
                                                         const GroundCoordinate& ground,
                                                         const OptimizationLayout& layout,
                                                         const FrameOptimizationSelection& selection,
                                                         std::size_t column,
                                                         double delta)
        {
            std::vector<double> update(layout.parameterCount(), 0.0);
            update[column] = delta;
            OptimizationUpdate optimization_update{model.instanceId(), model.definitionId(), std::move(update)};
            const auto perturbed = model.withOptimizationUpdate(optimization_update, selection);
            if (!perturbed)
            {
                return EvaluationResult<Projection>::failure(perturbed.errorCode(), perturbed.message());
            }
            return perturbed.value()->groundToImage(ground);
        }

        EvaluationResult<BundleProjectionLinearization> linearizeFrame(const FramePinholeModel& model,
                                                                       const GroundCoordinate& ground,
                                                                       const FrameOptimizationSelection& selection)
        {
            if (ground.frame != model.groundFrame())
            {
                return EvaluationResult<BundleProjectionLinearization>::failure(
                    CameraErrorCode::FrameMismatch, "bundle point frame does not match the frame camera");
            }
            const auto base = project(model, ground);
            if (!base)
            {
                return EvaluationResult<BundleProjectionLinearization>::failure(base.errorCode(), base.message());
            }

            BundleProjectionLinearization output;
            output.projection = base.value();
            output.layout = model.optimizationLayout(selection);
            output.parameterJacobian.assign(2U * output.layout.parameterCount(), 0.0);

            for (std::size_t axis = 0; axis < ground.position.size(); ++axis)
            {
                GroundCoordinate plus_point = ground;
                GroundCoordinate minus_point = ground;
                const double step = kPointStepMeters * std::max(1.0, std::abs(ground.position[axis]));
                plus_point.position[axis] += step;
                minus_point.position[axis] -= step;
                const auto plus = project(model, plus_point);
                const auto minus = project(model, minus_point);
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
                    const auto plus = perturbedProjection(model, ground, output.layout, selection, column, step);
                    const auto minus = perturbedProjection(model, ground, output.layout, selection, column, -step);
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

                    if (block.name == "pose.rotation" && component < 3U)
                    {
                        output.poseJacobian[component] = output.parameterJacobian[column];
                        output.poseJacobian[6U + component] = output.parameterJacobian[parameter_count + column];
                    }
                    else if (block.name == "pose.translation" && component < 3U)
                    {
                        output.poseJacobian[3U + component] = output.parameterJacobian[column];
                        output.poseJacobian[9U + component] = output.parameterJacobian[parameter_count + column];
                    }
                }
            }
            return EvaluationResult<BundleProjectionLinearization>::success(std::move(output));
        }

    } // namespace

    EvaluationResult<BundleProjectionLinearization> linearizeForBundle(const FramePinholeModel& model,
                                                                       const GroundCoordinate& ground,
                                                                       const FrameOptimizationSelection& selection)
    {
        return linearizeFrame(model, ground, selection);
    }

} // namespace placamera
