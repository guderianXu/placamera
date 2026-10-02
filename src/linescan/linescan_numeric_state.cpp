#include "placamera/linescan_numeric_state.h"

#include <cmath>
#include <string>
#include <utility>

namespace placamera
{
    namespace
    {

        void addBlock(OptimizationLayout* layout,
                      std::string name,
                      OptimizationParameterKind kind,
                      std::string unit,
                      std::size_t size,
                      bool affectsDefinition)
        {
            const std::size_t offset = layout->parameterCount();
            layout->blocks.push_back({std::move(name), kind, std::move(unit), offset, size, affectsDefinition});
        }

        Result<OptimizationLayout> makeLayout(const std::vector<TrajectorySample>& samples,
                                              const LineScanOptics& optics,
                                              const LineScanOptimizationSelection& selection)
        {
            if (!std::isfinite(selection.positionSecondDifferenceWeight) ||
                !std::isfinite(selection.rotationSecondDifferenceWeight) ||
                selection.positionSecondDifferenceWeight < 0.0 || selection.rotationSecondDifferenceWeight < 0.0)
            {
                return Result<OptimizationLayout>::failure(
                    CameraErrorCode::InvalidArgument, "line-scan smoothness weights must be finite and non-negative");
            }
            if (selection.detector.anyDetectorGeometryParameter() && !optics.detectorGeometry)
            {
                return Result<OptimizationLayout>::failure(
                    CameraErrorCode::InvalidArgument,
                    "line-scan detector-geometry parameters require detector geometry in the definition");
            }
            if (selection.calibration.any() && !optics.completeCalibration)
            {
                return Result<OptimizationLayout>::failure(
                    CameraErrorCode::InvalidArgument,
                    "complete line-scan calibration parameters require complete calibration in the definition");
            }

            OptimizationLayout layout;
            for (std::size_t index = 0; index < samples.size(); ++index)
            {
                const auto& constraints = samples[index].constraints;
                const std::string prefix = "trajectory.knots." + std::to_string(index);
                if (selection.knotPositions && !constraints.positionFixed)
                {
                    addBlock(&layout, prefix + ".position", OptimizationParameterKind::Translation, "m", 3, false);
                }
                if (selection.knotRotations && !constraints.rotationFixed)
                {
                    addBlock(&layout, prefix + ".rotation", OptimizationParameterKind::RotationVector, "rad", 3, false);
                }
            }
            if (selection.globalTranslation)
            {
                addBlock(&layout, "trajectory.translation", OptimizationParameterKind::Translation, "m", 3, false);
            }
            if (selection.globalRotation)
            {
                addBlock(&layout, "trajectory.rotation", OptimizationParameterKind::RotationVector, "rad", 3, false);
            }
            if (selection.timeOffset)
            {
                addBlock(&layout, "trajectory.time", OptimizationParameterKind::TimeOffset, "s", 1, false);
            }

            const auto& calibration = selection.calibration;
            if (calibration.f)
            {
                addBlock(&layout, "calibration.f", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (calibration.cx)
            {
                addBlock(&layout, "calibration.cx", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (calibration.cy)
            {
                addBlock(&layout, "calibration.cy", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (calibration.b1)
            {
                addBlock(&layout, "calibration.b1", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (calibration.b2)
            {
                addBlock(&layout, "calibration.b2", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (calibration.k1)
            {
                addBlock(&layout, "calibration.k1", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.k2)
            {
                addBlock(&layout, "calibration.k2", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.k3)
            {
                addBlock(&layout, "calibration.k3", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.k4)
            {
                addBlock(&layout, "calibration.k4", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.p1)
            {
                addBlock(&layout, "calibration.p1", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.p2)
            {
                addBlock(&layout, "calibration.p2", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.p3)
            {
                addBlock(&layout, "calibration.p3", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }
            if (calibration.p4)
            {
                addBlock(&layout, "calibration.p4", OptimizationParameterKind::Distortion, "normalized", 1, true);
            }

            const auto& mask = selection.detector;
            if (mask.focalLength)
            {
                addBlock(&layout, "detector.focal_length", OptimizationParameterKind::Intrinsics, "mm", 1, true);
            }
            if (mask.distortionK1)
            {
                addBlock(&layout, "detector.distortion_k1", OptimizationParameterKind::Distortion, "1", 1, true);
            }
            if (mask.samplePitch)
            {
                addBlock(&layout, "detector.sample_pitch", OptimizationParameterKind::Intrinsics, "mm", 1, true);
            }
            if (mask.principalSample)
            {
                addBlock(&layout, "detector.principal_sample", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (mask.detectorSampleSumming)
            {
                addBlock(&layout, "detector.sample_summing", OptimizationParameterKind::Intrinsics, "1", 1, true);
            }
            if (mask.detectorLineSumming)
            {
                addBlock(&layout, "detector.line_summing", OptimizationParameterKind::Intrinsics, "1", 1, true);
            }
            if (mask.detectorSampleOrigin)
            {
                addBlock(&layout, "detector.sample_origin", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (mask.detectorLineOrigin)
            {
                addBlock(&layout, "detector.line_origin", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (mask.startingDetectorSample)
            {
                addBlock(&layout, "detector.starting_sample", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (mask.startingDetectorLine)
            {
                addBlock(&layout, "detector.starting_line", OptimizationParameterKind::Intrinsics, "px", 1, true);
            }
            if (mask.focalToPixelSamples)
            {
                addBlock(&layout,
                         "detector.focal_to_pixel_samples",
                         OptimizationParameterKind::Intrinsics,
                         "mixed",
                         3,
                         true);
            }
            if (mask.focalToPixelLines)
            {
                addBlock(
                    &layout, "detector.focal_to_pixel_lines", OptimizationParameterKind::Intrinsics, "mixed", 3, true);
            }
            return Result<OptimizationLayout>::success(std::move(layout));
        }

    } // namespace

    bool LineScanCalibrationOptimizationMask::any() const noexcept
    {
        return f || cx || cy || b1 || b2 || k1 || k2 || k3 || k4 || p1 || p2 || p3 || p4;
    }

    bool LineScanDetectorOptimizationMask::anyDetectorGeometryParameter() const noexcept
    {
        return detectorSampleSumming || detectorLineSumming || detectorSampleOrigin || detectorLineOrigin ||
               startingDetectorSample || startingDetectorLine || focalToPixelSamples || focalToPixelLines;
    }

    bool LineScanDetectorOptimizationMask::any() const noexcept
    {
        return focalLength || distortionK1 || samplePitch || principalSample || anyDetectorGeometryParameter();
    }

    Result<LineScanNumericState> LineScanNumericState::fromModel(const LineScanModel& model,
                                                                 LineScanOptimizationSelection selection)
    {
        if (model.trajectory().frameComposed())
        {
            return Result<LineScanNumericState>::failure(
                CameraErrorCode::UnsupportedModel,
                "per-knot line-scan optimization requires a direct-sample trajectory");
        }
        const auto layout = makeLayout(model.trajectory().samples(), model.lineScanDefinition().optics(), selection);
        if (!layout)
        {
            return Result<LineScanNumericState>::failure(layout.error());
        }
        return Result<LineScanNumericState>::success(LineScanNumericState(model.instanceId(),
                                                                          model.definitionId(),
                                                                          model.imageId(),
                                                                          model.groundFrame(),
                                                                          model.imageSize(),
                                                                          model.captureTime(),
                                                                          model.trajectory().samples(),
                                                                          model.lineTiming(),
                                                                          model.trajectoryBias(),
                                                                          model.timeOffsetPrior(),
                                                                          model.lineScanDefinition().optics(),
                                                                          model.lineScanDefinition().pixelConvention(),
                                                                          std::move(selection),
                                                                          layout.value()));
    }

    const CameraInstanceId& LineScanNumericState::instanceId() const noexcept
    {
        return _instanceId;
    }
    const CameraDefinitionId& LineScanNumericState::definitionId() const noexcept
    {
        return _definitionId;
    }
    const ImageId& LineScanNumericState::imageId() const noexcept
    {
        return _imageId;
    }
    const FrameId& LineScanNumericState::groundFrame() const noexcept
    {
        return _groundFrame;
    }
    const ImageSize& LineScanNumericState::imageSize() const noexcept
    {
        return _imageSize;
    }
    const std::optional<TimeReference>& LineScanNumericState::captureTime() const noexcept
    {
        return _captureTime;
    }
    const std::vector<TrajectorySample>& LineScanNumericState::trajectorySamples() const noexcept
    {
        return _samples;
    }
    const LineTiming& LineScanNumericState::lineTiming() const noexcept
    {
        return _timing;
    }
    const LineScanTrajectoryBias& LineScanNumericState::trajectoryBias() const noexcept
    {
        return _bias;
    }
    const std::optional<LineScanTimeOffsetPrior>& LineScanNumericState::timeOffsetPrior() const noexcept
    {
        return _timeOffsetPrior;
    }
    const LineScanOptics& LineScanNumericState::optics() const noexcept
    {
        return _optics;
    }
    LineScanPixelConvention LineScanNumericState::pixelConvention() const noexcept
    {
        return _pixelConvention;
    }
    const LineScanOptimizationSelection& LineScanNumericState::selection() const noexcept
    {
        return _selection;
    }
    const OptimizationLayout& LineScanNumericState::optimizationLayout() const noexcept
    {
        return _layout;
    }
    bool LineScanNumericState::definitionDirty() const noexcept
    {
        return _definitionDirty;
    }

    LineScanNumericState::LineScanNumericState(CameraInstanceId instanceId,
                                               CameraDefinitionId definitionId,
                                               ImageId imageId,
                                               FrameId groundFrame,
                                               ImageSize imageSize,
                                               std::optional<TimeReference> captureTime,
                                               std::vector<TrajectorySample> samples,
                                               LineTiming timing,
                                               LineScanTrajectoryBias bias,
                                               std::optional<LineScanTimeOffsetPrior> timeOffsetPrior,
                                               LineScanOptics optics,
                                               LineScanPixelConvention pixelConvention,
                                               LineScanOptimizationSelection selection,
                                               OptimizationLayout layout)
        : _instanceId(std::move(instanceId)), _definitionId(std::move(definitionId)), _imageId(std::move(imageId)),
          _groundFrame(std::move(groundFrame)), _imageSize(imageSize), _captureTime(captureTime),
          _nominalSamples(samples), _samples(std::move(samples)), _positionDeltas(_samples.size()),
          _rotationDeltas(_samples.size()), _timing(std::move(timing)), _bias(bias), _timeOffsetPrior(timeOffsetPrior),
          _optics(std::move(optics)), _pixelConvention(pixelConvention), _selection(std::move(selection)),
          _layout(std::move(layout))
    {
    }

} // namespace placamera
