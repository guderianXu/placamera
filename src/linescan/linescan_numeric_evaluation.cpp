#include "placamera/linescan_numeric_state.h"

#include "internal/linescan_numeric_projection.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <utility>

namespace placamera
{
    EvaluationResult<TimeReference> LineScanNumericState::timeForLine(double line) const
    {
        if (!std::isfinite(line))
        {
            return EvaluationResult<TimeReference>::failure(CameraErrorCode::InvalidArgument,
                                                            "line coordinate must be finite");
        }
        const double minimum = _pixelConvention == LineScanPixelConvention::PixelCenter ? 0.5 : 0.0;
        if (line < minimum || line > minimum + static_cast<double>(_imageSize.lines - 1))
        {
            return EvaluationResult<TimeReference>::failure(CameraErrorCode::OutsideModelDomain,
                                                            "line coordinate lies outside the image");
        }
        double seconds = 0.0;
        if (_timing.segments.empty())
        {
            seconds = _timing.startTimeSeconds + (line - _timing.lineZero) * _timing.secondsPerLine;
        }
        else
        {
            auto segment = std::upper_bound(_timing.segments.begin(),
                                            _timing.segments.end(),
                                            line,
                                            [](double candidate, const LineRateSegment& value)
                                            { return candidate < value.startLine; });
            if (segment != _timing.segments.begin())
            {
                --segment;
            }
            seconds = segment->startTimeSeconds + (line - segment->startLine) * segment->secondsPerLine;
        }
        if (!std::isfinite(seconds))
        {
            return EvaluationResult<TimeReference>::failure(CameraErrorCode::OutsideModelDomain,
                                                            "line timing produced a non-finite time");
        }
        return EvaluationResult<TimeReference>::success(TimeReference{_timing.timeScale, seconds});
    }

    EvaluationResult<LineScanProjectionDetails> LineScanNumericState::projectAtLine(
        const GroundCoordinate& ground, double line, const EvaluationOptions& options) const
    {
        return internal::projectLineScanNumericAtLine(*this, ground, line, options);
    }

    EvaluationResult<Projection> LineScanNumericState::groundToImage(const GroundCoordinate& ground,
                                                                     const EvaluationOptions& options) const
    {
        return internal::projectLineScanNumeric(*this, ground, options);
    }

    EvaluationResult<ImagingLocus> LineScanNumericState::imageToImagingLocus(const ImageCoordinate& image,
                                                                             const EvaluationOptions& options) const
    {
        return internal::lineScanNumericImagingLocus(*this, image, options);
    }

    Result<CameraModelPtr<LineScanModel>>
    LineScanNumericState::toModel(CameraInstanceId resultInstanceId,
                                  std::optional<CameraDefinitionId> resultDefinitionId) const
    {
        if (_definitionDirty && !resultDefinitionId)
        {
            return Result<CameraModelPtr<LineScanModel>>::failure(
                CameraErrorCode::InvalidArgument,
                "modified line-scan detector optics require an explicit result definition identifier");
        }
        try
        {
            const CameraDefinitionId definition_id = resultDefinitionId ? *resultDefinitionId : _definitionId;
            const auto definition = LineScanDefinition::create(definition_id, _groundFrame, _optics, _pixelConvention);
            auto model = LineScanModel::create(std::move(resultInstanceId),
                                               _imageId,
                                               definition,
                                               _imageSize,
                                               LineScanTrajectory::create(_samples),
                                               _timing,
                                               _bias,
                                               _captureTime,
                                               _timeOffsetPrior);
            return Result<CameraModelPtr<LineScanModel>>::success(
                std::make_shared<const LineScanModel>(std::move(model)));
        }
        catch (const CameraValidationError& error)
        {
            return Result<CameraModelPtr<LineScanModel>>::failure(error.code(), error.what());
        }
        catch (const std::exception& error)
        {
            return Result<CameraModelPtr<LineScanModel>>::failure(CameraErrorCode::InvalidModelState, error.what());
        }
    }

} // namespace placamera
