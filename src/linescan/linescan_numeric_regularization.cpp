#include "placamera/linescan_numeric_state.h"

#include <cmath>

namespace placamera
{
    std::size_t LineScanNumericState::regularizationResidualCount() const noexcept
    {
        std::size_t count = 0;
        for (const auto& sample : _samples)
        {
            if (_selection.knotPositions && !sample.constraints.positionFixed && sample.constraints.positionSigmaMeters)
            {
                count += 3;
            }
            if (_selection.knotRotations && !sample.constraints.rotationFixed &&
                sample.constraints.rotationSigmaRadians)
            {
                count += 3;
            }
        }
        if (_selection.timeOffset && _timeOffsetPrior)
        {
            ++count;
        }
        for (std::size_t index = 1; index + 1 < _samples.size(); ++index)
        {
            if (_selection.positionSecondDifferenceWeight > 0.0 && _selection.knotPositions &&
                (!_samples[index - 1].constraints.positionFixed || !_samples[index].constraints.positionFixed ||
                 !_samples[index + 1].constraints.positionFixed))
            {
                count += 3;
            }
            if (_selection.rotationSecondDifferenceWeight > 0.0 && _selection.knotRotations &&
                (!_samples[index - 1].constraints.rotationFixed || !_samples[index].constraints.rotationFixed ||
                 !_samples[index + 1].constraints.rotationFixed))
            {
                count += 3;
            }
        }
        return count;
    }

    Result<void> LineScanNumericState::writeRegularizationResiduals(std::span<double> residuals) const
    {
        if (residuals.size() != regularizationResidualCount())
        {
            return Result<void>::failure(CameraErrorCode::InvalidArgument,
                                         "line-scan regularization buffer size is incorrect");
        }
        std::size_t cursor = 0;
        for (std::size_t index = 0; index < _samples.size(); ++index)
        {
            const auto& constraints = _samples[index].constraints;
            if (_selection.knotPositions && !constraints.positionFixed && constraints.positionSigmaMeters)
            {
                for (std::size_t axis = 0; axis < 3; ++axis)
                {
                    residuals[cursor++] = _positionDeltas[index][axis] / (*constraints.positionSigmaMeters)[axis];
                }
            }
            if (_selection.knotRotations && !constraints.rotationFixed && constraints.rotationSigmaRadians)
            {
                for (std::size_t axis = 0; axis < 3; ++axis)
                {
                    residuals[cursor++] = _rotationDeltas[index][axis] / (*constraints.rotationSigmaRadians)[axis];
                }
            }
        }
        if (_selection.timeOffset && _timeOffsetPrior)
        {
            residuals[cursor++] =
                (_bias.timeOffsetSeconds - _timeOffsetPrior->meanSeconds) / _timeOffsetPrior->sigmaSeconds;
        }

        const double position_scale = std::sqrt(_selection.positionSecondDifferenceWeight);
        const double rotation_scale = std::sqrt(_selection.rotationSecondDifferenceWeight);
        for (std::size_t index = 1; index + 1 < _samples.size(); ++index)
        {
            if (position_scale > 0.0 && _selection.knotPositions &&
                (!_samples[index - 1].constraints.positionFixed || !_samples[index].constraints.positionFixed ||
                 !_samples[index + 1].constraints.positionFixed))
            {
                for (std::size_t axis = 0; axis < 3; ++axis)
                {
                    residuals[cursor++] =
                        position_scale * (_positionDeltas[index - 1][axis] - 2.0 * _positionDeltas[index][axis] +
                                          _positionDeltas[index + 1][axis]);
                }
            }
            if (rotation_scale > 0.0 && _selection.knotRotations &&
                (!_samples[index - 1].constraints.rotationFixed || !_samples[index].constraints.rotationFixed ||
                 !_samples[index + 1].constraints.rotationFixed))
            {
                for (std::size_t axis = 0; axis < 3; ++axis)
                {
                    residuals[cursor++] =
                        rotation_scale * (_rotationDeltas[index - 1][axis] - 2.0 * _rotationDeltas[index][axis] +
                                          _rotationDeltas[index + 1][axis]);
                }
            }
        }
        return Result<void>::success();
    }

} // namespace placamera
