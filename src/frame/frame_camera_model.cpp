#include "placamera/frame_camera.h"

#include "internal/frame_camera_math.h"

#include <cmath>
#include <exception>
#include <memory>
#include <utility>

namespace placamera
{
    FramePinholeModel FramePinholeModel::create(CameraInstanceId instanceId,
                                                ImageId imageId,
                                                std::shared_ptr<const FramePinholeDefinition> definition,
                                                ImageSize imageSize,
                                                Pose pose,
                                                std::optional<TimeReference> captureTime,
                                                CameraAcquisitionState acquisition)
    {
        if (!definition)
        {
            throw CameraValidationError(CameraErrorCode::InvalidFrame, "pinhole model requires a definition");
        }
        if (!imageSize.isValid())
        {
            throw CameraValidationError(CameraErrorCode::InvalidImageSize, "pinhole model image size must be positive");
        }
        if (pose.frame != definition->groundFrame())
        {
            throw CameraValidationError(CameraErrorCode::FrameMismatch,
                                        "pinhole pose frame must match the definition ground frame");
        }
        pose = Pose::create(pose.frame, pose.center, pose.cameraToWorldRotation);
        if (captureTime && !std::isfinite(captureTime->seconds))
        {
            throw CameraValidationError(CameraErrorCode::InvalidTime,
                                        "pinhole capture time must contain finite seconds");
        }
        const auto valid_acquisition = validateCameraAcquisition(acquisition, &instanceId);
        if (!valid_acquisition)
        {
            throw CameraValidationError(valid_acquisition.errorCode(), valid_acquisition.message());
        }
        return FramePinholeModel(std::move(instanceId),
                                 std::move(imageId),
                                 std::move(definition),
                                 imageSize,
                                 std::move(pose),
                                 captureTime,
                                 std::move(acquisition));
    }

    const CameraInstanceId& FramePinholeModel::instanceId() const noexcept
    {
        return _instanceId;
    }

    const CameraDefinition& FramePinholeModel::definition() const noexcept
    {
        return *_definition;
    }

    const CameraDefinitionId& FramePinholeModel::definitionId() const noexcept
    {
        return _definition->definitionId();
    }

    const ImageId& FramePinholeModel::imageId() const noexcept
    {
        return _imageId;
    }

    std::string_view FramePinholeModel::modelType() const noexcept
    {
        return _definition->modelType();
    }

    int FramePinholeModel::parameterSchemaVersion() const noexcept
    {
        return _definition->parameterSchemaVersion();
    }

    const FrameId& FramePinholeModel::groundFrame() const noexcept
    {
        return _definition->groundFrame();
    }

    const ImageSize& FramePinholeModel::imageSize() const noexcept
    {
        return _imageSize;
    }

    const std::optional<TimeReference>& FramePinholeModel::captureTime() const noexcept
    {
        return _captureTime;
    }

    CapabilitySet FramePinholeModel::capabilities() const noexcept
    {
        return _definition->capabilities();
    }

    const FramePinholeDefinition& FramePinholeModel::pinholeDefinition() const noexcept
    {
        return *_definition;
    }

    const Pose& FramePinholeModel::pose() const noexcept
    {
        return _pose;
    }

    const CameraAcquisitionState& FramePinholeModel::acquisition() const noexcept
    {
        return _acquisition;
    }

    EvaluationResult<double> FramePinholeModel::signedDepth(const GroundCoordinate& ground) const
    {
        return internal::signedFrameDepth(groundFrame(), _definition->depthAxisFlipped(), _pose, ground);
    }

    EvaluationResult<Projection> FramePinholeModel::groundToImageSigned(const GroundCoordinate& ground,
                                                                        const EvaluationOptions& options) const
    {
        return internal::projectFrameSigned(groundFrame(),
                                            _imageSize,
                                            _captureTime,
                                            _definition->intrinsics(),
                                            _definition->distortion(),
                                            _definition->projectionModel(),
                                            _definition->depthAxisFlipped(),
                                            _pose,
                                            _acquisition,
                                            ground,
                                            options);
    }

    EvaluationResult<Projection> FramePinholeModel::groundToImage(const GroundCoordinate& ground,
                                                                  const EvaluationOptions& options) const
    {
        return internal::projectFrame(groundFrame(),
                                      _imageSize,
                                      _captureTime,
                                      _definition->intrinsics(),
                                      _definition->distortion(),
                                      _definition->projectionModel(),
                                      _definition->depthAxisFlipped(),
                                      _pose,
                                      _acquisition,
                                      ground,
                                      options);
    }

    EvaluationResult<Projection> FramePinholeModel::groundToImageAtLine(const GroundCoordinate& ground,
                                                                        double observationLine,
                                                                        const EvaluationOptions& options) const
    {
        return internal::projectFrameAtLine(groundFrame(),
                                            _imageSize,
                                            _captureTime,
                                            _definition->intrinsics(),
                                            _definition->distortion(),
                                            _definition->projectionModel(),
                                            _definition->depthAxisFlipped(),
                                            _pose,
                                            _acquisition,
                                            ground,
                                            observationLine,
                                            options);
    }

    EvaluationResult<ImagingLocus> FramePinholeModel::imageToImagingLocus(const ImageCoordinate& image,
                                                                          const EvaluationOptions& options) const
    {
        return internal::frameImagingLocus(groundFrame(),
                                           _imageSize,
                                           _captureTime,
                                           _definition->intrinsics(),
                                           _definition->distortion(),
                                           _definition->principalPointDecomposition(),
                                           _definition->projectionModel(),
                                           _definition->depthAxisFlipped(),
                                           _pose,
                                           _acquisition,
                                           image,
                                           options);
    }

    EvaluationResult<GroundCoordinate> FramePinholeModel::imageToGroundAtDepth(const ImageCoordinate& image,
                                                                               double positiveDepth,
                                                                               const EvaluationOptions& options) const
    {
        return internal::frameGroundAtDepth(groundFrame(),
                                            _imageSize,
                                            _definition->intrinsics(),
                                            _definition->distortion(),
                                            _definition->principalPointDecomposition(),
                                            _definition->projectionModel(),
                                            _definition->depthAxisFlipped(),
                                            _pose,
                                            _acquisition,
                                            image,
                                            positiveDepth,
                                            options);
    }

    FramePinholeModel FramePinholeModel::withPose(CameraInstanceId instanceId, Pose pose) const
    {
        CameraAcquisitionState acquisition = _acquisition;
        if (acquisition.masterCameraId && *acquisition.masterCameraId == instanceId)
        {
            throw CameraValidationError(CameraErrorCode::InvalidArgument,
                                        "rebound frame instance cannot equal its master camera identity");
        }
        return create(std::move(instanceId),
                      _imageId,
                      _definition,
                      _imageSize,
                      std::move(pose),
                      _captureTime,
                      std::move(acquisition));
    }

    FramePinholeModel FramePinholeModel::withImageSize(ImageSize imageSize) const
    {
        return create(_instanceId, _imageId, _definition, imageSize, _pose, _captureTime, _acquisition);
    }

    FramePinholeModel FramePinholeModel::normalizedForPositiveDepth(CameraDefinitionId definitionId,
                                                                    CameraInstanceId instanceId) const
    {
        const std::array<double, 3> axis_signs = _definition->positiveDepthAxisSigns();
        const Pose normalized_pose =
            Pose::create(_pose.frame,
                         _pose.center,
                         internal::normalizedCameraToWorldRotation(_pose.cameraToWorldRotation, axis_signs));
        return create(std::move(instanceId),
                      _imageId,
                      _definition->normalizedForPositiveDepth(std::move(definitionId)),
                      _imageSize,
                      normalized_pose,
                      _captureTime,
                      _acquisition);
    }

    FramePinholeModel::FramePinholeModel(CameraInstanceId instanceId,
                                         ImageId imageId,
                                         std::shared_ptr<const FramePinholeDefinition> definition,
                                         ImageSize imageSize,
                                         Pose pose,
                                         std::optional<TimeReference> captureTime,
                                         CameraAcquisitionState acquisition)
        : _instanceId(std::move(instanceId)), _imageId(std::move(imageId)), _definition(std::move(definition)),
          _imageSize(imageSize), _pose(std::move(pose)), _captureTime(captureTime), _acquisition(std::move(acquisition))
    {
    }

    Result<CameraModelPtr<CentralCameraModel>> bindCentralCamera(CentralCameraGeometry geometry,
                                                                 CentralCameraBinding binding)
    {
        try
        {
            auto model = FramePinholeModel::create(std::move(binding.instanceId),
                                                   std::move(binding.imageId),
                                                   std::move(geometry.definition),
                                                   binding.imageSize,
                                                   std::move(geometry.pose),
                                                   std::move(binding.captureTime),
                                                   std::move(binding.acquisition));
            return Result<CameraModelPtr<FramePinholeModel>>::success(
                std::make_shared<const FramePinholeModel>(std::move(model)));
        }
        catch (const CameraValidationError& error)
        {
            return Result<CameraModelPtr<FramePinholeModel>>::failure(error.code(), error.what());
        }
        catch (const std::exception& error)
        {
            return Result<CameraModelPtr<FramePinholeModel>>::failure(CameraErrorCode::InvalidModelState, error.what());
        }
    }

    Result<CameraModelPtr<FramePinholeModel>> bindFramePinhole(FramePinholeGeometry geometry,
                                                               FramePinholeBinding binding)
    {
        return bindCentralCamera(std::move(geometry), std::move(binding));
    }

} // namespace placamera
