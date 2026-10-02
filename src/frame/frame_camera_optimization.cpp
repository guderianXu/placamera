#include "placamera/frame_camera.h"

#include "../internal/plamatrix_rotation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace placamera
{

    namespace
    {

        void appendSelectedBlock(OptimizationLayout* layout,
                                 const FrameOptimizationSelection& selection,
                                 FrameOptimizationParameter parameter,
                                 std::string name,
                                 OptimizationParameterKind kind,
                                 std::string unit,
                                 std::size_t size,
                                 bool affectsDefinition)
        {
            if (!selection.enabled(parameter))
            {
                return;
            }
            layout->blocks.push_back(
                {std::move(name), kind, std::move(unit), layout->parameterCount(), size, affectsDefinition});
        }

        FrameOptimizationSelection effectiveSelection(const FramePinholeModel& model,
                                                      FrameOptimizationSelection selection) noexcept
        {
            const FrameProjectionModel projection = model.pinholeDefinition().projectionModel();
            if (projection == FrameProjectionModel::Spherical || projection == FrameProjectionModel::Cylindrical)
            {
                for (const FrameOptimizationParameter parameter : {FrameOptimizationParameter::F,
                                                                   FrameOptimizationParameter::Cx,
                                                                   FrameOptimizationParameter::Cy,
                                                                   FrameOptimizationParameter::B1,
                                                                   FrameOptimizationParameter::B2,
                                                                   FrameOptimizationParameter::K1,
                                                                   FrameOptimizationParameter::K2,
                                                                   FrameOptimizationParameter::K3,
                                                                   FrameOptimizationParameter::K4,
                                                                   FrameOptimizationParameter::P1,
                                                                   FrameOptimizationParameter::P2,
                                                                   FrameOptimizationParameter::P3,
                                                                   FrameOptimizationParameter::P4})
                {
                    selection.set(parameter, false);
                }
            }
            const CameraAcquisitionState& acquisition = model.acquisition();
            if (acquisition.rollingShutterMode == RollingShutterMode::Disabled)
            {
                selection.set(FrameOptimizationParameter::RollingRotation, false)
                    .set(FrameOptimizationParameter::RollingTranslation, false);
            }
            else if (acquisition.rollingShutterMode == RollingShutterMode::Regularized)
            {
                selection.set(FrameOptimizationParameter::RollingRotation, false);
            }
            return selection;
        }

        void applyPrincipalPointDelta(FrameCalibration* calibration, double cxDelta, double cyDelta) noexcept
        {
            if (calibration->principalPointDecomposition)
            {
                PrincipalPointDecomposition& principal = *calibration->principalPointDecomposition;
                if (cxDelta != 0.0)
                {
                    principal.cxOffset += cxDelta;
                    calibration->cx = principal.imageCenterX + principal.cxOffset;
                }
                if (cyDelta != 0.0)
                {
                    principal.cyOffset += cyDelta;
                    calibration->cy = principal.imageCenterY + principal.cyOffset;
                }
                return;
            }
            calibration->cx += cxDelta;
            calibration->cy += cyDelta;
        }

    } // namespace

    FrameOptimizationSelection FrameOptimizationSelection::none() noexcept
    {
        return {};
    }

    FrameOptimizationSelection FrameOptimizationSelection::all() noexcept
    {
        FrameOptimizationSelection selection;
        selection._enabled.fill(true);
        return selection;
    }

    bool FrameOptimizationSelection::enabled(FrameOptimizationParameter parameter) const noexcept
    {
        const std::size_t index = static_cast<std::size_t>(parameter);
        return index < _enabled.size() && _enabled[index];
    }

    FrameOptimizationSelection& FrameOptimizationSelection::set(FrameOptimizationParameter parameter,
                                                                bool enabledValue) noexcept
    {
        const std::size_t index = static_cast<std::size_t>(parameter);
        if (index < _enabled.size())
        {
            _enabled[index] = enabledValue;
        }
        return *this;
    }

    OptimizationLayout FramePinholeModel::optimizationLayout() const
    {
        if (_definition->projectionModel() != FrameProjectionModel::Perspective ||
            _acquisition.rollingShutterMode != RollingShutterMode::Disabled)
        {
            return optimizationLayout(FrameOptimizationSelection::all());
        }
        return {{{"pose.rotation", OptimizationParameterKind::RotationVector, "rad", 0, 3, false},
                 {"pose.translation", OptimizationParameterKind::Translation, "m", 3, 3, false},
                 {"intrinsics", OptimizationParameterKind::Intrinsics, "px", 6, 4, true},
                 {"distortion", OptimizationParameterKind::Distortion, "normalized", 10, 5, true}}};
    }

    OptimizationLayout FramePinholeModel::optimizationLayout(const FrameOptimizationSelection& selection) const
    {
        const FrameOptimizationSelection effective = effectiveSelection(*this, selection);
        OptimizationLayout layout;
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::PoseRotation,
                            "pose.rotation",
                            OptimizationParameterKind::RotationVector,
                            "rad",
                            3,
                            false);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::PoseTranslation,
                            "pose.translation",
                            OptimizationParameterKind::Translation,
                            "m",
                            3,
                            false);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::F,
                            "calibration.f",
                            OptimizationParameterKind::Intrinsics,
                            "px",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::Cx,
                            "calibration.cx",
                            OptimizationParameterKind::Intrinsics,
                            "px",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::Cy,
                            "calibration.cy",
                            OptimizationParameterKind::Intrinsics,
                            "px",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::B1,
                            "calibration.b1",
                            OptimizationParameterKind::Intrinsics,
                            "px",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::B2,
                            "calibration.b2",
                            OptimizationParameterKind::Intrinsics,
                            "px",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::K1,
                            "calibration.k1",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::K2,
                            "calibration.k2",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::K3,
                            "calibration.k3",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::K4,
                            "calibration.k4",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::P1,
                            "calibration.p1",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::P2,
                            "calibration.p2",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::P3,
                            "calibration.p3",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::P4,
                            "calibration.p4",
                            OptimizationParameterKind::Distortion,
                            "normalized",
                            1,
                            true);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::RollingRotation,
                            "rolling_shutter.rotation",
                            OptimizationParameterKind::RotationVector,
                            "rad",
                            3,
                            false);
        appendSelectedBlock(&layout,
                            effective,
                            FrameOptimizationParameter::RollingTranslation,
                            "rolling_shutter.translation",
                            OptimizationParameterKind::Translation,
                            "camera",
                            _acquisition.rollingShutterMode == RollingShutterMode::Regularized ? 2 : 3,
                            false);
        return layout;
    }

    Result<RasterModelPtr> FramePinholeModel::withOptimizationUpdate(const OptimizationUpdate& update) const
    {
        if (_definition->projectionModel() != FrameProjectionModel::Perspective ||
            _acquisition.rollingShutterMode != RollingShutterMode::Disabled)
        {
            return withOptimizationUpdate(update, FrameOptimizationSelection::all());
        }
        const OptimizationLayout layout = optimizationLayout();
        if (update.delta.size() != layout.parameterCount())
        {
            return Result<RasterModelPtr>::failure(CameraErrorCode::InvalidArgument,
                                                   "frame-pinhole optimization update must contain 15 values");
        }
        for (const double value : update.delta)
        {
            if (!std::isfinite(value))
            {
                return Result<RasterModelPtr>::failure(CameraErrorCode::InvalidArgument,
                                                       "frame-pinhole optimization update must be finite");
            }
        }

        bool changes_definition = false;
        for (std::size_t index = 6; index < update.delta.size(); ++index)
        {
            changes_definition = changes_definition || update.delta[index] != 0.0;
        }
        if (changes_definition && !update.definitionId)
        {
            return Result<RasterModelPtr>::failure(
                CameraErrorCode::InvalidArgument,
                "frame-pinhole calibration updates require an explicit result definition identifier");
        }

        try
        {
            const RotationMatrix rotation =
                internal::multiply(internal::angleAxisRotation(
                                       {update.delta[0], update.delta[1], update.delta[2]}),
                                   _pose.cameraToWorldRotation);
            const Vector3 center{{_pose.center[0] + update.delta[3],
                                  _pose.center[1] + update.delta[4],
                                  _pose.center[2] + update.delta[5]}};

            std::shared_ptr<const FramePinholeDefinition> definition = _definition;
            if (update.definitionId)
            {
                FrameCalibration calibration = _definition->calibration();
                calibration.f += update.delta[7];
                calibration.b1 += update.delta[6] - update.delta[7];
                applyPrincipalPointDelta(&calibration, update.delta[8], update.delta[9]);
                calibration.k1 += update.delta[10];
                calibration.k2 += update.delta[11];
                calibration.k3 += update.delta[12];
                if (_definition->distortion().tangentialConvention == BrownTangentialConvention::Metashape)
                {
                    calibration.p1 += update.delta[13];
                    calibration.p2 += update.delta[14];
                }
                else
                {
                    calibration.p2 += update.delta[13];
                    calibration.p1 += update.delta[14];
                }
                definition = FramePinholeDefinition::create(*update.definitionId,
                                                            calibration,
                                                            _definition->pixelConvention(),
                                                            _definition->groundFrame(),
                                                            _definition->depthAxisFlipped(),
                                                            _definition->intrinsics().pixelPitch,
                                                            _definition->intrinsics().uAxisSign,
                                                            _definition->intrinsics().vAxisSign,
                                                            _definition->projectionModel(),
                                                            _definition->sensorMount());
            }
            auto model = FramePinholeModel::create(update.instanceId,
                                                   _imageId,
                                                   std::move(definition),
                                                   _imageSize,
                                                   Pose::create(groundFrame(), center, rotation),
                                                   _captureTime,
                                                   _acquisition);
            return Result<RasterModelPtr>::success(std::make_shared<const FramePinholeModel>(std::move(model)));
        }
        catch (const CameraValidationError& error)
        {
            return Result<RasterModelPtr>::failure(error.code(), error.what());
        }
        catch (const std::exception& error)
        {
            return Result<RasterModelPtr>::failure(CameraErrorCode::InvalidModelState, error.what());
        }
    }

    Result<RasterModelPtr> FramePinholeModel::withOptimizationUpdate(const OptimizationUpdate& update,
                                                                     const FrameOptimizationSelection& selection) const
    {
        const FrameOptimizationSelection effective = effectiveSelection(*this, selection);
        const OptimizationLayout layout = optimizationLayout(effective);
        if (update.delta.size() != layout.parameterCount())
        {
            return Result<RasterModelPtr>::failure(
                CameraErrorCode::InvalidArgument,
                "selected frame optimization update size does not match its parameter layout");
        }
        for (const double value : update.delta)
        {
            if (!std::isfinite(value))
            {
                return Result<RasterModelPtr>::failure(CameraErrorCode::InvalidArgument,
                                                       "selected frame optimization update must be finite");
            }
        }

        std::size_t offset = 0;
        auto scalar = [&](FrameOptimizationParameter parameter)
        {
            if (!effective.enabled(parameter))
            {
                return 0.0;
            }
            return update.delta[offset++];
        };

        double rotation_x = 0.0;
        double rotation_y = 0.0;
        double rotation_z = 0.0;
        if (effective.enabled(FrameOptimizationParameter::PoseRotation))
        {
            rotation_x = update.delta[offset++];
            rotation_y = update.delta[offset++];
            rotation_z = update.delta[offset++];
        }
        Vector3 center = _pose.center;
        if (effective.enabled(FrameOptimizationParameter::PoseTranslation))
        {
            center[0] += update.delta[offset++];
            center[1] += update.delta[offset++];
            center[2] += update.delta[offset++];
        }

        FrameCalibration calibration = _definition->calibration();
        calibration.f += scalar(FrameOptimizationParameter::F);
        const double cx_delta = scalar(FrameOptimizationParameter::Cx);
        const double cy_delta = scalar(FrameOptimizationParameter::Cy);
        applyPrincipalPointDelta(&calibration, cx_delta, cy_delta);
        calibration.b1 += scalar(FrameOptimizationParameter::B1);
        calibration.b2 += scalar(FrameOptimizationParameter::B2);
        calibration.k1 += scalar(FrameOptimizationParameter::K1);
        calibration.k2 += scalar(FrameOptimizationParameter::K2);
        calibration.k3 += scalar(FrameOptimizationParameter::K3);
        calibration.k4 += scalar(FrameOptimizationParameter::K4);
        calibration.p1 += scalar(FrameOptimizationParameter::P1);
        calibration.p2 += scalar(FrameOptimizationParameter::P2);
        calibration.p3 += scalar(FrameOptimizationParameter::P3);
        calibration.p4 += scalar(FrameOptimizationParameter::P4);

        CameraAcquisitionState acquisition = _acquisition;
        if (effective.enabled(FrameOptimizationParameter::RollingRotation))
        {
            acquisition.rollingShutter.rotationVector[0] += update.delta[offset++];
            acquisition.rollingShutter.rotationVector[1] += update.delta[offset++];
            acquisition.rollingShutter.rotationVector[2] += update.delta[offset++];
        }
        if (effective.enabled(FrameOptimizationParameter::RollingTranslation))
        {
            acquisition.rollingShutter.translation[0] += update.delta[offset++];
            acquisition.rollingShutter.translation[1] += update.delta[offset++];
            if (acquisition.rollingShutterMode == RollingShutterMode::Full)
            {
                acquisition.rollingShutter.translation[2] += update.delta[offset++];
            }
        }

        const bool changes_definition =
            std::any_of(layout.blocks.begin(),
                        layout.blocks.end(),
                        [&](const OptimizationParameterBlock& block)
                        {
                            if (!block.affectsDefinition)
                            {
                                return false;
                            }
                            return std::any_of(update.delta.begin() + block.offset,
                                               update.delta.begin() + block.offset + block.size,
                                               [](double value) { return value != 0.0; });
                        });
        if (changes_definition && !update.definitionId)
        {
            return Result<RasterModelPtr>::failure(
                CameraErrorCode::InvalidArgument,
                "selected frame calibration updates require an explicit result definition identifier");
        }

        try
        {
            std::shared_ptr<const FramePinholeDefinition> definition = _definition;
            if (update.definitionId)
            {
                definition = FramePinholeDefinition::create(*update.definitionId,
                                                            calibration,
                                                            _definition->pixelConvention(),
                                                            _definition->groundFrame(),
                                                            _definition->depthAxisFlipped(),
                                                            _definition->intrinsics().pixelPitch,
                                                            _definition->intrinsics().uAxisSign,
                                                            _definition->intrinsics().vAxisSign,
                                                            _definition->projectionModel(),
                                                            _definition->sensorMount());
            }
            const RotationMatrix rotation =
                internal::multiply(internal::angleAxisRotation({rotation_x, rotation_y, rotation_z}),
                                   _pose.cameraToWorldRotation);
            auto model = FramePinholeModel::create(update.instanceId,
                                                   _imageId,
                                                   std::move(definition),
                                                   _imageSize,
                                                   Pose::create(groundFrame(), center, rotation),
                                                   _captureTime,
                                                   std::move(acquisition));
            return Result<RasterModelPtr>::success(std::make_shared<const FramePinholeModel>(std::move(model)));
        }
        catch (const CameraValidationError& error)
        {
            return Result<RasterModelPtr>::failure(error.code(), error.what());
        }
        catch (const std::exception& error)
        {
            return Result<RasterModelPtr>::failure(CameraErrorCode::InvalidModelState, error.what());
        }
    }

} // namespace placamera
