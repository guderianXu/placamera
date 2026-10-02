#include "placamera/dataset_formats.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>

namespace placamera
{
    namespace
    {
        bool matchesMetashapeInspectionView(const ImportedPixelCalibration& imported, double tolerance) noexcept
        {
            if (!imported.metashapeCalibration)
            {
                return true;
            }
            const MetashapeCalibration& exact = *imported.metashapeCalibration;
            const auto& matrix = imported.intrinsicMatrix;
            const BrownConradyDistortion& distortion = imported.distortion;
            const auto matches = [tolerance](double left, double right)
            { return std::isfinite(left) && std::isfinite(right) && std::abs(left - right) <= tolerance; };
            return matches(matrix[0], exact.f + exact.b1) && matches(matrix[1], exact.b2) &&
                   matches(matrix[2], exact.cx) && matches(matrix[4], exact.f) && matches(matrix[5], exact.cy) &&
                   matches(distortion.radialK1, exact.k1) && matches(distortion.radialK2, exact.k2) &&
                   matches(distortion.radialK3, exact.k3) && matches(distortion.radialK4, exact.k4) &&
                   matches(distortion.tangentialP1, exact.p1) && matches(distortion.tangentialP2, exact.p2) &&
                   matches(distortion.tangentialP3, exact.p3) && matches(distortion.tangentialP4, exact.p4) &&
                   distortion.tangentialConvention == BrownTangentialConvention::Metashape;
        }

    } // namespace

    bool CameraImportCompatibility::isExactlyRepresentable(double tolerance) const noexcept
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0 || unsupportedReason.has_value())
        {
            return false;
        }
        return std::all_of(sourceOnlyTerms.begin(),
                           sourceOnlyTerms.end(),
                           [tolerance](const CameraImportCalibrationTerm& term)
                           { return std::isfinite(term.value) && std::abs(term.value) <= tolerance; });
    }

    Result<CentralCameraGeometry> makeCentralCameraGeometry(const ImportedCamera& camera,
                                                            CameraDefinitionId definitionId,
                                                            FrameId worldFrame,
                                                            double tolerance)
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0)
        {
            return Result<FramePinholeGeometry>::failure(CameraErrorCode::InvalidArgument,
                                                         "dataset camera tolerance must be finite and non-negative");
        }
        if (!camera.compatibility.isExactlyRepresentable(tolerance))
        {
            const std::string reason = camera.compatibility.unsupportedReason.value_or(
                "dataset camera contains source calibration terms that the central-camera model cannot express");
            return Result<FramePinholeGeometry>::failure(CameraErrorCode::UnsupportedModel, reason);
        }
        if (!std::all_of(camera.calibration.intrinsicMatrix.begin(),
                         camera.calibration.intrinsicMatrix.end(),
                         [](double value) { return std::isfinite(value); }))
        {
            return Result<FramePinholeGeometry>::failure(CameraErrorCode::InvalidIntrinsics,
                                                         "dataset camera intrinsics must be finite");
        }
        if (!matchesMetashapeInspectionView(camera.calibration, tolerance))
        {
            return Result<FramePinholeGeometry>::failure(
                CameraErrorCode::InvalidModelState,
                "dataset camera's lossless Metashape calibration disagrees with its normalized inspection view");
        }

        const auto& intrinsics = camera.calibration.intrinsicMatrix;
        if (std::abs(intrinsics[3]) > tolerance || std::abs(intrinsics[6]) > tolerance ||
            std::abs(intrinsics[7]) > tolerance || std::abs(intrinsics[8] - 1.0) > tolerance)
        {
            return Result<FramePinholeGeometry>::failure(
                CameraErrorCode::UnsupportedModel, "dataset camera has a nonstandard homogeneous calibration row");
        }

        try
        {
            std::shared_ptr<const FramePinholeDefinition> definition;
            if (camera.calibration.metashapeCalibration)
            {
                definition = FramePinholeDefinition::create(std::move(definitionId),
                                                            *camera.calibration.metashapeCalibration,
                                                            PixelConvention::PixelCenter,
                                                            worldFrame,
                                                            false,
                                                            1.0,
                                                            1,
                                                            1,
                                                            camera.calibration.projectionModel);
            }
            else
            {
                definition = FramePinholeDefinition::create(
                    std::move(definitionId),
                    {intrinsics[0], intrinsics[4], intrinsics[2], intrinsics[5], 1.0, 1, 1, intrinsics[1]},
                    camera.calibration.distortion,
                    PixelConvention::PixelCenter,
                    worldFrame,
                    false,
                    camera.calibration.projectionModel);
            }
            auto pose = Pose::create(std::move(worldFrame), camera.center, camera.cameraToWorldRotation);
            return Result<FramePinholeGeometry>::success({std::move(definition), std::move(pose)});
        }
        catch (const CameraValidationError& error)
        {
            return Result<FramePinholeGeometry>::failure(error.code(), error.what());
        }
        catch (const std::exception& error)
        {
            return Result<FramePinholeGeometry>::failure(CameraErrorCode::InvalidModelState, error.what());
        }
    }

    Result<CentralCameraGeometry> makeDatasetCentralCamera(const DatasetFrameCamera& camera,
                                                           CameraDefinitionId definitionId,
                                                           FrameId worldFrame,
                                                           double tolerance)
    {
        return makeCentralCameraGeometry(camera, std::move(definitionId), std::move(worldFrame), tolerance);
    }

    Result<FramePinholeGeometry> makeDatasetFramePinhole(const ImportedCamera& camera,
                                                         CameraDefinitionId definitionId,
                                                         FrameId worldFrame,
                                                         double tolerance)
    {
        return makeCentralCameraGeometry(camera, std::move(definitionId), std::move(worldFrame), tolerance);
    }

} // namespace placamera
