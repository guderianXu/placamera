#include "placamera/reference/CameraReferenceComparator.h"

#include <placamera/frame_camera.h>

#include "../internal/plamatrix_rotation.h"

#include <algorithm>
#include <cmath>

namespace placamera::reference
{
    namespace
    {

        double norm(const std::array<double, 3>& vector)
        {
            return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
        }

    } // namespace

    CameraReferenceComparison CameraReferenceComparator::compare(const placamera::RasterModel& instance,
                                                                 const ResolvedCameraReference& reference)
    {
        if (reference.status != ReferenceResolutionStatus::Resolved || !reference.pose)
        {
            return {CameraReferenceComparisonStatus::UnresolvedReference,
                    std::nullopt,
                    std::nullopt,
                    reference.reason.empty() ? "reference is unresolved" : reference.reason};
        }

        const auto* frameCamera = dynamic_cast<const placamera::FramePinholeModel*>(&instance);
        if (!frameCamera || !instance.capabilities().contains(placamera::CapabilityKind::StaticPose))
        {
            return {CameraReferenceComparisonStatus::RegionModelWithoutStaticPose,
                    std::nullopt,
                    std::nullopt,
                    "camera model does not provide a static pose"};
        }

        const placamera::Pose& estimated = frameCamera->pose();
        if (estimated.frame != reference.pose->frame)
        {
            return {CameraReferenceComparisonStatus::FrameMismatch,
                    std::nullopt,
                    std::nullopt,
                    "estimated and reference poses use different coordinate frames"};
        }

        const std::array<double, 3> difference{estimated.center[0] - reference.pose->center[0],
                                               estimated.center[1] - reference.pose->center[1],
                                               estimated.center[2] - reference.pose->center[2]};
        const auto relative = placamera::internal::toMatrix(reference.pose->cameraToWorldRotation).transpose() *
                              placamera::internal::toMatrix(estimated.cameraToWorldRotation);
        const double trace = relative(0, 0) + relative(1, 1) + relative(2, 2);
        const double cosine = std::clamp((trace - 1.0) * 0.5, -1.0, 1.0);
        return {CameraReferenceComparisonStatus::Compared, norm(difference), std::acos(cosine), {}};
    }

} // namespace placamera::reference
