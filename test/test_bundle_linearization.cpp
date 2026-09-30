#include <placamera/bundle_linearization.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace
{

    using namespace placamera;

    FramePinholeModel makeFrameModel()
    {
        FrameIntrinsics intrinsics;
        intrinsics.focalX = 120.0;
        intrinsics.focalY = 118.0;
        intrinsics.principalX = 50.0;
        intrinsics.principalY = 45.0;
        intrinsics.skew = 0.3;
        BrownConradyDistortion distortion;
        distortion.radialK1 = 0.001;
        distortion.radialK2 = -0.0001;
        const FrameId frame("world");
        const auto definition = FramePinholeDefinition::create(
            CameraDefinitionId("definition"), intrinsics, distortion, PixelConvention::PixelCenter, frame);
        return FramePinholeModel::create(CameraInstanceId("instance"),
                                         ImageId("image"),
                                         definition,
                                         ImageSize{100, 100},
                                         Pose::create(frame,
                                                      Vector3{{0.2, -0.1, 0.3}},
                                                      RotationMatrix{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}));
    }

    TEST(BundleLinearizationTest, FrameReturnsPointPoseAndSelectedParameterJacobians)
    {
        const FramePinholeModel model = makeFrameModel();
        const GroundCoordinate point{FrameId("world"), Vector3{{1.0, 0.5, 8.0}}};
        FrameOptimizationSelection selection = FrameOptimizationSelection::none();
        selection.set(FrameOptimizationParameter::PoseRotation)
            .set(FrameOptimizationParameter::PoseTranslation)
            .set(FrameOptimizationParameter::F);

        const auto result = linearizeForBundle(model, point, selection);
        ASSERT_TRUE(result) << result.message();
        EXPECT_EQ(result.value().layout.parameterCount(), 7U);
        EXPECT_EQ(result.value().parameterJacobian.size(), 14U);
        GroundCoordinate plus_point = point;
        GroundCoordinate minus_point = point;
        plus_point.position[0] += 1.0e-6;
        minus_point.position[0] -= 1.0e-6;
        const auto plus = model.groundToImage(plus_point);
        const auto minus = model.groundToImage(minus_point);
        ASSERT_TRUE(plus);
        ASSERT_TRUE(minus);
        EXPECT_NEAR(
            result.value().pointJacobian[0], (plus.value().image.sample - minus.value().image.sample) / 2.0e-6, 1.0e-5);
        EXPECT_NEAR(
            result.value().pointJacobian[3], (plus.value().image.line - minus.value().image.line) / 2.0e-6, 1.0e-5);
        EXPECT_TRUE(std::isfinite(result.value().poseJacobian[0]));
        EXPECT_TRUE(std::isfinite(result.value().parameterJacobian[0]));
        EXPECT_TRUE(result.value().projection.positiveDepth.has_value());

        const auto full_result = linearizeForBundle(model, point, FrameOptimizationSelection::all());
        ASSERT_TRUE(full_result) << full_result.message();
        EXPECT_EQ(full_result.value().layout.parameterCount(), 19U);
        EXPECT_TRUE(std::all_of(full_result.value().parameterJacobian.begin(),
                                full_result.value().parameterJacobian.end(),
                                [](double value) { return std::isfinite(value); }));
    }

    TEST(BundleLinearizationTest, FrameRejectsPointBehindCamera)
    {
        const FramePinholeModel model = makeFrameModel();
        const GroundCoordinate point{FrameId("world"), Vector3{{0.0, 0.0, -2.0}}};
        const auto result = linearizeForBundle(model, point);
        EXPECT_FALSE(result);
        EXPECT_EQ(result.errorCode(), CameraErrorCode::OutsideModelDomain);
    }

    TEST(BundleLinearizationTest, LineScanReturnsObservationLineJacobian)
    {
        const FrameId frame("world");
        LineScanOptics optics;
        optics.focalLengthMillimeters = 10.0;
        optics.samplePitchMillimeters = 0.01;
        const auto definition = LineScanDefinition::create(
            CameraDefinitionId("line-definition"), frame, optics, LineScanPixelConvention::PixelCenter);
        const LineScanTrajectory trajectory = LineScanTrajectory::create(
            {TrajectorySample(TimeReference{TimeScale::Tdb, 0.0}, Vector3{{0.0, 0.0, 0.0}}),
             TrajectorySample(TimeReference{TimeScale::Tdb, 1.0}, Vector3{{0.0, 0.0, 0.0}})});
        LineTiming timing;
        timing.lineZero = 0.5;
        timing.startTimeSeconds = 0.0;
        timing.secondsPerLine = 0.01;
        const LineScanModel model = LineScanModel::create(CameraInstanceId("line-instance"),
                                                          ImageId("line-image"),
                                                          definition,
                                                          ImageSize{100, 100},
                                                          trajectory,
                                                          timing);
        const GroundCoordinate point{frame, Vector3{{0.5, 0.2, 10.0}}};
        const auto result = linearizeForBundle(model, point, 50.5);
        ASSERT_TRUE(result) << result.message();
        EXPECT_EQ(result.value().layout.parameterCount(), 1U);
        EXPECT_EQ(result.value().parameterJacobian.size(), 2U);
        EXPECT_TRUE(std::isfinite(result.value().pointJacobian[0]));
        EXPECT_TRUE(std::isfinite(result.value().poseJacobian[0]));
    }

} // namespace
