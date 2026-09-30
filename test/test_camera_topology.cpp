#include <placamera/camera_topology.h>
#include <placamera/rig_topology.h>

#include <gtest/gtest.h>

#include <vector>

namespace
{

    using namespace placamera;

    TEST(CameraTopologyTest, ValidatesIndependentSensorAndAcquisitionIdentity)
    {
        SensorMountState mount;
        mount.masterSensorId = CameraDefinitionId("master-sensor");
        EXPECT_TRUE(validateSensorMount(mount, nullptr));
        const CameraDefinitionId own_sensor("master-sensor");
        EXPECT_FALSE(validateSensorMount(mount, &own_sensor));

        CameraAcquisitionState acquisition;
        acquisition.captureGroupId = CaptureGroupId("capture-1");
        acquisition.masterCameraId = CameraInstanceId("master-camera");
        acquisition.layerIndex = 2;
        acquisition.role = CameraRole::Keyframe;
        EXPECT_TRUE(validateCameraAcquisition(acquisition, nullptr));
        const CameraInstanceId own_camera("master-camera");
        EXPECT_FALSE(validateCameraAcquisition(acquisition, &own_camera));
    }

    TEST(CameraTopologyTest, EnforcesRollingShutterModeDomain)
    {
        CameraAcquisitionState regularized;
        regularized.rollingShutterMode = RollingShutterMode::Regularized;
        regularized.rollingShutter.translation = {0.1, -0.2, 0.0};
        EXPECT_TRUE(validateCameraAcquisition(regularized));
        regularized.rollingShutter.rotationVector[0] = 0.01;
        EXPECT_FALSE(validateCameraAcquisition(regularized));

        CameraAcquisitionState full = regularized;
        full.rollingShutterMode = RollingShutterMode::Full;
        full.rollingShutter.translation[2] = 0.3;
        EXPECT_TRUE(validateCameraAcquisition(full));
    }

    TEST(CameraTopologyTest, ComposesInstalledNumericFrameStates)
    {
        const auto definition = FramePinholeDefinition::create(
            CameraDefinitionId("rig-definition"),
            FrameIntrinsics{100.0, 100.0, 50.0, 50.0},
            {},
            PixelConvention::PixelCenter,
            FrameId("world"));
        const auto model = FramePinholeModel::create(
            CameraInstanceId("rig-camera"),
            ImageId("rig-image"),
            definition,
            ImageSize{100, 100},
            Pose::create(FrameId("world"), {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}));
        const std::vector<FramePinholeNumericState> cameras{FramePinholeNumericState::fromModel(model)};

        RigTopology rig;
        rig.captures.push_back(RigCapture{});
        rig.captures.front().rigId = 0;
        rig.captures.front().captureId = 0;
        rig.captures.front().rigCenterInWorld = {2.0, 0.0, 0.0};
        rig.sensors.push_back(RigSensor{});
        rig.sensors.front().rigId = 0;
        rig.sensors.front().sensorId = 0;
        rig.cameraBindings.push_back({0, 0, 0, 0});

        const auto composed = composeRigCameras(cameras, rig);
        ASSERT_TRUE(composed) << composed.message();
        ASSERT_EQ(composed.value().size(), 1U);
        EXPECT_EQ(composed.value().front().pose().center, (Vector3{2.0, 0.0, 0.0}));
    }

} // namespace
