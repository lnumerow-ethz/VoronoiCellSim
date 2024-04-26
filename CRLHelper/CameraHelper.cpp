#include "CameraHelper.h"
#include "GeometryHelper.h"

void CameraHelper::rotateCameraFromDrag(CRLCamera &camera, const Vector2F &drag, F sensitivity) {
    Vector3F up_direction = camera.get_up_direction();
    Vector3F right_direction = camera.get_right_direction();

    camera.eye = GeometryHelper::rotateVectorAroundAxisAndCenter(camera.eye, up_direction, -drag.x() * sensitivity,
                                                                 camera.center);

    camera.eye = GeometryHelper::rotateVectorAroundAxisAndCenter(camera.eye, right_direction, -drag.y() * sensitivity,
                                                                 camera.center);

    Vector3F new_up_direction = GeometryHelper::rotateVectorAroundAxisAndCenter(up_direction, right_direction,
                                                                                -drag.y() * sensitivity, camera.center);
    camera.set_up_direction(new_up_direction);
}

void CameraHelper::panCameraFromDrag(CRLCamera &camera, const Vector2F &drag, F sensitivity) {
    Vector3F up_direction = camera.get_up_direction();
    Vector3F right_direction = camera.get_right_direction();

    Vector3F pan_delta = -drag.x() * right_direction * sensitivity + drag.y() * up_direction * sensitivity;

    camera.eye += pan_delta;
    camera.center += pan_delta;
}
