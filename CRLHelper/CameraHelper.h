#pragma once

#include "CRLSubApp.h"

namespace CameraHelper {
void rotateCameraFromDrag(CRLCamera &camera, const Vector2F &drag, F sensitivity);

void panCameraFromDrag(CRLCamera &camera, const Vector2F &drag, F sensitivity);
}  // namespace CameraHelper
