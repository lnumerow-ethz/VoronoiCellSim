#include <Eigen/Geometry>
#include "CRLHelper/GeometryHelper.h"

Vector3F GeometryHelper::rotateVectorAroundAxisAndCenter(const Vector3F &vec, const Vector3F &axis, F angle_rad,
                                                         const Vector3F &center) {
    // Normalize the rotation axis
    Vector3F unit_axis = axis.normalized();

    // Translate the vector such that the center of rotation becomes the origin
    Vector3F translated_vec = vec - center;

    // Create a quaternion representing the rotation
    Eigen::AngleAxisd rotation(angle_rad, unit_axis);

    // Rotate the vector
    Vector3F rotatedVec = rotation * translated_vec;

    // Translate back
    Vector3F result = rotatedVec + center;

    return result;
}

Vector2F GeometryHelper::rotateVector2D(const Vector2F &vec, F angle_rad) {
    F x = vec.x();
    F y = vec.y();
    F cost = cos(angle_rad);
    F sint = sin(angle_rad);
    return {x * cost - y * sint, x * sint + y * cost};
}

bool GeometryHelper::intersectLinesPointDirection2D(const Vector2F &p0, const Vector2F &d0, const Vector2F &p1,
                                                    const Vector2F &d1, Vector2F &t) {
    // Line 1: P = p0 + t0 * d0
    // Line 2: Q = p1 + t1 * d1
    // For intersection: p0 + t0 * d0 = p1 + t1 * d1
    // => t0 * d0 - t1 * d1 = p1 - p0

    Matrix<F, 2, 2> A;
    A << d0, -d1;
    Vector2F b = p1 - p0;

    // Check if lines are parallel (determinant is zero)
    F det = A.determinant();
    if (std::abs(det) < 1e-8) {
        // Lines are parallel or coincident, no unique intersection
        return false;
    }

    // Solve for t0 and t1
    t = A.inverse() * b;
    return true;
}

F GeometryHelper::windingNumber2D(const Vector2F &point, const Vector2F &edge0, const Vector2F &edge1) {
    F x0 = edge0(0);
    F y0 = edge0(1);
    F x1 = edge1(0);
    F y1 = edge1(1);

    F winding_number = atan2(y1 - point.y(), x1 - point.x()) - atan2(y0 - point.y(), x0 - point.x());
    if (winding_number > M_PI) winding_number -= 2 * M_PI;
    if (winding_number < -M_PI) winding_number += 2 * M_PI;
    return winding_number;
}

F GeometryHelper::windingNumber3D(const Vector3F &point, const Vector3F &face0, const Vector3F &face1,
                                  const Vector3F &face2) {
    Vector3F a = face0 - point;
    Vector3F b = face1 - point;
    Vector3F c = face2 - point;

    F al = a.norm();
    F bl = b.norm();
    F cl = c.norm();

    Vector3F t = b.cross(c);
    F numerator = a.dot(t);
    F denominator = al * bl * cl + a.dot(b) * cl + b.dot(c) * al + c.dot(a) * bl;

    F omega = 2.0 * atan2(numerator, denominator);
    if (omega > 2 * M_PI) omega -= 4 * M_PI;
    if (omega < -2 * M_PI) omega += 4 * M_PI;
    return omega;
}

F GeometryHelper::triangleAreaWithOrigin2D(const Vector2F &p0, const Vector2F &p1) {
    F x0 = p0.x();
    F y0 = p0.y();
    F x1 = p1.x();
    F y1 = p1.y();

    return 0.5e0 * x0 * y1 - 0.5e0 * x1 * y0;
}

F GeometryHelper::tetVolumeWithOrigin3D(const Vector3F &p0, const Vector3F &p1, const Vector3F &p2) {
    F x0 = p0.x();
    F y0 = p0.y();
    F z0 = p0.z();
    F x1 = p1.x();
    F y1 = p1.y();
    F z1 = p1.z();
    F x2 = p2.x();
    F y2 = p2.y();
    F z2 = p2.z();

    return (y1 * z2 - y2 * z1) * x0 / 0.6e1 + (-y0 * z2 + y2 * z0) * x1 / 0.6e1 + x2 * (y0 * z1 - y1 * z0) / 0.6e1;
}
