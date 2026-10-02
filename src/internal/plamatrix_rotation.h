#pragma once

#include "placamera/types.h"

#include <plamatrix/geometry/rotation.h>

#include <cmath>

namespace placamera::internal
{

    using Matrix3 = plamatrix::Matrix3d;
    using Vector3d = plamatrix::Vector3d;
    using Quaternion = plamatrix::Quaterniond;

    inline Matrix3 toMatrix(const RotationMatrix& values)
    {
        Matrix3 matrix;
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                matrix(row, column) = values[static_cast<std::size_t>(row * 3 + column)];
            }
        }
        return matrix;
    }

    inline RotationMatrix toRotationMatrix(const Matrix3& matrix)
    {
        RotationMatrix values{};
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                values[static_cast<std::size_t>(row * 3 + column)] = matrix(row, column);
            }
        }
        return values;
    }

    inline Matrix3 transpose(const RotationMatrix& values)
    {
        return toMatrix(values).transpose().eval();
    }

    inline RotationMatrix transposeArray(const RotationMatrix& values)
    {
        return toRotationMatrix(toMatrix(values).transpose().eval());
    }

    inline RotationMatrix multiply(const RotationMatrix& left, const RotationMatrix& right)
    {
        return toRotationMatrix((toMatrix(left) * toMatrix(right)).eval());
    }

    inline Vector3 multiply(const RotationMatrix& rotation, const Vector3& vector)
    {
        const Vector3d value(vector[0], vector[1], vector[2]);
        const Vector3d result = toMatrix(rotation) * value;
        return {result(0), result(1), result(2)};
    }

    inline Vector3 transposeMultiply(const RotationMatrix& rotation, const Vector3& vector)
    {
        const Vector3d value(vector[0], vector[1], vector[2]);
        const Vector3d result = toMatrix(rotation).transpose() * value;
        return {result(0), result(1), result(2)};
    }

    inline bool validRotation(const RotationMatrix& values) noexcept
    {
        const Matrix3 matrix = toMatrix(values);
        if (!matrix.allFinite())
        {
            return false;
        }
        constexpr double tolerance = 1.0e-8;
        Matrix3 inverse{};
        double determinant = 0.0;
        bool invertible = false;
        matrix.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.0);
        return invertible && (matrix.transpose() * matrix).isApprox(Matrix3::Identity(), tolerance) &&
               std::abs(determinant - 1.0) <= tolerance;
    }

    inline RotationMatrix angleAxisRotation(const Vector3& rotation) noexcept
    {
        const Vector3d vector(rotation[0], rotation[1], rotation[2]);
        const double angle = vector.norm();
        if (!(angle > 0.0) || !std::isfinite(angle))
        {
            return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        }
        const Vector3d axis(vector(0) / angle, vector(1) / angle, vector(2) / angle);
        return toRotationMatrix(plamatrix::AngleAxisd(angle, axis).toRotationMatrix());
    }

    inline bool validQuaternion(const std::array<double, 4>& values) noexcept
    {
        const double norm = std::hypot(std::hypot(values[0], values[1]), std::hypot(values[2], values[3]));
        return std::isfinite(norm) && norm > 0.0;
    }

    inline Quaternion quaternionFromRotation(const RotationMatrix& values)
    {
        return Quaternion(toMatrix(values)).normalized();
    }

    inline Quaternion quaternionFromScalarFirst(const std::array<double, 4>& values)
    {
        return Quaternion(values[0], values[1], values[2], values[3]).normalized();
    }

    inline RotationMatrix rotationFromQuaternion(const Quaternion& quaternion)
    {
        return toRotationMatrix(quaternion.normalized().toRotationMatrix());
    }

} // namespace placamera::internal
