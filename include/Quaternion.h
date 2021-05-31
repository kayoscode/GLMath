#ifndef INCLUDE_QUATERNION_H
#define INCLUDE_QUATERNION_H

#include "Vector3.h"
#include "Vector4.h"
#include "Matrix44.h"

#include "math.h"

#define MAX(a, b) ((a) > (b))? (a) : (b)

template<typename T>
class Quaternion
{
    public: 
        Quaternion(T x = 0, T y = 0, T z = 0, T w = 1) 
            :x(x), y(y), z(z), w(w)
        {
        }

        Quaternion(const Vector3<T>& axis, T angle) 
            :x(0), y(0), z(0), w(1)
        {
            setToAxisAngle(axis, angle);
        }

        Quaternion(const Quaternion<T>& q) 
            :x(q.x), y(q.y), z(q.z), w(q.w)
        {

        }

        Quaternion<T>& setIdentity() {
            x = 0;
            y = 0;
            z = 0;
            w = 1;

            return *this;
        }

        T length() const {
            return sqrt(x * x + y * y + z * z + w * w);
        }

        Quaternion<T>& normalize() {
            T mag = length();

            if(mag == 0) return *this;

            x /= mag;
            y /= mag;
            z /= mag;
            w /= mag;

            return *this;
        }

        Quaternion<T>& negate() {
            x = -x;
            y = -y;
            z = -z;
            w = -w;

            return *this;
        }

        Quaternion<T> operator*(const Quaternion<T>& right) {
            return Quaternion<T>(
                this->x * right.w + this->w * right.x + this->y * right.z - this->z * right.y, this->y * right.w + this->w * right.y + this->z * right.x - this->x * right.z, this->z * right.w + this->w * right.z + this->x * right.y - this->y * right.x,
                this->w * right.w - this->x * right.x - this->y * right.y - this->z * right.z);
        }

        Quaternion<T>& operator*=(const Quaternion<T>& right) {
            this->x = this->x * right.w + this->w * right.x + this->y * right.z - this->z * right.y;
            this->y = this->y * right.w + this->w * right.y + this->z * right.x - this->x * right.z;
            this->z = this->z * right.w + this->w * right.z + this->x * right.y - this->y * right.x;
            this->w = this->w * right.w - this->x * right.x - this->y * right.y - this->z * right.z;

            return *this;
        }

        Vector3<T> operator*(const Vector3<T>& right) {
            Matrix44 trans = toMatrix();
            Vector4 mod = trans * Vector4(right.x, right.y, right.z, 0);
            return Vector3(mod.x, mod.y, mod.z);
        }

        Vector4<T> operator*(const Vector4<T>& right) {
            Matrix44 trans = toMatrix();
            return trans * right;
        }

        Matrix44<T> toMatrix() const {
            Matrix44<T> matrix;

            T xy = x * y;
            T xz = x * z;
            T xw = x * w;
            T yz = y * z;
            T yw = y * w;
            T zw = z * w;
            T xSquared = x * x;
            T ySquared = y * y;
            T zSquared = z * z;

            matrix.data[0][0] = 1 - 2 * (ySquared + zSquared);
            matrix.data[0][1] = 2 * (xy - zw);
            matrix.data[0][2] = 2 * (xz + yw);
            matrix.data[0][3] = 0;
            matrix.data[1][0] = 2 * (xy + zw);
            matrix.data[1][1] = 1 - 2 * (xSquared + zSquared);
            matrix.data[1][2] = 2 * (yz - xw);
            matrix.data[1][3] = 0;
            matrix.data[2][0] = 2 * (xz - yw);
            matrix.data[2][1] = 2 * (yz + xw);
            matrix.data[2][2] = 1 - 2 * (xSquared + ySquared);
            matrix.data[2][3] = 0;
            matrix.data[3][0] = 0;
            matrix.data[3][1] = 0;
            matrix.data[3][2] = 0;
            matrix.data[3][3] = 1;

            return matrix;
        }

        Quaternion<T>& setToAxisAngle(const Vector3<T>& axis, T angle) {
            Matrix44 rot;
            rot.rotate(axis, angle);

            setMatrix(rot);
            normalize();

            return *this;
        }

        static Quaternion<T> slerp(const Quaternion<T>& a, const Quaternion<T>& b, float blend) {
            Quaternion<T> result;

            T dot = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
            T blendI = 1.0f - blend;

            if (dot < 0)
            {
                result.w = blendI * a.w + blend * -b.w;
                result.x = blendI * a.x + blend * -b.x;
                result.y = blendI * a.y + blend * -b.y;
                result.z = blendI * a.z + blend * -b.z;
            }
            else
            {
                result.w = blendI * a.w + blend * b.w;
                result.x = blendI * a.x + blend * b.x;
                result.y = blendI * a.y + blend * b.y;
                result.z = blendI * a.z + blend * b.z;
            }

            result.normalize();
            return result;
        }

        Quaternion<T>& lookRotation(const Vector3<T>& f, const Vector3<T>& u) {
            Vector3<T> forward(f);
            Vector3<T> up(u);

            forward.normalize();
            up.normalize();
            Vector3<T> right = forward % up;

            Matrix44<T> rot;
            rot.data[0][0] = right.x;
            rot.data[1][0] = right.y;
            rot.data[2][0] = right.z;

            rot.data[0][1] = up.x;
            rot.data[1][1] = up.y;
            rot.data[2][1] = up.z;

            rot.data[0][2] = forward.x;
            rot.data[1][2] = forward.y;
            rot.data[2][2] = forward.z;

            setMatrix(rot);
            normalize();

            return *this;
        }

        Quaternion<T>& setMatrix(const Matrix44<T>& matrix) {
            T m00, m01, m02;
            T m10, m11, m12;
            T m20, m21, m22;

            m00 = matrix.data[0][0];
            m01 = matrix.data[0][1];
            m02 = matrix.data[0][2];

            m10 = matrix.data[1][0];
            m11 = matrix.data[1][1];
            m12 = matrix.data[1][2];

            m20 = matrix.data[2][0];
            m21 = matrix.data[2][1];
            m22 = matrix.data[2][2];

            T s;
            T tr = m00 + m11 + m22;
            if (tr >= 0.0)
            {
                s = (T)sqrt(tr + 1.0);
                w = s * 0.5f;
                s = 0.5f / s;
                x = (m21 - m12) * s;
                y = (m02 - m20) * s;
                z = (m10 - m01) * s;
            }
            else
            {
                T max = MAX(MAX(m00, m11), m22);
                if (max == m00)
                {
                    s = (T)sqrt(m00 - (m11 + m22) + 1.0);
                    x = s * 0.5f;
                    s = 0.5f / s;
                    y = (m01 + m10) * s;
                    z = (m20 + m02) * s;
                    w = (m21 - m12) * s;
                }
                else if (max == m11)
                {
                    s = (T)sqrt(m11 - (m22 + m00) + 1.0);
                    y = s * 0.5f;
                    s = 0.5f / s;
                    z = (m12 + m21) * s;
                    x = (m01 + m10) * s;
                    w = (m02 - m20) * s;
                }
                else
                {
                    s = (T)sqrt(m22 - (m00 + m11) + 1.0);
                    z = s * 0.5f;
                    s = 0.5f / s;
                    x = (m20 + m02) * s;
                    y = (m12 + m21) * s;
                    w = (m10 - m01) * s;
                }
            }

            return *this;
        }

        T x, y, z, w;
};

typedef Quaternion<float> Quaternionf;
typedef Quaternion<double> Quaterniond;
typedef Quaternion<int> Quaternioni;

#endif