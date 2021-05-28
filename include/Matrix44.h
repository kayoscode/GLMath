#ifndef INCLUDE_MATRIX44_H
#define INCLUDE_MATRIX44_H

#include "Vector2.h"
#include "Vector3.h"
#include "Vector4.h"
#include <math.h>

template <typename T = float>
struct Matrix44{
    Matrix44(T m00 = 1, T m01 = 0, T m02 = 0, T m03 = 0,
                T m10 = 0, T m11 = 1, T m12 = 0, T m13 = 0,
                T m20 = 0, T m21 = 0, T m22 = 1, T m23 = 0,
                T m30 = 0, T m31 = 0, T m32 = 0, T m33 = 1)
        :m00(m00), m01(m01), m02(m02), m03(m03),
        m10(m10), m11(m11), m12(m12), m13(m12),
        m20(m20), m21(m21), m22(m22), m23(m23),
        m30(m30), m31(m31), m32(m32), m33(m33)
    {}

    Matrix44(const Matrix44<T>& mat)
        :m00(mat.m00), m01(mat.m01), m02(mat.m02), m03(mat.m03),
        m10(mat.m10), m11(mat.m11), m12(mat.m12), m13(mat.m12),
        m20(mat.m20), m21(mat.m21), m22(mat.m22), m23(mat.m23),
        m30(mat.m30), m31(mat.m31), m32(mat.m32), m33(mat.m33)
    {}

    Matrix44<T>& operator=(const Matrix44<T>& mat){
        m00 = mat.m00; 
        m01 = mat.m01; 
        m02 = mat.m02; 
        m03 = mat.m03;

        m10 = mat.m10; 
        m11 = mat.m11; 
        m12 = mat.m12; 
        m13 = mat.m13;

        m20 = mat.m20; 
        m21 = mat.m21; 
        m22 = mat.m22; 
        m23 = mat.m23;

        m30 = mat.m30; 
        m31 = mat.m31; 
        m32 = mat.m32; 
        m33 = mat.m33;

        return *this;
    }

    Matrix44<T>& negate(){
        m00 = -m00; 
        m01 = -m01; 
        m02 = -m02; 
        m03 = -m03;

        m10 = -m10; 
        m11 = -m11; 
        m12 = -m12; 
        m13 = -m13;

        m20 = -m20; 
        m21 = -m21; 
        m22 = -m22; 
        m23 = -m23;

        m30 = -m30; 
        m31 = -m31; 
        m32 = -m32; 
        m33 = -m33;

        return *this;
    }


    Matrix44<T>& setIdentity(){
        m00 = 1; 
        m01 = 0; 
        m02 = 0; 
        m03 = 0;

        m10 = 0; 
        m11 = 1; 
        m12 = 0; 
        m13 = 0;

        m20 = 0; 
        m21 = 0; 
        m22 = 1; 
        m23 = 0;

        m30 = 0; 
        m31 = 0; 
        m32 = 0; 
        m33 = 1;

        return *this;
    }

    Matrix44<T>& setZero(){
        m00 = 0; 
        m01 = 0; 
        m02 = 0; 
        m03 = 0;

        m10 = 0; 
        m11 = 0; 
        m12 = 0; 
        m13 = 0;

        m20 = 0; 
        m21 = 0; 
        m22 = 0; 
        m23 = 0;

        m30 = 0; 
        m31 = 0; 
        m32 = 0; 
        m33 = 0;

        return *this;
    }

    Matrix44<T> operator+(const Matrix44<T>& mat) const {
        Matrix44<T> ret;

        ret.m00 = this->m00 + mat.m00;
        ret.m01 = this->m01 + mat.m01;
        ret.m02 = this->m02 + mat.m02;
        ret.m03 = this->m03 + mat.m03;
        ret.m10 = this->m10 + mat.m10;
        ret.m11 = this->m11 + mat.m11;
        ret.m12 = this->m12 + mat.m12;
        ret.m13 = this->m13 + mat.m13;
        ret.m20 = this->m20 + mat.m20;
        ret.m21 = this->m21 + mat.m21;
        ret.m22 = this->m22 + mat.m22;
        ret.m23 = this->m23 + mat.m23;
        ret.m30 = this->m30 + mat.m30;
        ret.m31 = this->m31 + mat.m31;
        ret.m32 = this->m32 + mat.m32;
        ret.m33 = this->m33 + mat.m33;

        return ret;
    }

    Matrix44<T>& operator+=(const Matrix44<T>& mat){
        this->m00 = this->m00 + mat.m00;
        this->m01 = this->m01 + mat.m01;
        this->m02 = this->m02 + mat.m02;
        this->m03 = this->m03 + mat.m03;
        this->m10 = this->m10 + mat.m10;
        this->m11 = this->m11 + mat.m11;
        this->m12 = this->m12 + mat.m12;
        this->m13 = this->m13 + mat.m13;
        this->m20 = this->m20 + mat.m20;
        this->m21 = this->m21 + mat.m21;
        this->m22 = this->m22 + mat.m22;
        this->m23 = this->m23 + mat.m23;
        this->m30 = this->m30 + mat.m30;
        this->m31 = this->m31 + mat.m31;
        this->m32 = this->m32 + mat.m32;
        this->m33 = this->m33 + mat.m33;
    }

    Matrix44<T> operator-(const Matrix44<T>& mat) const {
        Matrix44<T> ret;

        ret.m00 = this->m00 - mat.m00;
        ret.m01 = this->m01 - mat.m01;
        ret.m02 = this->m02 - mat.m02;
        ret.m03 = this->m03 - mat.m03;
        ret.m10 = this->m10 - mat.m10;
        ret.m11 = this->m11 - mat.m11;
        ret.m12 = this->m12 - mat.m12;
        ret.m13 = this->m13 - mat.m13;
        ret.m20 = this->m20 - mat.m20;
        ret.m21 = this->m21 - mat.m21;
        ret.m22 = this->m22 - mat.m22;
        ret.m23 = this->m23 - mat.m23;
        ret.m30 = this->m30 - mat.m30;
        ret.m31 = this->m31 - mat.m31;
        ret.m32 = this->m32 - mat.m32;
        ret.m33 = this->m33 - mat.m33;

        return ret;
    }

    Matrix44<T>& operator-=(const Matrix44<T>& mat){
        this->m00 = this->m00 - mat.m00;
        this->m01 = this->m01 - mat.m01;
        this->m02 = this->m02 - mat.m02;
        this->m03 = this->m03 - mat.m03;
        this->m10 = this->m10 - mat.m10;
        this->m11 = this->m11 - mat.m11;
        this->m12 = this->m12 - mat.m12;
        this->m13 = this->m13 - mat.m13;
        this->m20 = this->m20 - mat.m20;
        this->m21 = this->m21 - mat.m21;
        this->m22 = this->m22 - mat.m22;
        this->m23 = this->m23 - mat.m23;
        this->m30 = this->m30 - mat.m30;
        this->m31 = this->m31 - mat.m31;
        this->m32 = this->m32 - mat.m32;
        this->m33 = this->m33 - mat.m33;
    }

    Matrix44<T> operator*(const Matrix44<T>& mat) const {
        Matrix44<T> ret;

        T m00 = this->m00 * mat.m00 + this->m10 * mat.m01 + this->m20 * mat.m02 + this->m30 * mat.m03;
        T m01 = this->m01 * mat.m00 + this->m11 * mat.m01 + this->m21 * mat.m02 + this->m31 * mat.m03;
        T m02 = this->m02 * mat.m00 + this->m12 * mat.m01 + this->m22 * mat.m02 + this->m32 * mat.m03;
        T m03 = this->m03 * mat.m00 + this->m13 * mat.m01 + this->m23 * mat.m02 + this->m33 * mat.m03;
        T m10 = this->m00 * mat.m10 + this->m10 * mat.m11 + this->m20 * mat.m12 + this->m30 * mat.m13;
        T m11 = this->m01 * mat.m10 + this->m11 * mat.m11 + this->m21 * mat.m12 + this->m31 * mat.m13;
        T m12 = this->m02 * mat.m10 + this->m12 * mat.m11 + this->m22 * mat.m12 + this->m32 * mat.m13;
        T m13 = this->m03 * mat.m10 + this->m13 * mat.m11 + this->m23 * mat.m12 + this->m33 * mat.m13;
        T m20 = this->m00 * mat.m20 + this->m10 * mat.m21 + this->m20 * mat.m22 + this->m30 * mat.m23;
        T m21 = this->m01 * mat.m20 + this->m11 * mat.m21 + this->m21 * mat.m22 + this->m31 * mat.m23;
        T m22 = this->m02 * mat.m20 + this->m12 * mat.m21 + this->m22 * mat.m22 + this->m32 * mat.m23;
        T m23 = this->m03 * mat.m20 + this->m13 * mat.m21 + this->m23 * mat.m22 + this->m33 * mat.m23;
        T m30 = this->m00 * mat.m30 + this->m10 * mat.m31 + this->m20 * mat.m32 + this->m30 * mat.m33;
        T m31 = this->m01 * mat.m30 + this->m11 * mat.m31 + this->m21 * mat.m32 + this->m31 * mat.m33;
        T m32 = this->m02 * mat.m30 + this->m12 * mat.m31 + this->m22 * mat.m32 + this->m32 * mat.m33;
        T m33 = this->m03 * mat.m30 + this->m13 * mat.m31 + this->m23 * mat.m32 + this->m33 * mat.m33;

        ret.m00 = m00;
        ret.m01 = m01;
        ret.m02 = m02;
        ret.m03 = m03;
        ret.m10 = m10;
        ret.m11 = m11;
        ret.m12 = m12;
        ret.m13 = m13;
        ret.m20 = m20;
        ret.m21 = m21;
        ret.m22 = m22;
        ret.m23 = m23;
        ret.m30 = m30;
        ret.m31 = m31;
        ret.m32 = m32;
        ret.m33 = m33;

        return ret;
    }

    Matrix44<T>& operator*=(const Matrix44<T>& mat){
        T m00 = this->m00 * mat.m00 + this->m10 * mat.m01 + this->m20 * mat.m02 + this->m30 * mat.m03;
        T m01 = this->m01 * mat.m00 + this->m11 * mat.m01 + this->m21 * mat.m02 + this->m31 * mat.m03;
        T m02 = this->m02 * mat.m00 + this->m12 * mat.m01 + this->m22 * mat.m02 + this->m32 * mat.m03;
        T m03 = this->m03 * mat.m00 + this->m13 * mat.m01 + this->m23 * mat.m02 + this->m33 * mat.m03;
        T m10 = this->m00 * mat.m10 + this->m10 * mat.m11 + this->m20 * mat.m12 + this->m30 * mat.m13;
        T m11 = this->m01 * mat.m10 + this->m11 * mat.m11 + this->m21 * mat.m12 + this->m31 * mat.m13;
        T m12 = this->m02 * mat.m10 + this->m12 * mat.m11 + this->m22 * mat.m12 + this->m32 * mat.m13;
        T m13 = this->m03 * mat.m10 + this->m13 * mat.m11 + this->m23 * mat.m12 + this->m33 * mat.m13;
        T m20 = this->m00 * mat.m20 + this->m10 * mat.m21 + this->m20 * mat.m22 + this->m30 * mat.m23;
        T m21 = this->m01 * mat.m20 + this->m11 * mat.m21 + this->m21 * mat.m22 + this->m31 * mat.m23;
        T m22 = this->m02 * mat.m20 + this->m12 * mat.m21 + this->m22 * mat.m22 + this->m32 * mat.m23;
        T m23 = this->m03 * mat.m20 + this->m13 * mat.m21 + this->m23 * mat.m22 + this->m33 * mat.m23;
        T m30 = this->m00 * mat.m30 + this->m10 * mat.m31 + this->m20 * mat.m32 + this->m30 * mat.m33;
        T m31 = this->m01 * mat.m30 + this->m11 * mat.m31 + this->m21 * mat.m32 + this->m31 * mat.m33;
        T m32 = this->m02 * mat.m30 + this->m12 * mat.m31 + this->m22 * mat.m32 + this->m32 * mat.m33;
        T m33 = this->m03 * mat.m30 + this->m13 * mat.m31 + this->m23 * mat.m32 + this->m33 * mat.m33;

        this->m00 = m00;
        this->m01 = m01;
        this->m02 = m02;
        this->m03 = m03;
        this->m10 = m10;
        this->m11 = m11;
        this->m12 = m12;
        this->m13 = m13;
        this->m20 = m20;
        this->m21 = m21;
        this->m22 = m22;
        this->m23 = m23;
        this->m30 = m30;
        this->m31 = m31;
        this->m32 = m32;
        this->m33 = m33;

        return *this;
    }

    void transform(const Vector4<T>& right, Vector4<T>& dest) {
        T x = m00 * right.x + m10 * right.y + m20 * right.z + m30 * right.w;
        T y = m01 * right.x + m11 * right.y + m21 * right.z + m31 * right.w;
        T z = m02 * right.x + m12 * right.y + m22 * right.z + m32 * right.w;
        T w = m03 * right.x + m13 * right.y + m23 * right.z + m33 * right.w;

        dest.x = x;
        dest.y = y;
        dest.z = z;
        dest.w = w;
    }

    Vector4<T> operator*(const Vector4<T>& right) {
        Vector4<T> ret;
        
        T x = m00 * right.x + m10 * right.y + m20 * right.z + m30 * right.w;
        T y = m01 * right.x + m11 * right.y + m21 * right.z + m31 * right.w;
        T z = m02 * right.x + m12 * right.y + m22 * right.z + m32 * right.w;
        T w = m03 * right.x + m13 * right.y + m23 * right.z + m33 * right.w;

        ret.x = x;
        ret.y = y;
        ret.z = z;
        ret.w = w;

        return ret;
    }

    void scale(const Vector3<T>& scale) {
        this->m00 = this->m00 * scale.x;
        this->m01 = this->m01 * scale.x;
        this->m02 = this->m02 * scale.x;
        this->m03 = this->m03 * scale.x;
        this->m10 = this->m10 * scale.y;
        this->m11 = this->m11 * scale.y;
        this->m12 = this->m12 * scale.y;
        this->m13 = this->m13 * scale.y;
        this->m20 = this->m20 * scale.z;
        this->m21 = this->m21 * scale.z;
        this->m22 = this->m22 * scale.z;
        this->m23 = this->m23 * scale.z;
    }

    void translate(const Vector3<T>& vec) {
        this->m30 += this->m00 * vec.x + this->m10 * vec.y + this->m20 * vec.z;
        this->m31 += this->m01 * vec.x + this->m11 * vec.y + this->m21 * vec.z;
        this->m32 += this->m02 * vec.x + this->m12 * vec.y + this->m22 * vec.z;
        this->m33 += this->m03 * vec.x + this->m13 * vec.y + this->m23 * vec.z;
    }

    void rotate(const Vector3<T>& axis, T angle) {
        T c = (T)cos(angle);
        T s = (T)sin(angle);

        T oneminusc = (T)1 - c;
        T xy = axis.x * axis.y;
        T yz = axis.y * axis.z;
        T xz = axis.x * axis.z;
        T xs = axis.x * s;
        T ys = axis.y * s;
        T zs = axis.z * s;

        T f00 = axis.x * axis.x * oneminusc + c;
        T f01 = xy * oneminusc + zs;
        T f02 = xz * oneminusc - ys;

        T f10 = xy * oneminusc - zs;
        T f11 = axis.y * axis.y * oneminusc + c;
        T f12 = yz * oneminusc + xs;

        T f20 = xz * oneminusc + ys;
        T f21 = yz * oneminusc - xs;
        T f22 = axis.z * axis.z * oneminusc + c;

        T t00 = this->m00 * f00 + this->m10 * f01 + this->m20 * f02;
        T t01 = this->m01 * f00 + this->m11 * f01 + this->m21 * f02;
        T t02 = this->m02 * f00 + this->m12 * f01 + this->m22 * f02;
        T t03 = this->m03 * f00 + this->m13 * f01 + this->m23 * f02;
        T t10 = this->m00 * f10 + this->m10 * f11 + this->m20 * f12;
        T t11 = this->m01 * f10 + this->m11 * f11 + this->m21 * f12;
        T t12 = this->m02 * f10 + this->m12 * f11 + this->m22 * f12;
        T t13 = this->m03 * f10 + this->m13 * f11 + this->m23 * f12;

        this->m20 = this->m00 * f20 + this->m10 * f21 + this->m20 * f22;
        this->m21 = this->m01 * f20 + this->m11 * f21 + this->m21 * f22;
        this->m22 = this->m02 * f20 + this->m12 * f21 + this->m22 * f22;
        this->m23 = this->m03 * f20 + this->m13 * f21 + this->m23 * f22;
        this->m00 = t00;
        this->m01 = t01;
        this->m02 = t02;
        this->m03 = t03;
        this->m10 = t10;
        this->m11 = t11;
        this->m12 = t12;
        this->m13 = t13;
    }

    Matrix44<T>& transpose(){
        T m00 = this->m00;
        T m01 = this->m10;
        T m02 = this->m20;
        T m03 = this->m30;
        T m10 = this->m01;
        T m11 = this->m11;
        T m12 = this->m21;
        T m13 = this->m31;
        T m20 = this->m02;
        T m21 = this->m12;
        T m22 = this->m22;
        T m23 = this->m32;
        T m30 = this->m03;
        T m31 = this->m13;
        T m32 = this->m23;
        T m33 = this->m33;

        this->m00 = m00;
        this->m01 = m01;
        this->m02 = m02;
        this->m03 = m03;
        this->m10 = m10;
        this->m11 = m11;
        this->m12 = m12;
        this->m13 = m13;
        this->m20 = m20;
        this->m21 = m21;
        this->m22 = m22;
        this->m23 = m23;
        this->m30 = m30;
        this->m31 = m31;
        this->m32 = m32;
        this->m33 = m33;

        return *this;
    }

    T det() {
        T f = m00 * ((m11 * m22 * m33 + m12 * m23 * m31 + m13 * m21 * m32)
                - m13 * m22 * m31
                - m11 * m23 * m32
                - m12 * m21 * m33);

        f -= m01 * ((m10 * m22 * m33 + m12 * m23 * m30 + m13 * m20 * m32)
            - m13 * m22 * m30
            - m10 * m23 * m32
            - m12 * m20 * m33);

        f += m02 * ((m10 * m21 * m33 + m11 * m23 * m30 + m13 * m20 * m31)
            - m13 * m21 * m30
            - m10 * m23 * m31
            - m11 * m20 * m33);

        f -= m03 * ((m10 * m21 * m32 + m11 * m22 * m30 + m12 * m20 * m31)
            - m12 * m21 * m30
            - m10 * m22 * m31
            - m11 * m20 * m32);

        return f;
    }

    static T det33(T t00, T t01, T t02,
                            T t10, T t11, T t12,
                            T t20, T t21, T t22)
    {
        return   t00 * (t11 * t22 - t12 * t21)
            + t01 * (t12 * t20 - t10 * t22)
            + t02 * (t10 * t21 - t11 * t20);
    }

    Matrix44<T> inverse() {
        T determinant = this->det();
        Matrix44<T> ret;

        if (determinant != 0) {
            T determinant_inv = (T)1/determinant;

            //first row
            T t00 =  det33(this->m11, this->m12, this->m13, this->m21, this->m22, this->m23, this->m31, this->m32, this->m33);
            T t01 = -det33(this->m10, this->m12, this->m13, this->m20, this->m22, this->m23, this->m30, this->m32, this->m33);
            T t02 =  det33(this->m10, this->m11, this->m13, this->m20, this->m21, this->m23, this->m30, this->m31, this->m33);
            T t03 = -det33(this->m10, this->m11, this->m12, this->m20, this->m21, this->m22, this->m30, this->m31, this->m32);

            // second row
            T t10 = -det33(this->m01, this->m02, this->m03, this->m21, this->m22, this->m23, this->m31, this->m32, this->m33);
            T t11 =  det33(this->m00, this->m02, this->m03, this->m20, this->m22, this->m23, this->m30, this->m32, this->m33);
            T t12 = -det33(this->m00, this->m01, this->m03, this->m20, this->m21, this->m23, this->m30, this->m31, this->m33);
            T t13 =  det33(this->m00, this->m01, this->m02, this->m20, this->m21, this->m22, this->m30, this->m31, this->m32);

            // third row
            T t20 =  det33(this->m01, this->m02, this->m03, this->m11, this->m12, this->m13, this->m31, this->m32, this->m33);
            T t21 = -det33(this->m00, this->m02, this->m03, this->m10, this->m12, this->m13, this->m30, this->m32, this->m33);
            T t22 =  det33(this->m00, this->m01, this->m03, this->m10, this->m11, this->m13, this->m30, this->m31, this->m33);
            T t23 = -det33(this->m00, this->m01, this->m02, this->m10, this->m11, this->m12, this->m30, this->m31, this->m32);

            // fourth row
            T t30 = -det33(this->m01, this->m02, this->m03, this->m11, this->m12, this->m13, this->m21, this->m22, this->m23);
            T t31 =  det33(this->m00, this->m02, this->m03, this->m10, this->m12, this->m13, this->m20, this->m22, this->m23);
            T t32 = -det33(this->m00, this->m01, this->m03, this->m10, this->m11, this->m13, this->m20, this->m21, this->m23);
            T t33 =  det33(this->m00, this->m01, this->m02, this->m10, this->m11, this->m12, this->m20, this->m21, this->m22);

            //transpose and divide by the determinant
            ret.m00 = t00*determinant_inv;
            ret.m11 = t11*determinant_inv;
            ret.m22 = t22*determinant_inv;
            ret.m33 = t33*determinant_inv;
            ret.m01 = t10*determinant_inv;
            ret.m10 = t01*determinant_inv;
            ret.m20 = t02*determinant_inv;
            ret.m02 = t20*determinant_inv;
            ret.m12 = t21*determinant_inv;
            ret.m21 = t12*determinant_inv;
            ret.m03 = t30*determinant_inv;
            ret.m30 = t03*determinant_inv;
            ret.m13 = t31*determinant_inv;
            ret.m31 = t13*determinant_inv;
            ret.m32 = t23*determinant_inv;
            ret.m23 = t32*determinant_inv;
        }

        return ret;
    }

    Matrix44<T>& invert(){
        T determinant = this->det();

        if (determinant != 0) {
            T determinant_inv = (T)1/determinant;

            //first row
            T t00 =  det33(this->m11, this->m12, this->m13, this->m21, this->m22, this->m23, this->m31, this->m32, this->m33);
            T t01 = -det33(this->m10, this->m12, this->m13, this->m20, this->m22, this->m23, this->m30, this->m32, this->m33);
            T t02 =  det33(this->m10, this->m11, this->m13, this->m20, this->m21, this->m23, this->m30, this->m31, this->m33);
            T t03 = -det33(this->m10, this->m11, this->m12, this->m20, this->m21, this->m22, this->m30, this->m31, this->m32);

            // second row
            T t10 = -det33(this->m01, this->m02, this->m03, this->m21, this->m22, this->m23, this->m31, this->m32, this->m33);
            T t11 =  det33(this->m00, this->m02, this->m03, this->m20, this->m22, this->m23, this->m30, this->m32, this->m33);
            T t12 = -det33(this->m00, this->m01, this->m03, this->m20, this->m21, this->m23, this->m30, this->m31, this->m33);
            T t13 =  det33(this->m00, this->m01, this->m02, this->m20, this->m21, this->m22, this->m30, this->m31, this->m32);

            // third row
            T t20 =  det33(this->m01, this->m02, this->m03, this->m11, this->m12, this->m13, this->m31, this->m32, this->m33);
            T t21 = -det33(this->m00, this->m02, this->m03, this->m10, this->m12, this->m13, this->m30, this->m32, this->m33);
            T t22 =  det33(this->m00, this->m01, this->m03, this->m10, this->m11, this->m13, this->m30, this->m31, this->m33);
            T t23 = -det33(this->m00, this->m01, this->m02, this->m10, this->m11, this->m12, this->m30, this->m31, this->m32);

            // fourth row
            T t30 = -det33(this->m01, this->m02, this->m03, this->m11, this->m12, this->m13, this->m21, this->m22, this->m23);
            T t31 =  det33(this->m00, this->m02, this->m03, this->m10, this->m12, this->m13, this->m20, this->m22, this->m23);
            T t32 = -det33(this->m00, this->m01, this->m03, this->m10, this->m11, this->m13, this->m20, this->m21, this->m23);
            T t33 =  det33(this->m00, this->m01, this->m02, this->m10, this->m11, this->m12, this->m20, this->m21, this->m22);

            //transpose and divide by the determinant
            this->m00 = t00*determinant_inv;
            this->m11 = t11*determinant_inv;
            this->m22 = t22*determinant_inv;
            this->m33 = t33*determinant_inv;
            this->m01 = t10*determinant_inv;
            this->m10 = t01*determinant_inv;
            this->m20 = t02*determinant_inv;
            this->m02 = t20*determinant_inv;
            this->m12 = t21*determinant_inv;
            this->m21 = t12*determinant_inv;
            this->m03 = t30*determinant_inv;
            this->m30 = t03*determinant_inv;
            this->m13 = t31*determinant_inv;
            this->m31 = t13*determinant_inv;
            this->m32 = t23*determinant_inv;
            this->m23 = t32*determinant_inv;
        }

        return *this;
    }

    operator std::string() const {
        char buff[200];
        snprintf(buff, sizeof(buff), "\n"
                                    "[%4.4f] [%4.4f] [%4.4f] [%4.4f]\n"
                                    "[%4.4f] [%4.4f] [%4.4f] [%4.4f]\n"
                                    "[%4.4f] [%4.4f] [%4.4f] [%4.4f]\n"
                                    "[%4.4f] [%4.4f] [%4.4f] [%4.4f]",
                                    m00, m01, m02, m03, 
                                    m10, m11, m12, m13,
                                    m20, m21, m22, m23,
                                    m30, m31, m32, m33);

        return std::string(buff);
    }

    T m00, m01, m02, m03;
    T m10, m11, m12, m13;
    T m20, m21, m22, m23;
    T m30, m31, m32, m33;
};

typedef Matrix44<float> Matrix44f;
typedef Matrix44<double> Matrix44d;

#endif