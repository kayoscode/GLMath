#ifndef INCLUDE_MATRIX33_H
#define INCLUDE_MATRIX33_H

#include <string>

template<typename T = float>

struct Matrix33 {
	Matrix33(T m00 = 1, T m01 = 0, T m02 = 0, T m10 = 0, T m11 = 1, T m12 = 0,
				T m20 = 0, T m21 = 0, T m22 = 1)
		:m00(m00), m01(m01), m02(m02),
		m10(m10), m11(m11), m12(m12),
		m20(m20), m21(m21), m22(m22)
	{}

	Matrix33(const Matrix33<T>& mat)
		:m00(mat.m00), m01(mat.m01), m02(mat.m02),
		m10(mat.m10), m11(mat.m11), m12(mat.m12),
		m20(mat.m20), m21(mat.m21), m22(mat.m22)
	{}

	Matrix33<T>& operator=(const Matrix33<T>& mat){
		m00 = mat.m00;
		m01 = mat.m01;
		m02 = mat.m02;

		m10 = mat.m10;
		m11 = mat.m11;
		m12 = mat.m12;

		m20 = mat.m20;
		m21 = mat.m21;
		m22 = mat.m22;

		return *this;
	}

	Matrix33<T>& setIdentity(){
		m00 = 1;
		m01 = 0;
		m02 = 0;

		m10 = 0;
		m11 = 1;
		m12 = 0;

		m20 = 0;
		m21 = 0;
		m22 = 0;

		return *this;
	}

	Matrix33<T>& setZero(){
		m00 = 0;
		m01 = 0;
		m02 = 0;

		m10 = 0;
		m11 = 0;
		m12 = 0;

		m20 = 0;
		m21 = 0;
		m22 = 0;

		return *this;
	}

	Matrix33<T> operator+(const Matrix33<T>& right){
		Matrix33<T> ret;

		ret.m00 = this->m00 + right.m00;
		ret.m01 = this->m01 + right.m01;
		ret.m02 = this->m02 + right.m02;
		ret.m10 = this->m10 + right.m10;
		ret.m11 = this->m11 + right.m11;
		ret.m12 = this->m12 + right.m12;
		ret.m20 = this->m20 + right.m20;
		ret.m21 = this->m21 + right.m21;
		ret.m22 = this->m22 + right.m22;

		return ret;
	}

	Matrix33<T>& operator+=(const Matrix33<T>& right){
		this->m00 = this->m00 + right.m00;
		this->m01 = this->m01 + right.m01;
		this->m02 = this->m02 + right.m02;
		this->m10 = this->m10 + right.m10;
		this->m11 = this->m11 + right.m11;
		this->m12 = this->m12 + right.m12;
		this->m20 = this->m20 + right.m20;
		this->m21 = this->m21 + right.m21;
		this->m22 = this->m22 + right.m22;

		return *this;
	}

	Matrix33<T> operator-(const Matrix33<T>& right){
		Matrix33<T> ret;

		ret.m00 = this->m00 - right.m00;
		ret.m01 = this->m01 - right.m01;
		ret.m02 = this->m02 - right.m02;
		ret.m10 = this->m10 - right.m10;
		ret.m11 = this->m11 - right.m11;
		ret.m12 = this->m12 - right.m12;
		ret.m20 = this->m20 - right.m20;
		ret.m21 = this->m21 - right.m21;
		ret.m22 = this->m22 - right.m22;

		return ret;
	}

	Matrix33<T>& operator-=(const Matrix33<T>& right){
		this->m00 = this->m00 - right.m00;
		this->m01 = this->m01 - right.m01;
		this->m02 = this->m02 - right.m02;
		this->m10 = this->m10 - right.m10;
		this->m11 = this->m11 - right.m11;
		this->m12 = this->m12 - right.m12;
		this->m20 = this->m20 - right.m20;
		this->m21 = this->m21 - right.m21;
		this->m22 = this->m22 - right.m22;

		return *this;
	}

	Matrix33<T> operator*(const Matrix33<T>& right){
		Matrix33<T> ret;

		T m00 = this->m00 * right.m00 + this->m10 * right.m01 + this->m20 * right.m02;
		T m01 = this->m01 * right.m00 + this->m11 * right.m01 + this->m21 * right.m02;
		T m02 = this->m02 * right.m00 + this->m12 * right.m01 + this->m22 * right.m02;
		T m10 = this->m00 * right.m10 + this->m10 * right.m11 + this->m20 * right.m12;
		T m11 = this->m01 * right.m10 + this->m11 * right.m11 + this->m21 * right.m12;
		T m12 = this->m02 * right.m10 + this->m12 * right.m11 + this->m22 * right.m12;
		T m20 = this->m00 * right.m20 + this->m10 * right.m21 + this->m20 * right.m22;
		T m21 = this->m01 * right.m20 + this->m11 * right.m21 + this->m21 * right.m22;
		T m22 = this->m02 * right.m20 + this->m12 * right.m21 + this->m22 * right.m22;

		ret.m00 = m00;
		ret.m01 = m01;
		ret.m02 = m02;
		ret.m10 = m10;
		ret.m11 = m11;
		ret.m12 = m12;
		ret.m20 = m20;
		ret.m21 = m21;
		ret.m22 = m22;

		return ret;
	}

	Matrix33<T>& operator*=(const Matrix33<T>& right){
		T m00 = this->m00 * right.m00 + this->m10 * right.m01 + this->m20 * right.m02;
		T m01 = this->m01 * right.m00 + this->m11 * right.m01 + this->m21 * right.m02;
		T m02 = this->m02 * right.m00 + this->m12 * right.m01 + this->m22 * right.m02;
		T m10 = this->m00 * right.m10 + this->m10 * right.m11 + this->m20 * right.m12;
		T m11 = this->m01 * right.m10 + this->m11 * right.m11 + this->m21 * right.m12;
		T m12 = this->m02 * right.m10 + this->m12 * right.m11 + this->m22 * right.m12;
		T m20 = this->m00 * right.m20 + this->m10 * right.m21 + this->m20 * right.m22;
		T m21 = this->m01 * right.m20 + this->m11 * right.m21 + this->m21 * right.m22;
		T m22 = this->m02 * right.m20 + this->m12 * right.m21 + this->m22 * right.m22;

		this->m00 = m00;
		this->m01 = m01;
		this->m02 = m02;
		this->m10 = m10;
		this->m11 = m11;
		this->m12 = m12;
		this->m20 = m20;
		this->m21 = m21;
		this->m22 = m22;

		return *this;
	}

	Matrix33<T>& transpose(){
		T m00 = this->m00;
		T m01 = this->m10;
		T m02 = this->m20;
		T m10 = this->m01;
		T m11 = this->m11;
		T m12 = this->m21;
		T m20 = this->m02;
		T m21 = this->m12;
		T m22 = this->m22;

		this->m00 = m00;
		this->m01 = m01;
		this->m02 = m02;
		this->m10 = m10;
		this->m11 = m11;
		this->m12 = m12;
		this->m20 = m20;
		this->m21 = m21;
		this->m22 = m22;

		return *this;
	}

	T det(){
		return m00 * (m11 * m22 - m12 * m21)
		+ m01 * (m12 * m20 - m10 * m22)
		+ m02 * (m10 * m21 - m11 * m20);
	}

	operator std::string(){
		char buff[100];
		snprintf(buff, sizeof(buff), "\n"
									"[%4.4f] [%4.4f] [%4.4f]\n"
									"[%4.4f] [%4.4f] [%4.4f]\n"
									"[%4.4f] [%4.4f] [%4.4f]",
									m00, m01, m02, m10, m11, m12, m20, m21, m22);

		return std::string(buff);
	}

	Matrix33<T>& invert(){
		T determinant = this->det();

		if(determinant != 0){
			T determinant_inv = (T)1/determinant;

			T t00 = this->m11 * this->m22 - this->m12* this->m21;
			T t01 = - this->m10 * this->m22 + this->m12 * this->m20;
			T t02 = this->m10 * this->m21 - this->m11 * this->m20;
			T t10 = - this->m01 * this->m22 + this->m02 * this->m21;
			T t11 = this->m00 * this->m22 - this->m02 * this->m20;
			T t12 = - this->m00 * this->m21 + this->m01 * this->m20;
			T t20 = this->m01 * this->m12 - this->m02 * this->m11;
			T t21 = -this->m00 * this->m12 + this->m02 * this->m10;
			T t22 = this->m00 * this->m11 - this->m01 * this->m10;

			this->m00 = t00*determinant_inv;
			this->m11 = t11*determinant_inv;
			this->m22 = t22*determinant_inv;
			this->m01 = t10*determinant_inv;
			this->m10 = t01*determinant_inv;
			this->m20 = t02*determinant_inv;
			this->m02 = t20*determinant_inv;
			this->m12 = t21*determinant_inv;
			this->m21 = t12*determinant_inv;
		}

		return *this;
	}

	Matrix33<T>& negate(){
		m00 = -m00;
		m01 = -m01;
		m02 = -m02;

		m10 = -m10;
		m11 = -m11;
		m12 = -m12;

		m20 = -m20;
		m21 = -m21;
		m22 = -m22;

		return *this;
	}

	/**
	 * TODO: create vector3 transform
	 * */
	
	
	T m00, m01, m02;
	T m10, m11, m12;
	T m20, m21, m22;
};

typedef Matrix33<float> Matrix33f;
typedef Matrix33<double> Matrix33d;

#endif