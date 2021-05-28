#ifndef INCLUDE_MATRIX22_H
#define INCLUDE_MATRIX22_H

#include <string>
#include <sstream>

template<class T = float>
struct Matrix22 {
    public:
        /**
         * constructor
         * */
        Matrix22(T m00 = 1, T m01 = 0, T m10 = 0, T m11 = 1)
            :m00(m00), m01(m01),
            m10(m10), m11(m11)
        {
        }

        /**
         * copy constructor
         * */
        Matrix22(const Matrix22<T>& mat)
            :m00(mat.m00), m01(mat.m01),
            m10(mat.m10), m11(mat.m11)
        {
        }

        /**
         * assignment operator
         * */
        Matrix22<T>& operator=(const Matrix22<T>& mat){
            m00 = mat.m00;
            m01 = mat.m01;
            m10 = mat.m10;
            m11 = mat.m11;

            return *this;
        }

        /**
         * addition operator
         * */
        Matrix22<T> operator+(const Matrix22<T>& mat){
            Matrix22<T> ret;

            ret.m00 = m00 + mat.m00;
            ret.m01 = m01 + mat.m01;
            ret.m10 = m10 + mat.m10;
            ret.m11 = m11 + mat.m11;

            return ret;
        }

        /**
         * addition operator
         * */
        Matrix22<T>& operator+=(const Matrix22<T>& mat){
            m00 += mat.m00;
            m01 += mat.m01;
            m01 += mat.m10;
            m11 += mat.m11;

            return *this;
        }

        /**
         * subtraction operator
         * */
        Matrix22<T> operator-(const Matrix22<T>& mat){
            Matrix22<T> ret;

            ret.m00 = m00 - mat.m00;
            ret.m01 = m01 - mat.m01;
            ret.m10 = m10 - mat.m10;
            ret.m11 = m11 - mat.m11;

            return ret;
        }

        /**
         * subtracts a matrix
         * */
        Matrix22<T>& operator-=(const Matrix22<T>& mat){
            m00 -= mat.m00;
            m01 -= mat.m01;
            m01 -= mat.m10;
            m11 -= mat.m11;

            return *this;
        }

        /**
         * multiply by a matrix
         * */
        Matrix22<T> operator*(const Matrix22<T>& mat){
            Matrix22<T> ret;

            T m00 = this->m00 * mat.m00 + this->m10 * mat.m01;
            T m01 = this->m01 * mat.m00 + this->m11 * mat.m01;
            T m10 = this->m00 * mat.m10 + this->m10 * mat.m11;
            T m11 = this->m01 * mat.m10 + this->m11 * mat.m11;

            ret.m00 = m00;
            ret.m01 = m01;
            ret.m10 = m10;
            ret.m11 = m11;
        
            return ret;
        }

        /**
         * multiply by another 2x2 matrix
         * */
        Matrix22<T> operator*=(const Matrix22<T>& mat){
            T m00 = this->m00 * mat.m00 + this->m10 * mat.m01;
            T m01 = this->m01 * mat.m00 + this->m11 * mat.m01;
            T m10 = this->m00 * mat.m10 + this->m10 * mat.m11;
            T m11 = this->m01 * mat.m10 + this->m11 * mat.m11;

            this->m00 = m00;
            this->m01 = m01;
            this->m10 = m10;
            this->m11 = m11;

            return *this;
        }

        /**
         * TODO: transform to vector2
         * */

        /**
         * gets the matrix's transpose
         * */
        Matrix22<T>& transpose(){
            T m01 = this->m10;
            T m10 = this->m01;

            this->m01 = m01;
            this->m10 = m10;

            return *this;
        }

        /**
         * inverts the matrix
         * */
        Matrix22<T>& invert(){
            //set to the determinant
            T determinant = det();

            if(determinant != 0) {
                T determinant_inv = ((T)1)/determinant;

                T t00 =  this->m11*determinant_inv;
                T t01 = -this->m01*determinant_inv;
                T t11 =  this->m00*determinant_inv;
                T t10 = -this->m10*determinant_inv;

                this->m00 = t00;
                this->m01 = t01;
                this->m10 = t10;
                this->m11 = t11;
            }

            return *this;
        }

        /**
         * sets the matrix to identity
         * */
        Matrix22<T>& setIdentity(){
            m00 = 1;
            m01 = 0;
            m10 = 0;
            m11 = 1;

            return *this;
        }

        /**
         * negates the matrix
         * */
        Matrix22<T>& negate(){
            m00 = -m00;
            m01 = -m01;
            m10 = -m10;
            m11 = -m11;

            return *this;
        }

        /**
         * sets the matrix to 0
         * */
        Matrix22<T>& setZero(){
            m00 = 0;
            m01 = 0;
            m10 = 0;
            m11 = 0;

            return *this;
        }

        /**
         * gets the determinant of the 2x2 matrix
         * */
        T det(){
            return m00 * m11 - m01 * m10;
        }

        /**
         * converts the matrix to a printable string
         * */
        operator std::string(){
            char buff[100];
            snprintf(buff, sizeof(buff), "\n[%4.4f] [%4.4f]\n[%4.4f] [%4.4f]", m00, m01, m10, m11);

            return std::string(buff);
        }

    /**
     * TODO: use the vector2 class to scale it
     * */

    T m00, m01;
    T m10, m11;
};

typedef Matrix22<float> Matrix22f;
typedef Matrix22<double> Matrix22d;

#endif