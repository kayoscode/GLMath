#ifndef INCLUDE_VECTOR2_H
#define INCLUDE_VECTOR2_H

#include <cmath>
#include <string>

template <typename T = float>
struct Vector2
{
    /**
     * Standard constructor.
     */
    Vector2(T x) : x(x), y(x)
    {
    }

    Vector2(T x = 0, T y = 0) : x(x), y(y)
    {
    }

    // copy constructor
    Vector2(const Vector2<T> &other) : x(other.x), y(other.y)
    {
    }

    /**
     * Implicit cast
     */
    operator std::string() const
    {
        std::stringstream builder;
        builder << "X: " << x << " |Y: " << y;
        return builder.str();
    }

    /**
     * Assignment operator
     */
    Vector2<T> &operator=(const Vector2<T> &other)
    {
        x = other.x;
        y = other.y;
        return *this;
    }

    Vector2<T> operator+(const Vector2<T> &other) const
    {
        return Vector2<T>(x + other.x, y + other.y);
    }

    Vector2<T> &operator+=(const Vector2<T> &other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2<T> operator-(const Vector2<T> &other) const
    {
        return Vector2<T>(x - other.x, y - other.y);
    }

    Vector2<T> &operator-=(Vector2<T> &other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2<T> operator*(const T &scale) const
    {
        return Vector2<T>(x * scale, y * scale);
    }

    Vector2<T> &operator*=(const T &scale)
    {
        x *= scale;
        y *= scale;
        return *this;
    }

    /**
     * Compute the dot product.
     */
    T operator*(const Vector2<T> &other) const
    {
        return x * other.x + y * other.y;
    }

    /**
     * Normalize to get length 1.
     */
    Vector2<T> &normalize()
    {
        T mag = magnitude();
        x /= mag;
        y /= mag;
        return *this;
    }

    Vector2<T> retNormalized() const
    {
        T mag = magnitude();
        return Vector2<T>(x / mag, y / mag);
    }

    T magnitude() const
    {
        return sqrtl(x * x + y * y);
    }

    T magnitudeSquared() const
    {
        return x * x + y * y;
    }

    Vector2<T> negate()
    {
        return Vector2<T>(x * -1, y * -1);
    }

    // Swizzling operators.
    Vector2<T> xy() const { return Vector2<T>(x, y); }
    Vector2<T> yx() const { return Vector2<T>(y, x); }

    // data
    T x, y;
};

typedef Vector2<> Vector2f;
typedef Vector2<double> Vector2d;
typedef Vector2<int> Vector2i;

#endif