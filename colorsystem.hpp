
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
//
// policy: currently single-headered, simple, low overhead, no ICC support.
//

#include <stdio.h>
#include <tuple>
#include <array>

#pragma once

namespace ColorSystem
{

class Vector3
{
  public:
    typedef std::array<float, 3> vec3;
    vec3 v_;
    constexpr Vector3(const float a, const float b, const float c) : v_({a, b, c}) { ; }
    constexpr float x() const { return v_[0]; }
    constexpr float y() const { return v_[1]; }
    constexpr float z() const { return v_[2]; }
    constexpr float operator[](const int &i) const { return v_[i]; }
};

class Matrix3
{
  public:
    typedef std::array<float, 9> matrix;
    matrix m_;

    constexpr int M(const int x, const int y) const { return x + y * 3; }
    constexpr int I(const int y, const int x) const { return ((x - 1) + (y - 1) * 3); }
    constexpr Matrix3(
        const float &a00, const float &a01, const float &a02,
        const float &a10, const float &a11, const float &a12,
        const float &a20, const float &a21, const float &a22) : m_({a00, a01, a02, a10, a11, a12, a20, a21, a22}) { ; }
    constexpr Matrix3(void) : m_({1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f}) { ; }

    constexpr Matrix3 diag(const Vector3 &v) const
    {
        return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]);
    }

    constexpr float operator[](const int &i) const
    {
        return m_[i];
    }
    constexpr Matrix3 mul(const Matrix3 &a, const Matrix3 &b) const
    {
        return Matrix3(
            a[M(0, 0)] * b[M(0, 0)] + a[M(1, 0)] * b[M(0, 1)] + a[M(2, 0)] * b[M(0, 2)],
            a[M(0, 0)] * b[M(1, 0)] + a[M(1, 0)] * b[M(1, 1)] + a[M(2, 0)] * b[M(1, 2)],
            a[M(0, 0)] * b[M(2, 0)] + a[M(1, 0)] * b[M(2, 1)] + a[M(2, 0)] * b[M(2, 2)],

            a[M(0, 1)] * b[M(0, 0)] + a[M(1, 1)] * b[M(0, 1)] + a[M(2, 1)] * b[M(0, 2)],
            a[M(0, 1)] * b[M(1, 0)] + a[M(1, 1)] * b[M(1, 1)] + a[M(2, 1)] * b[M(1, 2)],
            a[M(0, 1)] * b[M(2, 0)] + a[M(1, 1)] * b[M(2, 1)] + a[M(2, 1)] * b[M(2, 2)],

            a[M(0, 2)] * b[M(0, 0)] + a[M(1, 2)] * b[M(0, 1)] + a[M(2, 2)] * b[M(0, 2)],
            a[M(0, 2)] * b[M(1, 0)] + a[M(1, 2)] * b[M(1, 1)] + a[M(2, 2)] * b[M(1, 2)],
            a[M(0, 2)] * b[M(2, 0)] + a[M(1, 2)] * b[M(2, 1)] + a[M(2, 2)] * b[M(2, 2)]);
    }
    constexpr Matrix3 mul(const Matrix3 &b) const
    {
        return mul(*this, b);
    }

    constexpr Vector3 apply(const Matrix3 &m, const Vector3 &v) const
    {
        return Vector3(
            m[M(0, 0)] * v[0] + m[M(1, 0)] * v[1] + m[M(2, 0)] * v[2],
            m[M(0, 1)] * v[0] + m[M(1, 1)] * v[1] + m[M(2, 1)] * v[2],
            m[M(0, 2)] * v[0] + m[M(1, 2)] * v[1] + m[M(2, 2)] * v[2]);
    }
    constexpr Vector3 apply(const Vector3 &v) const
    {
        return apply(*this, v);
    }

    constexpr float det(const Matrix3 &m) const
    {
        return m[I(1, 1)] * m[I(2, 2)] * m[I(3, 3)] +
               m[I(2, 1)] * m[I(3, 2)] * m[I(1, 3)] +
               m[I(3, 1)] * m[I(1, 2)] * m[I(2, 3)] -
               m[I(1, 1)] * m[I(3, 2)] * m[I(2, 3)] -
               m[I(3, 1)] * m[I(2, 2)] * m[I(1, 3)] -
               m[I(2, 1)] * m[I(1, 2)] * m[I(3, 3)];
    }
    constexpr float det(void) const
    {
        return det(*this);
    }

    constexpr Matrix3 mul(const Matrix3 &m, const float &a) const
    {
        return Matrix3(
            m[0] * a, m[1] * a, m[2] * a,
            m[3] * a, m[4] * a, m[5] * a,
            m[6] * a, m[7] * a, m[8] * a);
    }
    constexpr Matrix3 mul(const float a) const
    {
        return mul(*this, a);
    }
    constexpr Matrix3 div(const Matrix3 &m, const float &a) const
    {
        return mul(m, 1.f / a);
    }
    constexpr Matrix3 div(const float a) const
    {
        return div(*this, a);
    }

    constexpr Matrix3 add(const Matrix3 &a, const Matrix3 &b) const
    {
        return Matrix3(
            a[0] + b[0], a[1] + b[1], a[2] + b[2],
            a[3] + b[3], a[4] + b[4], a[5] + b[5],
            a[6] + b[6], a[7] + b[7], a[8] + b[8]);
    }
    constexpr Matrix3 add(const Matrix3 &b) const
    {
        return add(*this, b);
    }

    constexpr Matrix3 invert(const Matrix3 &a) const
    {
        return mul(
            Matrix3(
                (m_[I(2, 2)] * m_[I(3, 3)] - m_[I(2, 3)] * m_[I(3, 2)]),
                (m_[I(1, 3)] * m_[I(3, 2)] - m_[I(1, 2)] * m_[I(3, 3)]),
                (m_[I(1, 2)] * m_[I(2, 3)] - m_[I(1, 3)] * m_[I(2, 2)]),
                (m_[I(2, 3)] * m_[I(3, 1)] - m_[I(2, 1)] * m_[I(3, 3)]),
                (m_[I(1, 1)] * m_[I(3, 3)] - m_[I(1, 3)] * m_[I(3, 1)]),
                (m_[I(1, 3)] * m_[I(2, 1)] - m_[I(1, 1)] * m_[I(2, 3)]),
                (m_[I(2, 1)] * m_[I(3, 2)] - m_[I(2, 2)] * m_[I(3, 1)]),
                (m_[I(1, 2)] * m_[I(3, 1)] - m_[I(1, 1)] * m_[I(3, 2)]),
                (m_[I(1, 1)] * m_[I(2, 2)] - m_[I(1, 2)] * m_[I(2, 1)])),
            1.f / det(a));
    }
    constexpr Matrix3 invert(void) const
    {
        return invert(*this);
    }
};

class Tristimulus
{
  public:
    Vector3 v_;
    constexpr Tristimulus(const float &a, const float &b, const float &c) : v_(a, b, c)
    {
        ;
    }
    constexpr Tristimulus() : v_(0, 0, 0) { ; }
    constexpr Tristimulus(const Vector3 &v) : v_(v) { ; }

    constexpr float operator[](const int &i) const
    {
        return v_[i];
    }
    constexpr const Vector3 &vec3(void) const { return v_; }
};

class Gamut
{
  public:
    Matrix3 toXYZ_;
    Matrix3 fromXYZ_;

    constexpr float z_from_xy(const float &x, const float &y) const { return 1 - x - y; }
    constexpr float X_from_xyY(const float &x, const float &y, const float &Y) const { return x * Y / y; }
    constexpr float Y_from_xyY(const float &x, const float &y, const float &Y) const { return Y; }
    constexpr float Z_from_xyY(const float &x, const float &y, const float &Y) const { return z_from_xy(x, y) * Y / y; }
    constexpr Matrix3 primMat(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB) const
    {
        return Matrix3(xR, xG, xB, yR, yG, yB, z_from_xy(xR, yR), z_from_xy(xG, yG), z_from_xy(xB, yB));
    }
    constexpr Matrix3 diag(const Vector3 &v) const
    {
        return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]);
    }
    constexpr Matrix3 mulDiag(const Matrix3 &m, const Vector3 &v) const
    {
        return m.mul(diag(m.invert().apply(v)));
    }
    constexpr Matrix3 fromPrimaries(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB, const float &xW, const float &yW) const
    {
        return mulDiag(primMat(xR, yR, xG, yG, xB, yB), Vector3(X_from_xyY(xW, yW, 1.f), Y_from_xyY(xW, yW, 1.f), Z_from_xyY(xW, yW, 1.f)));
    }

    constexpr Gamut(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB, const float &xW, const float &yW)
        : toXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW)), fromXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW).invert())
    {
        ;
    }
    const Matrix3 &toXYZ(void) const { return toXYZ_; }
    const Matrix3 &fromXYZ(void) const { return fromXYZ_; }

    constexpr Tristimulus toXYZ(const Tristimulus &tri) const
    {
        return Tristimulus(toXYZ_.apply(tri.vec3()));
    }
    constexpr Tristimulus fromXYZ(const Tristimulus &tri) const
    {
        return Tristimulus(fromXYZ_.apply(tri.vec3()));
    }
};

class OTF
{
  public:
    typedef enum {
        OTF_LINEAR,
        OTF_GAMMA, // simplest gamma
        OTF_SRGB,
        OTF_ST2084,
        OTF_HLG // Hybrid-log-gamma
    } TYPE;
    // TODO
    OTF()
    {
        ;
    }
};

class Spectrum
{
  public:
    typedef std::array<float, 400> spectrum; // 380-780, 1nm, fixed.
    spectrum s_;

    Spectrum()
    {
        ;
    }

    Spectrum blackbody(const float temp) const
    {
    }

    Tristimulus toXYZ(const Spectrum &) const
    {
        ;
    }
};

class Delta
{
  public:
    Delta()
    {
        ;
    }
    static float deltaAB()
    {
    }
};

void
dump(const Matrix3 &m)
{
    printf("-------\n");
    printf("%f,%f,%f\n", m[0], m[1], m[2]);
    printf("%f,%f,%f\n", m[3], m[4], m[5]);
    printf("%f,%f,%f\n", m[6], m[7], m[8]);
}
void dump(const Vector3 &v)
{
    printf("v-------\n");
    printf("%f,%f,%f\n", v[0], v[1], v[2]);
}
//
} // namespace ColorSystem
