//
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
//
// policy: currently single-headered, simple, low overhead, no ICC support.
//
#pragma once
#ifndef colorsystem_hpp__abf47c16efbc4a80838738ff9b8a0eea
#define colorsystem_hpp__abf47c16efbc4a80838738ff9b8a0eea 1

#include <algorithm>
#include <array>
#include <assert.h>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <vector>

namespace ColorSystem
{
namespace util
{
    namespace Detail
    {
        float constexpr sqrtNewtonRaphsonF(float x, float curr, float prev)
        {
            return curr == prev ? curr : sqrtNewtonRaphsonF(x, 0.5f * (curr + x / curr), curr);
        }
    } // namespace Detail

    /*
     * Constexpr version of the square root
     * Return value:
     *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
     *   - Otherwise, returns NaN
     */
    // https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
    float constexpr sqrtf(float x)
    {
        return x >= 0.f && x < std::numeric_limits<float>::infinity() ? Detail::sqrtNewtonRaphsonF(x, x, 0.f)
                                                                      : std::numeric_limits<float>::quiet_NaN();
    }
} // namespace util

static const float PI = 3.14159265358979323846f;

class Vector3
{
  public:
    typedef std::array<float, 3> vec3;
    vec3                         v_;
    constexpr Vector3(const float a, const float b, const float c) : v_({a, b, c}) { ; }
    constexpr float        operator[](const int &i) const { return v_[i]; }
    constexpr float        x() const { return v_[0]; }
    constexpr float        y() const { return v_[1]; }
    constexpr float        z() const { return v_[2]; }
    constexpr auto         size() const { return v_.size(); }
    auto                   begin() const { return v_.begin(); }
    auto                   end() const { return v_.end(); }
    constexpr float        dot(const Vector3 &a) const { return v_[0] * a[0] + v_[1] * a[1] + v_[2] * a[2]; }
    static constexpr float dot(const Vector3 &a, const Vector3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
};

class Matrix3
{
  public:
    typedef std::array<float, 9> matrix;
    matrix                       m_;

  private:
    static constexpr int M(const int x, const int y) { return x + y * 3; }
    static constexpr int I(const int y, const int x) { return ((x - 1) + (y - 1) * 3); }

  public:
    constexpr Matrix3(const float &a00, const float &a01, const float &a02, const float &a10, const float &a11,
        const float &a12, const float &a20, const float &a21, const float &a22)
        : m_({a00, a01, a02, a10, a11, a12, a20, a21, a22})
    {
        ;
    }
    constexpr Matrix3(void) : m_({1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f}) { ; }
    static constexpr Matrix3 fromArray(const float *p)
    {
        return Matrix3(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
    }
    static constexpr Matrix3 diag(const Vector3 &v) { return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]); }

    constexpr float operator[](const int &i) const { return m_[i]; }

    constexpr Vector3        row(const int i) const { return Vector3(m_[M(0, i)], m_[M(1, i)], m_[M(2, i)]); }
    constexpr Vector3        col(const int i) const { return Vector3(m_[M(i, 0)], m_[M(i, 1)], m_[M(i, 2)]); }
    static constexpr Matrix3 transpose(const Matrix3 &a)
    {
        return Matrix3(
            a[M(0, 0)], a[M(0, 1)], a[M(0, 2)], a[M(1, 0)], a[M(1, 1)], a[M(1, 2)], a[M(2, 0)], a[M(2, 1)], a[M(2, 2)]);
    }
    constexpr Matrix3 transpose(void) const { return transpose(*this); }

    static constexpr Matrix3 mul(const Matrix3 &a, const Matrix3 &b)
    {
        return Matrix3(a[M(0, 0)] * b[M(0, 0)] + a[M(1, 0)] * b[M(0, 1)] + a[M(2, 0)] * b[M(0, 2)],
            a[M(0, 0)] * b[M(1, 0)] + a[M(1, 0)] * b[M(1, 1)] + a[M(2, 0)] * b[M(1, 2)],
            a[M(0, 0)] * b[M(2, 0)] + a[M(1, 0)] * b[M(2, 1)] + a[M(2, 0)] * b[M(2, 2)],

            a[M(0, 1)] * b[M(0, 0)] + a[M(1, 1)] * b[M(0, 1)] + a[M(2, 1)] * b[M(0, 2)],
            a[M(0, 1)] * b[M(1, 0)] + a[M(1, 1)] * b[M(1, 1)] + a[M(2, 1)] * b[M(1, 2)],
            a[M(0, 1)] * b[M(2, 0)] + a[M(1, 1)] * b[M(2, 1)] + a[M(2, 1)] * b[M(2, 2)],

            a[M(0, 2)] * b[M(0, 0)] + a[M(1, 2)] * b[M(0, 1)] + a[M(2, 2)] * b[M(0, 2)],
            a[M(0, 2)] * b[M(1, 0)] + a[M(1, 2)] * b[M(1, 1)] + a[M(2, 2)] * b[M(1, 2)],
            a[M(0, 2)] * b[M(2, 0)] + a[M(1, 2)] * b[M(2, 1)] + a[M(2, 2)] * b[M(2, 2)]);
    }
    constexpr Matrix3 mul(const Matrix3 &b) const { return mul(*this, b); }

    static constexpr Vector3 apply(const Matrix3 &m, const Vector3 &v)
    {
        return Vector3(m[M(0, 0)] * v[0] + m[M(1, 0)] * v[1] + m[M(2, 0)] * v[2],
            m[M(0, 1)] * v[0] + m[M(1, 1)] * v[1] + m[M(2, 1)] * v[2],
            m[M(0, 2)] * v[0] + m[M(1, 2)] * v[1] + m[M(2, 2)] * v[2]);
    }
    constexpr Vector3 apply(const Vector3 &v) const { return apply(*this, v); }

    static constexpr float det(const Matrix3 &m)
    {
        return m[I(1, 1)] * m[I(2, 2)] * m[I(3, 3)] + m[I(2, 1)] * m[I(3, 2)] * m[I(1, 3)] +
               m[I(3, 1)] * m[I(1, 2)] * m[I(2, 3)] - m[I(1, 1)] * m[I(3, 2)] * m[I(2, 3)] -
               m[I(3, 1)] * m[I(2, 2)] * m[I(1, 3)] - m[I(2, 1)] * m[I(1, 2)] * m[I(3, 3)];
    }
    constexpr float det(void) const { return det(*this); }

    static constexpr Matrix3 mul(const Matrix3 &m, const float &a)
    {
        return Matrix3(m[0] * a, m[1] * a, m[2] * a, m[3] * a, m[4] * a, m[5] * a, m[6] * a, m[7] * a, m[8] * a);
    }
    constexpr Matrix3        mul(const float a) const { return mul(*this, a); }
    static constexpr Matrix3 div(const Matrix3 &m, const float &a) { return mul(m, 1.f / a); }
    constexpr Matrix3        div(const float a) const { return div(*this, a); }

    static constexpr Matrix3 add(const Matrix3 &a, const Matrix3 &b)
    {
        return Matrix3(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4], a[5] + b[5], a[6] + b[6],
            a[7] + b[7], a[8] + b[8]);
    }
    constexpr Matrix3 add(const Matrix3 &b) const { return add(*this, b); }

    static constexpr Matrix3 invert(const Matrix3 &a)
    {
        return mul(Matrix3((a.m_[I(2, 2)] * a.m_[I(3, 3)] - a.m_[I(2, 3)] * a.m_[I(3, 2)]),
                       (a.m_[I(1, 3)] * a.m_[I(3, 2)] - a.m_[I(1, 2)] * a.m_[I(3, 3)]),
                       (a.m_[I(1, 2)] * a.m_[I(2, 3)] - a.m_[I(1, 3)] * a.m_[I(2, 2)]),
                       (a.m_[I(2, 3)] * a.m_[I(3, 1)] - a.m_[I(2, 1)] * a.m_[I(3, 3)]),
                       (a.m_[I(1, 1)] * a.m_[I(3, 3)] - a.m_[I(1, 3)] * a.m_[I(3, 1)]),
                       (a.m_[I(1, 3)] * a.m_[I(2, 1)] - a.m_[I(1, 1)] * a.m_[I(2, 3)]),
                       (a.m_[I(2, 1)] * a.m_[I(3, 2)] - a.m_[I(2, 2)] * a.m_[I(3, 1)]),
                       (a.m_[I(1, 2)] * a.m_[I(3, 1)] - a.m_[I(1, 1)] * a.m_[I(3, 2)]),
                       (a.m_[I(1, 1)] * a.m_[I(2, 2)] - a.m_[I(1, 2)] * a.m_[I(2, 1)])),
            1.f / det(a));
    }
    constexpr Matrix3 invert(void) const { return invert(*this); }
    constexpr auto    size() const { return m_.size(); }
    auto              begin() const { return m_.begin(); }
    auto              end() const { return m_.end(); }
};

class Tristimulus
{
  public:
    Vector3 v_;

    constexpr Tristimulus(const float &a, const float &b, const float &c) : v_(a, b, c) { ; }
    constexpr Tristimulus() : v_(0, 0, 0) { ; }
    constexpr Tristimulus(const Vector3 &v) : v_(v) { ; }
    constexpr Tristimulus(const float &v) : v_(v, v, v) { ; }

    constexpr float operator[](const int &i) const { return v_[i]; }
    constexpr auto  size() const { return v_.size(); }
    auto            begin() const { return v_.begin(); }
    auto            end() const { return v_.end(); }

    constexpr const Vector3 &vec3(void) const { return v_; }
    constexpr float          a() const { return v_[0]; }
    constexpr float          b() const { return v_[1]; }
    constexpr float          c() const { return v_[2]; }

    static constexpr Tristimulus scale(const Tristimulus &t, const float &s)
    {
        return Tristimulus(t[0] * s, t[1] * s, t[2] * s);
    }
    constexpr Tristimulus scale(const float &s) const { return scale(*this, s); }
    constexpr Tristimulus operator*(const float &s) const { return scale(s); }
    constexpr Tristimulus operator/(const float &s) const { return scale(1.f / s); }

    // apply color transform matrix
    static constexpr Tristimulus apply(const Tristimulus &t, const Matrix3 &m)
    {
        return Tristimulus(m.apply(t.vec3()));
    }
    constexpr Tristimulus apply(const Matrix3 &m) const { return apply(*this, m); }

    static constexpr Tristimulus add(const Tristimulus &a, const Tristimulus &b)
    {
        return Tristimulus(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
    }
    constexpr const Tristimulus add(const Tristimulus &b) const { return add(*this, b); }
    constexpr const Tristimulus operator+(const Tristimulus &b) const { return add(*this, b); }
    // per-element re-lighting
    static constexpr Tristimulus mul(const Tristimulus &a, const Tristimulus &b)
    {
        return Tristimulus(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
    }
    constexpr const Tristimulus mul(const Tristimulus &b) const { return mul(*this, b); }
    constexpr const Tristimulus operator*(const Tristimulus &b) const { return mul(*this, b); }
    static constexpr float      dot(const Tristimulus &a, const Tristimulus &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    constexpr float              dot(const Tristimulus &b) const { return dot(*this, b); }
    static constexpr const float mini(const float &a, const float &b) { return (a < b) ? a : b; }
    static constexpr const float maxi(const float &a, const float &b) { return (a > b) ? a : b; }
    static constexpr Tristimulus min(const Tristimulus &a, const Tristimulus &b)
    {
        return Tristimulus(mini(a[0], b[0]), mini(a[1], b[1]), mini(a[2], b[2]));
    }
    static constexpr Tristimulus max(const Tristimulus &a, const Tristimulus &b)
    {
        return Tristimulus(maxi(a[0], b[0]), maxi(a[1], b[1]), maxi(a[2], b[2]));
    }
    constexpr Tristimulus min(const Tristimulus &a) const { return min(*this, a); }
    constexpr Tristimulus max(const Tristimulus &a) const { return max(*this, a); }

    constexpr float min3(void) const { return mini(mini(a(), b()), c()); }
    constexpr float max3(void) const { return maxi(maxi(a(), b()), c()); }

    constexpr Tristimulus clip(const float &l, const float &h) const
    {
        return max(min(*this, Tristimulus(h)), Tristimulus(l));
    }
    static constexpr Tristimulus clip(const Tristimulus &t, const float &l, const float &h) { return t.clip(l, h); }
    constexpr Tristimulus        positive() const { return max(*this, Tristimulus(0.f)); }
    static constexpr Tristimulus positive(const Tristimulus &t) { return t.positive(); }
    constexpr bool               isNegative(const float &a) const { return (a < 0.f); }
    constexpr bool         hasNegative() const { return isNegative(v_[0]) || isNegative(v_[1]) || isNegative(v_[2]); }
    static constexpr bool  hasNegative(const Tristimulus &t) { return t.hasNegative(); }
    static constexpr float abs_f(const float f) { return (f < 0.f) ? -f : f; }

    static constexpr float z_from_xy(const float &x, const float &y) { return 1 - x - y; }
    static constexpr float X_from_Yxy(const float &Y, const float &x, const float &y)
    {
        return (abs_f(y) < 1e-8f) ? 0.f : x * Y / y;
    }
    static constexpr float Y_from_Yxy(const float &Y, const float &x, const float &y)
    {
        (void)x;
        (void)y;
        return Y;
    }
    static constexpr float Z_from_Yxy(const float &Y, const float &x, const float &y)
    {
        return (abs_f(y) < 1e-8f) ? 0.f : z_from_xy(x, y) * Y / y;
    }
    static constexpr float Y_from_XYZ(const float &x, const float &y, const float &z)
    {
        (void)x;
        (void)z;
        return y;
    }
    static constexpr float x_from_XYZ(const float &x, const float &y, const float &z)
    {
        return (abs_f(x + y + z) < 1e-8f) ? 0.3127f : x / (x + y + z);
    }
    static constexpr float y_from_XYZ(const float &x, const float &y, const float &z)
    {
        return (abs_f(x + y + z) < 1e-8f) ? 0.3290f : y / (x + y + z);
    }
    static constexpr Tristimulus fromYxy(const float &Y, const float &x, const float &y)
    {
        return Tristimulus(X_from_Yxy(Y, x, y), Y_from_Yxy(Y, x, y), Z_from_Yxy(Y, x, y));
    }
    static constexpr Tristimulus toYxy(const float &X, const float &Y, const float &Z)
    {
        return Tristimulus(Y_from_XYZ(X, Y, Z), x_from_XYZ(X, Y, Z), y_from_XYZ(X, Y, Z));
    }
    static constexpr Tristimulus fromYxy(const Tristimulus &Yxy) { return fromYxy(Yxy[0], Yxy[1], Yxy[2]); }
    static constexpr Tristimulus toYxy(const Tristimulus &XYZ) { return toYxy(XYZ[0], XYZ[1], XYZ[2]); }
    //
    constexpr Tristimulus fromYxy(void) const { return fromYxy(*this); }
    constexpr Tristimulus toYxy(void) const { return toYxy(*this); }
    // support u'v'
    //  u = 4x / (-2x + 12y + 3)
    //  v = 6y / (-2x + 12y + 3)
    //  x = 3u / (2u - 8v + 4)
    //  y = 2v / (2u - 8v + 4)
    //  u' = 4x / (-2x + 12y + 3)    [ = u ]
    //  v' = 9y / (-2x + 12y + 3)    [ = 1.5v ]
    //  x = 9u' / (6u' - 16v' + 12)
    //  y = 4v' / (6u' - 16v' + 12)
    static constexpr float u_from_xy(const float &x, const float &y) { return 4.f * x / (-2.f * x + 12.f * y + 3.f); }
    static constexpr float v_from_xy(const float &x, const float &y) { return 6.f * y / (-2.f * x + 12.f * y + 3.f); }
    static constexpr float x_from_uv(const float &u, const float &v) { return 3.f * u / (2.f * u - 8.f * v + 4.f); }
    static constexpr float y_from_uv(const float &u, const float &v) { return 2.f * v / (2.f * u - 8.f * v + 4.f); }
    static constexpr Tristimulus YxyToYuv(const Tristimulus &Yxy)
    {
        return Tristimulus(Yxy[0], u_from_xy(Yxy[1], Yxy[2]), v_from_xy(Yxy[1], Yxy[2]));
    }
    static constexpr Tristimulus YuvToYxy(const Tristimulus &Yuv)
    {
        return Tristimulus(Yuv[0], x_from_uv(Yuv[1], Yuv[2]), y_from_uv(Yuv[1], Yuv[2]));
    }
    static constexpr Tristimulus toYuv(const Tristimulus &XYZ) { return YxyToYuv(toYxy(XYZ)); }
    static constexpr Tristimulus toYuv(const float &X, const float &Y, const float &Z)
    {
        return toYuv(Tristimulus(X, Y, Z));
    }
    static constexpr Tristimulus fromYuv(const Tristimulus &Yuv) { return fromYxy(YuvToYxy(Yuv)); }
    static constexpr Tristimulus fromYuv(const float &X, const float &Y, const float &Z)
    {
        return fromYuv(Tristimulus(X, Y, Z));
    }
    constexpr Tristimulus toYuv(void) const { return toYuv(*this); }
    constexpr Tristimulus fromYuv(void) const { return fromYuv(*this); }

    // uv only used in CCT.
    static constexpr float blackbody_x_approx(const float &T)
    {
        return (float)(0.811973208 +
                       T * (-0.00016747 + T * (0.00000000331067 +
                                                  T * (0.00000000000690001 +
                                                          T * (-1.39671E-15 + T * (1.11951E-19 + T * -3.31672E-24))))));
    }
    static constexpr float blackbody_y_approx(const float &T)
    {
        return (float)(0.136977825 +
                       T * (0.000317653 + T * (-0.000000127946 +
                                                  T * (0.0000000000222415 + T * (-1.81438E-15 + T * 5.67564E-20)))));
    }
    // XYZ from Color Temperature with Y.
    // cubic spline approx https://en.wikipedia.org/wiki/Planckian_locus#Approximation
    static constexpr Tristimulus fromPlanckianLocus(const float &T, const float &Y = 1.f)
    {
        const float x =
            (T < 4000.f) ? ((-0.2661239e9f / (T * T * T)) - (0.2343580e6f / (T * T)) + (0.8776956e3f / T) + 0.179910f)
                         : ((-3.0258469e9f / (T * T * T)) + (2.1070379e6f / (T * T)) + (0.2226347e3f / T) + 0.240390f);
        const float y =
            (T < 2222.f)
                ? ((-1.1063814f * x * x * x) - (1.34811020f * x * x) + (2.18555832f * x) - 0.20219683f)
                : ((T < 4000.f) ? ((-1.1063814f * x * x * x - 1.34811020f * x * x + 2.18555832f * x - 0.20219683f))
                                : ((3.0817580f * x * x * x - 5.87338670f * x * x + 3.75112997f * x - 0.37001483f)));
        return Tristimulus(Y, x, y).fromYxy();
    }
    static constexpr Tristimulus fromCT(const float &T, const float Y = 1.f)
    {
        return Tristimulus(Y, blackbody_x_approx(T), blackbody_y_approx(T)).fromYxy();
    }
    static constexpr float CCT_x_approx(const float &T, const float &dUV)
    {
        const float x1 = blackbody_x_approx(T);
        const float y1 = blackbody_y_approx(T);
        const float u1 = u_from_xy(x1, y1);
        const float v1 = v_from_xy(x1, y1);
        const float x2 = blackbody_x_approx(T - 1.f);
        const float y2 = blackbody_y_approx(T - 1.f);
        const float u2 = u_from_xy(x2, y2);
        const float v2 = v_from_xy(x2, y2);
        const float du = u2 - u1;
        const float dv = v2 - v1;
        const float l  = util::sqrtf(du * du + dv * dv);
        return x_from_uv(u1 - dUV * du / l, v1 + dUV * dv / l);
    }
    static constexpr float CCT_y_approx(const float &T, const float &dUV)
    {
        const float x1 = blackbody_x_approx(T);
        const float y1 = blackbody_y_approx(T);
        const float u1 = u_from_xy(x1, y1);
        const float v1 = v_from_xy(x1, y1);
        const float x2 = blackbody_x_approx(T - 1.f);
        const float y2 = blackbody_y_approx(T - 1.f);
        const float u2 = u_from_xy(x2, y2);
        const float v2 = v_from_xy(x2, y2);
        const float du = u2 - u1;
        const float dv = v2 - v1;
        const float l  = util::sqrtf(du * du + dv * dv);
        return y_from_uv(u1 - dUV * du / l, v1 + dUV * dv / l);
    }
    static constexpr Tristimulus fromCCT(const float &T, const float &dUV, const float Y = 1.f)
    {
        return Tristimulus(Y, CCT_x_approx(T, dUV), CCT_y_approx(T, dUV)).fromYuv();
    }

    // Lab
    static constexpr float CIELAB_curve(const float &f)
    {
        const float threshold = 216.f / 24389.0f;
        const float K         = 24389.0f / 27.0f;
        return (f > 1.f) ? 1.f : ((f > threshold) ? powf(f, 1.f / 3.f) : ((K * f + 16.f) / 116.f));
    }
    static constexpr float CIELAB_decurve(const float &f)
    {
        const float K = (3.f / 29.f) * (3.f / 29.f) * (3.f / 29.f);
        return (f > 1.f) ? 1.f : ((f > 6.f / 29.f) ? powf(f, 3.f) : (116.f * f - 16.f) * K);
    }
    static constexpr Tristimulus toCIELAB(const Tristimulus &t, const Tristimulus &white)
    {
        const float x0 = white[0];
        const float y0 = white[1];
        const float z0 = white[2];
        const float x1 = CIELAB_curve(t[0] / x0);
        const float y1 = CIELAB_curve(t[1] / y0);
        const float z1 = CIELAB_curve(t[2] / z0);
        return Tristimulus(116.f * y1 - 16.f, 500.f * (x1 - y1), 200.f * (y1 - z1));
    }
    static constexpr Tristimulus fromCIELAB(const Tristimulus &t, const Tristimulus &white)
    {
        const float x0 = white[0];
        const float y0 = white[1];
        const float z0 = white[2];
        const float fy = (t[0] + 16.f) / 116.f;
        const float fx = fy + (t[1] / 500.f);
        const float fz = fy - (t[2] / 200.f);
        return Tristimulus(CIELAB_decurve(fx) * x0, CIELAB_decurve(fy) * y0, CIELAB_decurve(fz) * z0);
    }
    constexpr Tristimulus toCIELAB(const Tristimulus &white) const { return toCIELAB(*this, white); }
    constexpr Tristimulus fromCIELAB(const Tristimulus &white) const { return fromCIELAB(*this, white); }
    // CIELAB uses D50 by default.
    constexpr Tristimulus toCIELAB(void) const { return toCIELAB(*this, Tristimulus(0.9642f, 1.0f, 0.8249f)); }
    constexpr Tristimulus fromCIELAB(void) const { return fromCIELAB(*this, Tristimulus(0.9642f, 1.0f, 0.8249f)); }
    // HSV
    static constexpr float mod360(const float &r)
    {
        return (r < 0.f) ? mod360(r + 360.f) : ((r > 360.f) ? mod360(r - 360.f) : r);
    }
    static Tristimulus toHSV_atan(const Tristimulus &t)
    {
        const float max = maxi(maxi(t[0], t[1]), t[2]);
        const float min = mini(mini(t[0], t[1]), t[2]);
        return Tristimulus(mod360((180.f / PI) * atan2f(sqrtf(3.f) * (t[1] - t[2]), 2.f * t[0] - t[1] - t[2])),
            (max == 0.f) ? 0.f : (max - min) / max, max);
    }
    Tristimulus                  toHSV_atan(void) const { return toHSV_atan(*this); }
    static constexpr Tristimulus toHSV(const Tristimulus &t)
    {
        const float max = maxi(maxi(t[0], t[1]), t[2]);
        const float min = mini(mini(t[0], t[1]), t[2]);
        const float d   = max - min;
        return Tristimulus(
            mod360(((max == min) ? 0.f
                                 : ((max == t[0]) ? (60.f * (t[1] - t[2]) / d)
                                                  : ((max == t[1]) ? (60.f * (t[2] - t[0]) / d + 120.f)
                                                                   : (60.f * (t[0] - t[1]) / d + 240.f))))),
            (max == 0.f) ? 0.f : d / max, max);
    }

    static constexpr Tristimulus fromHSV(const Tristimulus &t)
    {
        const float h = mod360(t[0]);
        const int   d = (int)(h / 60.f);
        const int   i = d % 6;
        const float r = (h / 60.f) - d;

        const float s = t[1];
        const float v = t[2];
        const float m = v * (1.0f - s);
        const float n = v * (1.0f - s * r);
        const float p = v * (1.0f - s * (1.0f - r));

        return (i == 0) ? Tristimulus(v, p, m)
                        : ((i == 1) ? Tristimulus(n, v, m)
                                    : ((i == 2) ? Tristimulus(m, v, p)
                                                : ((i == 3) ? Tristimulus(m, n, v)
                                                            : ((i == 4) ? Tristimulus(p, m, v)
                                                                        : ((i == 5) ? Tristimulus(v, m, n)
                                                                                    : Tristimulus(0, 0, 0))))));
    }

    constexpr Tristimulus toHSV(void) const { return toHSV(*this); }
    constexpr Tristimulus fromHSV(void) const { return fromHSV(*this); }
};

class Gamut
{
  public:
    const char *name_;
    Matrix3     toXYZ_;
    Matrix3     fromXYZ_;

    static constexpr Matrix3 primMat(
        const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB)
    {
        return Matrix3(xR, xG, xB, yR, yG, yB, Tristimulus::z_from_xy(xR, yR), Tristimulus::z_from_xy(xG, yG),
            Tristimulus::z_from_xy(xB, yB));
    }
    static constexpr Matrix3 diag(const Vector3 &v) { return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]); }
    static constexpr Matrix3 mulDiag(const Matrix3 &m, const Vector3 &v) { return m.mul(diag(m.invert().apply(v))); }
    static constexpr Matrix3 fromPrimaries(const float &xR, const float &yR, const float &xG, const float &yG,
        const float &xB, const float &yB, const float &xW, const float &yW)
    {
        return mulDiag(primMat(xR, yR, xG, yG, xB, yB), Tristimulus::fromYxy(1.f, xW, yW).vec3());
    }
    constexpr Gamut(const char *name, const Matrix3 &fromXYZ) : name_(name), toXYZ_(fromXYZ.invert()), fromXYZ_(fromXYZ)
    {
        ;
    }
    constexpr Gamut(const char *name, const float &xR, const float &yR, const float &xG, const float &yG,
        const float &xB, const float &yB, const float &xW, const float &yW)
        : name_(name),
          toXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW)),
          fromXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW).invert())
    {
        ;
    }
    const char *name(void) const { return name_; }

    constexpr Matrix3     toXYZ(void) const { return toXYZ_; }
    constexpr Matrix3     fromXYZ(void) const { return fromXYZ_; }
    constexpr Tristimulus toXYZ(const Tristimulus &tri) const { return Tristimulus(toXYZ_.apply(tri.vec3())); }
    constexpr Tristimulus fromXYZ(const Tristimulus &tri) const { return Tristimulus(fromXYZ_.apply(tri.vec3())); }

    constexpr Vector3 primaryVector(void) const
    {
        return Vector3(
            toXYZ_[0] + toXYZ_[3] + toXYZ_[6], toXYZ_[1] + toXYZ_[4] + toXYZ_[7], toXYZ_[2] + toXYZ_[5] + toXYZ_[8]);
    }
    constexpr Matrix3 primaryMatrix(void) const
    {
        const Vector3 t(primaryVector());
        return Matrix3(toXYZ_[0] / t[0], toXYZ_[1] / t[1], toXYZ_[2] / t[2], toXYZ_[3] / t[0], toXYZ_[4] / t[1],
            toXYZ_[5] / t[2], toXYZ_[6] / t[0], toXYZ_[7] / t[1], toXYZ_[8] / t[2]);
    }
    constexpr Tristimulus primaryWhite() const { return Tristimulus(primaryMatrix().apply(primaryVector())); }
    constexpr Tristimulus primaryRed() const
    {
        const Matrix3 n(primaryMatrix());
        return Tristimulus(n[0], n[3], n[6]);
    }
    constexpr Tristimulus primaryGreen() const
    {
        const Matrix3 n(primaryMatrix());
        return Tristimulus(n[1], n[4], n[7]);
    }
    constexpr Tristimulus primaryBlue() const
    {
        const Matrix3 n(primaryMatrix());
        return Tristimulus(n[2], n[5], n[8]);
    }
};
class OTF
{
  public:
    typedef enum
    {
        LINEAR,
        GAMMA, // simplest gamma
        SRGB,
        BT709,
        ST2084,
        SLOG2,
        HLG // Hybrid-log-gamma
    } TYPE;

    static float       gamma(const float &v, const float &g) { return powf(v, 1.f / g); }
    static float       degamma(const float &v, const float &g) { return powf(v, g); }
    static const float ST2084_to_Y(const float &pixel) // pixel should be 0-1
    {
        const float pq_m1 = 0.1593017578125f; // ( 2610.0 / 4096.0 ) / 4.0;
        const float pq_m2 = 78.84375f;        // ( 2523.0 / 4096.0 ) * 128.0;
        const float pq_c1 = 0.8359375f;       // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
        const float pq_c2 = 18.8515625f;      // ( 2413.0 / 4096.0 ) * 32.0;
        const float pq_c3 = 18.6875f;         // ( 2392.0 / 4096.0 ) * 32.0;
        const float pq_C  = 100.0f;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this assumes full range (0 - 1)
        float Np = powf(pixel, 1.0f / pq_m2);
        float L  = Np - pq_c1;
        if (L < 0.0)
            L = 0.0;
        L = L / (pq_c2 - pq_c3 * Np);
        L = powf(L, 1.0f / pq_m1);
        return L * pq_C; // returns 0-100, 1=100cd/m^2
    }

    static const float Y_to_ST2084(const float &C) // C should be 0-100, 1=100cd/m^2
    {
        if (C <= 0.f)
            return 0.f;
        if (C >= 100.f)
            return 1.f;
        const float pq_m1 = 0.1593017578125f; // ( 2610.0 / 4096.0 ) / 4.0;
        const float pq_m2 = 78.84375f;        // ( 2523.0 / 4096.0 ) * 128.0;
        const float pq_c1 = 0.8359375f;       // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
        const float pq_c2 = 18.8515625f;      // ( 2413.0 / 4096.0 ) * 32.0;
        const float pq_c3 = 18.6875f;         // ( 2392.0 / 4096.0 ) * 32.0;
        const float pq_C  = 100.0f;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this returns full range (0 - 1)
        float L  = C / pq_C;
        float Lm = powf(L, pq_m1);
        float N  = (pq_c1 + pq_c2 * Lm) / (1.0f + pq_c3 * Lm);
        N        = powf(N, pq_m2);
        return N;
    }
    static const float Y_to_sRGB(const float &C) // returns signal, 0-1, input 0-1
    {
        return (C < 0.f)
                   ? 0.f
                   : ((C > 1.f) ? 1.f : ((C < 0.0031308f) ? C * 12.92f : (1.055f * powf(C, 1.0f / 2.4f) - 0.055f)));
    }
    static const float sRGB_to_Y(const float &C) // returns 0-1, 1=100 nits
    {
        return (C < 0.f) ? 0.f : ((C > 1.f) ? 1.f : ((C < 0.04045f) ? C / 12.92f : powf((C + 0.055f) / 1.055f, 2.4f)));
    }
    static const float Y_to_BT709(const float &C) // returns signal, 0-1, input 0-1
    {
        return (C < 0.f) ? 0.f : ((C > 1.f) ? 1.f : ((C < 0.018f) ? C * 4.50f : (1.099f * powf(C, 0.45f) - 0.099f)));
    }
    static const float BT709_to_Y(const float &C) // returns nits, 0-100[cd/m^2]
    {
        return (C < 0.f) ? 0.f
                         : ((C > 1.f) ? 1.f : ((C < 0.081f) ? C / 4.50f : powf((C + 0.099f) / 1.099f, 1.f / 0.45f)));
    }

    static const float Y_to_HLG(const float &C)
    {
        const float a = 0.17883277f;
        const float b = 0.28466892f;
        const float c = 0.55991073f;
        return (C < 0.f) ? 0.f : ((C < 1.f) ? (0.5f * sqrtf(C)) : (a * logf(C - b) + c));
    }

    static const float HLG_to_Y(const float &C)
    {
        const float a = 0.17883277f;
        const float b = 0.28466892f;
        const float c = 0.55991073f;
        return (C < 0.f) ? 0.f : ((C <= 0.5f) ? (4.f * C * C) : exp((C - c) / a) + b);
    }

    static const float CV_to_IRE_SLog2(const float &cv)
    {
        const float BLACK = 64.f / 1024.f;
        const float WV    = 876.f / 1024.f; // 940-64
        return (cv - BLACK) / WV;
    }
    static const float IRE_to_CV_SLog2(const float &ire)
    {
        const float BLACK = 64.f / 1024.f;
        const float WV    = 876.f / 1024.f; // 940-64
        return (ire * WV) + BLACK;
    }
    static const float Y_to_SLog2(const float &x) // returns signal, 0-1, input 0-1
    {
        const float y = (x < 0.f) ? x * 3.53881278538813f + 0.030001222851889303f
                                  : (0.432699f * log10f(155.0f * x / 219.0f + 0.037584f) + 0.616596f) + 0.03f;
        return IRE_to_CV_SLog2(y);
    }
    static const float SLog2_to_Y(const float &C) // returns 0-1, 1=100cd/m^2
    {
        const float x = CV_to_IRE_SLog2(C);
        const float y = (x >= 0.030001222851889303f)
                            ? 219.0f * (powf(10.0f, ((x - 0.616596f - 0.03f) / 0.432699f)) - 0.037584f) / 155.0f
                            : (x - 0.030001222851889303f) / 3.53881278538813f;
        return (y > 0.f) ? y : 0.f;
    }
    static const Tristimulus toScreen(TYPE type, const Tristimulus &scene, const float g = 1.f)
    {
        switch (type)
        {
        case GAMMA:
        {
            return Tristimulus(gamma(scene[0], g), gamma(scene[1], g), gamma(scene[2], g));
        }
        break;
        case SRGB:
        {
            return Tristimulus(Y_to_sRGB(scene[0]), Y_to_sRGB(scene[1]), Y_to_sRGB(scene[2]));
        }
        break;
        case BT709:
        {
            return Tristimulus(Y_to_BT709(scene[0]), Y_to_BT709(scene[1]), Y_to_BT709(scene[2]));
        }
        break;
        case ST2084:
        {
            return Tristimulus(Y_to_ST2084(scene[0]), Y_to_ST2084(scene[1]), Y_to_ST2084(scene[2]));
        }
        break;
        case SLOG2:
        {
            return Tristimulus(Y_to_SLog2(scene[0]), Y_to_SLog2(scene[1]), Y_to_SLog2(scene[2]));
        }
        break;
        case HLG:
        {
            return Tristimulus(Y_to_HLG(scene[0]), Y_to_HLG(scene[1]), Y_to_HLG(scene[2]));
        }
        case LINEAR:
        default:
            return scene;
        }
    }
    static const Tristimulus toScene(TYPE type, const Tristimulus &screen, const float g = 1.f)
    {
        switch (type)
        {
        case GAMMA:
        {
            return Tristimulus(degamma(screen[0], g), degamma(screen[1], g), degamma(screen[2], g));
        }
        break;
        case SRGB:
        {
            return Tristimulus(sRGB_to_Y(screen[0]), sRGB_to_Y(screen[1]), sRGB_to_Y(screen[2]));
        }
        break;
        case BT709:
        {
            return Tristimulus(BT709_to_Y(screen[0]), BT709_to_Y(screen[1]), BT709_to_Y(screen[2]));
        }
        break;
        case ST2084:
        {
            return Tristimulus(ST2084_to_Y(screen[0]), ST2084_to_Y(screen[1]), ST2084_to_Y(screen[2]));
        }
        case SLOG2:
        {
            return Tristimulus(SLog2_to_Y(screen[0]), SLog2_to_Y(screen[1]), SLog2_to_Y(screen[2]));
        }
        case HLG:
        {
            return Tristimulus(HLG_to_Y(screen[0]), HLG_to_Y(screen[1]), HLG_to_Y(screen[2]));
        }
        break;
        case LINEAR:
        default:
            return screen;
        }
    }
};

class MemoryStream
{
  public:
    typedef size_t Pointer;
    typedef enum
    {
        BUFFERBYTES = 1,
        BUFFERBITS  = 8,
    } PARAM;
    typedef enum
    {
        TOP,
        CURRENT,
        END
    } CURSOR;
    typedef enum
    {
        LITTLE,
        BIG,
    } ENDIAN;

  public:
    MemoryStream() { init(); }
    MemoryStream(const void *data,size_t sz){init(data,sz);}
    std::vector<uint8_t>  buffer_;
    std::vector<uint8_t> &buffer() { return buffer_; }

    int     ptr_;
    ENDIAN  readEndian_;
    ENDIAN  writeEndian_;

    void init(const void *data = NULL,size_t sz = 0)
    {
        if (data)
        {
            buffer_.resize(sz);
            memcpy(&(buffer_.front()), data, sz);
        }
        else
        {
            buffer_.resize(sz);
        }
        ptr_         = 0;
        readEndian_  = LITTLE;
        writeEndian_ = LITTLE;
    }

    size_t  ftell(void) { return (size_t)ptr_; }
    Pointer mark(void) { return ptr_; }
    int     fseek(int offs, CURSOR cursor)
    {
        int ptr_p = ptr_;
        switch (cursor)
        {
        case TOP:
        {
            ptr_ = offs;
        }
        break;
        case CURRENT:
        {
            ptr_ += offs;
        }
        break;
        case END:
        {
            ptr_ = (int)buffer_.size() - offs - 1;
        }
        break;
        }
        assert(ptr_ >= 0);
        assert(ptr_ <= buffer_.size());
        return ptr_;
    }

    void seekTo(Pointer p) { fseek((int)p, TOP); }

    void *ptr(void) { return (void *)&(buffer_[ptr_]); }
    void *ptr(Pointer p) { return (void *)&(buffer_[p]); }
    bool  feof(void) { return ptr_ >= buffer_.size(); }

    // read absolute.
    uint8_t get_at(size_t abs)
    {
        if (abs < buffer_.size())
            return buffer_[abs];
        return 0;
    }
    void put_at(size_t abs, uint8_t val)
    {
        if (abs < buffer_.size())
            buffer_.resize(abs);
        buffer_[abs] = val;
    }
    size_t read_at(void *buffer, size_t s, Pointer p)
    {
        int r      = (int)buffer_.size() - (int)p;
        if( r < 0 )
            return 0;
        size_t remain = (s < r) ? s : r;
        memcpy(buffer, ptr(p), remain);
        return remain;
    }
    size_t write_at(const void *buffer, size_t s, Pointer p)
    {
        size_t s1 = buffer_.size();
        size_t s2 = p + s;
        if (s2 > s1)
            buffer_.resize(s2);
        memcpy(ptr(p), buffer, s);
        return s;
    }
    size_t read(void *buffer, size_t size)
    {
        size_t sz = read_at(buffer, size, ptr_);
        ptr_ += (int)sz;
        return sz;
    }
    size_t write(const void *buffer, size_t size)
    {
        size_t s = size;
        write_at(buffer, s, ptr_);
        ptr_ += (int)s;
        return s;
    }
    uint8_t getUint8(void)
    {
        if (ptr_ < buffer_.size())
        {
            uint8_t c = buffer_[ptr_];
            ptr_++;
            return c;
        }
        return 0;
    }
    void putUint8(const uint8_t uc)
    {
        size_t s = buffer_.size();
        if (ptr_ < s)
        {
            buffer_[ptr_] = uc;
            return;
        }
        else if (ptr_ == s)
        {
            buffer_.push_back(uc);
            ptr_++;
            return;
        }
        assert(false);
    }
    void setReadEndian(const ENDIAN e) { readEndian_ = e; }
    void setWriteEndian(const ENDIAN e) { writeEndian_ = e; }
    const int8_t getInt8(void)
    {
        uint8_t uc = getUint8();
        return *(int8_t *)&uc;
    }
    const uint16_t getUint16(void)
    {
        uint8_t a = getUint8();
        uint8_t b = getUint8();
        return (readEndian_ == BIG) ? ((a << 8) | b) : ((b << 8) | a);
    }
    const int16_t getInt16(void)
    {
        uint8_t a = getUint8();
        uint8_t b = getUint8();
        return (readEndian_ == BIG) ? ((*(int8_t *)&a << 8) | b) : ((*(int8_t *)&b << 8) | a);
    }
    const uint32_t getUint32(void)
    {
        uint16_t a = getUint16();
        uint16_t b = getUint16();
        return (readEndian_ == BIG) ? ((a << 16) | b) : ((b << 16) | a);
    }
    const int32_t getInt32(void)
    {
        uint16_t a = getUint16();
        uint16_t b = getUint16();
        return (readEndian_ == BIG) ? ((*(int16_t *)&a << 16) | b) : ((*(int16_t *)&b << 16) | a);
    }
    const uint64_t getUint64(void)
    {
        uint64_t lo, hi;
        if (readEndian_ == BIG)
        {
            hi = getUint32();
            lo = getUint32();
        }
        else
        {
            lo = getUint32();
            hi = getUint32();
        }
        return (lo | hi << 32);
    }
    const int64_t getInt64(void)
    {
        uint64_t lo;
        int64_t  hi;
        if (readEndian_ == BIG)
        {
            hi = getInt32();
            lo = getUint32();
        }
        else
        {
            lo = getUint32();
            hi = getInt32();
        }
        return (lo | hi << 32);
    }
    const float getFloat(void)
    {
        union
        {
            uint32_t ui;
            float     f;
        } x;
        x.ui = getUint32();
        return x.f;
    }
    const double getDouble(void)
    {
        union
        {
            uint64_t ull;
            double   f;
        } x;
        x.ull = getUint64();
        return x.f;
    }
    const std::string getSubstr(int len)
    {
        std::string str;
        for (int i = 0; i < len - 1; i++)
        {
            uint8_t uc = getUint8();
            str.push_back(uc);
        }
        if (len > 0)
            getUint8();
        return str;
    }
    const std::string get_string(void)
    {
        int len = getUint32();
        if (len > 0)
            return getSubstr(len);
        return "";
    }

    void align(int a)
    {
        size_t s = ftell();
        int m = (a - (s % a)) % a;
        fseek(m, CURRENT);
    }
    Pointer getPointer32(void)
    {
        Pointer p = getUint32();
        return p;
    }
    Pointer putPointer32(void)
    {
        Pointer p = mark();
        putUint32(0);
        return p;
    }
    Pointer getPointer64(void)
    {
        Pointer p = getUint64();
        return p;
    }
    Pointer putPointer64(void)
    {
        Pointer p = mark();
        putUint64(0);
        return p;
    }

    void adjustPointer32(Pointer p)
    {
        uint32_t ui  = (uint32_t)buffer_.size();
        uint8_t  uc0 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc1 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc2 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc3 = ui & 0xFF;
        if (writeEndian_ == BIG)
        {
            buffer_[p]     = uc3;
            buffer_[p + 1] = uc2;
            buffer_[p + 2] = uc1;
            buffer_[p + 3] = uc0;
        }
        else
        {
            buffer_[p]     = uc0;
            buffer_[p + 1] = uc1;
            buffer_[p + 2] = uc2;
            buffer_[p + 3] = uc3;
        }
    }
    void adjustPointer64(Pointer p)
    {
        uint64_t ui  = buffer_.size();
        uint8_t  uc0 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc1 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc2 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc3 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc4 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc5 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc6 = ui & 0xFF;
        ui >>= 8;
        uint8_t uc7 = ui & 0xFF;
        if (writeEndian_ == BIG)
        {
            buffer_[p]     = uc7;
            buffer_[p + 1] = uc6;
            buffer_[p + 2] = uc5;
            buffer_[p + 3] = uc4;
            buffer_[p + 4] = uc3;
            buffer_[p + 5] = uc2;
            buffer_[p + 6] = uc1;
            buffer_[p + 7] = uc0;
        }
        else
        {
            buffer_[p]     = uc0;
            buffer_[p + 1] = uc1;
            buffer_[p + 2] = uc2;
            buffer_[p + 3] = uc3;
            buffer_[p + 4] = uc4;
            buffer_[p + 5] = uc5;
            buffer_[p + 6] = uc6;
            buffer_[p + 7] = uc7;
        }
    }

    void putInt8(const int8_t v) { putUint8(*(uint8_t *)&v); }
    void putUint16(const uint16_t u)
    {
        if (writeEndian_ == BIG)
        {
            putUint8((u >> 8) & 0xFF);
            putUint8(u & 0xFF);
        }
        else
        {
            // LE
            putUint8(u & 0xFF);
            putUint8((u >> 8) & 0xFF);
        }
    }
    void putInt16(const int16_t u)
    {
        if (writeEndian_ == BIG)
        {
            putInt8((u >> 8) & 0xFF);
            putUint8(u & 0xFF);
        }
        else
        {
            // LE
            putUint8(u & 0xFF);
            putInt8((u >> 8) & 0xFF);
        }
    }
    void putUint32(const uint32_t u)
    {
        if (writeEndian_ == BIG)
        {
            putUint16((u >> 16) & 0xFFFF);
            putUint16(u & 0xFFFF);
        }
        else
        {
            putUint16(u & 0xFFFF);
            putUint16((u >> 16) & 0xFFFF);
        }
    }
    void putInt32(const int32_t u)
    {
        if (writeEndian_ == BIG)
        {
            putInt16((u >> 16) & 0xFFFF);
            putUint16(u & 0xFFFF);
        }
        else
        {
            putUint16(u & 0xFFFF);
            putInt16((u >> 16) & 0xFFFF);
        }
    }

    void putUint64(const uint64_t v)
    {
        if (writeEndian_ == BIG)
        {
            putUint32(v >> 32);
            putUint32(v & 0xFFFFFFFF);
        }
        else
        {
            putUint32(v & 0xFFFFFFFF);
            putUint32(v >> 32);
        }
    }
    void putInt64(const int64_t v)
    {
        if (writeEndian_ == BIG)
        {
            putInt32(v >> 32);
            putUint32(v & 0xFFFFFFFF);
        }
        else
        {
            putUint32(v & 0xFFFFFFFF);
            putInt32(v >> 32);
        }
    }
    void putFloat(const float v) { putUint32(*(uint32_t *)&v); } // endian henkan is done by uint.
    void putDouble(const double v) { putUint64(*(uint64_t *)&v); }
    void putSubstr(const std::string &str, size_t len)
    {
        for (int i = 0; i < len - 1; i++)
        {
            putUint8(str[i]);
        }
        if (len > 0)
            putUint8(0);
    }
    void putString(const std::string &str) { putSubstr(str, str.size()); }
};

static Gamut loadGamutFromICCProfileMemory(const void *mem, size_t size) 
{
    MemoryStream stream(mem, size);
    return Gamut("",Matrix3(1,0,0,0,1,0,0,0,1));
}


class Spectrum
{
  public:
    typedef std::array<float, 400> spectrum; // 380-780, 1nm, fixed.
    spectrum                       s_;

    constexpr Spectrum() : s_{} { ; }
    constexpr Spectrum(const spectrum &s) : s_(s) { ; }
    static constexpr float lerp(const float a, const float b, const float r) { return a * (1 - r) + b * r; }
    static constexpr float fetch(
        const float *sample, const int samples, const float &lo, const float &hi, const float &t)
    {
        const int   count = samples - 1;
        const int   index = (int)(count * (t - lo) / (hi - lo));
        const float l     = index * (hi - lo) / count;
        const float h     = (index + 1) * (hi - lo) / count;
        return (index >= count)
                   ? sample[count]
                   : (index <= 0) ? sample[0] : lerp(sample[index], sample[index + 1], (t - l - lo) / (h - l));
    }
    constexpr float fetch(const float &lambda) const { return fetch(s_.data(), 400, 380.f, 779.f, lambda); }

    Spectrum(const float *sample, const int samples = 400, const float &lo = 380.f, const float &hi = 779.f)
    {
        for (int i = 0; i < 400; i++)
        {
            s_[i] = fetch(sample, samples, lo, hi, 380.f + i);
        }
    }
    constexpr float operator[](const int i) const { return s_[i]; }

    const spectrum &s(void) const { return s_; }

    static const double planck(const double &T, // temperature (Kelvin)
        const double &                       l)                        // wavelength (meter)
    {
        static const double h    = 6.62606896e-34; // Plank constant         J s
        static const double c    = 2.99792458e+8;  // Speed of light         m/s
        static const double k    = 1.38064880e-23; // Boltzmann constant     J/K
        static const double hc   = 1.984832809e-25;
        static const double hcc  = 5.950379064e-17;
        static const double hc_k = 1.438776827e-2;
        static const double arg1 = 2 * hcc;                                     // J*s*m/s*m/s                 = J*m^2/s
        static const double arg2 = (h * c) / k;                                 // J*s*m/s / (J/K) = J*m/(J/K) = m * K
        return (float)(arg1 * pow(l, -5) / (exp(hc_k / (l * T)) - 1.0)) / 1e9f; // in W/m^3 sr
    }

    static const Spectrum blackbody(const float temp)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = (float)planck(temp, (double)(380.f + i) * 1e-9);
        }
        return Spectrum(s);
    }
    static const Spectrum E(const float e = 1.f)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = e;
        }
        return Spectrum(s);
    }
    static const Spectrum mul(const Spectrum &a, const Spectrum &b)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = a[i] * b[i];
        }
        return Spectrum(s);
    }
    static const Spectrum add(const Spectrum &a, const Spectrum &b)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = a[i] + b[i];
        }
        return Spectrum(s);
    }
    const Spectrum         operator*(const Spectrum &b) const { return mul(*this, b); }
    const Spectrum         operator+(const Spectrum &b) const { return add(*this, b); }
    static constexpr float sumHelper(const Spectrum &a, const int i)
    {
        return (i > 0) ? sumHelper(a, i - 1) + a[i] : a[0];
    }
    static constexpr float sum(const Spectrum &a) { return sumHelper(a, 399); }
    constexpr float        sum() const { return sum(*this); }
    static constexpr float dotHelper(const Spectrum &a, const Spectrum &b, const int i)
    {
        return (i > 0) ? a[i] * b[i] + dotHelper(a, b, i - 1) : a[0] * b[0];
    }
    static constexpr float dot(const Spectrum &a, const Spectrum &b) { return dotHelper(a, b, 399); }
    constexpr float        dot(const Spectrum &s) const { return dotHelper(*this, s, 399); }
    static constexpr float dotHelper3(const Spectrum &a, const Spectrum &b, const Spectrum &c, const int i)
    {
        return (i > 0) ? a[i] * b[i] * c[i] + dotHelper3(a, b, c, i - 1) : a[0] * b[0] * c[0];
    }
    static constexpr float dot3(const Spectrum &a, const Spectrum &b, const Spectrum &c)
    {
        return dotHelper3(a, b, c, 399);
    }
};

// standard illuminants
static constexpr Tristimulus Illuminant_A(1.09850f, 1.f, 0.35585f);
static constexpr Tristimulus Illuminant_B(1.99072f, 1.f, 0.85223f);
static constexpr Tristimulus Illuminant_C(0.98074f, 1.0f, 1.18232f);
static constexpr Tristimulus Illuminant_D50(0.96422f, 1.0f, 0.82521f);
static constexpr Tristimulus Illuminant_D55(0.95682f, 1.0f, 0.92149f);
static constexpr Tristimulus Illuminant_D60(0.95265f, 1.0f, 1.00883f); // by ACES TB-2014-004.pdf
static constexpr Tristimulus Illuminant_D65(0.95047f, 1.0f, 1.08883f);
static constexpr Tristimulus Illuminant_D75(0.94972f, 1.0f, 1.22638f);
static constexpr Tristimulus Illuminant_E(1.f, 1.f, 1.f);
static constexpr Tristimulus Illuminant_F2(0.99186f, 1.f, 0.67393f);
static constexpr Tristimulus Illuminant_F7(0.95041f, 1.f, 1.08747f);
static constexpr Tristimulus Illuminant_F11(1.00962f, 1.f, 0.64350f);

// xR,yR,xG,yG,xB,yB,xW,yW
static constexpr Gamut AdobeRGB("AdobeRGB", 0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
static constexpr Gamut Rec709("Rec.709", 0.64f, 0.33f, 0.30f, 0.60f, 0.15f, 0.06f, 0.3127f, 0.3290f);
static constexpr Gamut Rec2020("Rec.2020", 0.708f, 0.292f, 0.17f, 0.797f, 0.131f, 0.046f, 0.3127f, 0.3290f);
static constexpr Gamut DCI_P3("DCI P3", 0.68f, 0.32f, 0.265f, 0.69f, 0.15f, 0.06f, 0.314f, 0.351f);
static constexpr Gamut S_Gamut("S-Gamut", 0.73f, 0.28f, 0.14f, 0.855f, 0.10f, -0.05f, 0.3127f, 0.3290f);
static constexpr Gamut S_Gamut3_Cine(
    "S-Gamut3.Cine", 0.766f, 0.275f, 0.225f, 0.800f, 0.089f, -0.087f, 0.3127f, 0.3290f);
static constexpr Gamut ACEScg("ACES cg", 0.713f, 0.293f, 0.165f, 0.830f, 0.128f, 0.044f, 0.32168f, 0.33767f);     // AP1
static constexpr Gamut ACES2065("ACES 2065", 0.73470f, 0.26530f, 0.f, 1.f, 0.0001f, -0.077f, 0.32168f, 0.33767f); // AP0
static constexpr Gamut LMS("LMS",
    Matrix3(0.8951f, 0.2664f, -0.1614f, -0.7502f, 1.7135f, 0.0367f, 0.0389f, -0.0685f, 1.0296f)); // fromXYZ matrix.
static constexpr Gamut XYZ("XYZ", Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1));

// returns Gamut convert matrix
static constexpr Matrix3 GamutConvert(const Gamut &src, const Gamut &dst) { return dst.fromXYZ().mul(src.toXYZ()); }

// returns Bradford adaptation matrix
static constexpr Matrix3 Bradford(const Tristimulus &white_src, const Tristimulus &white_dst)
{
    const Tristimulus &lms_src(LMS.fromXYZ(white_src));
    const Tristimulus &lms_dst(LMS.fromXYZ(white_dst));
    const Matrix3      scale(
        Matrix3::diag(Vector3(lms_dst[0] / lms_src[0], lms_dst[1] / lms_src[1], lms_dst[2] / lms_src[2])));
    return LMS.toXYZ().mul(scale).mul(LMS.fromXYZ());
}

static constexpr Tristimulus XYZ_to_ICtCp(const Tristimulus& xyz)
{
    const Tristimulus lms = LMS.fromXYZ(xyz);
    const Tristimulus pq = Tristimulus(
        OTF::Y_to_ST2084(lms[0]/100.f),
        OTF::Y_to_ST2084(lms[1]/100.f),
        OTF::Y_to_ST2084(lms[2]/100.f)
    ); // map 10000 nits to 0-100
    return Tristimulus(
        (pq[0]+pq[1])*0.5f,
        (6610.f * pq[0] - 13613.f * pq[1] + 7003.f * pq[2])/4096.f,
        (17933.f * pq[0] - 17390.f * pq[1] + 543.f * pq[2])/4096.f
    );
}

// Color Differences
class Delta
{
  public:
    static const float UV(const Tristimulus &a_Yuv, const Tristimulus &b_Yuv) // a, b both are XYZ
    {
        return sqrtf((a_Yuv[1] - b_Yuv[1]) * (a_Yuv[1] - b_Yuv[1]) + (a_Yuv[2] - b_Yuv[2]) * (a_Yuv[2] - b_Yuv[2]));
    }
    static const float E76(const Tristimulus &a_LAB, const Tristimulus &b_LAB)
    {
        return sqrtf((a_LAB[0] - b_LAB[0]) * (a_LAB[0] - b_LAB[0]) + (a_LAB[1] - b_LAB[1]) * (a_LAB[1] - b_LAB[1]) +
                     (a_LAB[2] - b_LAB[2]) * (a_LAB[2] - b_LAB[2]));
    }
    static const float E00(const Tristimulus &lab1, const Tristimulus &lab2, const float &Kl = 1.f,
        const float &Kc = 1.f, const float &Kh = 1.f)
    {
        const float PI      = 3.14159265358979323846264338327950288f;
        const float L1      = lab1[0];
        const float a1      = lab1[1];
        const float b1      = lab1[2];
        const float L2      = lab2[0];
        const float a2      = lab2[1];
        const float b2      = lab2[2];
        const float Lbar    = (L1 + L2) / 2.f;
        const float C1      = sqrtf(a1 * a1 + b1 * b1);
        const float C2      = sqrtf(a2 * a2 + b2 * b2);
        const float Cbar    = (C1 + C2) / 2.f;
        const float C7      = powf(Cbar, 7.f);
        const float pow25_7 = 25.f * 25.f * 25.f * 25.f * 25.f * 25.f * 25.f;
        const float G       = (1.f - sqrtf(C7 / (C7 + pow25_7))) / 2.f;
        const float ad1     = a1 * (1.f + G);
        const float ad2     = a2 * (1.f + G);
        const float Cd1     = sqrtf(ad1 * ad1 + b1 * b1);
        const float Cd2     = sqrtf(ad2 * ad2 + b2 * b2);
        const float CdBar   = (Cd1 + Cd2) / 2.f;
        const float h1      = fmodf(360.f + atan2f(b1, ad1) * 180.0f / PI, 360.f);
        const float h2      = fmodf(360.f + atan2f(b2, ad2) * 180.0f / PI, 360.f);
        const float HdBar   = (fabs(h1 - h2) > 180.f ? (h1 + h2 + 360.f) : (h1 + h2)) / 2.f;
        const float T1      = 1.f - 0.17f * cosf(PI * (1.f * HdBar - 30.f) / 180.f);
        const float T2      = 0.24f * cosf(PI * (2.f * HdBar) / 180.f) + 0.32f * cosf(PI * (3.f * HdBar + 6.f) / 180.f);
        const float T3      = 0.20f * cosf(PI * (4.f * HdBar - 63.f) / 180.f);
        const float T       = T1 + T2 - T3;
        const float deltah  = (fabs(h2 - h1) <= 180.f) ? h2 - h1 : ((h2 <= h1) ? h2 - h1 + 360.f : h2 - h1 - 360.f);
        const float deltaL  = L2 - L1;
        const float deltaC  = Cd2 - Cd1;
        const float deltaH  = 2.f * sqrtf(Cd1 * Cd2) * sinf(PI * deltah / (180.f * 2.f));
        const float Lbar2   = (Lbar - 50.f) * (Lbar - 50.f);
        const float Sl      = 1.f + 0.015f * Lbar2 / sqrtf(20.f + Lbar2);
        const float Sc      = 1.f + 0.045f * CdBar;
        const float Sh      = 1.f + 0.015f * CdBar * T;
        const float HdBar2  = (HdBar - 275.f) * (HdBar - 275.f) / (25.f * 25.f);
        const float deltaTheta = 30.f * expf(-HdBar2);
        const float CdBar7     = powf(CdBar, 7.f);
        const float Rc         = 2.f * sqrtf(CdBar7 / (CdBar7 + pow25_7));
        const float Rt         = -Rc * sinf(2.f * deltaTheta * PI / 180.f);
        const float dl         = deltaL / (Kl * Sl);
        const float dc         = deltaC / (Kc * Sc);
        const float dh         = deltaH / (Kh * Sh);

        return sqrtf(dl * dl + dc * dc + dh * dh + Rt * dc * dh);
    }

    //https://calman.spectracal.com/delta-ictcp-color-difference-metric.html
    static const float ICtCp(const Tristimulus& a_xyz, const Tristimulus& b_xyz)
    {
        const Tristimulus a_itp = XYZ_to_ICtCp(a_xyz);
        const Tristimulus b_itp = XYZ_to_ICtCp(b_xyz);
        const float dI = a_itp[0] - b_itp[0];
        const float dT = a_itp[1] - b_itp[1];
        const float dP = a_itp[2] - b_itp[2];
        return sqrtf(dI*dI + dT*dT*0.25f + dP*dP);
    }
};


class Observer
{
  public:
    Spectrum    X_, Y_, Z_;
    Tristimulus normalize_;
    constexpr Observer(const Spectrum &X, const Spectrum &Y, const Spectrum &Z)
        : X_(X), Y_(Y), Z_(Z), normalize_(1.f / X.sum(), 1.f / Y.sum(), 1.f / Z.sum())
    {
        ;
    }
    static constexpr Tristimulus SpectrumIntegrate(
        const Spectrum &s, const Spectrum &x, const Spectrum &y, const Spectrum &z)
    {
        return Tristimulus(Spectrum::dot(s, x), Spectrum::dot(s, y), Spectrum::dot(s, z));
    }
    static constexpr Tristimulus SpectrumIntegrate3(
        const Spectrum &r, const Spectrum &l, const Spectrum &x, const Spectrum &y, const Spectrum &z)
    {
        return Tristimulus(Spectrum::dot3(r, l, x), Spectrum::dot3(r, l, y), Spectrum::dot3(r, l, z));
    }
    constexpr Tristimulus fromSpectrum(const Spectrum &s) const
    {
        return SpectrumIntegrate(s, X_, Y_, Z_) * normalize_;
    }
    constexpr Tristimulus fromReflectanceAndLight(const Spectrum &r, const Spectrum &l) const // r:reflectance, l:light
    {
        return SpectrumIntegrate3(r, l, X_, Y_, Z_) * normalize_;
    }
    constexpr float lumensFromMonochromaticFlux(const float &lambda, const float &watt) const // return lm
    {
        return Y_.fetch(lambda) * watt * 683.0f; // photoptic luminance efficiency
    }
    constexpr Tristimulus candellasFromMonochromaticRadiance(
        const float &lambda, const float &watt_per_sr_m2) const // return cd/m^2
    {
        return Tristimulus(X_.fetch(lambda), Y_.fetch(lambda), Z_.fetch(lambda)) * 683.0f * watt_per_sr_m2;
    }
    constexpr Tristimulus xyz(const float lambda) const
    {
        return Tristimulus(X_.fetch(lambda), Y_.fetch(lambda), Z_.fetch(lambda)).mul(normalize_);
    }
};

static constexpr Spectrum CIE1931_X({0.00136800000f, 0.00150205000f, 0.00164232800f, 0.00180238200f, 0.00199575700f,
    0.00223600000f, 0.00253538500f, 0.00289260300f, 0.00330082900f, 0.00375323600f, 0.00424300000f, 0.00476238900f,
    0.00533004800f, 0.00597871200f, 0.00674111700f, 0.00765000000f, 0.00875137300f, 0.01002888000f, 0.01142170000f,
    0.01286901000f, 0.01431000000f, 0.01570443000f, 0.01714744000f, 0.01878122000f, 0.02074801000f, 0.02319000000f,
    0.02620736000f, 0.02978248000f, 0.03388092000f, 0.03846824000f, 0.04351000000f, 0.04899560000f, 0.05502260000f,
    0.06171880000f, 0.06921200000f, 0.07763000000f, 0.08695811000f, 0.09717672000f, 0.10840630000f, 0.12076720000f,
    0.13438000000f, 0.14935820000f, 0.16539570000f, 0.18198310000f, 0.19861100000f, 0.21477000000f, 0.23018680000f,
    0.24487970000f, 0.25877730000f, 0.27180790000f, 0.28390000000f, 0.29494380000f, 0.30489650000f, 0.31378730000f,
    0.32164540000f, 0.32850000000f, 0.33435130000f, 0.33921010000f, 0.34312130000f, 0.34612960000f, 0.34828000000f,
    0.34959990000f, 0.35014740000f, 0.35001300000f, 0.34928700000f, 0.34806000000f, 0.34637330000f, 0.34426240000f,
    0.34180880000f, 0.33909410000f, 0.33620000000f, 0.33319770000f, 0.33004110000f, 0.32663570000f, 0.32288680000f,
    0.31870000000f, 0.31402510000f, 0.30888400000f, 0.30329040000f, 0.29725790000f, 0.29080000000f, 0.28397010000f,
    0.27672140000f, 0.26891780000f, 0.26042270000f, 0.25110000000f, 0.24084750000f, 0.22985120000f, 0.21840720000f,
    0.20681150000f, 0.19536000000f, 0.18421360000f, 0.17332730000f, 0.16268810000f, 0.15228330000f, 0.14210000000f,
    0.13217860000f, 0.12256960000f, 0.11327520000f, 0.10429790000f, 0.09564000000f, 0.08729955000f, 0.07930804000f,
    0.07171776000f, 0.06458099000f, 0.05795001000f, 0.05186211000f, 0.04628152000f, 0.04115088000f, 0.03641283000f,
    0.03201000000f, 0.02791720000f, 0.02414440000f, 0.02068700000f, 0.01754040000f, 0.01470000000f, 0.01216179000f,
    0.00991996000f, 0.00796724000f, 0.00629634600f, 0.00490000000f, 0.00377717300f, 0.00294532000f, 0.00242488000f,
    0.00223629300f, 0.00240000000f, 0.00292552000f, 0.00383656000f, 0.00517484000f, 0.00698208000f, 0.00930000000f,
    0.01214949000f, 0.01553588000f, 0.01947752000f, 0.02399277000f, 0.02910000000f, 0.03481485000f, 0.04112016000f,
    0.04798504000f, 0.05537861000f, 0.06327000000f, 0.07163501000f, 0.08046224000f, 0.08973996000f, 0.09945645000f,
    0.10960000000f, 0.12016740000f, 0.13111450000f, 0.14236790000f, 0.15385420000f, 0.16550000000f, 0.17725710000f,
    0.18914000000f, 0.20116940000f, 0.21336580000f, 0.22574990000f, 0.23832090000f, 0.25106680000f, 0.26399220000f,
    0.27710170000f, 0.29040000000f, 0.30389120000f, 0.31757260000f, 0.33143840000f, 0.34548280000f, 0.35970000000f,
    0.37408390000f, 0.38863960000f, 0.40337840000f, 0.41831150000f, 0.43344990000f, 0.44879530000f, 0.46433600000f,
    0.48006400000f, 0.49597130000f, 0.51205010000f, 0.52829590000f, 0.54469160000f, 0.56120940000f, 0.57782150000f,
    0.59450000000f, 0.61122090000f, 0.62797580000f, 0.64476020000f, 0.66156970000f, 0.67840000000f, 0.69523920000f,
    0.71205860000f, 0.72882840000f, 0.74551880000f, 0.76210000000f, 0.77854320000f, 0.79482560000f, 0.81092640000f,
    0.82682480000f, 0.84250000000f, 0.85793250000f, 0.87308160000f, 0.88789440000f, 0.90231810000f, 0.91630000000f,
    0.92979950000f, 0.94279840000f, 0.95527760000f, 0.96721790000f, 0.97860000000f, 0.98938560000f, 0.99954880000f,
    1.00908920000f, 1.01800640000f, 1.02630000000f, 1.03398270000f, 1.04098600000f, 1.04718800000f, 1.05246670000f,
    1.05670000000f, 1.05979440000f, 1.06179920000f, 1.06280680000f, 1.06290960000f, 1.06220000000f, 1.06073520000f,
    1.05844360000f, 1.05522440000f, 1.05097680000f, 1.04560000000f, 1.03903690000f, 1.03136080000f, 1.02266620000f,
    1.01304770000f, 1.00260000000f, 0.99136750000f, 0.97933140000f, 0.96649160000f, 0.95284790000f, 0.93840000000f,
    0.92319400000f, 0.90724400000f, 0.89050200000f, 0.87292000000f, 0.85444990000f, 0.83508400000f, 0.81494600000f,
    0.79418600000f, 0.77295400000f, 0.75140000000f, 0.72958360000f, 0.70758880000f, 0.68560220000f, 0.66381040000f,
    0.64240000000f, 0.62151490000f, 0.60111380000f, 0.58110520000f, 0.56139770000f, 0.54190000000f, 0.52259950000f,
    0.50354640000f, 0.48474360000f, 0.46619390000f, 0.44790000000f, 0.42986130000f, 0.41209800000f, 0.39464400000f,
    0.37753330000f, 0.36080000000f, 0.34445630000f, 0.32851680000f, 0.31301920000f, 0.29800110000f, 0.28350000000f,
    0.26954480000f, 0.25611840000f, 0.24318960000f, 0.23072720000f, 0.21870000000f, 0.20709710000f, 0.19592320000f,
    0.18517080000f, 0.17483230000f, 0.16490000000f, 0.15536670000f, 0.14623000000f, 0.13749000000f, 0.12914670000f,
    0.12120000000f, 0.11363970000f, 0.10646500000f, 0.09969044000f, 0.09333061000f, 0.08740000000f, 0.08190096000f,
    0.07680428000f, 0.07207712000f, 0.06768664000f, 0.06360000000f, 0.05980685000f, 0.05628216000f, 0.05297104000f,
    0.04981861000f, 0.04677000000f, 0.04378405000f, 0.04087536000f, 0.03807264000f, 0.03540461000f, 0.03290000000f,
    0.03056419000f, 0.02838056000f, 0.02634484000f, 0.02445275000f, 0.02270000000f, 0.02108429000f, 0.01959988000f,
    0.01823732000f, 0.01698717000f, 0.01584000000f, 0.01479064000f, 0.01383132000f, 0.01294868000f, 0.01212920000f,
    0.01135916000f, 0.01062935000f, 0.00993884600f, 0.00928842200f, 0.00867885400f, 0.00811091600f, 0.00758238800f,
    0.00708874600f, 0.00662731300f, 0.00619540800f, 0.00579034600f, 0.00540982600f, 0.00505258300f, 0.00471751200f,
    0.00440350700f, 0.00410945700f, 0.00383391300f, 0.00357574800f, 0.00333434200f, 0.00310907500f, 0.00289932700f,
    0.00270434800f, 0.00252302000f, 0.00235416800f, 0.00219661600f, 0.00204919000f, 0.00191096000f, 0.00178143800f,
    0.00166011000f, 0.00154645900f, 0.00143997100f, 0.00134004200f, 0.00124627500f, 0.00115847100f, 0.00107643000f,
    0.00099994930f, 0.00092873580f, 0.00086243320f, 0.00080075030f, 0.00074339600f, 0.00069007860f, 0.00064051560f,
    0.00059450210f, 0.00055186460f, 0.00051242900f, 0.00047602130f, 0.00044245360f, 0.00041151170f, 0.00038298140f,
    0.00035664910f, 0.00033230110f, 0.00030975860f, 0.00028888710f, 0.00026953940f, 0.00025156820f, 0.00023482610f,
    0.00021917100f, 0.00020452580f, 0.00019084050f, 0.00017806540f, 0.00016615050f, 0.00015502360f, 0.00014462190f,
    0.00013490980f, 0.00012585200f, 0.00011741300f, 0.00010955150f, 0.00010222450f, 0.00009539445f, 0.00008902390f,
    0.00008307527f, 0.00007751269f, 0.00007231304f, 0.00006745778f, 0.00006292844f, 0.00005870652f, 0.00005477028f,
    0.00005109918f, 0.00004767654f, 0.00004448567f});
static constexpr Spectrum CIE1931_Y({0.00003900000f, 0.00004282640f, 0.00004691460f, 0.00005158960f, 0.00005717640f,
    0.00006400000f, 0.00007234421f, 0.00008221224f, 0.00009350816f, 0.00010613610f, 0.00012000000f, 0.00013498400f,
    0.00015149200f, 0.00017020800f, 0.00019181600f, 0.00021700000f, 0.00024690670f, 0.00028124000f, 0.00031852000f,
    0.00035726670f, 0.00039600000f, 0.00043371470f, 0.00047302400f, 0.00051787600f, 0.00057221870f, 0.00064000000f,
    0.00072456000f, 0.00082550000f, 0.00094116000f, 0.00106988000f, 0.00121000000f, 0.00136209100f, 0.00153075200f,
    0.00172036800f, 0.00193532300f, 0.00218000000f, 0.00245480000f, 0.00276400000f, 0.00311780000f, 0.00352640000f,
    0.00400000000f, 0.00454624000f, 0.00515932000f, 0.00582928000f, 0.00654616000f, 0.00730000000f, 0.00808650700f,
    0.00890872000f, 0.00976768000f, 0.01066443000f, 0.01160000000f, 0.01257317000f, 0.01358272000f, 0.01462968000f,
    0.01571509000f, 0.01684000000f, 0.01800736000f, 0.01921448000f, 0.02045392000f, 0.02171824000f, 0.02300000000f,
    0.02429461000f, 0.02561024000f, 0.02695857000f, 0.02835125000f, 0.02980000000f, 0.03131083000f, 0.03288368000f,
    0.03452112000f, 0.03622571000f, 0.03800000000f, 0.03984667000f, 0.04176800000f, 0.04376600000f, 0.04584267000f,
    0.04800000000f, 0.05024368000f, 0.05257304000f, 0.05498056000f, 0.05745872000f, 0.06000000000f, 0.06260197000f,
    0.06527752000f, 0.06804208000f, 0.07091109000f, 0.07390000000f, 0.07701600000f, 0.08026640000f, 0.08366680000f,
    0.08723280000f, 0.09098000000f, 0.09491755000f, 0.09904584000f, 0.10336740000f, 0.10788460000f, 0.11260000000f,
    0.11753200000f, 0.12267440000f, 0.12799280000f, 0.13345280000f, 0.13902000000f, 0.14467640000f, 0.15046930000f,
    0.15646190000f, 0.16271770000f, 0.16930000000f, 0.17624310000f, 0.18355810000f, 0.19127350000f, 0.19941800000f,
    0.20802000000f, 0.21711990000f, 0.22673450000f, 0.23685710000f, 0.24748120000f, 0.25860000000f, 0.27018490000f,
    0.28229390000f, 0.29505050000f, 0.30857800000f, 0.32300000000f, 0.33840210000f, 0.35468580000f, 0.37169860000f,
    0.38928750000f, 0.40730000000f, 0.42562990000f, 0.44430960000f, 0.46339440000f, 0.48293950000f, 0.50300000000f,
    0.52356930000f, 0.54451200000f, 0.56569000000f, 0.58696530000f, 0.60820000000f, 0.62934560000f, 0.65030680000f,
    0.67087520000f, 0.69084240000f, 0.71000000000f, 0.72818520000f, 0.74546360000f, 0.76196940000f, 0.77783680000f,
    0.79320000000f, 0.80811040000f, 0.82249620000f, 0.83630680000f, 0.84949160000f, 0.86200000000f, 0.87381080000f,
    0.88496240000f, 0.89549360000f, 0.90544320000f, 0.91485010000f, 0.92373480000f, 0.93209240000f, 0.93992260000f,
    0.94722520000f, 0.95400000000f, 0.96025610000f, 0.96600740000f, 0.97126060000f, 0.97602250000f, 0.98030000000f,
    0.98409240000f, 0.98741820000f, 0.99031280000f, 0.99281160000f, 0.99495010000f, 0.99671080000f, 0.99809830000f,
    0.99911200000f, 0.99974820000f, 1.00000000000f, 0.99985670000f, 0.99930460000f, 0.99832550000f, 0.99689870000f,
    0.99500000000f, 0.99260050000f, 0.98974260000f, 0.98644440000f, 0.98272410000f, 0.97860000000f, 0.97408370000f,
    0.96917120000f, 0.96385680000f, 0.95813490000f, 0.95200000000f, 0.94545040000f, 0.93849920000f, 0.93116280000f,
    0.92345760000f, 0.91540000000f, 0.90700640000f, 0.89827720000f, 0.88920480000f, 0.87978160000f, 0.87000000000f,
    0.85986130000f, 0.84939200000f, 0.83862200000f, 0.82758130000f, 0.81630000000f, 0.80479470000f, 0.79308200000f,
    0.78119200000f, 0.76915470000f, 0.75700000000f, 0.74475410000f, 0.73242240000f, 0.72000360000f, 0.70749650000f,
    0.69490000000f, 0.68221920000f, 0.66947160000f, 0.65667440000f, 0.64384480000f, 0.63100000000f, 0.61815550000f,
    0.60531440000f, 0.59247560000f, 0.57963790000f, 0.56680000000f, 0.55396110000f, 0.54113720000f, 0.52835280000f,
    0.51563230000f, 0.50300000000f, 0.49046880000f, 0.47803040000f, 0.46567760000f, 0.45340320000f, 0.44120000000f,
    0.42908000000f, 0.41703600000f, 0.40503200000f, 0.39303200000f, 0.38100000000f, 0.36891840000f, 0.35682720000f,
    0.34477680000f, 0.33281760000f, 0.32100000000f, 0.30933810000f, 0.29785040000f, 0.28659360000f, 0.27562450000f,
    0.26500000000f, 0.25476320000f, 0.24488960000f, 0.23533440000f, 0.22605280000f, 0.21700000000f, 0.20816160000f,
    0.19954880000f, 0.19115520000f, 0.18297440000f, 0.17500000000f, 0.16722350000f, 0.15964640000f, 0.15227760000f,
    0.14512590000f, 0.13820000000f, 0.13150030000f, 0.12502480000f, 0.11877920000f, 0.11276910000f, 0.10700000000f,
    0.10147620000f, 0.09618864000f, 0.09112296000f, 0.08626485000f, 0.08160000000f, 0.07712064000f, 0.07282552000f,
    0.06871008000f, 0.06476976000f, 0.06100000000f, 0.05739621000f, 0.05395504000f, 0.05067376000f, 0.04754965000f,
    0.04458000000f, 0.04175872000f, 0.03908496000f, 0.03656384000f, 0.03420048000f, 0.03200000000f, 0.02996261000f,
    0.02807664000f, 0.02632936000f, 0.02470805000f, 0.02320000000f, 0.02180077000f, 0.02050112000f, 0.01928108000f,
    0.01812069000f, 0.01700000000f, 0.01590379000f, 0.01483718000f, 0.01381068000f, 0.01283478000f, 0.01192000000f,
    0.01106831000f, 0.01027339000f, 0.00953331100f, 0.00884615700f, 0.00821000000f, 0.00762378100f, 0.00708542400f,
    0.00659147600f, 0.00613848500f, 0.00572300000f, 0.00534305900f, 0.00499579600f, 0.00467640400f, 0.00438007500f,
    0.00410200000f, 0.00383845300f, 0.00358909900f, 0.00335421900f, 0.00313409300f, 0.00292900000f, 0.00273813900f,
    0.00255987600f, 0.00239324400f, 0.00223727500f, 0.00209100000f, 0.00195358700f, 0.00182458000f, 0.00170358000f,
    0.00159018700f, 0.00148400000f, 0.00138449600f, 0.00129126800f, 0.00120409200f, 0.00112274400f, 0.00104700000f,
    0.00097658960f, 0.00091110880f, 0.00085013320f, 0.00079323840f, 0.00074000000f, 0.00069008270f, 0.00064331000f,
    0.00059949600f, 0.00055845470f, 0.00052000000f, 0.00048391360f, 0.00045005280f, 0.00041834520f, 0.00038871840f,
    0.00036110000f, 0.00033538350f, 0.00031144040f, 0.00028916560f, 0.00026845390f, 0.00024920000f, 0.00023130190f,
    0.00021468560f, 0.00019928840f, 0.00018504750f, 0.00017190000f, 0.00015977810f, 0.00014860440f, 0.00013830160f,
    0.00012879250f, 0.00012000000f, 0.00011185950f, 0.00010432240f, 0.00009733560f, 0.00009084587f, 0.00008480000f,
    0.00007914667f, 0.00007385800f, 0.00006891600f, 0.00006430267f, 0.00006000000f, 0.00005598187f, 0.00005222560f,
    0.00004871840f, 0.00004544747f, 0.00004240000f, 0.00003956104f, 0.00003691512f, 0.00003444868f, 0.00003214816f,
    0.00003000000f, 0.00002799125f, 0.00002611356f, 0.00002436024f, 0.00002272461f, 0.00002120000f, 0.00001977855f,
    0.00001845285f, 0.00001721687f, 0.00001606459f});
static constexpr Spectrum CIE1931_Z({0.00645000100f, 0.00708321600f, 0.00774548800f, 0.00850115200f, 0.00941454400f,
    0.01054999000f, 0.01196580000f, 0.01365587000f, 0.01558805000f, 0.01773015000f, 0.02005001000f, 0.02251136000f,
    0.02520288000f, 0.02827972000f, 0.03189704000f, 0.03621000000f, 0.04143771000f, 0.04750372000f, 0.05411988000f,
    0.06099803000f, 0.06785001000f, 0.07448632000f, 0.08136156000f, 0.08915364000f, 0.09854048000f, 0.11020000000f,
    0.12461330000f, 0.14170170000f, 0.16130350000f, 0.18325680000f, 0.20740000000f, 0.23369210000f, 0.26261140000f,
    0.29477460000f, 0.33079850000f, 0.37130000000f, 0.41620910000f, 0.46546420000f, 0.51969480000f, 0.57953030000f,
    0.64560000000f, 0.71848380000f, 0.79671330000f, 0.87784590000f, 0.95943900000f, 1.03905010000f, 1.11536730000f,
    1.18849710000f, 1.25812330000f, 1.32392960000f, 1.38560000000f, 1.44263520000f, 1.49480350000f, 1.54219030000f,
    1.58488070000f, 1.62296000000f, 1.65640480000f, 1.68529590000f, 1.70987450000f, 1.73038210000f, 1.74706000000f,
    1.76004460000f, 1.76962330000f, 1.77626370000f, 1.78043340000f, 1.78260000000f, 1.78296820000f, 1.78169980000f,
    1.77919820000f, 1.77586710000f, 1.77211000000f, 1.76825890000f, 1.76403900000f, 1.75894380000f, 1.75246630000f,
    1.74410000000f, 1.73355950000f, 1.72085810000f, 1.70593690000f, 1.68873720000f, 1.66920000000f, 1.64752870000f,
    1.62341270000f, 1.59602230000f, 1.56452800000f, 1.52810000000f, 1.48611140000f, 1.43952150000f, 1.38987990000f,
    1.33873620000f, 1.28764000000f, 1.23742230000f, 1.18782430000f, 1.13876110000f, 1.09014800000f, 1.04190000000f,
    0.99419760000f, 0.94734730000f, 0.90145310000f, 0.85661930000f, 0.81295010000f, 0.77051730000f, 0.72944480000f,
    0.68991360000f, 0.65210490000f, 0.61620000000f, 0.58232860000f, 0.55041620000f, 0.52033760000f, 0.49196730000f,
    0.46518000000f, 0.43992460000f, 0.41618360000f, 0.39388220000f, 0.37294590000f, 0.35330000000f, 0.33485780000f,
    0.31755210000f, 0.30133750000f, 0.28616860000f, 0.27200000000f, 0.25881710000f, 0.24648380000f, 0.23477180000f,
    0.22345330000f, 0.21230000000f, 0.20116920000f, 0.19011960000f, 0.17922540000f, 0.16856080000f, 0.15820000000f,
    0.14813830000f, 0.13837580000f, 0.12899420000f, 0.12007510000f, 0.11170000000f, 0.10390480000f, 0.09666748000f,
    0.08998272000f, 0.08384531000f, 0.07824999000f, 0.07320899000f, 0.06867816000f, 0.06456784000f, 0.06078835000f,
    0.05725001000f, 0.05390435000f, 0.05074664000f, 0.04775276000f, 0.04489859000f, 0.04216000000f, 0.03950728000f,
    0.03693564000f, 0.03445836000f, 0.03208872000f, 0.02984000000f, 0.02771181000f, 0.02569444000f, 0.02378716000f,
    0.02198925000f, 0.02030000000f, 0.01871805000f, 0.01724036000f, 0.01586364000f, 0.01458461000f, 0.01340000000f,
    0.01230723000f, 0.01130188000f, 0.01037792000f, 0.00952930600f, 0.00874999900f, 0.00803520000f, 0.00738160000f,
    0.00678540000f, 0.00624280000f, 0.00574999900f, 0.00530360000f, 0.00489980000f, 0.00453420000f, 0.00420240000f,
    0.00390000000f, 0.00362320000f, 0.00337060000f, 0.00314140000f, 0.00293480000f, 0.00274999900f, 0.00258520000f,
    0.00243860000f, 0.00230940000f, 0.00219680000f, 0.00210000000f, 0.00201773300f, 0.00194820000f, 0.00188980000f,
    0.00184093300f, 0.00180000000f, 0.00176626700f, 0.00173780000f, 0.00171120000f, 0.00168306700f, 0.00165000100f,
    0.00161013300f, 0.00156440000f, 0.00151360000f, 0.00145853300f, 0.00140000000f, 0.00133666700f, 0.00127000000f,
    0.00120500000f, 0.00114666700f, 0.00110000000f, 0.00106880000f, 0.00104940000f, 0.00103560000f, 0.00102120000f,
    0.00100000000f, 0.00096864000f, 0.00092992000f, 0.00088688000f, 0.00084256000f, 0.00080000000f, 0.00076096000f,
    0.00072368000f, 0.00068592000f, 0.00064544000f, 0.00060000000f, 0.00054786670f, 0.00049160000f, 0.00043540000f,
    0.00038346670f, 0.00034000000f, 0.00030725330f, 0.00028316000f, 0.00026544000f, 0.00025181330f, 0.00024000000f,
    0.00022954670f, 0.00022064000f, 0.00021196000f, 0.00020218670f, 0.00019000000f, 0.00017421330f, 0.00015564000f,
    0.00013596000f, 0.00011685330f, 0.00010000000f, 0.00008613333f, 0.00007460000f, 0.00006500000f, 0.00005693333f,
    0.00004999999f, 0.00004416000f, 0.00003948000f, 0.00003572000f, 0.00003264000f, 0.00003000000f, 0.00002765333f,
    0.00002556000f, 0.00002364000f, 0.00002181333f, 0.00002000000f, 0.00001813333f, 0.00001620000f, 0.00001420000f,
    0.00001213333f, 0.00001000000f, 0.00000773333f, 0.00000540000f, 0.00000320000f, 0.00000133333f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f});

static constexpr Spectrum CIE2012_X({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0037696470f, 0.0045324160f,
    0.0054465530f, 0.0065388680f, 0.0078396990f, 0.0093829670f, 0.0112060800f, 0.0133496500f, 0.0158569000f,
    0.0187728600f, 0.0221430200f, 0.0260128500f, 0.0304303600f, 0.0354432500f, 0.0410964000f, 0.0474298600f,
    0.0544739400f, 0.0622361200f, 0.0707004800f, 0.0798251300f, 0.0895380300f, 0.0997484800f, 0.1104019000f,
    0.1214566000f, 0.1328741000f, 0.1446214000f, 0.1566468000f, 0.1687901000f, 0.1808328000f, 0.1925216000f,
    0.2035729000f, 0.2137531000f, 0.2231348000f, 0.2319245000f, 0.2403892000f, 0.2488523000f, 0.2575896000f,
    0.2664991000f, 0.2753532000f, 0.2838921000f, 0.2918246000f, 0.2989200000f, 0.3052993000f, 0.3112031000f,
    0.3169047000f, 0.3227087000f, 0.3288194000f, 0.3349242000f, 0.3405452000f, 0.3451688000f, 0.3482554000f,
    0.3494153000f, 0.3489075000f, 0.3471746000f, 0.3446705000f, 0.3418483000f, 0.3390240000f, 0.3359926000f,
    0.3324276000f, 0.3280157000f, 0.3224637000f, 0.3156225000f, 0.3078201000f, 0.2994771000f, 0.2909776000f,
    0.2826646000f, 0.2747962000f, 0.2674312000f, 0.2605847000f, 0.2542749000f, 0.2485254000f, 0.2433039000f,
    0.2383414000f, 0.2333253000f, 0.2279619000f, 0.2219781000f, 0.2151735000f, 0.2075619000f, 0.1992183000f,
    0.1902290000f, 0.1806905000f, 0.1707154000f, 0.1604471000f, 0.1500244000f, 0.1395705000f, 0.1291920000f,
    0.1189859000f, 0.1090615000f, 0.0995142400f, 0.0904185000f, 0.0818289500f, 0.0737681700f, 0.0661947700f,
    0.0590638000f, 0.0523424200f, 0.0460086500f, 0.0400615400f, 0.0345437300f, 0.0294909100f, 0.0249214000f,
    0.0208398100f, 0.0172359100f, 0.0140792400f, 0.0113451600f, 0.0090196580f, 0.0070977310f, 0.0055711450f,
    0.0043945660f, 0.0035163030f, 0.0028876380f, 0.0024615880f, 0.0022063480f, 0.0021495590f, 0.0023370910f,
    0.0028189310f, 0.0036491780f, 0.0048913590f, 0.0066293640f, 0.0089429020f, 0.0119022400f, 0.0155698900f,
    0.0199766800f, 0.0250469800f, 0.0306753000f, 0.0367499900f, 0.0431517100f, 0.0497858400f, 0.0566855400f,
    0.0639165100f, 0.0715435200f, 0.0796291700f, 0.0882147300f, 0.0972697800f, 0.1067504000f, 0.1166192000f,
    0.1268468000f, 0.1374060000f, 0.1482471000f, 0.1593076000f, 0.1705181000f, 0.1818026000f, 0.1931090000f,
    0.2045085000f, 0.2161166000f, 0.2280650000f, 0.2405015000f, 0.2535441000f, 0.2671300000f, 0.2811351000f,
    0.2954164000f, 0.3098117000f, 0.3241678000f, 0.3384319000f, 0.3525786000f, 0.3665839000f, 0.3804244000f,
    0.3940988000f, 0.4076972000f, 0.4213484000f, 0.4352003000f, 0.4494206000f, 0.4641616000f, 0.4794395000f,
    0.4952180000f, 0.5114395000f, 0.5280233000f, 0.5448696000f, 0.5618898000f, 0.5790137000f, 0.5961882000f,
    0.6133784000f, 0.6305897000f, 0.6479223000f, 0.6654866000f, 0.6833782000f, 0.7016774000f, 0.7204110000f,
    0.7394495000f, 0.7586285000f, 0.7777885000f, 0.7967750000f, 0.8154530000f, 0.8337389000f, 0.8515493000f,
    0.8687862000f, 0.8853376000f, 0.9011588000f, 0.9165278000f, 0.9318245000f, 0.9474524000f, 0.9638388000f,
    0.9812596000f, 0.9992953000f, 1.0173430000f, 1.0347900000f, 1.0510110000f, 1.0655220000f, 1.0784210000f,
    1.0899440000f, 1.1003200000f, 1.1097670000f, 1.1184380000f, 1.1262660000f, 1.1331380000f, 1.1389520000f,
    1.1436200000f, 1.1470950000f, 1.1494640000f, 1.1508380000f, 1.1513260000f, 1.1510330000f, 1.1500020000f,
    1.1480610000f, 1.1449980000f, 1.1406220000f, 1.1347570000f, 1.1272980000f, 1.1183420000f, 1.1080330000f,
    1.0965150000f, 1.0839280000f, 1.0703870000f, 1.0559340000f, 1.0405920000f, 1.0243850000f, 1.0073440000f,
    0.9895268000f, 0.9711213000f, 0.9523257000f, 0.9333248000f, 0.9142877000f, 0.8952798000f, 0.8760157000f,
    0.8561607000f, 0.8354235000f, 0.8135565000f, 0.7904565000f, 0.7664364000f, 0.7418777000f, 0.7171219000f,
    0.6924717000f, 0.6681600000f, 0.6442697000f, 0.6208450000f, 0.5979243000f, 0.5755410000f, 0.5537296000f,
    0.5325412000f, 0.5120218000f, 0.4922070000f, 0.4731224000f, 0.4547417000f, 0.4368719000f, 0.4193121000f,
    0.4018980000f, 0.3844986000f, 0.3670592000f, 0.3497167000f, 0.3326305000f, 0.3159341000f, 0.2997374000f,
    0.2841189000f, 0.2691053000f, 0.2547077000f, 0.2409319000f, 0.2277792000f, 0.2152431000f, 0.2033010000f,
    0.1919276000f, 0.1810987000f, 0.1707914000f, 0.1609842000f, 0.1516577000f, 0.1427936000f, 0.1343737000f,
    0.1263808000f, 0.1187979000f, 0.1116088000f, 0.1047975000f, 0.0983483500f, 0.0922459700f, 0.0864750600f,
    0.0810198600f, 0.0758651400f, 0.0709963300f, 0.0663996000f, 0.0620622500f, 0.0579740900f, 0.0541253300f,
    0.0505060000f, 0.0471060600f, 0.0439141100f, 0.0409141100f, 0.0380906700f, 0.0354303400f, 0.0329213800f,
    0.0305567200f, 0.0283414600f, 0.0262803300f, 0.0243746500f, 0.0226230600f, 0.0210193500f, 0.0195464700f,
    0.0181872700f, 0.0169272700f, 0.0157541700f, 0.0146585400f, 0.0136357100f, 0.0126820500f, 0.0117939400f,
    0.0109677800f, 0.0101996400f, 0.0094843170f, 0.0088168510f, 0.0081929210f, 0.0076087500f, 0.0070613910f,
    0.0065495090f, 0.0060719700f, 0.0056274760f, 0.0052146080f, 0.0048318480f, 0.0044775790f, 0.0041501660f,
    0.0038479880f, 0.0035694520f, 0.0033128570f, 0.0030760220f, 0.0028568940f, 0.0026536810f, 0.0024648210f,
    0.0022890600f, 0.0021256940f, 0.0019741210f, 0.0018337230f, 0.0017038760f, 0.0015839040f, 0.0014729390f,
    0.0013701510f, 0.0012748030f, 0.0011862380f, 0.0011038710f, 0.0010271940f, 0.0009557493f, 0.0008891262f,
    0.0008269535f, 0.0007689351f, 0.0007149425f, 0.0006648590f, 0.0006185421f, 0.0005758303f, 0.0005365046f,
    0.0005001842f, 0.0004665005f, 0.0004351386f, 0.0004058303f, 0.0003783733f, 0.0003526892f, 0.0003287199f,
    0.0003063998f, 0.0002856577f, 0.0002664108f, 0.0002485462f, 0.0002319529f, 0.0002165300f, 0.0002021853f,
    0.0001888338f, 0.0001763935f, 0.0001647895f, 0.0001539542f, 0.0001438270f, 0.0001343572f, 0.0001255141f,
    0.0001172706f, 0.0001095983f, 0.0001024685f, 0.0000958472f, 0.0000896832f, 0.0000839273f, 0.0000785371f,
    0.0000734755f, 0.0000687158f, 0.0000642526f, 0.0000600829f, 0.0000562010f, 0.0000525987f, 0.0000492628f,
    0.0000461662f, 0.0000432821f, 0.0000405872f});

static constexpr Spectrum CIE2012_Y({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0004146161f, 0.0005028333f,
    0.0006084991f, 0.0007344436f, 0.0008837389f, 0.0010596460f, 0.0012655320f, 0.0015047530f, 0.0017804930f,
    0.0020955720f, 0.0024521940f, 0.0028522160f, 0.0032991150f, 0.0037974660f, 0.0043527680f, 0.0049717170f,
    0.0056610140f, 0.0064216150f, 0.0072503120f, 0.0081401730f, 0.0090798600f, 0.0100560800f, 0.0110645600f,
    0.0121052200f, 0.0131801400f, 0.0142937700f, 0.0154500400f, 0.0166409300f, 0.0178530200f, 0.0190701800f,
    0.0202736900f, 0.0214480500f, 0.0226004100f, 0.0237478900f, 0.0249124700f, 0.0261210600f, 0.0273992300f,
    0.0287499300f, 0.0301690900f, 0.0316514500f, 0.0331903800f, 0.0347791200f, 0.0364149500f, 0.0380956900f,
    0.0398184300f, 0.0415794000f, 0.0433709800f, 0.0451718000f, 0.0469542000f, 0.0486871800f, 0.0503365700f,
    0.0518761100f, 0.0533221800f, 0.0547060300f, 0.0560633500f, 0.0574339300f, 0.0588510700f, 0.0603080900f,
    0.0617864400f, 0.0632657000f, 0.0647235200f, 0.0661474900f, 0.0675725600f, 0.0690492800f, 0.0706328000f,
    0.0723833900f, 0.0743596000f, 0.0765938300f, 0.0791143600f, 0.0819534500f, 0.0851481600f, 0.0887265700f,
    0.0926600800f, 0.0968972300f, 0.1013746000f, 0.1060145000f, 0.1107377000f, 0.1155111000f, 0.1203122000f,
    0.1251161000f, 0.1298957000f, 0.1346299000f, 0.1393309000f, 0.1440235000f, 0.1487372000f, 0.1535066000f,
    0.1583644000f, 0.1633199000f, 0.1683761000f, 0.1735365000f, 0.1788048000f, 0.1841819000f, 0.1896559000f,
    0.1952101000f, 0.2008259000f, 0.2064828000f, 0.2121826000f, 0.2180279000f, 0.2241586000f, 0.2307302000f,
    0.2379160000f, 0.2458706000f, 0.2546023000f, 0.2640760000f, 0.2742490000f, 0.2850680000f, 0.2964837000f,
    0.3085010000f, 0.3211393000f, 0.3344175000f, 0.3483536000f, 0.3629601000f, 0.3782275000f, 0.3941359000f,
    0.4106582000f, 0.4277595000f, 0.4453993000f, 0.4635396000f, 0.4821376000f, 0.5011430000f, 0.5204972000f,
    0.5401387000f, 0.5600208000f, 0.5800972000f, 0.6003172000f, 0.6206256000f, 0.6409398000f, 0.6610772000f,
    0.6808134000f, 0.6999044000f, 0.7180890000f, 0.7351593000f, 0.7511821000f, 0.7663143000f, 0.7807352000f,
    0.7946448000f, 0.8082074000f, 0.8213817000f, 0.8340701000f, 0.8461711000f, 0.8575799000f, 0.8682408000f,
    0.8783061000f, 0.8879907000f, 0.8975211000f, 0.9071347000f, 0.9169947000f, 0.9269295000f, 0.9366731000f,
    0.9459482000f, 0.9544675000f, 0.9619834000f, 0.9684390000f, 0.9738289000f, 0.9781519000f, 0.9814106000f,
    0.9836669000f, 0.9852081000f, 0.9863813000f, 0.9875357000f, 0.9890228000f, 0.9910811000f, 0.9934913000f,
    0.9959172000f, 0.9980205000f, 0.9994608000f, 0.9999930000f, 0.9997557000f, 0.9989839000f, 0.9979123000f,
    0.9967737000f, 0.9957356000f, 0.9947115000f, 0.9935534000f, 0.9921156000f, 0.9902549000f, 0.9878596000f,
    0.9849324000f, 0.9815036000f, 0.9776035000f, 0.9732611000f, 0.9684764000f, 0.9631369000f, 0.9571062000f,
    0.9502540000f, 0.9424569000f, 0.9336897000f, 0.9242893000f, 0.9146707000f, 0.9052333000f, 0.8963613000f,
    0.8883069000f, 0.8808462000f, 0.8736445000f, 0.8663755000f, 0.8587203000f, 0.8504295000f, 0.8415047000f,
    0.8320109000f, 0.8220154000f, 0.8115868000f, 0.8007874000f, 0.7896515000f, 0.7782053000f, 0.7664733000f,
    0.7544785000f, 0.7422473000f, 0.7298229000f, 0.7172525000f, 0.7045818000f, 0.6918553000f, 0.6791009000f,
    0.6662846000f, 0.6533595000f, 0.6402807000f, 0.6270066000f, 0.6135148000f, 0.5998494000f, 0.5860682000f,
    0.5722261000f, 0.5583746000f, 0.5445535000f, 0.5307673000f, 0.5170130000f, 0.5032889000f, 0.4895950000f,
    0.4759442000f, 0.4623958000f, 0.4490154000f, 0.4358622000f, 0.4229897000f, 0.4104152000f, 0.3980356000f,
    0.3857300000f, 0.3733907000f, 0.3609245000f, 0.3482860000f, 0.3355702000f, 0.3228963000f, 0.3103704000f,
    0.2980865000f, 0.2861160000f, 0.2744822000f, 0.2631953000f, 0.2522628000f, 0.2416902000f, 0.2314809000f,
    0.2216378000f, 0.2121622000f, 0.2030542000f, 0.1943124000f, 0.1859227000f, 0.1778274000f, 0.1699654000f,
    0.1622841000f, 0.1547397000f, 0.1473081000f, 0.1400169000f, 0.1329013000f, 0.1259913000f, 0.1193120000f,
    0.1128820000f, 0.1067113000f, 0.1008052000f, 0.0951665300f, 0.0897959400f, 0.0846904400f, 0.0798400900f,
    0.0752337200f, 0.0708606100f, 0.0667104500f, 0.0627736000f, 0.0590417900f, 0.0555070300f, 0.0521613900f,
    0.0489969900f, 0.0460057800f, 0.0431788500f, 0.0405075500f, 0.0379837600f, 0.0355998200f, 0.0333485600f,
    0.0312233200f, 0.0292178000f, 0.0273260100f, 0.0255422300f, 0.0238612100f, 0.0222785900f, 0.0207902000f,
    0.0193918500f, 0.0180793900f, 0.0168481700f, 0.0156918800f, 0.0146044600f, 0.0135806200f, 0.0126157300f,
    0.0117069600f, 0.0108560800f, 0.0100647600f, 0.0093333760f, 0.0086612840f, 0.0080460480f, 0.0074811300f,
    0.0069599870f, 0.0064770700f, 0.0060276770f, 0.0056081690f, 0.0052166910f, 0.0048517850f, 0.0045120080f,
    0.0041959410f, 0.0039020570f, 0.0036283710f, 0.0033730050f, 0.0031343150f, 0.0029108640f, 0.0027015280f,
    0.0025057960f, 0.0023232310f, 0.0021533330f, 0.0019955570f, 0.0018493160f, 0.0017139760f, 0.0015888990f,
    0.0014734530f, 0.0013670220f, 0.0012689540f, 0.0011784210f, 0.0010946440f, 0.0010169430f, 0.0009447269f,
    0.0008775171f, 0.0008150438f, 0.0007570755f, 0.0007033755f, 0.0006537050f, 0.0006078048f, 0.0005653435f,
    0.0005260046f, 0.0004895061f, 0.0004555970f, 0.0004240548f, 0.0003946860f, 0.0003673178f, 0.0003417941f,
    0.0003179738f, 0.0002957441f, 0.0002750558f, 0.0002558640f, 0.0002381142f, 0.0002217445f, 0.0002066711f,
    0.0001927474f, 0.0001798315f, 0.0001678023f, 0.0001565566f, 0.0001460168f, 0.0001361535f, 0.0001269451f,
    0.0001183671f, 0.0001103928f, 0.0001029908f, 0.0000961184f, 0.0000897332f, 0.0000837969f, 0.0000782744f,
    0.0000731331f, 0.0000683414f, 0.0000638704f, 0.0000596939f, 0.0000557886f, 0.0000521351f, 0.0000487218f,
    0.0000455385f, 0.0000425744f, 0.0000398188f, 0.0000372588f, 0.0000348747f, 0.0000326477f, 0.0000305614f,
    0.0000286018f, 0.0000267584f, 0.0000250294f, 0.0000234137f, 0.0000219091f, 0.0000205126f, 0.0000192190f,
    0.0000180180f, 0.0000168990f, 0.0000158531f});
static constexpr Spectrum CIE2012_Z({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0184726000f, 0.0222110100f,
    0.0266981900f, 0.0320693700f, 0.0384783200f, 0.0460978400f, 0.0551195300f, 0.0657525700f, 0.0782211300f,
    0.0927601300f, 0.1096090000f, 0.1290077000f, 0.1512047000f, 0.1764441000f, 0.2049517000f, 0.2369246000f,
    0.2725123000f, 0.3117820000f, 0.3547064000f, 0.4011473000f, 0.4508369000f, 0.5034164000f, 0.5586361000f,
    0.6162734000f, 0.6760982000f, 0.7378822000f, 0.8013019000f, 0.8655573000f, 0.9295791000f, 0.9921293000f,
    1.0518210000f, 1.1075090000f, 1.1595270000f, 1.2088690000f, 1.2568340000f, 1.3050080000f, 1.3547580000f,
    1.4055940000f, 1.4564140000f, 1.5059600000f, 1.5528260000f, 1.5959020000f, 1.6357680000f, 1.6735730000f,
    1.7106040000f, 1.7482800000f, 1.7875040000f, 1.8266090000f, 1.8631080000f, 1.8943320000f, 1.9174790000f,
    1.9305290000f, 1.9348190000f, 1.9326500000f, 1.9263950000f, 1.9184370000f, 1.9104300000f, 1.9012240000f,
    1.8890000000f, 1.8719960000f, 1.8485450000f, 1.8177920000f, 1.7816270000f, 1.7425140000f, 1.7027490000f,
    1.6644390000f, 1.6292070000f, 1.5973600000f, 1.5688960000f, 1.5438230000f, 1.5221570000f, 1.5036110000f,
    1.4866730000f, 1.4695950000f, 1.4507090000f, 1.4284400000f, 1.4015870000f, 1.3700940000f, 1.3342200000f,
    1.2942750000f, 1.2506100000f, 1.2036960000f, 1.1543160000f, 1.1032840000f, 1.0513470000f, 0.9991789000f,
    0.9473958000f, 0.8966222000f, 0.8473981000f, 0.8001576000f, 0.7552379000f, 0.7127879000f, 0.6725198000f,
    0.6340976000f, 0.5972433000f, 0.5617313000f, 0.5274921000f, 0.4948809000f, 0.4642586000f, 0.4358841000f,
    0.4099313000f, 0.3864261000f, 0.3650566000f, 0.3454812000f, 0.3274095000f, 0.3105939000f, 0.2948102000f,
    0.2798194000f, 0.2654100000f, 0.2514084000f, 0.2376753000f, 0.2241211000f, 0.2107484000f, 0.1975839000f,
    0.1846574000f, 0.1720018000f, 0.1596918000f, 0.1479415000f, 0.1369428000f, 0.1268279000f, 0.1176796000f,
    0.1094970000f, 0.1020943000f, 0.0952799300f, 0.0889007500f, 0.0828354800f, 0.0770098200f, 0.0714400100f,
    0.0661543600f, 0.0611719900f, 0.0565040700f, 0.0521512100f, 0.0480956600f, 0.0443172000f, 0.0407973400f,
    0.0375191200f, 0.0344684600f, 0.0316376400f, 0.0290190100f, 0.0266036400f, 0.0243816400f, 0.0223409700f,
    0.0204641500f, 0.0187345600f, 0.0171378800f, 0.0156617400f, 0.0142964400f, 0.0130370200f, 0.0118789700f,
    0.0108172500f, 0.0098464700f, 0.0089606870f, 0.0081528110f, 0.0074160250f, 0.0067441150f, 0.0061314210f,
    0.0055727780f, 0.0050634630f, 0.0045991690f, 0.0041759710f, 0.0037902910f, 0.0034389520f, 0.0031193410f,
    0.0028290380f, 0.0025657220f, 0.0023271860f, 0.0021112800f, 0.0019157660f, 0.0017385890f, 0.0015779200f,
    0.0014321280f, 0.0012997810f, 0.0011796670f, 0.0010706940f, 0.0009718623f, 0.0008822531f, 0.0008010231f,
    0.0007273884f, 0.0006606347f, 0.0006001146f, 0.0005452416f, 0.0004954847f, 0.0004503642f, 0.0004094455f,
    0.0003723345f, 0.0003386739f, 0.0003081396f, 0.0002804370f, 0.0002552996f, 0.0002324859f, 0.0002117772f,
    0.0001929758f, 0.0001759024f, 0.0001603947f, 0.0001463059f, 0.0001335031f, 0.0001218660f, 0.0001112857f,
    0.0001016634f, 0.0000929100f, 0.0000849447f, 0.0000776943f, 0.0000710925f, 0.0000650794f, 0.0000596006f,
    0.0000546071f, 0.0000500542f, 0.0000459016f, 0.0000421127f, 0.0000386544f, 0.0000354966f, 0.0000326122f,
    0.0000299764f, 0.0000275669f, 0.0000253634f, 0.0000233474f, 0.0000215022f, 0.0000198127f, 0.0000182650f,
    0.0000168467f, 0.0000155463f, 0.0000143536f, 0.0000132592f, 0.0000122544f, 0.0000113317f, 0.0000104839f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f});

// taken from http://cvrl.ioo.ucl.ac.uk/, lerped by 1nm
static constexpr Spectrum CIE1931_JuddVos_X({0.0026899f, 0.00321402f, 0.00373814f, 0.00426226f, 0.00478638f, 0.0053105f,
    0.0064046f, 0.0074987f, 0.0085928f, 0.0096869f, 0.010781f, 0.0127832f, 0.0147854f, 0.0167876f, 0.0187898f,
    0.020792f, 0.0242298f, 0.0276676f, 0.0311054f, 0.0345432f, 0.037981f, 0.0430162f, 0.0480514f, 0.0530866f,
    0.0581218f, 0.063157f, 0.0705138f, 0.0778706f, 0.0852274f, 0.0925842f, 0.099941f, 0.1116008f, 0.1232606f,
    0.1349204f, 0.1465802f, 0.15824f, 0.172488f, 0.186736f, 0.200984f, 0.215232f, 0.22948f, 0.2398f, 0.25012f, 0.26044f,
    0.27076f, 0.28108f, 0.287054f, 0.293028f, 0.299002f, 0.304976f, 0.31095f, 0.314904f, 0.318858f, 0.322812f,
    0.326766f, 0.33072f, 0.331248f, 0.331776f, 0.332304f, 0.332832f, 0.33336f, 0.330032f, 0.326704f, 0.323376f,
    0.320048f, 0.31672f, 0.31114f, 0.30556f, 0.29998f, 0.2944f, 0.28882f, 0.282994f, 0.277168f, 0.271342f, 0.265516f,
    0.25969f, 0.254304f, 0.248918f, 0.243532f, 0.238146f, 0.23276f, 0.228206f, 0.223652f, 0.219098f, 0.214544f,
    0.20999f, 0.202944f, 0.195898f, 0.188852f, 0.181806f, 0.17476f, 0.166382f, 0.158004f, 0.149626f, 0.141248f,
    0.13287f, 0.1246848f, 0.1164996f, 0.1083144f, 0.1001292f, 0.091944f, 0.0849522f, 0.0779604f, 0.0709686f, 0.0639768f,
    0.056985f, 0.0519342f, 0.0468834f, 0.0418326f, 0.0367818f, 0.031731f, 0.0283074f, 0.0248838f, 0.0214602f,
    0.0180366f, 0.014613f, 0.01266022f, 0.01070744f, 0.00875466f, 0.00680188f, 0.0048491f, 0.00434358f, 0.00383806f,
    0.00333254f, 0.00282702f, 0.0023215f, 0.00371518f, 0.00510886f, 0.00650254f, 0.00789622f, 0.0092899f, 0.01328752f,
    0.01728514f, 0.02128276f, 0.02528038f, 0.029278f, 0.0361806f, 0.0430832f, 0.0499858f, 0.0568884f, 0.063791f,
    0.0731948f, 0.0825986f, 0.0920024f, 0.1014062f, 0.11081f, 0.122032f, 0.133254f, 0.144476f, 0.155698f, 0.16692f,
    0.179072f, 0.191224f, 0.203376f, 0.215528f, 0.22768f, 0.240682f, 0.253684f, 0.266686f, 0.279688f, 0.29269f,
    0.306602f, 0.320514f, 0.334426f, 0.348338f, 0.36225f, 0.37707f, 0.39189f, 0.40671f, 0.42153f, 0.43635f, 0.452106f,
    0.467862f, 0.483618f, 0.499374f, 0.51513f, 0.5316f, 0.54807f, 0.56454f, 0.58101f, 0.59748f, 0.614226f, 0.630972f,
    0.647718f, 0.664464f, 0.68121f, 0.697818f, 0.714426f, 0.731034f, 0.747642f, 0.76425f, 0.780188f, 0.796126f,
    0.812064f, 0.828002f, 0.84394f, 0.858422f, 0.872904f, 0.887386f, 0.901868f, 0.91635f, 0.928486f, 0.940622f,
    0.952758f, 0.964894f, 0.97703f, 0.986224f, 0.995418f, 1.004612f, 1.013806f, 1.023f, 1.02866f, 1.03432f, 1.03998f,
    1.04564f, 1.0513f, 1.05204f, 1.05278f, 1.05352f, 1.05426f, 1.055f, 1.05124f, 1.04748f, 1.04372f, 1.03996f, 1.0362f,
    1.027438f, 1.018676f, 1.009914f, 1.001152f, 0.99239f, 0.979634f, 0.966878f, 0.954122f, 0.941366f, 0.92861f,
    0.91158f, 0.89455f, 0.87752f, 0.86049f, 0.84346f, 0.822734f, 0.802008f, 0.781282f, 0.760556f, 0.73983f, 0.718442f,
    0.697054f, 0.675666f, 0.654278f, 0.63289f, 0.613014f, 0.593138f, 0.573262f, 0.553386f, 0.53351f, 0.514932f,
    0.496354f, 0.477776f, 0.459198f, 0.44062f, 0.423402f, 0.406184f, 0.388966f, 0.371748f, 0.35453f, 0.339348f,
    0.324166f, 0.308984f, 0.293802f, 0.27862f, 0.265866f, 0.253112f, 0.240358f, 0.227604f, 0.21485f, 0.204202f,
    0.193554f, 0.182906f, 0.172258f, 0.16161f, 0.152928f, 0.144246f, 0.135564f, 0.126882f, 0.1182f, 0.1117106f,
    0.1052212f, 0.0987318f, 0.0922424f, 0.085753f, 0.0812178f, 0.0766826f, 0.0721474f, 0.0676122f, 0.063077f,
    0.0596284f, 0.0561798f, 0.0527312f, 0.0492826f, 0.045834f, 0.0430786f, 0.0403232f, 0.0375678f, 0.0348124f,
    0.032057f, 0.030083f, 0.028109f, 0.026135f, 0.024161f, 0.022187f, 0.020872f, 0.019557f, 0.018242f, 0.016927f,
    0.015612f, 0.0147092f, 0.0138064f, 0.0129036f, 0.0120008f, 0.011098f, 0.01046306f, 0.00982812f, 0.00919318f,
    0.00855824f, 0.0079233f, 0.00746926f, 0.00701522f, 0.00656118f, 0.00610714f, 0.0056531f, 0.00532326f, 0.00499342f,
    0.00466358f, 0.00433374f, 0.0040039f, 0.00376818f, 0.00353246f, 0.00329674f, 0.00306102f, 0.0028253f, 0.00265918f,
    0.00249306f, 0.00232694f, 0.00216082f, 0.0019947f, 0.00187564f, 0.00175658f, 0.00163752f, 0.00151846f, 0.0013994f,
    0.00131348f, 0.00122756f, 0.00114164f, 0.00105572f, 0.0009698f, 0.000909534f, 0.000849268f, 0.000789002f,
    0.000728736f, 0.00066847f, 0.000627058f, 0.000585646f, 0.000544234f, 0.000502822f, 0.00046141f, 0.000433274f,
    0.000405138f, 0.000377002f, 0.000348866f, 0.00032073f, 0.00030173f, 0.00028273f, 0.00026373f, 0.00024473f,
    0.00022573f, 0.00021253f, 0.00019933f, 0.00018613f, 0.00017293f, 0.00015973f, 0.000150334f, 0.000140938f,
    0.000131542f, 0.000122146f, 0.00011275f, 0.000106103f, 9.94552E-05f, 9.28078E-05f, 8.61604E-05f, 0.000079513f,
    7.48278E-05f, 7.01426E-05f, 6.54574E-05f, 6.07722E-05f, 0.000056087f, 5.27778E-05f, 4.94686E-05f, 4.61594E-05f,
    4.28502E-05f});

static constexpr Spectrum CIE1931_JuddVos_Y({0.0002f, 0.000239112f, 0.000278224f, 0.000317336f, 0.000356448f,
    0.00039556f, 0.000476448f, 0.000557336f, 0.000638224f, 0.000719112f, 0.0008f, 0.00094914f, 0.00109828f, 0.00124742f,
    0.00139656f, 0.0015457f, 0.00179656f, 0.00204742f, 0.00229828f, 0.00254914f, 0.0028f, 0.00317124f, 0.00354248f,
    0.00391372f, 0.00428496f, 0.0046562f, 0.00520496f, 0.00575372f, 0.00630248f, 0.00685124f, 0.0074f, 0.0082758f,
    0.0091516f, 0.0100274f, 0.0109032f, 0.011779f, 0.0129232f, 0.0140674f, 0.0152116f, 0.0163558f, 0.0175f, 0.0185356f,
    0.0195712f, 0.0206068f, 0.0216424f, 0.022678f, 0.0236024f, 0.0245268f, 0.0254512f, 0.0263756f, 0.0273f, 0.0283568f,
    0.0294136f, 0.0304704f, 0.0315272f, 0.032584f, 0.0336472f, 0.0347104f, 0.0357736f, 0.0368368f, 0.0379f, 0.0387982f,
    0.0396964f, 0.0405946f, 0.0414928f, 0.042391f, 0.0432728f, 0.0441546f, 0.0450364f, 0.0459182f, 0.0468f, 0.0478644f,
    0.0489288f, 0.0499932f, 0.0510576f, 0.052122f, 0.0536976f, 0.0552732f, 0.0568488f, 0.0584244f, 0.06f, 0.0625884f,
    0.0651768f, 0.0677652f, 0.0703536f, 0.072942f, 0.0765496f, 0.0801572f, 0.0837648f, 0.0873724f, 0.09098f, 0.095352f,
    0.099724f, 0.104096f, 0.108468f, 0.11284f, 0.118076f, 0.123312f, 0.128548f, 0.133784f, 0.13902f, 0.14519f, 0.15136f,
    0.15753f, 0.1637f, 0.16987f, 0.1775f, 0.18513f, 0.19276f, 0.20039f, 0.20802f, 0.218032f, 0.228044f, 0.238056f,
    0.248068f, 0.25808f, 0.271064f, 0.284048f, 0.297032f, 0.310016f, 0.323f, 0.33948f, 0.35596f, 0.37244f, 0.38892f,
    0.4054f, 0.42492f, 0.44444f, 0.46396f, 0.48348f, 0.503f, 0.524022f, 0.545044f, 0.566066f, 0.587088f, 0.60811f,
    0.628488f, 0.648866f, 0.669244f, 0.689622f, 0.71f, 0.72702f, 0.74404f, 0.76106f, 0.77808f, 0.7951f, 0.80848f,
    0.82186f, 0.83524f, 0.84862f, 0.862f, 0.87261f, 0.88322f, 0.89383f, 0.90444f, 0.91505f, 0.92284f, 0.93063f,
    0.93842f, 0.94621f, 0.954f, 0.959208f, 0.964416f, 0.969624f, 0.974832f, 0.98004f, 0.983022f, 0.986004f, 0.988986f,
    0.991968f, 0.99495f, 0.99598f, 0.99701f, 0.99804f, 0.99907f, 1.0001f, 0.99908f, 0.99806f, 0.99704f, 0.99602f,
    0.995f, 0.99175f, 0.9885f, 0.98525f, 0.982f, 0.97875f, 0.9734f, 0.96805f, 0.9627f, 0.95735f, 0.952f, 0.944716f,
    0.937432f, 0.930148f, 0.922864f, 0.91558f, 0.906464f, 0.897348f, 0.888232f, 0.879116f, 0.87f, 0.859246f, 0.848492f,
    0.837738f, 0.826984f, 0.81623f, 0.804384f, 0.792538f, 0.780692f, 0.768846f, 0.757f, 0.744566f, 0.732132f, 0.719698f,
    0.707264f, 0.69483f, 0.682064f, 0.669298f, 0.656532f, 0.643766f, 0.631f, 0.618108f, 0.605216f, 0.592324f, 0.579432f,
    0.56654f, 0.553832f, 0.541124f, 0.528416f, 0.515708f, 0.503f, 0.490744f, 0.478488f, 0.466232f, 0.453976f, 0.44172f,
    0.429576f, 0.417432f, 0.405288f, 0.393144f, 0.381f, 0.368904f, 0.356808f, 0.344712f, 0.332616f, 0.32052f, 0.309416f,
    0.298312f, 0.287208f, 0.276104f, 0.265f, 0.255404f, 0.245808f, 0.236212f, 0.226616f, 0.21702f, 0.208616f, 0.200212f,
    0.191808f, 0.183404f, 0.175f, 0.167624f, 0.160248f, 0.152872f, 0.145496f, 0.13812f, 0.131896f, 0.125672f, 0.119448f,
    0.113224f, 0.107f, 0.1019304f, 0.0968608f, 0.0917912f, 0.0867216f, 0.081652f, 0.0775216f, 0.0733912f, 0.0692608f,
    0.0651304f, 0.061f, 0.0576654f, 0.0543308f, 0.0509962f, 0.0476616f, 0.044327f, 0.0418616f, 0.0393962f, 0.0369308f,
    0.0344654f, 0.032f, 0.0302908f, 0.0285816f, 0.0268724f, 0.0251632f, 0.023454f, 0.0221632f, 0.0208724f, 0.0195816f,
    0.0182908f, 0.017f, 0.0159744f, 0.0149488f, 0.0139232f, 0.0128976f, 0.011872f, 0.0111396f, 0.0104072f, 0.0096748f,
    0.0089424f, 0.00821f, 0.00772246f, 0.00723492f, 0.00674738f, 0.00625984f, 0.0057723f, 0.00543824f, 0.00510418f,
    0.00477012f, 0.00443606f, 0.004102f, 0.00386742f, 0.00363284f, 0.00339826f, 0.00316368f, 0.0029291f, 0.00276148f,
    0.00259386f, 0.00242624f, 0.00225862f, 0.002091f, 0.00196924f, 0.00184748f, 0.00172572f, 0.00160396f, 0.0014822f,
    0.00139516f, 0.00130812f, 0.00122108f, 0.00113404f, 0.001047f, 0.00098563f, 0.00092426f, 0.00086289f, 0.00080152f,
    0.00074015f, 0.00069612f, 0.00065209f, 0.00060806f, 0.00056403f, 0.00052f, 0.000488186f, 0.000456372f, 0.000424558f,
    0.000392744f, 0.00036093f, 0.000338584f, 0.000316238f, 0.000293892f, 0.000271546f, 0.0002492f, 0.000233822f,
    0.000218444f, 0.000203066f, 0.000187688f, 0.00017231f, 0.000161848f, 0.000151386f, 0.000140924f, 0.000130462f,
    0.00012f, 0.000112924f, 0.000105848f, 0.000098772f, 0.000091696f, 0.00008462f, 0.000079696f, 0.000074772f,
    0.000069848f, 0.000064924f, 0.00006f, 5.64892E-05f, 5.29784E-05f, 4.94676E-05f, 4.59568E-05f, 0.000042446f,
    3.99568E-05f, 3.74676E-05f, 3.49784E-05f, 3.24892E-05f, 0.00003f, 0.000028242f, 0.000026484f, 0.000024726f,
    0.000022968f, 0.00002121f, 1.99658E-05f, 1.87216E-05f, 1.74774E-05f, 1.62332E-05f});

static constexpr Spectrum CIE1931_JuddVos_Z({0.01226f, 0.0146524f, 0.0170448f, 0.0194372f, 0.0218296f, 0.024222f,
    0.0292276f, 0.0342332f, 0.0392388f, 0.0442444f, 0.04925f, 0.058427f, 0.067604f, 0.076781f, 0.085958f, 0.095135f,
    0.110926f, 0.126717f, 0.142508f, 0.158299f, 0.17409f, 0.197298f, 0.220506f, 0.243714f, 0.266922f, 0.29013f,
    0.32421f, 0.35829f, 0.39237f, 0.42645f, 0.46053f, 0.514756f, 0.568982f, 0.623208f, 0.677434f, 0.73166f, 0.798488f,
    0.865316f, 0.932144f, 0.998972f, 1.0658f, 1.11556f, 1.16532f, 1.21508f, 1.26484f, 1.3146f, 1.34512f, 1.37564f,
    1.40616f, 1.43668f, 1.4672f, 1.48968f, 1.51216f, 1.53464f, 1.55712f, 1.5796f, 1.587f, 1.5944f, 1.6018f, 1.6092f,
    1.6166f, 1.60692f, 1.59724f, 1.58756f, 1.57788f, 1.5682f, 1.5489f, 1.5296f, 1.5103f, 1.491f, 1.4717f, 1.45216f,
    1.43262f, 1.41308f, 1.39354f, 1.374f, 1.35754f, 1.34108f, 1.32462f, 1.30816f, 1.2917f, 1.28048f, 1.26926f, 1.25804f,
    1.24682f, 1.2356f, 1.21124f, 1.18688f, 1.16252f, 1.13816f, 1.1138f, 1.07948f, 1.04516f, 1.01084f, 0.97652f, 0.9422f,
    0.904952f, 0.867704f, 0.830456f, 0.793208f, 0.75596f, 0.722048f, 0.688136f, 0.654224f, 0.620312f, 0.5864f,
    0.558458f, 0.530516f, 0.502574f, 0.474632f, 0.44669f, 0.425584f, 0.404478f, 0.383372f, 0.362266f, 0.34116f,
    0.325802f, 0.310444f, 0.295086f, 0.279728f, 0.26437f, 0.252684f, 0.240998f, 0.229312f, 0.217626f, 0.20594f,
    0.195642f, 0.185344f, 0.175046f, 0.164748f, 0.15445f, 0.145396f, 0.136342f, 0.127288f, 0.118234f, 0.10918f,
    0.102661f, 0.096142f, 0.089623f, 0.083104f, 0.076585f, 0.0725134f, 0.0684418f, 0.0643702f, 0.0602986f, 0.056227f,
    0.0532548f, 0.0502826f, 0.0473104f, 0.0443382f, 0.041366f, 0.0389634f, 0.0365608f, 0.0341582f, 0.0317556f,
    0.029353f, 0.0274908f, 0.0256286f, 0.0237664f, 0.0219042f, 0.020042f, 0.018696f, 0.01735f, 0.016004f, 0.014658f,
    0.013312f, 0.01240606f, 0.01150012f, 0.01059418f, 0.00968824f, 0.0087823f, 0.0081973f, 0.0076123f, 0.0070273f,
    0.0064423f, 0.0058573f, 0.0054957f, 0.0051341f, 0.0047725f, 0.0044109f, 0.0040493f, 0.00382378f, 0.00359826f,
    0.00337274f, 0.00314722f, 0.0029217f, 0.00279278f, 0.00266386f, 0.00253494f, 0.00240602f, 0.0022771f, 0.0022158f,
    0.0021545f, 0.0020932f, 0.0020319f, 0.0019706f, 0.0019378f, 0.001905f, 0.0018722f, 0.0018394f, 0.0018066f,
    0.00175426f, 0.00170192f, 0.00164958f, 0.00159724f, 0.0015449f, 0.00148288f, 0.00142086f, 0.00135884f, 0.00129682f,
    0.0012348f, 0.00121138f, 0.00118796f, 0.00116454f, 0.00114112f, 0.0011177f, 0.001075288f, 0.001032876f,
    0.000990464f, 0.000948052f, 0.00090564f, 0.000863446f, 0.000821252f, 0.000779058f, 0.000736864f, 0.00069467f,
    0.000641506f, 0.000588342f, 0.000535178f, 0.000482014f, 0.00042885f, 0.000406714f, 0.000384578f, 0.000362442f,
    0.000340306f, 0.00031817f, 0.000305732f, 0.000293294f, 0.000280856f, 0.000268418f, 0.00025598f, 0.000236142f,
    0.000216304f, 0.000196466f, 0.000176628f, 0.00015679f, 0.000144971f, 0.000133152f, 0.000121332f, 0.000109513f,
    0.000097694f, 0.000091944f, 0.000086194f, 0.000080444f, 0.000074694f, 0.000068944f, 6.53882E-05f, 6.18324E-05f,
    5.82766E-05f, 5.47208E-05f, 0.000051165f, 4.81352E-05f, 4.51054E-05f, 4.20756E-05f, 3.90458E-05f, 0.000036016f,
    3.36604E-05f, 3.13048E-05f, 2.89492E-05f, 2.65936E-05f, 0.000024238f, 2.27734E-05f, 2.13088E-05f, 1.98442E-05f,
    1.83796E-05f, 0.000016915f, 1.59132E-05f, 1.49114E-05f, 1.39096E-05f, 1.29078E-05f, 0.000011906f, 1.11546E-05f,
    1.04032E-05f, 9.65174E-06f, 8.90032E-06f, 8.1489E-06f, 7.63924E-06f, 7.12958E-06f, 6.61992E-06f, 6.11026E-06f,
    5.6006E-06f, 5.27136E-06f, 4.94212E-06f, 4.61288E-06f, 4.28364E-06f, 3.9544E-06f, 3.72176E-06f, 3.48912E-06f,
    3.25648E-06f, 3.02384E-06f, 2.7912E-06f, 2.61648E-06f, 2.44176E-06f, 2.26704E-06f, 2.09232E-06f, 1.9176E-06f,
    1.79678E-06f, 1.67596E-06f, 1.55514E-06f, 1.43432E-06f, 1.3135E-06f, 1.23384E-06f, 1.15418E-06f, 1.07451E-06f,
    9.94852E-07f, 9.1519E-07f, 8.61686E-07f, 8.08182E-07f, 7.54678E-07f, 7.01174E-07f, 6.4767E-07f, 6.1084E-07f,
    5.7401E-07f, 5.3718E-07f, 5.0035E-07f, 4.6352E-07f, 4.37424E-07f, 4.11328E-07f, 3.85232E-07f, 3.59136E-07f,
    3.3304E-07f, 3.14078E-07f, 2.95116E-07f, 2.76154E-07f, 2.57192E-07f, 2.3823E-07f, 2.24636E-07f, 2.11042E-07f,
    1.97448E-07f, 1.83854E-07f, 1.7026E-07f, 1.60622E-07f, 1.50984E-07f, 1.41346E-07f, 1.31708E-07f, 1.2207E-07f,
    1.15077E-07f, 1.08085E-07f, 1.01092E-07f, 9.40996E-08f, 8.7107E-08f, 8.19766E-08f, 7.68462E-08f, 7.17158E-08f,
    6.65854E-08f, 6.1455E-08f, 5.77964E-08f, 5.41378E-08f, 5.04792E-08f, 4.68206E-08f, 4.3162E-08f, 4.06054E-08f,
    3.80488E-08f, 3.54922E-08f, 3.29356E-08f, 3.0379E-08f, 2.8614E-08f, 2.6849E-08f, 2.5084E-08f, 2.3319E-08f,
    2.1554E-08f, 2.03418E-08f, 1.91296E-08f, 1.79174E-08f, 1.67052E-08f, 1.5493E-08f, 1.46352E-08f, 1.37774E-08f,
    1.29196E-08f, 1.20618E-08f, 1.1204E-08f, 1.05807E-08f, 9.95732E-09f, 9.33398E-09f, 8.71064E-09f, 8.0873E-09f,
    7.63664E-09f, 7.18598E-09f, 6.73532E-09f, 6.28466E-09f, 5.834E-09f, 5.5094E-09f, 5.1848E-09f, 4.8602E-09f,
    4.5356E-09f, 4.211E-09f, 3.97646E-09f, 3.74192E-09f, 3.50738E-09f, 3.27284E-09f});

// standard illuminant D65( linear interpolated to 1nm )
static constexpr Spectrum CIE_D65({49.97550f, 50.44276f, 50.91002f, 51.37728f, 51.84454f, 52.31180f, 52.77908f,
    53.24636f, 53.71364f, 54.18092f, 54.64820f, 57.45886f, 60.26952f, 63.08018f, 65.89084f, 68.70150f, 71.51218f,
    74.32286f, 77.13354f, 79.94422f, 82.75490f, 83.62800f, 84.50110f, 85.37420f, 86.24730f, 87.12040f, 87.99352f,
    88.86664f, 89.73976f, 90.61288f, 91.48600f, 91.68058f, 91.87516f, 92.06974f, 92.26432f, 92.45890f, 92.65348f,
    92.84806f, 93.04264f, 93.23722f, 93.43180f, 92.75684f, 92.08188f, 91.40692f, 90.73196f, 90.05700f, 89.38206f,
    88.70712f, 88.03218f, 87.35724f, 86.68230f, 88.50056f, 90.31882f, 92.13708f, 93.95534f, 95.77360f, 97.59188f,
    99.41016f, 101.22844f, 103.04672f, 104.86500f, 106.07920f, 107.29340f, 108.50760f, 109.72180f, 110.93600f,
    112.15040f, 113.36480f, 114.57920f, 115.79360f, 117.00800f, 117.08840f, 117.16880f, 117.24920f, 117.32960f,
    117.41000f, 117.49040f, 117.57080f, 117.65120f, 117.73160f, 117.81200f, 117.51680f, 117.22160f, 116.92640f,
    116.63120f, 116.33600f, 116.04100f, 115.74600f, 115.45100f, 115.15600f, 114.86100f, 114.96720f, 115.07340f,
    115.17960f, 115.28580f, 115.39200f, 115.49820f, 115.60440f, 115.71060f, 115.81680f, 115.92300f, 115.21180f,
    114.50060f, 113.78940f, 113.07820f, 112.36700f, 111.65580f, 110.94460f, 110.23340f, 109.52220f, 108.81100f,
    108.86520f, 108.91940f, 108.97360f, 109.02780f, 109.08200f, 109.13640f, 109.19080f, 109.24520f, 109.29960f,
    109.35400f, 109.19880f, 109.04360f, 108.88840f, 108.73320f, 108.57800f, 108.42280f, 108.26760f, 108.11240f,
    107.95720f, 107.80200f, 107.50080f, 107.19960f, 106.89840f, 106.59720f, 106.29600f, 105.99480f, 105.69360f,
    105.39240f, 105.09120f, 104.79000f, 105.07980f, 105.36960f, 105.65940f, 105.94920f, 106.23900f, 106.52900f,
    106.81900f, 107.10900f, 107.39900f, 107.68900f, 107.36060f, 107.03220f, 106.70380f, 106.37540f, 106.04700f,
    105.71860f, 105.39020f, 105.06180f, 104.73340f, 104.40500f, 104.36900f, 104.33300f, 104.29700f, 104.26100f,
    104.22500f, 104.18920f, 104.15340f, 104.11760f, 104.08180f, 104.04600f, 103.64140f, 103.23680f, 102.83220f,
    102.42760f, 102.02300f, 101.61840f, 101.21380f, 100.80920f, 100.40460f, 100.00000f, 99.63342f, 99.26684f, 98.90026f,
    98.53368f, 98.16710f, 97.80052f, 97.43394f, 97.06736f, 96.70078f, 96.33420f, 96.27958f, 96.22496f, 96.17034f,
    96.11572f, 96.06110f, 96.00648f, 95.95186f, 95.89724f, 95.84262f, 95.78800f, 95.07776f, 94.36752f, 93.65728f,
    92.94704f, 92.23680f, 91.52656f, 90.81632f, 90.10608f, 89.39584f, 88.68560f, 88.81766f, 88.94972f, 89.08178f,
    89.21384f, 89.34590f, 89.47796f, 89.61002f, 89.74208f, 89.87414f, 90.00620f, 89.96548f, 89.92476f, 89.88404f,
    89.84332f, 89.80260f, 89.76190f, 89.72120f, 89.68050f, 89.63980f, 89.59910f, 89.40906f, 89.21902f, 89.02898f,
    88.83894f, 88.64890f, 88.45886f, 88.26882f, 88.07878f, 87.88874f, 87.69870f, 87.25768f, 86.81666f, 86.37564f,
    85.93462f, 85.49360f, 85.05260f, 84.61160f, 84.17060f, 83.72960f, 83.28860f, 83.32966f, 83.37072f, 83.41178f,
    83.45284f, 83.49390f, 83.53496f, 83.57602f, 83.61708f, 83.65814f, 83.69920f, 83.33196f, 82.96472f, 82.59748f,
    82.23024f, 81.86300f, 81.49576f, 81.12852f, 80.76128f, 80.39404f, 80.02680f, 80.04558f, 80.06436f, 80.08314f,
    80.10192f, 80.12070f, 80.13948f, 80.15826f, 80.17704f, 80.19582f, 80.21460f, 80.42092f, 80.62724f, 80.83356f,
    81.03988f, 81.24620f, 81.45252f, 81.65884f, 81.86516f, 82.07148f, 82.27780f, 81.87844f, 81.47908f, 81.07972f,
    80.68036f, 80.28100f, 79.88164f, 79.48228f, 79.08292f, 78.68356f, 78.28420f, 77.42790f, 76.57160f, 75.71530f,
    74.85900f, 74.00270f, 73.14642f, 72.29014f, 71.43386f, 70.57758f, 69.72130f, 69.91008f, 70.09886f, 70.28764f,
    70.47642f, 70.66520f, 70.85398f, 71.04276f, 71.23154f, 71.42032f, 71.60910f, 71.88308f, 72.15706f, 72.43104f,
    72.70502f, 72.97900f, 73.25300f, 73.52700f, 73.80100f, 74.07500f, 74.34900f, 73.07450f, 71.80000f, 70.52550f,
    69.25100f, 67.97650f, 66.70200f, 65.42750f, 64.15300f, 62.87850f, 61.60400f, 62.43216f, 63.26032f, 64.08848f,
    64.91664f, 65.74480f, 66.57296f, 67.40112f, 68.22928f, 69.05744f, 69.88560f, 70.40574f, 70.92588f, 71.44602f,
    71.96616f, 72.48630f, 73.00644f, 73.52658f, 74.04672f, 74.56686f, 75.08700f, 73.93756f, 72.78812f, 71.63868f,
    70.48924f, 69.33980f, 68.19038f, 67.04096f, 65.89154f, 64.74212f, 63.59270f, 61.87524f, 60.15778f, 58.44032f,
    56.72286f, 55.00540f, 53.28796f, 51.57052f, 49.85308f, 48.13564f, 46.41820f, 48.45692f, 50.49564f, 52.53436f,
    54.57308f, 56.61180f, 58.65052f, 60.68924f, 62.72796f, 64.76668f, 66.80540f, 66.46314f, 66.12088f, 65.77862f,
    65.43636f, 65.09410f, 64.75184f, 64.40958f, 64.06732f, 63.72506f});

// --- constants.

// Standard observers
static constexpr Observer CIE1931(CIE1931_X, CIE1931_Y, CIE1931_Z);
static constexpr Observer CIE31JV(CIE1931_JuddVos_X, CIE1931_JuddVos_Y, CIE1931_JuddVos_Z);
static constexpr Observer CIE2012(CIE2012_X, CIE2012_Y, CIE2012_Z);

// ------------------ IES TM-30-15 spectrums.

namespace TM_30_15
{
    static constexpr Spectrum CES01({0.619470f, 0.620290f, 0.621110f, 0.621930f, 0.622750f, 0.623570f, 0.624390f,
        0.625210f, 0.626020f, 0.626840f, 0.627660f, 0.628470f, 0.629280f, 0.630100f, 0.630910f, 0.631720f, 0.632530f,
        0.633340f, 0.634150f, 0.634960f, 0.635900f, 0.636500f, 0.637260f, 0.638140f, 0.639110f, 0.640150f, 0.641230f,
        0.642330f, 0.643400f, 0.644440f, 0.645400f, 0.646270f, 0.647040f, 0.647730f, 0.648340f, 0.648880f, 0.649350f,
        0.649760f, 0.650110f, 0.650430f, 0.650700f, 0.650940f, 0.651140f, 0.651300f, 0.651410f, 0.651470f, 0.651470f,
        0.651410f, 0.651280f, 0.651080f, 0.650800f, 0.650440f, 0.650000f, 0.649480f, 0.648890f, 0.648230f, 0.647510f,
        0.646720f, 0.645870f, 0.644960f, 0.644000f, 0.642990f, 0.641910f, 0.640770f, 0.639550f, 0.638250f, 0.636850f,
        0.635340f, 0.633720f, 0.631980f, 0.630100f, 0.628090f, 0.625950f, 0.623710f, 0.621370f, 0.618950f, 0.616480f,
        0.613960f, 0.611410f, 0.608850f, 0.606300f, 0.603760f, 0.601230f, 0.598710f, 0.596170f, 0.593610f, 0.591030f,
        0.588410f, 0.585730f, 0.583000f, 0.580200f, 0.577330f, 0.574390f, 0.571430f, 0.568450f, 0.565480f, 0.562530f,
        0.559640f, 0.556830f, 0.554100f, 0.551500f, 0.549030f, 0.546670f, 0.544420f, 0.542260f, 0.540160f, 0.538100f,
        0.536080f, 0.534070f, 0.532050f, 0.530000f, 0.527920f, 0.525800f, 0.523670f, 0.521540f, 0.519420f, 0.517330f,
        0.515280f, 0.513280f, 0.511350f, 0.509500f, 0.507750f, 0.506090f, 0.504510f, 0.503030f, 0.501630f, 0.500320f,
        0.499080f, 0.497910f, 0.496820f, 0.495800f, 0.494850f, 0.493960f, 0.493140f, 0.492390f, 0.491710f, 0.491100f,
        0.490570f, 0.490100f, 0.489710f, 0.489400f, 0.489160f, 0.488990f, 0.488890f, 0.488840f, 0.488840f, 0.488890f,
        0.488970f, 0.489090f, 0.489240f, 0.489400f, 0.489580f, 0.489790f, 0.490030f, 0.490320f, 0.490680f, 0.491110f,
        0.491620f, 0.492230f, 0.492960f, 0.493800f, 0.494780f, 0.495890f, 0.497150f, 0.498550f, 0.500110f, 0.501820f,
        0.503690f, 0.505720f, 0.507920f, 0.510300f, 0.512850f, 0.515580f, 0.518500f, 0.521600f, 0.524880f, 0.528360f,
        0.532020f, 0.535880f, 0.539940f, 0.544200f, 0.548660f, 0.553310f, 0.558140f, 0.563150f, 0.568320f, 0.573640f,
        0.579110f, 0.584720f, 0.590450f, 0.596300f, 0.602250f, 0.608300f, 0.614410f, 0.620580f, 0.626790f, 0.633020f,
        0.639260f, 0.645470f, 0.651660f, 0.657800f, 0.663880f, 0.669880f, 0.675810f, 0.681660f, 0.687430f, 0.693120f,
        0.698710f, 0.704210f, 0.709600f, 0.714900f, 0.720090f, 0.725170f, 0.730140f, 0.734990f, 0.739740f, 0.744360f,
        0.748880f, 0.753270f, 0.757550f, 0.761700f, 0.765730f, 0.769650f, 0.773450f, 0.777140f, 0.780720f, 0.784210f,
        0.787590f, 0.790880f, 0.794080f, 0.797200f, 0.800230f, 0.803190f, 0.806060f, 0.808850f, 0.811560f, 0.814190f,
        0.816730f, 0.819200f, 0.821590f, 0.823900f, 0.826130f, 0.828280f, 0.830350f, 0.832350f, 0.834270f, 0.836120f,
        0.837900f, 0.839600f, 0.841240f, 0.842800f, 0.844300f, 0.845730f, 0.847110f, 0.848430f, 0.849710f, 0.850950f,
        0.852150f, 0.853320f, 0.854470f, 0.855600f, 0.856710f, 0.857810f, 0.858890f, 0.859950f, 0.860980f, 0.861990f,
        0.862970f, 0.863910f, 0.864830f, 0.865700f, 0.866540f, 0.867340f, 0.868100f, 0.868840f, 0.869550f, 0.870230f,
        0.870900f, 0.871540f, 0.872180f, 0.872800f, 0.873410f, 0.874020f, 0.874620f, 0.875200f, 0.875780f, 0.876350f,
        0.876910f, 0.877450f, 0.877980f, 0.878500f, 0.879010f, 0.879500f, 0.879990f, 0.880460f, 0.880940f, 0.881410f,
        0.881880f, 0.882350f, 0.882820f, 0.883300f, 0.883790f, 0.884280f, 0.884780f, 0.885280f, 0.885790f, 0.886310f,
        0.886830f, 0.887350f, 0.887870f, 0.888400f, 0.888930f, 0.889450f, 0.889960f, 0.890450f, 0.890910f, 0.891350f,
        0.891740f, 0.892080f, 0.892370f, 0.892600f, 0.892960f, 0.893270f, 0.893580f, 0.893890f, 0.894200f, 0.894510f,
        0.894810f, 0.895120f, 0.895420f, 0.895730f, 0.896030f, 0.896330f, 0.896630f, 0.896930f, 0.897230f, 0.897530f,
        0.897830f, 0.898130f, 0.898430f, 0.898720f, 0.899020f, 0.899310f, 0.899610f, 0.899900f, 0.900190f, 0.900480f,
        0.900770f, 0.901060f, 0.901350f, 0.901640f, 0.901930f, 0.902220f, 0.902500f, 0.902790f, 0.903070f, 0.903360f,
        0.903640f, 0.903920f, 0.904200f, 0.904490f, 0.904770f, 0.905050f, 0.905320f, 0.905600f, 0.905880f, 0.906160f,
        0.906430f, 0.906710f, 0.906980f, 0.907260f, 0.907530f, 0.907800f, 0.908070f, 0.908340f, 0.908610f, 0.908880f,
        0.909150f, 0.909420f, 0.909690f, 0.909950f, 0.910220f, 0.910480f, 0.910750f, 0.911010f, 0.911280f, 0.911540f,
        0.911800f, 0.912060f, 0.912320f, 0.912580f, 0.912840f, 0.913100f, 0.913350f, 0.913610f, 0.913870f, 0.914120f,
        0.914380f, 0.914630f, 0.914880f});
    static constexpr Spectrum CES02({0.251730f, 0.252210f, 0.252690f, 0.253170f, 0.253650f, 0.254130f, 0.254610f,
        0.255100f, 0.255580f, 0.256060f, 0.256550f, 0.257030f, 0.257520f, 0.258000f, 0.258490f, 0.258980f, 0.259470f,
        0.259960f, 0.260440f, 0.260930f, 0.261500f, 0.261870f, 0.262340f, 0.262870f, 0.263460f, 0.264090f, 0.264720f,
        0.265350f, 0.265950f, 0.266510f, 0.267000f, 0.267410f, 0.267740f, 0.268010f, 0.268200f, 0.268340f, 0.268430f,
        0.268470f, 0.268480f, 0.268450f, 0.268400f, 0.268330f, 0.268240f, 0.268120f, 0.267970f, 0.267790f, 0.267580f,
        0.267330f, 0.267030f, 0.266690f, 0.266300f, 0.265860f, 0.265360f, 0.264800f, 0.264180f, 0.263490f, 0.262740f,
        0.261920f, 0.261020f, 0.260050f, 0.259000f, 0.257870f, 0.256650f, 0.255340f, 0.253930f, 0.252410f, 0.250790f,
        0.249050f, 0.247200f, 0.245210f, 0.243100f, 0.240850f, 0.238490f, 0.236020f, 0.233470f, 0.230860f, 0.228210f,
        0.225520f, 0.222830f, 0.220150f, 0.217500f, 0.214890f, 0.212320f, 0.209770f, 0.207240f, 0.204710f, 0.202180f,
        0.199630f, 0.197060f, 0.194450f, 0.191800f, 0.189100f, 0.186360f, 0.183600f, 0.180850f, 0.178130f, 0.175450f,
        0.172830f, 0.170310f, 0.167890f, 0.165600f, 0.163450f, 0.161440f, 0.159530f, 0.157730f, 0.156020f, 0.154360f,
        0.152760f, 0.151200f, 0.149650f, 0.148100f, 0.146540f, 0.144970f, 0.143400f, 0.141830f, 0.140260f, 0.138710f,
        0.137170f, 0.135650f, 0.134160f, 0.132700f, 0.131280f, 0.129890f, 0.128560f, 0.127270f, 0.126040f, 0.124870f,
        0.123770f, 0.122740f, 0.121780f, 0.120900f, 0.120100f, 0.119390f, 0.118750f, 0.118180f, 0.117690f, 0.117250f,
        0.116880f, 0.116570f, 0.116310f, 0.116100f, 0.115940f, 0.115810f, 0.115730f, 0.115690f, 0.115690f, 0.115710f,
        0.115770f, 0.115850f, 0.115970f, 0.116100f, 0.116260f, 0.116430f, 0.116630f, 0.116850f, 0.117100f, 0.117370f,
        0.117660f, 0.117980f, 0.118330f, 0.118700f, 0.119100f, 0.119540f, 0.120020f, 0.120550f, 0.121140f, 0.121790f,
        0.122520f, 0.123320f, 0.124210f, 0.125200f, 0.126290f, 0.127490f, 0.128830f, 0.130320f, 0.131980f, 0.133810f,
        0.135850f, 0.138100f, 0.140570f, 0.143300f, 0.146290f, 0.149540f, 0.153060f, 0.156850f, 0.160920f, 0.165270f,
        0.169890f, 0.174810f, 0.180010f, 0.185500f, 0.191280f, 0.197350f, 0.203670f, 0.210250f, 0.217070f, 0.224110f,
        0.231360f, 0.238800f, 0.246420f, 0.254200f, 0.262130f, 0.270200f, 0.278390f, 0.286680f, 0.295050f, 0.303490f,
        0.311980f, 0.320510f, 0.329050f, 0.337600f, 0.346130f, 0.354640f, 0.363120f, 0.371540f, 0.379910f, 0.388210f,
        0.396430f, 0.404560f, 0.412580f, 0.420500f, 0.428290f, 0.435970f, 0.443520f, 0.450960f, 0.458280f, 0.465480f,
        0.472570f, 0.479560f, 0.486430f, 0.493200f, 0.499860f, 0.506410f, 0.512840f, 0.519150f, 0.525330f, 0.531360f,
        0.537260f, 0.543000f, 0.548580f, 0.554000f, 0.559250f, 0.564340f, 0.569270f, 0.574070f, 0.578720f, 0.583250f,
        0.587670f, 0.591970f, 0.596180f, 0.600300f, 0.604340f, 0.608290f, 0.612160f, 0.615940f, 0.619630f, 0.623240f,
        0.626750f, 0.630160f, 0.633480f, 0.636700f, 0.639820f, 0.642850f, 0.645800f, 0.648680f, 0.651510f, 0.654280f,
        0.657020f, 0.659720f, 0.662420f, 0.665100f, 0.667790f, 0.670470f, 0.673160f, 0.675830f, 0.678500f, 0.681160f,
        0.683800f, 0.686420f, 0.689020f, 0.691600f, 0.694150f, 0.696680f, 0.699200f, 0.701710f, 0.704220f, 0.706750f,
        0.709290f, 0.711860f, 0.714460f, 0.717100f, 0.719780f, 0.722490f, 0.725200f, 0.727880f, 0.730500f, 0.733060f,
        0.735510f, 0.737830f, 0.740000f, 0.742000f, 0.743810f, 0.745450f, 0.746960f, 0.748390f, 0.749760f, 0.751110f,
        0.752490f, 0.753920f, 0.755440f, 0.757100f, 0.758910f, 0.760840f, 0.762840f, 0.764850f, 0.766840f, 0.768750f,
        0.770530f, 0.772130f, 0.773500f, 0.774600f, 0.776280f, 0.777720f, 0.779160f, 0.780590f, 0.782010f, 0.783430f,
        0.784840f, 0.786240f, 0.787640f, 0.789030f, 0.790410f, 0.791780f, 0.793150f, 0.794520f, 0.795870f, 0.797220f,
        0.798570f, 0.799900f, 0.801230f, 0.802560f, 0.803870f, 0.805180f, 0.806480f, 0.807780f, 0.809070f, 0.810350f,
        0.811630f, 0.812900f, 0.814160f, 0.815420f, 0.816670f, 0.817910f, 0.819150f, 0.820380f, 0.821610f, 0.822820f,
        0.824030f, 0.825240f, 0.826440f, 0.827630f, 0.828810f, 0.829990f, 0.831160f, 0.832330f, 0.833490f, 0.834640f,
        0.835790f, 0.836930f, 0.838060f, 0.839190f, 0.840310f, 0.841420f, 0.842530f, 0.843630f, 0.844730f, 0.845820f,
        0.846900f, 0.847980f, 0.849050f, 0.850110f, 0.851170f, 0.852220f, 0.853270f, 0.854310f, 0.855340f, 0.856370f,
        0.857390f, 0.858400f, 0.859410f, 0.860420f, 0.861410f, 0.862400f, 0.863390f, 0.864370f, 0.865340f, 0.866310f,
        0.867270f, 0.868230f, 0.869180f});
    static constexpr Spectrum CES03({0.000008f, 0.000641f, 0.001980f, 0.003848f, 0.005605f, 0.006630f, 0.008000f,
        0.009467f, 0.011731f, 0.014499f, 0.016937f, 0.019039f, 0.020151f, 0.019922f, 0.018672f, 0.017729f, 0.017524f,
        0.016769f, 0.015711f, 0.013973f, 0.011674f, 0.009649f, 0.008111f, 0.006818f, 0.005726f, 0.005109f, 0.004350f,
        0.003574f, 0.003147f, 0.003359f, 0.003785f, 0.004424f, 0.005293f, 0.006272f, 0.007368f, 0.008737f, 0.010293f,
        0.012035f, 0.014144f, 0.016532f, 0.018883f, 0.021246f, 0.023596f, 0.025631f, 0.027417f, 0.029720f, 0.032465f,
        0.035204f, 0.037738f, 0.040150f, 0.041988f, 0.043189f, 0.044156f, 0.045524f, 0.047211f, 0.048964f, 0.050410f,
        0.051750f, 0.052855f, 0.053792f, 0.054566f, 0.055445f, 0.056440f, 0.057083f, 0.057070f, 0.056742f, 0.056205f,
        0.055547f, 0.054937f, 0.054413f, 0.054100f, 0.053902f, 0.053774f, 0.053765f, 0.054203f, 0.054871f, 0.055664f,
        0.056679f, 0.057649f, 0.058413f, 0.059081f, 0.059649f, 0.060026f, 0.060341f, 0.060549f, 0.060752f, 0.061032f,
        0.061311f, 0.061645f, 0.062179f, 0.062592f, 0.062958f, 0.063186f, 0.063270f, 0.063130f, 0.063044f, 0.062775f,
        0.062329f, 0.061926f, 0.061441f, 0.060687f, 0.060151f, 0.059726f, 0.059534f, 0.059488f, 0.059669f, 0.059799f,
        0.059946f, 0.059870f, 0.059691f, 0.059439f, 0.058867f, 0.058073f, 0.057351f, 0.056628f, 0.056008f, 0.055706f,
        0.055466f, 0.055030f, 0.054521f, 0.053919f, 0.053180f, 0.052790f, 0.052549f, 0.052383f, 0.052021f, 0.051481f,
        0.050688f, 0.049959f, 0.049342f, 0.048948f, 0.048681f, 0.048467f, 0.048208f, 0.047985f, 0.047624f, 0.047382f,
        0.047261f, 0.047091f, 0.046776f, 0.046471f, 0.046118f, 0.045696f, 0.045187f, 0.044646f, 0.043976f, 0.043331f,
        0.042642f, 0.041988f, 0.041440f, 0.041024f, 0.040629f, 0.040410f, 0.040326f, 0.040347f, 0.040464f, 0.040826f,
        0.041158f, 0.041509f, 0.041707f, 0.041868f, 0.041924f, 0.041975f, 0.041805f, 0.041803f, 0.041972f, 0.042084f,
        0.041975f, 0.041883f, 0.041602f, 0.041234f, 0.040886f, 0.040741f, 0.040890f, 0.041396f, 0.041849f, 0.042197f,
        0.042496f, 0.042671f, 0.042718f, 0.043051f, 0.043563f, 0.044035f, 0.044529f, 0.044840f, 0.044837f, 0.044615f,
        0.044444f, 0.044162f, 0.044104f, 0.044294f, 0.044551f, 0.044879f, 0.045463f, 0.045961f, 0.046315f, 0.046699f,
        0.047125f, 0.047595f, 0.048351f, 0.049141f, 0.049935f, 0.050704f, 0.051166f, 0.051322f, 0.051510f, 0.051805f,
        0.052008f, 0.052446f, 0.052802f, 0.053052f, 0.053202f, 0.053442f, 0.053637f, 0.054138f, 0.054632f, 0.055210f,
        0.055836f, 0.056751f, 0.057576f, 0.058628f, 0.059789f, 0.061282f, 0.062793f, 0.064374f, 0.066035f, 0.067686f,
        0.069226f, 0.070720f, 0.072334f, 0.074268f, 0.076501f, 0.078768f, 0.080920f, 0.083096f, 0.085163f, 0.087222f,
        0.089622f, 0.092294f, 0.095102f, 0.097902f, 0.100340f, 0.102350f, 0.104610f, 0.107020f, 0.109280f, 0.111990f,
        0.114760f, 0.117120f, 0.118960f, 0.120640f, 0.122030f, 0.123520f, 0.125020f, 0.126390f, 0.127490f, 0.128650f,
        0.129480f, 0.130400f, 0.131630f, 0.133220f, 0.134750f, 0.136320f, 0.137750f, 0.139260f, 0.140940f, 0.142940f,
        0.145140f, 0.147240f, 0.148850f, 0.150030f, 0.151240f, 0.152320f, 0.153510f, 0.155000f, 0.156560f, 0.157590f,
        0.158540f, 0.159530f, 0.160630f, 0.162070f, 0.163850f, 0.165440f, 0.167210f, 0.169080f, 0.170950f, 0.172770f,
        0.175020f, 0.176930f, 0.178480f, 0.180040f, 0.181510f, 0.182910f, 0.184390f, 0.185700f, 0.186640f, 0.187150f,
        0.187270f, 0.186890f, 0.186660f, 0.186110f, 0.185490f, 0.185230f, 0.185170f, 0.185460f, 0.186500f, 0.188170f,
        0.189730f, 0.191670f, 0.193570f, 0.195290f, 0.197300f, 0.199360f, 0.200970f, 0.202330f, 0.202100f, 0.202700f,
        0.203930f, 0.203930f, 0.203930f, 0.203930f, 0.208280f, 0.209520f, 0.210770f, 0.212020f, 0.213280f, 0.214550f,
        0.215820f, 0.217090f, 0.218370f, 0.219660f, 0.220950f, 0.222250f, 0.223550f, 0.224850f, 0.226170f, 0.227490f,
        0.228810f, 0.230140f, 0.231470f, 0.232810f, 0.234160f, 0.235510f, 0.236870f, 0.238230f, 0.239590f, 0.240970f,
        0.242340f, 0.243730f, 0.245110f, 0.246510f, 0.247910f, 0.249310f, 0.250720f, 0.252140f, 0.253560f, 0.254980f,
        0.256410f, 0.257850f, 0.259290f, 0.260730f, 0.262190f, 0.263640f, 0.265100f, 0.266570f, 0.268040f, 0.269520f,
        0.271000f, 0.272490f, 0.273980f, 0.275480f, 0.276980f, 0.278490f, 0.280010f, 0.281520f, 0.283050f, 0.284570f,
        0.286110f, 0.287640f, 0.289190f, 0.290740f, 0.292290f, 0.293840f, 0.295410f, 0.296970f, 0.298550f, 0.300120f,
        0.301700f, 0.303290f, 0.304880f, 0.306480f, 0.308080f, 0.309680f, 0.311290f, 0.312900f, 0.314520f, 0.316140f,
        0.317770f, 0.319400f, 0.321040f});
    static constexpr Spectrum CES04({0.473760f, 0.477183f, 0.478810f, 0.478934f, 0.478504f, 0.477931f, 0.476276f,
        0.473907f, 0.471000f, 0.467751f, 0.464059f, 0.459961f, 0.455314f, 0.450429f, 0.445506f, 0.441016f, 0.436857f,
        0.432851f, 0.428936f, 0.425136f, 0.421566f, 0.418564f, 0.415470f, 0.412248f, 0.408937f, 0.405590f, 0.402245f,
        0.398946f, 0.395730f, 0.392583f, 0.389496f, 0.386511f, 0.383638f, 0.380829f, 0.378146f, 0.375665f, 0.373356f,
        0.371153f, 0.369062f, 0.367027f, 0.365006f, 0.363034f, 0.361127f, 0.359287f, 0.357609f, 0.356108f, 0.354693f,
        0.353333f, 0.352038f, 0.350730f, 0.349358f, 0.347942f, 0.346500f, 0.345023f, 0.343520f, 0.342000f, 0.340440f,
        0.338850f, 0.337233f, 0.335564f, 0.333824f, 0.332017f, 0.330128f, 0.328158f, 0.326140f, 0.324106f, 0.322044f,
        0.319922f, 0.317701f, 0.315345f, 0.312822f, 0.310134f, 0.307317f, 0.304392f, 0.301351f, 0.298162f, 0.294796f,
        0.291244f, 0.287533f, 0.283708f, 0.279821f, 0.275935f, 0.272075f, 0.268214f, 0.264330f, 0.260421f, 0.256480f,
        0.252518f, 0.248569f, 0.244652f, 0.240781f, 0.236962f, 0.233183f, 0.229433f, 0.225723f, 0.222057f, 0.218439f,
        0.214887f, 0.211401f, 0.207976f, 0.204617f, 0.201324f, 0.198093f, 0.194935f, 0.191851f, 0.188837f, 0.185896f,
        0.183029f, 0.180233f, 0.177517f, 0.174885f, 0.172343f, 0.169907f, 0.167590f, 0.165405f, 0.163365f, 0.161482f,
        0.159759f, 0.158203f, 0.156816f, 0.155600f, 0.154554f, 0.153674f, 0.152946f, 0.152357f, 0.151895f, 0.151556f,
        0.151335f, 0.151241f, 0.151273f, 0.151434f, 0.151720f, 0.152130f, 0.152662f, 0.153323f, 0.154120f, 0.155063f,
        0.156156f, 0.157403f, 0.158803f, 0.160363f, 0.162086f, 0.163980f, 0.166052f, 0.168308f, 0.170750f, 0.173384f,
        0.176209f, 0.179227f, 0.182445f, 0.185858f, 0.189457f, 0.193242f, 0.197200f, 0.201315f, 0.205585f, 0.210003f,
        0.214554f, 0.219231f, 0.224024f, 0.228918f, 0.233907f, 0.238994f, 0.244167f, 0.249425f, 0.254772f, 0.260204f,
        0.265718f, 0.271317f, 0.276989f, 0.282721f, 0.288507f, 0.294333f, 0.300190f, 0.306085f, 0.312023f, 0.317996f,
        0.323994f, 0.330007f, 0.336013f, 0.342003f, 0.347965f, 0.353876f, 0.359718f, 0.365485f, 0.371163f, 0.376734f,
        0.382193f, 0.387538f, 0.392750f, 0.397823f, 0.402754f, 0.407546f, 0.412215f, 0.416779f, 0.421226f, 0.425540f,
        0.429721f, 0.433757f, 0.437633f, 0.441361f, 0.444956f, 0.448415f, 0.451736f, 0.454918f, 0.457955f, 0.460847f,
        0.463605f, 0.466222f, 0.468701f, 0.471052f, 0.473280f, 0.475385f, 0.477382f, 0.479278f, 0.481066f, 0.482744f,
        0.484313f, 0.485775f, 0.487139f, 0.488420f, 0.489630f, 0.490766f, 0.491845f, 0.492880f, 0.493864f, 0.494797f,
        0.495698f, 0.496569f, 0.497412f, 0.498244f, 0.499075f, 0.499901f, 0.500730f, 0.501532f, 0.502290f, 0.503011f,
        0.503712f, 0.504383f, 0.505033f, 0.505666f, 0.506266f, 0.506813f, 0.507310f, 0.507764f, 0.508194f, 0.508609f,
        0.509017f, 0.509406f, 0.509769f, 0.510089f, 0.510375f, 0.510639f, 0.510897f, 0.511155f, 0.511424f, 0.511687f,
        0.511936f, 0.512168f, 0.512390f, 0.512617f, 0.512870f, 0.513141f, 0.513432f, 0.513730f, 0.514023f, 0.514293f,
        0.514551f, 0.514792f, 0.515017f, 0.515219f, 0.515394f, 0.515528f, 0.515617f, 0.515658f, 0.515662f, 0.515634f,
        0.515581f, 0.515502f, 0.515400f, 0.515266f, 0.515089f, 0.514861f, 0.514610f, 0.514354f, 0.514104f, 0.513878f,
        0.513688f, 0.513511f, 0.513351f, 0.513222f, 0.513151f, 0.513144f, 0.513201f, 0.513300f, 0.513426f, 0.513472f,
        0.513449f, 0.513422f, 0.513398f, 0.513368f, 0.513380f, 0.513358f, 0.513302f, 0.513217f, 0.513114f, 0.513005f,
        0.512908f, 0.512815f, 0.512719f, 0.512615f, 0.512505f, 0.512387f, 0.512273f, 0.512171f, 0.512082f, 0.512008f,
        0.511955f, 0.511919f, 0.511902f, 0.511913f, 0.511965f, 0.512084f, 0.512288f, 0.512568f, 0.512910f, 0.513293f,
        0.513685f, 0.514067f, 0.514454f, 0.514850f, 0.515251f, 0.515648f, 0.516024f, 0.516361f, 0.516666f, 0.516938f,
        0.517270f, 0.517667f, 0.518073f, 0.518488f, 0.518912f, 0.519287f, 0.519661f, 0.520036f, 0.520410f, 0.520785f,
        0.521159f, 0.521534f, 0.521908f, 0.522282f, 0.522657f, 0.523031f, 0.523405f, 0.523780f, 0.524154f, 0.524528f,
        0.524902f, 0.525276f, 0.525650f, 0.526024f, 0.526398f, 0.526772f, 0.527146f, 0.527520f, 0.527894f, 0.528268f,
        0.528642f, 0.529016f, 0.529390f, 0.529763f, 0.530137f, 0.530511f, 0.530884f, 0.531258f, 0.531632f, 0.532005f,
        0.532379f, 0.532752f, 0.533126f, 0.533499f, 0.533873f, 0.534246f, 0.534620f, 0.534993f, 0.535367f, 0.535740f,
        0.536113f, 0.536487f, 0.536860f, 0.537233f, 0.537606f, 0.537979f, 0.538352f, 0.538725f, 0.539098f, 0.539471f,
        0.539844f, 0.540216f, 0.540590f});
    static constexpr Spectrum CES05({0.133510f, 0.134070f, 0.134530f, 0.134900f, 0.135190f, 0.135400f, 0.135540f,
        0.135630f, 0.135670f, 0.135680f, 0.135650f, 0.135600f, 0.135550f, 0.135490f, 0.135440f, 0.135400f, 0.135390f,
        0.135420f, 0.135490f, 0.135610f, 0.135790f, 0.136040f, 0.136350f, 0.136720f, 0.137120f, 0.137570f, 0.138040f,
        0.138530f, 0.139020f, 0.139520f, 0.140010f, 0.140480f, 0.140920f, 0.141320f, 0.141670f, 0.141950f, 0.142150f,
        0.142260f, 0.142270f, 0.142170f, 0.141940f, 0.141570f, 0.141080f, 0.140460f, 0.139720f, 0.138870f, 0.137930f,
        0.136890f, 0.135760f, 0.134550f, 0.133270f, 0.131920f, 0.130520f, 0.129050f, 0.127530f, 0.125950f, 0.124330f,
        0.122660f, 0.120940f, 0.119180f, 0.117380f, 0.115550f, 0.113680f, 0.111790f, 0.109870f, 0.107940f, 0.106000f,
        0.104040f, 0.102090f, 0.100140f, 0.098190f, 0.096254f, 0.094332f, 0.092427f, 0.090541f, 0.088677f, 0.086838f,
        0.085025f, 0.083243f, 0.081494f, 0.079780f, 0.078104f, 0.076468f, 0.074874f, 0.073325f, 0.071823f, 0.070368f,
        0.068965f, 0.067614f, 0.066319f, 0.065080f, 0.063900f, 0.062779f, 0.061715f, 0.060710f, 0.059762f, 0.058870f,
        0.058034f, 0.057254f, 0.056530f, 0.055860f, 0.055244f, 0.054677f, 0.054156f, 0.053675f, 0.053231f, 0.052819f,
        0.052434f, 0.052073f, 0.051729f, 0.051400f, 0.051081f, 0.050772f, 0.050473f, 0.050185f, 0.049907f, 0.049640f,
        0.049385f, 0.049141f, 0.048910f, 0.048690f, 0.048483f, 0.048292f, 0.048118f, 0.047964f, 0.047834f, 0.047728f,
        0.047650f, 0.047603f, 0.047589f, 0.047610f, 0.047669f, 0.047768f, 0.047906f, 0.048087f, 0.048310f, 0.048577f,
        0.048890f, 0.049249f, 0.049655f, 0.050110f, 0.050614f, 0.051162f, 0.051750f, 0.052373f, 0.053024f, 0.053700f,
        0.054395f, 0.055103f, 0.055820f, 0.056540f, 0.057259f, 0.057974f, 0.058685f, 0.059387f, 0.060081f, 0.060764f,
        0.061434f, 0.062090f, 0.062729f, 0.063350f, 0.063952f, 0.064541f, 0.065121f, 0.065699f, 0.066281f, 0.066873f,
        0.067481f, 0.068111f, 0.068769f, 0.069460f, 0.070192f, 0.070976f, 0.071825f, 0.072750f, 0.073763f, 0.074878f,
        0.076105f, 0.077458f, 0.078949f, 0.080590f, 0.082390f, 0.084350f, 0.086465f, 0.088734f, 0.091152f, 0.093718f,
        0.096428f, 0.099278f, 0.102270f, 0.105390f, 0.108640f, 0.112000f, 0.115440f, 0.118940f, 0.122470f, 0.126010f,
        0.129540f, 0.133020f, 0.136430f, 0.139760f, 0.142970f, 0.146080f, 0.149060f, 0.151940f, 0.154710f, 0.157370f,
        0.159930f, 0.162380f, 0.164730f, 0.166980f, 0.169140f, 0.171240f, 0.173330f, 0.175440f, 0.177630f, 0.179940f,
        0.182410f, 0.185070f, 0.187980f, 0.191180f, 0.194710f, 0.198590f, 0.202850f, 0.207520f, 0.212630f, 0.218190f,
        0.224240f, 0.230810f, 0.237910f, 0.245570f, 0.253810f, 0.262630f, 0.271980f, 0.281850f, 0.292210f, 0.303050f,
        0.314330f, 0.326040f, 0.338150f, 0.350630f, 0.363460f, 0.376590f, 0.389970f, 0.403540f, 0.417260f, 0.431060f,
        0.444900f, 0.458730f, 0.472490f, 0.486130f, 0.499600f, 0.512880f, 0.525960f, 0.538800f, 0.551400f, 0.563740f,
        0.575790f, 0.587540f, 0.598970f, 0.610070f, 0.620820f, 0.631220f, 0.641290f, 0.651050f, 0.660510f, 0.669670f,
        0.678560f, 0.687190f, 0.695560f, 0.703700f, 0.711610f, 0.719300f, 0.726770f, 0.734020f, 0.741050f, 0.747870f,
        0.754470f, 0.760870f, 0.767050f, 0.773030f, 0.778800f, 0.784370f, 0.789750f, 0.794940f, 0.799950f, 0.804790f,
        0.809450f, 0.813950f, 0.818290f, 0.822480f, 0.826520f, 0.830420f, 0.834190f, 0.837820f, 0.841330f, 0.844710f,
        0.847980f, 0.851130f, 0.854180f, 0.857120f, 0.859970f, 0.862720f, 0.865390f, 0.867980f, 0.870500f, 0.872950f,
        0.875340f, 0.877680f, 0.879970f, 0.882220f, 0.884440f, 0.886620f, 0.888760f, 0.890860f, 0.892920f, 0.894940f,
        0.896920f, 0.898850f, 0.900730f, 0.902560f, 0.904340f, 0.906070f, 0.907740f, 0.909360f, 0.910910f, 0.912410f,
        0.913850f, 0.915220f, 0.916520f, 0.917760f, 0.918930f, 0.920030f, 0.921080f, 0.922080f, 0.923040f, 0.923950f,
        0.924840f, 0.925700f, 0.926540f, 0.927370f, 0.928190f, 0.929020f, 0.929850f, 0.930690f, 0.931560f, 0.932450f,
        0.933370f, 0.934330f, 0.935330f, 0.936390f, 0.936740f, 0.937530f, 0.938310f, 0.939080f, 0.939850f, 0.940600f,
        0.941350f, 0.942080f, 0.942810f, 0.943530f, 0.944240f, 0.944940f, 0.945640f, 0.946330f, 0.947000f, 0.947670f,
        0.948330f, 0.948990f, 0.949630f, 0.950270f, 0.950900f, 0.951530f, 0.952140f, 0.952750f, 0.953350f, 0.953940f,
        0.954530f, 0.955110f, 0.955680f, 0.956250f, 0.956810f, 0.957360f, 0.957900f, 0.958440f, 0.958970f, 0.959500f,
        0.960010f, 0.960530f, 0.961030f, 0.961530f, 0.962030f, 0.962510f, 0.963000f, 0.963470f, 0.963940f, 0.964400f,
        0.964860f, 0.965310f, 0.965760f});
    static constexpr Spectrum CES06({0.045250f, 0.046774f, 0.048347f, 0.049969f, 0.051643f, 0.053370f, 0.055151f,
        0.056989f, 0.058883f, 0.060837f, 0.062851f, 0.064927f, 0.067066f, 0.069271f, 0.071543f, 0.073884f, 0.076295f,
        0.078778f, 0.081334f, 0.083966f, 0.082492f, 0.087222f, 0.091715f, 0.095975f, 0.100010f, 0.103820f, 0.107400f,
        0.110780f, 0.113940f, 0.116900f, 0.119650f, 0.122200f, 0.124560f, 0.126740f, 0.128720f, 0.130530f, 0.132160f,
        0.133610f, 0.134900f, 0.136030f, 0.136990f, 0.137810f, 0.138480f, 0.139010f, 0.139430f, 0.139750f, 0.139970f,
        0.140110f, 0.140170f, 0.140180f, 0.140150f, 0.140080f, 0.139970f, 0.139840f, 0.139660f, 0.139460f, 0.139210f,
        0.138940f, 0.138620f, 0.138270f, 0.137880f, 0.137460f, 0.136990f, 0.136500f, 0.135970f, 0.135410f, 0.134820f,
        0.134190f, 0.133540f, 0.132870f, 0.132170f, 0.131440f, 0.130690f, 0.129930f, 0.129140f, 0.128350f, 0.127540f,
        0.126730f, 0.125920f, 0.125100f, 0.124280f, 0.123470f, 0.122650f, 0.121840f, 0.121030f, 0.120210f, 0.119390f,
        0.118560f, 0.117720f, 0.116870f, 0.116000f, 0.115120f, 0.114230f, 0.113340f, 0.112450f, 0.111570f, 0.110690f,
        0.109840f, 0.109000f, 0.108200f, 0.107430f, 0.106690f, 0.106000f, 0.105340f, 0.104720f, 0.104130f, 0.103580f,
        0.103070f, 0.102580f, 0.102130f, 0.101710f, 0.101320f, 0.100960f, 0.100620f, 0.100290f, 0.099987f, 0.099693f,
        0.099406f, 0.099124f, 0.098842f, 0.098557f, 0.098267f, 0.097973f, 0.097681f, 0.097393f, 0.097115f, 0.096849f,
        0.096600f, 0.096372f, 0.096169f, 0.095995f, 0.095852f, 0.095744f, 0.095670f, 0.095631f, 0.095629f, 0.095663f,
        0.095736f, 0.095848f, 0.095999f, 0.096192f, 0.096425f, 0.096700f, 0.097016f, 0.097372f, 0.097768f, 0.098203f,
        0.098678f, 0.099191f, 0.099742f, 0.100330f, 0.100960f, 0.101620f, 0.102330f, 0.103080f, 0.103880f, 0.104720f,
        0.105620f, 0.106560f, 0.107560f, 0.108610f, 0.109720f, 0.110890f, 0.112110f, 0.113390f, 0.114720f, 0.116110f,
        0.117560f, 0.119050f, 0.120610f, 0.122210f, 0.123870f, 0.125580f, 0.127350f, 0.129160f, 0.131040f, 0.132960f,
        0.134950f, 0.136990f, 0.139080f, 0.141230f, 0.143440f, 0.145700f, 0.148010f, 0.150360f, 0.152750f, 0.155180f,
        0.157640f, 0.160130f, 0.162650f, 0.165180f, 0.167730f, 0.170290f, 0.172850f, 0.175410f, 0.177960f, 0.180490f,
        0.183000f, 0.185470f, 0.187920f, 0.190310f, 0.192660f, 0.194960f, 0.197220f, 0.199420f, 0.201570f, 0.203680f,
        0.205730f, 0.207740f, 0.209690f, 0.211600f, 0.213460f, 0.215270f, 0.217010f, 0.218690f, 0.220300f, 0.221830f,
        0.223280f, 0.224640f, 0.225910f, 0.227080f, 0.228140f, 0.229100f, 0.229960f, 0.230720f, 0.231390f, 0.231950f,
        0.232420f, 0.232800f, 0.233090f, 0.233280f, 0.233390f, 0.233430f, 0.233400f, 0.233310f, 0.233180f, 0.233020f,
        0.232840f, 0.232650f, 0.232470f, 0.232300f, 0.232160f, 0.232060f, 0.232010f, 0.232040f, 0.232150f, 0.232370f,
        0.232700f, 0.233170f, 0.233790f, 0.234570f, 0.235520f, 0.236670f, 0.238010f, 0.239560f, 0.241320f, 0.243310f,
        0.245530f, 0.247990f, 0.250710f, 0.253690f, 0.256930f, 0.260460f, 0.264290f, 0.268420f, 0.272870f, 0.277640f,
        0.282770f, 0.288240f, 0.294080f, 0.300300f, 0.306910f, 0.313890f, 0.321240f, 0.328920f, 0.336950f, 0.345280f,
        0.353930f, 0.362850f, 0.372060f, 0.381510f, 0.391210f, 0.401110f, 0.411150f, 0.421290f, 0.431490f, 0.441700f,
        0.451870f, 0.461940f, 0.471890f, 0.481650f, 0.491190f, 0.500510f, 0.509610f, 0.518510f, 0.527200f, 0.535700f,
        0.544010f, 0.552140f, 0.560100f, 0.567890f, 0.575520f, 0.582980f, 0.590270f, 0.597380f, 0.604300f, 0.611030f,
        0.617560f, 0.623890f, 0.630000f, 0.635890f, 0.641560f, 0.646990f, 0.652180f, 0.657120f, 0.661810f, 0.666240f,
        0.670410f, 0.674300f, 0.677910f, 0.681230f, 0.689250f, 0.693970f, 0.698640f, 0.703270f, 0.707860f, 0.712410f,
        0.716920f, 0.721380f, 0.725800f, 0.730170f, 0.734510f, 0.738790f, 0.743040f, 0.747230f, 0.751380f, 0.755490f,
        0.759550f, 0.763560f, 0.767530f, 0.771450f, 0.775320f, 0.779140f, 0.782920f, 0.786660f, 0.790340f, 0.793980f,
        0.797570f, 0.801110f, 0.804610f, 0.808060f, 0.811470f, 0.814820f, 0.818130f, 0.821400f, 0.824620f, 0.827790f,
        0.830920f, 0.834000f, 0.837040f, 0.840030f, 0.842980f, 0.845880f, 0.848740f, 0.851550f, 0.854330f, 0.857050f,
        0.859740f, 0.862380f, 0.864980f, 0.867540f, 0.870060f, 0.872540f, 0.874980f, 0.877380f, 0.879730f, 0.882050f,
        0.884330f, 0.886570f, 0.888780f, 0.890940f, 0.893070f, 0.895160f, 0.897220f, 0.899240f, 0.901220f, 0.903170f,
        0.905090f, 0.906970f, 0.908820f, 0.910630f, 0.912410f, 0.914160f, 0.915880f, 0.917570f, 0.919230f, 0.920850f,
        0.922450f, 0.924010f, 0.925550f});
    static constexpr Spectrum CES07({0.060013f, 0.059927f, 0.059841f, 0.059755f, 0.059670f, 0.059584f, 0.059499f,
        0.059413f, 0.059328f, 0.059243f, 0.059158f, 0.059073f, 0.058988f, 0.058903f, 0.058819f, 0.058734f, 0.058650f,
        0.058566f, 0.058482f, 0.058398f, 0.058300f, 0.058238f, 0.058159f, 0.058068f, 0.057967f, 0.057858f, 0.057745f,
        0.057630f, 0.057515f, 0.057404f, 0.057300f, 0.057204f, 0.057118f, 0.057042f, 0.056976f, 0.056919f, 0.056874f,
        0.056839f, 0.056814f, 0.056802f, 0.056800f, 0.056810f, 0.056830f, 0.056859f, 0.056896f, 0.056939f, 0.056987f,
        0.057038f, 0.057092f, 0.057146f, 0.057200f, 0.057252f, 0.057302f, 0.057350f, 0.057395f, 0.057437f, 0.057477f,
        0.057513f, 0.057545f, 0.057575f, 0.057600f, 0.057621f, 0.057636f, 0.057643f, 0.057640f, 0.057625f, 0.057596f,
        0.057551f, 0.057488f, 0.057405f, 0.057300f, 0.057172f, 0.057022f, 0.056852f, 0.056665f, 0.056463f, 0.056248f,
        0.056022f, 0.055787f, 0.055546f, 0.055300f, 0.055052f, 0.054803f, 0.054554f, 0.054306f, 0.054060f, 0.053817f,
        0.053579f, 0.053346f, 0.053119f, 0.052900f, 0.052689f, 0.052487f, 0.052292f, 0.052104f, 0.051923f, 0.051748f,
        0.051579f, 0.051415f, 0.051255f, 0.051100f, 0.050949f, 0.050801f, 0.050658f, 0.050519f, 0.050385f, 0.050257f,
        0.050134f, 0.050016f, 0.049905f, 0.049800f, 0.049702f, 0.049609f, 0.049522f, 0.049440f, 0.049361f, 0.049285f,
        0.049212f, 0.049141f, 0.049070f, 0.049000f, 0.048930f, 0.048860f, 0.048791f, 0.048724f, 0.048659f, 0.048598f,
        0.048541f, 0.048488f, 0.048441f, 0.048400f, 0.048366f, 0.048337f, 0.048314f, 0.048295f, 0.048278f, 0.048263f,
        0.048249f, 0.048234f, 0.048218f, 0.048200f, 0.048179f, 0.048154f, 0.048127f, 0.048098f, 0.048067f, 0.048035f,
        0.048002f, 0.047968f, 0.047934f, 0.047900f, 0.047867f, 0.047835f, 0.047805f, 0.047777f, 0.047753f, 0.047732f,
        0.047716f, 0.047705f, 0.047699f, 0.047700f, 0.047707f, 0.047721f, 0.047741f, 0.047766f, 0.047796f, 0.047830f,
        0.047868f, 0.047910f, 0.047954f, 0.048000f, 0.048048f, 0.048100f, 0.048156f, 0.048218f, 0.048288f, 0.048366f,
        0.048455f, 0.048556f, 0.048671f, 0.048800f, 0.048946f, 0.049109f, 0.049293f, 0.049499f, 0.049729f, 0.049984f,
        0.050267f, 0.050579f, 0.050923f, 0.051300f, 0.051713f, 0.052170f, 0.052681f, 0.053254f, 0.053897f, 0.054621f,
        0.055435f, 0.056346f, 0.057365f, 0.058500f, 0.059761f, 0.061163f, 0.062720f, 0.064446f, 0.066357f, 0.068467f,
        0.070792f, 0.073346f, 0.076144f, 0.079200f, 0.082530f, 0.086152f, 0.090082f, 0.094338f, 0.098938f, 0.103900f,
        0.109240f, 0.114970f, 0.121120f, 0.127700f, 0.134720f, 0.142170f, 0.150020f, 0.158250f, 0.166830f, 0.175750f,
        0.184980f, 0.194500f, 0.204280f, 0.214300f, 0.224540f, 0.234950f, 0.245490f, 0.256110f, 0.266770f, 0.277410f,
        0.288010f, 0.298500f, 0.308850f, 0.319000f, 0.328920f, 0.338580f, 0.347970f, 0.357050f, 0.365800f, 0.374220f,
        0.382260f, 0.389930f, 0.397180f, 0.404000f, 0.410380f, 0.416320f, 0.421850f, 0.426980f, 0.431720f, 0.436100f,
        0.440130f, 0.443830f, 0.447210f, 0.450300f, 0.453110f, 0.455650f, 0.457950f, 0.460020f, 0.461880f, 0.463540f,
        0.465030f, 0.466360f, 0.467540f, 0.468600f, 0.469550f, 0.470390f, 0.471130f, 0.471770f, 0.472300f, 0.472750f,
        0.473090f, 0.473350f, 0.473520f, 0.473600f, 0.473600f, 0.473520f, 0.473380f, 0.473190f, 0.472960f, 0.472690f,
        0.472400f, 0.472100f, 0.471790f, 0.471500f, 0.471220f, 0.470970f, 0.470730f, 0.470500f, 0.470300f, 0.470110f,
        0.469940f, 0.469780f, 0.469630f, 0.469500f, 0.469380f, 0.469270f, 0.469160f, 0.469040f, 0.468910f, 0.468770f,
        0.468600f, 0.468400f, 0.468170f, 0.467900f, 0.467580f, 0.467230f, 0.466850f, 0.466460f, 0.466070f, 0.465690f,
        0.465330f, 0.465010f, 0.464730f, 0.464500f, 0.464160f, 0.463860f, 0.463560f, 0.463260f, 0.462970f, 0.462670f,
        0.462370f, 0.462070f, 0.461780f, 0.461480f, 0.461180f, 0.460880f, 0.460580f, 0.460290f, 0.459990f, 0.459690f,
        0.459390f, 0.459100f, 0.458800f, 0.458500f, 0.458210f, 0.457910f, 0.457610f, 0.457310f, 0.457020f, 0.456720f,
        0.456420f, 0.456120f, 0.455830f, 0.455530f, 0.455230f, 0.454940f, 0.454640f, 0.454340f, 0.454040f, 0.453750f,
        0.453450f, 0.453150f, 0.452860f, 0.452560f, 0.452260f, 0.451970f, 0.451670f, 0.451370f, 0.451080f, 0.450780f,
        0.450480f, 0.450190f, 0.449890f, 0.449590f, 0.449300f, 0.449000f, 0.448710f, 0.448410f, 0.448110f, 0.447820f,
        0.447520f, 0.447220f, 0.446930f, 0.446630f, 0.446340f, 0.446040f, 0.445740f, 0.445450f, 0.445150f, 0.444860f,
        0.444560f, 0.444270f, 0.443970f, 0.443670f, 0.443380f, 0.443080f, 0.442790f, 0.442490f, 0.442200f, 0.441900f,
        0.441610f, 0.441310f, 0.441010f});
    static constexpr Spectrum CES08({0.058180f, 0.058588f, 0.058930f, 0.059212f, 0.059437f, 0.059610f, 0.059734f,
        0.059815f, 0.059857f, 0.059864f, 0.059840f, 0.059790f, 0.059717f, 0.059627f, 0.059523f, 0.059410f, 0.059293f,
        0.059175f, 0.059060f, 0.058954f, 0.058860f, 0.058782f, 0.058717f, 0.058664f, 0.058619f, 0.058579f, 0.058541f,
        0.058504f, 0.058463f, 0.058416f, 0.058360f, 0.058293f, 0.058218f, 0.058138f, 0.058055f, 0.057973f, 0.057896f,
        0.057826f, 0.057766f, 0.057719f, 0.057690f, 0.057679f, 0.057686f, 0.057705f, 0.057735f, 0.057772f, 0.057813f,
        0.057855f, 0.057893f, 0.057926f, 0.057950f, 0.057962f, 0.057963f, 0.057952f, 0.057932f, 0.057902f, 0.057864f,
        0.057818f, 0.057764f, 0.057705f, 0.057640f, 0.057570f, 0.057498f, 0.057425f, 0.057354f, 0.057286f, 0.057223f,
        0.057169f, 0.057123f, 0.057090f, 0.057070f, 0.057065f, 0.057073f, 0.057091f, 0.057116f, 0.057145f, 0.057176f,
        0.057204f, 0.057228f, 0.057244f, 0.057250f, 0.057243f, 0.057224f, 0.057195f, 0.057158f, 0.057113f, 0.057064f,
        0.057012f, 0.056957f, 0.056903f, 0.056850f, 0.056800f, 0.056754f, 0.056710f, 0.056669f, 0.056631f, 0.056595f,
        0.056561f, 0.056529f, 0.056499f, 0.056470f, 0.056443f, 0.056415f, 0.056388f, 0.056358f, 0.056326f, 0.056290f,
        0.056250f, 0.056204f, 0.056151f, 0.056090f, 0.056021f, 0.055945f, 0.055864f, 0.055779f, 0.055692f, 0.055605f,
        0.055520f, 0.055438f, 0.055361f, 0.055290f, 0.055227f, 0.055173f, 0.055125f, 0.055085f, 0.055051f, 0.055023f,
        0.055000f, 0.054982f, 0.054969f, 0.054960f, 0.054954f, 0.054951f, 0.054950f, 0.054949f, 0.054948f, 0.054946f,
        0.054942f, 0.054936f, 0.054925f, 0.054910f, 0.054889f, 0.054864f, 0.054835f, 0.054803f, 0.054768f, 0.054732f,
        0.054695f, 0.054659f, 0.054623f, 0.054590f, 0.054559f, 0.054531f, 0.054506f, 0.054483f, 0.054464f, 0.054447f,
        0.054433f, 0.054423f, 0.054415f, 0.054410f, 0.054408f, 0.054409f, 0.054412f, 0.054416f, 0.054421f, 0.054426f,
        0.054431f, 0.054436f, 0.054439f, 0.054440f, 0.054439f, 0.054437f, 0.054436f, 0.054438f, 0.054443f, 0.054453f,
        0.054471f, 0.054497f, 0.054532f, 0.054580f, 0.054640f, 0.054714f, 0.054799f, 0.054896f, 0.055005f, 0.055126f,
        0.055257f, 0.055398f, 0.055549f, 0.055710f, 0.055881f, 0.056070f, 0.056283f, 0.056530f, 0.056818f, 0.057154f,
        0.057548f, 0.058006f, 0.058538f, 0.059150f, 0.059852f, 0.060654f, 0.061570f, 0.062611f, 0.063788f, 0.065116f,
        0.066604f, 0.068267f, 0.070114f, 0.072160f, 0.074415f, 0.076892f, 0.079602f, 0.082557f, 0.085769f, 0.089248f,
        0.093007f, 0.097058f, 0.101410f, 0.106080f, 0.111070f, 0.116380f, 0.122010f, 0.127940f, 0.134170f, 0.140700f,
        0.147510f, 0.154610f, 0.161980f, 0.169620f, 0.177520f, 0.185660f, 0.194010f, 0.202570f, 0.211300f, 0.220180f,
        0.229200f, 0.238330f, 0.247550f, 0.256840f, 0.266180f, 0.275540f, 0.284920f, 0.294300f, 0.303660f, 0.312980f,
        0.322250f, 0.331460f, 0.340570f, 0.349590f, 0.358490f, 0.367260f, 0.375880f, 0.384350f, 0.392660f, 0.400780f,
        0.408720f, 0.416450f, 0.423960f, 0.431250f, 0.438300f, 0.445110f, 0.451690f, 0.458050f, 0.464170f, 0.470080f,
        0.475780f, 0.481260f, 0.486530f, 0.491600f, 0.496470f, 0.501150f, 0.505630f, 0.509930f, 0.514060f, 0.518010f,
        0.521800f, 0.525420f, 0.528880f, 0.532200f, 0.535370f, 0.538400f, 0.541310f, 0.544090f, 0.546760f, 0.549320f,
        0.551790f, 0.554170f, 0.556470f, 0.558700f, 0.560860f, 0.562970f, 0.565020f, 0.567020f, 0.568970f, 0.570870f,
        0.572730f, 0.574560f, 0.576360f, 0.578120f, 0.579860f, 0.581570f, 0.583270f, 0.584950f, 0.586620f, 0.588270f,
        0.589910f, 0.591550f, 0.593190f, 0.594830f, 0.596470f, 0.598110f, 0.599750f, 0.601390f, 0.603030f, 0.604660f,
        0.606280f, 0.607890f, 0.609500f, 0.611090f, 0.612670f, 0.614230f, 0.615770f, 0.617290f, 0.618780f, 0.620240f,
        0.621670f, 0.623070f, 0.624420f, 0.625740f, 0.627010f, 0.628240f, 0.629440f, 0.630600f, 0.631730f, 0.632840f,
        0.633930f, 0.634990f, 0.636050f, 0.637090f, 0.638130f, 0.639160f, 0.640200f, 0.641240f, 0.642280f, 0.643350f,
        0.644420f, 0.645520f, 0.646640f, 0.647790f, 0.648670f, 0.649720f, 0.650770f, 0.651810f, 0.652850f, 0.653890f,
        0.654930f, 0.655970f, 0.657000f, 0.658040f, 0.659070f, 0.660100f, 0.661130f, 0.662160f, 0.663190f, 0.664210f,
        0.665240f, 0.666260f, 0.667280f, 0.668300f, 0.669320f, 0.670330f, 0.671350f, 0.672360f, 0.673370f, 0.674380f,
        0.675390f, 0.676390f, 0.677400f, 0.678400f, 0.679400f, 0.680400f, 0.681400f, 0.682400f, 0.683390f, 0.684390f,
        0.685380f, 0.686370f, 0.687350f, 0.688340f, 0.689330f, 0.690310f, 0.691290f, 0.692270f, 0.693250f, 0.694220f,
        0.695200f, 0.696170f, 0.697140f});
    static constexpr Spectrum CES09({0.041573f, 0.041603f, 0.041633f, 0.041663f, 0.041693f, 0.041723f, 0.041753f,
        0.041783f, 0.041813f, 0.041843f, 0.041873f, 0.041903f, 0.041933f, 0.041964f, 0.041994f, 0.042024f, 0.042054f,
        0.042085f, 0.042115f, 0.042145f, 0.042249f, 0.042248f, 0.042252f, 0.042261f, 0.042274f, 0.042292f, 0.042314f,
        0.042342f, 0.042375f, 0.042413f, 0.042455f, 0.042503f, 0.042556f, 0.042615f, 0.042679f, 0.042748f, 0.042822f,
        0.042902f, 0.042988f, 0.043080f, 0.043177f, 0.043279f, 0.043388f, 0.043501f, 0.043619f, 0.043742f, 0.043869f,
        0.044000f, 0.044134f, 0.044272f, 0.044413f, 0.044557f, 0.044703f, 0.044852f, 0.045003f, 0.045157f, 0.045313f,
        0.045471f, 0.045631f, 0.045794f, 0.045959f, 0.046126f, 0.046294f, 0.046465f, 0.046638f, 0.046812f, 0.046989f,
        0.047167f, 0.047346f, 0.047528f, 0.047711f, 0.047895f, 0.048079f, 0.048262f, 0.048441f, 0.048616f, 0.048784f,
        0.048944f, 0.049094f, 0.049233f, 0.049359f, 0.049472f, 0.049571f, 0.049659f, 0.049737f, 0.049806f, 0.049869f,
        0.049926f, 0.049980f, 0.050031f, 0.050081f, 0.050132f, 0.050189f, 0.050256f, 0.050339f, 0.050441f, 0.050568f,
        0.050725f, 0.050916f, 0.051146f, 0.051420f, 0.051742f, 0.052106f, 0.052507f, 0.052939f, 0.053397f, 0.053873f,
        0.054363f, 0.054860f, 0.055358f, 0.055851f, 0.056334f, 0.056801f, 0.057245f, 0.057662f, 0.058045f, 0.058388f,
        0.058686f, 0.058934f, 0.059124f, 0.059252f, 0.059313f, 0.059312f, 0.059253f, 0.059143f, 0.058987f, 0.058790f,
        0.058557f, 0.058296f, 0.058010f, 0.057706f, 0.057389f, 0.057060f, 0.056721f, 0.056374f, 0.056020f, 0.055661f,
        0.055298f, 0.054933f, 0.054567f, 0.054203f, 0.053840f, 0.053479f, 0.053119f, 0.052760f, 0.052402f, 0.052043f,
        0.051684f, 0.051323f, 0.050961f, 0.050596f, 0.050229f, 0.049858f, 0.049483f, 0.049103f, 0.048717f, 0.048325f,
        0.047926f, 0.047520f, 0.047104f, 0.046680f, 0.046246f, 0.045805f, 0.045358f, 0.044907f, 0.044455f, 0.044004f,
        0.043556f, 0.043112f, 0.042676f, 0.042249f, 0.041834f, 0.041435f, 0.041056f, 0.040702f, 0.040378f, 0.040089f,
        0.039838f, 0.039630f, 0.039471f, 0.039364f, 0.039313f, 0.039315f, 0.039367f, 0.039465f, 0.039604f, 0.039783f,
        0.039996f, 0.040240f, 0.040512f, 0.040807f, 0.041122f, 0.041457f, 0.041811f, 0.042185f, 0.042578f, 0.042990f,
        0.043421f, 0.043870f, 0.044339f, 0.044825f, 0.045330f, 0.045851f, 0.046387f, 0.046936f, 0.047495f, 0.048063f,
        0.048639f, 0.049220f, 0.049804f, 0.050390f, 0.050975f, 0.051556f, 0.052128f, 0.052687f, 0.053230f, 0.053751f,
        0.054247f, 0.054714f, 0.055147f, 0.055542f, 0.055897f, 0.056212f, 0.056492f, 0.056739f, 0.056957f, 0.057148f,
        0.057315f, 0.057462f, 0.057591f, 0.057706f, 0.057810f, 0.057903f, 0.057987f, 0.058064f, 0.058134f, 0.058199f,
        0.058260f, 0.058317f, 0.058373f, 0.058428f, 0.058483f, 0.058538f, 0.058593f, 0.058648f, 0.058701f, 0.058754f,
        0.058804f, 0.058853f, 0.058899f, 0.058943f, 0.058984f, 0.059025f, 0.059068f, 0.059116f, 0.059172f, 0.059238f,
        0.059318f, 0.059414f, 0.059528f, 0.059664f, 0.059823f, 0.060004f, 0.060204f, 0.060421f, 0.060651f, 0.060894f,
        0.061146f, 0.061404f, 0.061667f, 0.061931f, 0.062195f, 0.062456f, 0.062712f, 0.062961f, 0.063201f, 0.063431f,
        0.063647f, 0.063849f, 0.064033f, 0.064198f, 0.064343f, 0.064468f, 0.064574f, 0.064664f, 0.064738f, 0.064797f,
        0.064844f, 0.064879f, 0.064904f, 0.064919f, 0.064928f, 0.064929f, 0.064923f, 0.064911f, 0.064892f, 0.064867f,
        0.064836f, 0.064800f, 0.064759f, 0.064713f, 0.064663f, 0.064608f, 0.064548f, 0.064483f, 0.064413f, 0.064339f,
        0.064260f, 0.064176f, 0.064086f, 0.063992f, 0.063893f, 0.063788f, 0.063678f, 0.063563f, 0.063442f, 0.063316f,
        0.063184f, 0.063047f, 0.062904f, 0.062756f, 0.062728f, 0.062617f, 0.062505f, 0.062394f, 0.062283f, 0.062173f,
        0.062062f, 0.061952f, 0.061842f, 0.061732f, 0.061622f, 0.061512f, 0.061403f, 0.061294f, 0.061184f, 0.061076f,
        0.060967f, 0.060858f, 0.060750f, 0.060642f, 0.060534f, 0.060426f, 0.060318f, 0.060211f, 0.060104f, 0.059997f,
        0.059890f, 0.059783f, 0.059676f, 0.059570f, 0.059464f, 0.059358f, 0.059252f, 0.059146f, 0.059041f, 0.058936f,
        0.058830f, 0.058725f, 0.058621f, 0.058516f, 0.058412f, 0.058307f, 0.058203f, 0.058099f, 0.057996f, 0.057892f,
        0.057789f, 0.057686f, 0.057582f, 0.057480f, 0.057377f, 0.057274f, 0.057172f, 0.057070f, 0.056968f, 0.056866f,
        0.056764f, 0.056663f, 0.056562f, 0.056460f, 0.056359f, 0.056259f, 0.056158f, 0.056057f, 0.055957f, 0.055857f,
        0.055757f, 0.055657f, 0.055558f, 0.055458f, 0.055359f, 0.055260f, 0.055161f, 0.055062f, 0.054963f, 0.054865f,
        0.054767f, 0.054668f, 0.054570f});
    static constexpr Spectrum CES10({0.176020f, 0.177150f, 0.178280f, 0.179420f, 0.180570f, 0.181720f, 0.182870f,
        0.184040f, 0.185200f, 0.186380f, 0.187560f, 0.188740f, 0.189930f, 0.191130f, 0.192330f, 0.193540f, 0.194760f,
        0.195970f, 0.197200f, 0.198430f, 0.199830f, 0.200810f, 0.202010f, 0.203370f, 0.204810f, 0.206280f, 0.207690f,
        0.208990f, 0.210110f, 0.210970f, 0.211520f, 0.211690f, 0.211520f, 0.211060f, 0.210350f, 0.209440f, 0.208380f,
        0.207200f, 0.205970f, 0.204730f, 0.203520f, 0.202380f, 0.201310f, 0.200320f, 0.199380f, 0.198500f, 0.197670f,
        0.196890f, 0.196140f, 0.195420f, 0.194730f, 0.194050f, 0.193390f, 0.192750f, 0.192120f, 0.191490f, 0.190880f,
        0.190260f, 0.189650f, 0.189030f, 0.188420f, 0.187790f, 0.187170f, 0.186540f, 0.185910f, 0.185290f, 0.184670f,
        0.184060f, 0.183460f, 0.182880f, 0.182310f, 0.181750f, 0.181220f, 0.180710f, 0.180230f, 0.179780f, 0.179350f,
        0.178970f, 0.178620f, 0.178300f, 0.178040f, 0.177810f, 0.177620f, 0.177460f, 0.177320f, 0.177180f, 0.177050f,
        0.176910f, 0.176760f, 0.176580f, 0.176370f, 0.176110f, 0.175810f, 0.175470f, 0.175090f, 0.174670f, 0.174200f,
        0.173690f, 0.173140f, 0.172550f, 0.171920f, 0.171250f, 0.170550f, 0.169820f, 0.169080f, 0.168340f, 0.167600f,
        0.166880f, 0.166170f, 0.165510f, 0.164880f, 0.164300f, 0.163780f, 0.163300f, 0.162880f, 0.162500f, 0.162170f,
        0.161890f, 0.161660f, 0.161480f, 0.161350f, 0.161260f, 0.161210f, 0.161190f, 0.161190f, 0.161200f, 0.161210f,
        0.161220f, 0.161210f, 0.161180f, 0.161110f, 0.161000f, 0.160870f, 0.160730f, 0.160590f, 0.160470f, 0.160390f,
        0.160360f, 0.160400f, 0.160530f, 0.160750f, 0.161090f, 0.161550f, 0.162110f, 0.162790f, 0.163570f, 0.164460f,
        0.165460f, 0.166550f, 0.167750f, 0.169040f, 0.170430f, 0.171920f, 0.173490f, 0.175150f, 0.176890f, 0.178710f,
        0.180610f, 0.182580f, 0.184630f, 0.186730f, 0.188910f, 0.191160f, 0.193500f, 0.195940f, 0.198500f, 0.201200f,
        0.204030f, 0.207030f, 0.210190f, 0.213550f, 0.217100f, 0.220850f, 0.224820f, 0.229020f, 0.233450f, 0.238120f,
        0.243040f, 0.248220f, 0.253660f, 0.259380f, 0.265390f, 0.271740f, 0.278490f, 0.285700f, 0.293410f, 0.301680f,
        0.310580f, 0.320160f, 0.330470f, 0.341570f, 0.353500f, 0.366210f, 0.379620f, 0.393680f, 0.408300f, 0.423420f,
        0.438960f, 0.454870f, 0.471070f, 0.487490f, 0.504060f, 0.520710f, 0.537390f, 0.554010f, 0.570530f, 0.586860f,
        0.602950f, 0.618720f, 0.634120f, 0.649080f, 0.663530f, 0.677460f, 0.690850f, 0.703690f, 0.715960f, 0.727650f,
        0.738740f, 0.749220f, 0.759080f, 0.768300f, 0.776870f, 0.784800f, 0.792140f, 0.798890f, 0.805080f, 0.810730f,
        0.815880f, 0.820540f, 0.824730f, 0.828490f, 0.831830f, 0.834800f, 0.837430f, 0.839770f, 0.841850f, 0.843710f,
        0.845400f, 0.846950f, 0.848400f, 0.849790f, 0.851150f, 0.852490f, 0.853810f, 0.855100f, 0.856370f, 0.857600f,
        0.858800f, 0.859960f, 0.861090f, 0.862170f, 0.863210f, 0.864210f, 0.865170f, 0.866080f, 0.866950f, 0.867770f,
        0.868550f, 0.869290f, 0.869990f, 0.870640f, 0.871250f, 0.871820f, 0.872360f, 0.872870f, 0.873370f, 0.873850f,
        0.874320f, 0.874790f, 0.875270f, 0.875750f, 0.876240f, 0.876750f, 0.877250f, 0.877750f, 0.878220f, 0.878680f,
        0.879100f, 0.879490f, 0.879820f, 0.880110f, 0.880330f, 0.880490f, 0.880610f, 0.880690f, 0.880740f, 0.880760f,
        0.880760f, 0.880760f, 0.880750f, 0.880750f, 0.880760f, 0.880780f, 0.880810f, 0.880840f, 0.880860f, 0.880880f,
        0.880890f, 0.880890f, 0.880870f, 0.880830f, 0.880770f, 0.880690f, 0.880600f, 0.880510f, 0.880440f, 0.880380f,
        0.880350f, 0.880350f, 0.880400f, 0.880500f, 0.880660f, 0.880870f, 0.881110f, 0.881380f, 0.881670f, 0.881960f,
        0.882240f, 0.882500f, 0.882720f, 0.882910f, 0.883180f, 0.883410f, 0.883650f, 0.883890f, 0.884120f, 0.884360f,
        0.884590f, 0.884830f, 0.885060f, 0.885290f, 0.885520f, 0.885760f, 0.885990f, 0.886220f, 0.886450f, 0.886680f,
        0.886910f, 0.887140f, 0.887370f, 0.887600f, 0.887830f, 0.888060f, 0.888290f, 0.888510f, 0.888740f, 0.888970f,
        0.889190f, 0.889420f, 0.889640f, 0.889870f, 0.890090f, 0.890320f, 0.890540f, 0.890770f, 0.890990f, 0.891210f,
        0.891430f, 0.891660f, 0.891880f, 0.892100f, 0.892320f, 0.892540f, 0.892760f, 0.892980f, 0.893200f, 0.893420f,
        0.893630f, 0.893850f, 0.894070f, 0.894290f, 0.894500f, 0.894720f, 0.894940f, 0.895150f, 0.895370f, 0.895580f,
        0.895800f, 0.896010f, 0.896220f, 0.896440f, 0.896650f, 0.896860f, 0.897070f, 0.897290f, 0.897500f, 0.897710f,
        0.897920f, 0.898130f, 0.898340f, 0.898550f, 0.898760f, 0.898960f, 0.899170f, 0.899380f, 0.899590f, 0.899800f,
        0.900000f, 0.900210f, 0.900410f});
    static constexpr Spectrum CES11({0.039511f, 0.040462f, 0.041434f, 0.042429f, 0.043447f, 0.044489f, 0.045554f,
        0.046643f, 0.047757f, 0.048896f, 0.050061f, 0.051253f, 0.052471f, 0.053716f, 0.054989f, 0.056291f, 0.057621f,
        0.058981f, 0.060371f, 0.061792f, 0.060810f, 0.063415f, 0.065866f, 0.068166f, 0.070318f, 0.072327f, 0.074195f,
        0.075925f, 0.077521f, 0.078986f, 0.080324f, 0.081538f, 0.082630f, 0.083606f, 0.084466f, 0.085217f, 0.085859f,
        0.086397f, 0.086835f, 0.087175f, 0.087420f, 0.087576f, 0.087650f, 0.087652f, 0.087591f, 0.087477f, 0.087319f,
        0.087126f, 0.086908f, 0.086675f, 0.086435f, 0.086195f, 0.085956f, 0.085711f, 0.085459f, 0.085196f, 0.084917f,
        0.084619f, 0.084299f, 0.083953f, 0.083576f, 0.083168f, 0.082726f, 0.082251f, 0.081745f, 0.081206f, 0.080637f,
        0.080037f, 0.079407f, 0.078747f, 0.078057f, 0.077339f, 0.076595f, 0.075826f, 0.075036f, 0.074227f, 0.073400f,
        0.072560f, 0.071707f, 0.070845f, 0.069976f, 0.069101f, 0.068221f, 0.067335f, 0.066442f, 0.065542f, 0.064633f,
        0.063717f, 0.062791f, 0.061855f, 0.060908f, 0.059952f, 0.058989f, 0.058027f, 0.057069f, 0.056123f, 0.055192f,
        0.054284f, 0.053402f, 0.052553f, 0.051742f, 0.050974f, 0.050249f, 0.049566f, 0.048925f, 0.048326f, 0.047767f,
        0.047248f, 0.046769f, 0.046329f, 0.045928f, 0.045563f, 0.045233f, 0.044933f, 0.044660f, 0.044410f, 0.044178f,
        0.043961f, 0.043756f, 0.043559f, 0.043365f, 0.043172f, 0.042981f, 0.042793f, 0.042609f, 0.042432f, 0.042262f,
        0.042101f, 0.041951f, 0.041814f, 0.041690f, 0.041581f, 0.041488f, 0.041412f, 0.041352f, 0.041311f, 0.041288f,
        0.041284f, 0.041300f, 0.041337f, 0.041394f, 0.041473f, 0.041577f, 0.041709f, 0.041871f, 0.042067f, 0.042300f,
        0.042573f, 0.042889f, 0.043250f, 0.043661f, 0.044125f, 0.044656f, 0.045268f, 0.045974f, 0.046791f, 0.047732f,
        0.048811f, 0.050043f, 0.051442f, 0.053024f, 0.054801f, 0.056791f, 0.059008f, 0.061468f, 0.064186f, 0.067178f,
        0.070459f, 0.074046f, 0.077953f, 0.082197f, 0.086789f, 0.091729f, 0.097014f, 0.102640f, 0.108610f, 0.114910f,
        0.121540f, 0.128500f, 0.135790f, 0.143400f, 0.151320f, 0.159530f, 0.167990f, 0.176660f, 0.185520f, 0.194520f,
        0.203640f, 0.212830f, 0.222070f, 0.231310f, 0.240540f, 0.249700f, 0.258770f, 0.267720f, 0.276510f, 0.285120f,
        0.293490f, 0.301620f, 0.309450f, 0.316960f, 0.324120f, 0.330930f, 0.337400f, 0.343530f, 0.349330f, 0.354800f,
        0.359950f, 0.364780f, 0.369310f, 0.373530f, 0.377460f, 0.381100f, 0.384480f, 0.387610f, 0.390510f, 0.393180f,
        0.395660f, 0.397940f, 0.400060f, 0.402010f, 0.403830f, 0.405520f, 0.407100f, 0.408570f, 0.409960f, 0.411270f,
        0.412510f, 0.413710f, 0.414870f, 0.416010f, 0.417140f, 0.418260f, 0.419390f, 0.420540f, 0.421720f, 0.422930f,
        0.424190f, 0.425500f, 0.426880f, 0.428330f, 0.429860f, 0.431480f, 0.433180f, 0.434960f, 0.436840f, 0.438800f,
        0.440850f, 0.442990f, 0.445220f, 0.447550f, 0.449960f, 0.452470f, 0.455050f, 0.457710f, 0.460450f, 0.463240f,
        0.466090f, 0.469000f, 0.471950f, 0.474950f, 0.477980f, 0.481060f, 0.484190f, 0.487390f, 0.490660f, 0.494030f,
        0.497490f, 0.501050f, 0.504740f, 0.508550f, 0.512510f, 0.516590f, 0.520810f, 0.525150f, 0.529620f, 0.534210f,
        0.538920f, 0.543740f, 0.548660f, 0.553690f, 0.558820f, 0.564030f, 0.569280f, 0.574550f, 0.579820f, 0.585050f,
        0.590230f, 0.595320f, 0.600300f, 0.605140f, 0.609820f, 0.614360f, 0.618760f, 0.623060f, 0.627260f, 0.631390f,
        0.635460f, 0.639490f, 0.643510f, 0.647520f, 0.651550f, 0.655580f, 0.659600f, 0.663620f, 0.667600f, 0.671550f,
        0.675450f, 0.679290f, 0.683060f, 0.686750f, 0.690340f, 0.693840f, 0.697220f, 0.700480f, 0.703600f, 0.706580f,
        0.709400f, 0.712050f, 0.714520f, 0.716810f, 0.721660f, 0.724740f, 0.727800f, 0.730850f, 0.733860f, 0.736860f,
        0.739840f, 0.742790f, 0.745720f, 0.748630f, 0.751520f, 0.754390f, 0.757230f, 0.760050f, 0.762850f, 0.765620f,
        0.768370f, 0.771100f, 0.773810f, 0.776500f, 0.779160f, 0.781800f, 0.784410f, 0.787010f, 0.789580f, 0.792120f,
        0.794650f, 0.797150f, 0.799630f, 0.802090f, 0.804520f, 0.806930f, 0.809320f, 0.811690f, 0.814030f, 0.816350f,
        0.818650f, 0.820920f, 0.823180f, 0.825410f, 0.827620f, 0.829800f, 0.831970f, 0.834110f, 0.836230f, 0.838330f,
        0.840410f, 0.842460f, 0.844490f, 0.846510f, 0.848500f, 0.850470f, 0.852410f, 0.854340f, 0.856250f, 0.858130f,
        0.860000f, 0.861840f, 0.863670f, 0.865470f, 0.867250f, 0.869020f, 0.870760f, 0.872480f, 0.874190f, 0.875870f,
        0.877540f, 0.879180f, 0.880810f, 0.882420f, 0.884010f, 0.885580f, 0.887130f, 0.888660f, 0.890180f, 0.891670f,
        0.893150f, 0.894610f, 0.896060f});
    static constexpr Spectrum CES12({0.241670f, 0.237160f, 0.232670f, 0.227810f, 0.222130f, 0.216220f, 0.210720f,
        0.204830f, 0.198910f, 0.192980f, 0.186870f, 0.180780f, 0.175420f, 0.170620f, 0.166080f, 0.162040f, 0.158330f,
        0.154020f, 0.150850f, 0.148150f, 0.145140f, 0.142150f, 0.139490f, 0.136460f, 0.133080f, 0.129810f, 0.126550f,
        0.123400f, 0.120610f, 0.118020f, 0.115360f, 0.113190f, 0.111050f, 0.108840f, 0.107070f, 0.105840f, 0.105190f,
        0.104970f, 0.104900f, 0.104720f, 0.104300f, 0.103830f, 0.102810f, 0.101450f, 0.099815f, 0.098708f, 0.097967f,
        0.096978f, 0.095877f, 0.095158f, 0.094833f, 0.094687f, 0.094808f, 0.095574f, 0.096729f, 0.097671f, 0.098137f,
        0.098167f, 0.098439f, 0.098831f, 0.098292f, 0.097482f, 0.096940f, 0.096416f, 0.095543f, 0.094735f, 0.095001f,
        0.095409f, 0.095288f, 0.094672f, 0.093947f, 0.093808f, 0.093406f, 0.092917f, 0.092431f, 0.091913f, 0.091230f,
        0.090197f, 0.088798f, 0.087422f, 0.086025f, 0.084641f, 0.083285f, 0.082259f, 0.082097f, 0.081328f, 0.080320f,
        0.079021f, 0.077459f, 0.076197f, 0.075411f, 0.074464f, 0.073308f, 0.072032f, 0.070652f, 0.068276f, 0.066607f,
        0.065116f, 0.063658f, 0.062596f, 0.061592f, 0.060002f, 0.058833f, 0.058186f, 0.057778f, 0.057188f, 0.056242f,
        0.055110f, 0.054417f, 0.054342f, 0.054162f, 0.053723f, 0.053394f, 0.052825f, 0.051979f, 0.050865f, 0.050007f,
        0.049561f, 0.049321f, 0.049053f, 0.048672f, 0.048363f, 0.048595f, 0.049059f, 0.049520f, 0.049948f, 0.050598f,
        0.050997f, 0.051167f, 0.051604f, 0.052351f, 0.053125f, 0.053939f, 0.054441f, 0.054939f, 0.055646f, 0.056306f,
        0.056772f, 0.057617f, 0.058880f, 0.059795f, 0.060268f, 0.060745f, 0.061345f, 0.062326f, 0.063376f, 0.064193f,
        0.065223f, 0.067008f, 0.069284f, 0.071189f, 0.073067f, 0.075259f, 0.077662f, 0.080018f, 0.082032f, 0.084366f,
        0.087440f, 0.090541f, 0.093574f, 0.096584f, 0.099907f, 0.103400f, 0.106830f, 0.110170f, 0.113740f, 0.117600f,
        0.121370f, 0.125090f, 0.128740f, 0.132500f, 0.136330f, 0.140610f, 0.145390f, 0.150190f, 0.155160f, 0.160280f,
        0.165590f, 0.171190f, 0.176960f, 0.183310f, 0.189980f, 0.196770f, 0.203500f, 0.210000f, 0.216580f, 0.223300f,
        0.230350f, 0.237750f, 0.245500f, 0.253340f, 0.260960f, 0.268430f, 0.275780f, 0.283420f, 0.291350f, 0.299790f,
        0.308750f, 0.316890f, 0.325510f, 0.333920f, 0.342430f, 0.351640f, 0.361310f, 0.370680f, 0.379660f, 0.388760f,
        0.398000f, 0.407260f, 0.417390f, 0.426910f, 0.436060f, 0.445220f, 0.454250f, 0.462800f, 0.471370f, 0.480490f,
        0.489300f, 0.497370f, 0.505020f, 0.512600f, 0.520250f, 0.528610f, 0.535610f, 0.542050f, 0.548740f, 0.555880f,
        0.562180f, 0.567910f, 0.574020f, 0.579600f, 0.584770f, 0.589490f, 0.593210f, 0.598320f, 0.603020f, 0.607080f,
        0.610670f, 0.614310f, 0.618260f, 0.622300f, 0.625890f, 0.629400f, 0.632680f, 0.635720f, 0.637690f, 0.639890f,
        0.642050f, 0.643840f, 0.645530f, 0.647590f, 0.649760f, 0.651870f, 0.653790f, 0.655390f, 0.657050f, 0.659020f,
        0.660820f, 0.661820f, 0.663880f, 0.666060f, 0.667590f, 0.668480f, 0.670130f, 0.671120f, 0.672390f, 0.674080f,
        0.675780f, 0.677160f, 0.678650f, 0.678830f, 0.678630f, 0.678060f, 0.677150f, 0.675990f, 0.675620f, 0.675860f,
        0.675900f, 0.675440f, 0.675220f, 0.675550f, 0.676060f, 0.676750f, 0.677620f, 0.678220f, 0.678850f, 0.679440f,
        0.678950f, 0.678020f, 0.677900f, 0.678210f, 0.678340f, 0.678110f, 0.677300f, 0.677170f, 0.678430f, 0.679340f,
        0.679820f, 0.680560f, 0.681540f, 0.681610f, 0.681270f, 0.681350f, 0.682300f, 0.683680f, 0.684120f, 0.683410f,
        0.683070f, 0.682890f, 0.682770f, 0.682730f, 0.683040f, 0.683700f, 0.683670f, 0.683630f, 0.683250f, 0.682010f,
        0.686530f, 0.683100f, 0.686760f, 0.686760f, 0.685530f, 0.685750f, 0.685960f, 0.686170f, 0.686390f, 0.686600f,
        0.686810f, 0.687030f, 0.687240f, 0.687450f, 0.687660f, 0.687880f, 0.688090f, 0.688300f, 0.688510f, 0.688730f,
        0.688940f, 0.689150f, 0.689360f, 0.689570f, 0.689790f, 0.690000f, 0.690210f, 0.690420f, 0.690630f, 0.690850f,
        0.691060f, 0.691270f, 0.691480f, 0.691690f, 0.691900f, 0.692110f, 0.692320f, 0.692530f, 0.692750f, 0.692960f,
        0.693170f, 0.693380f, 0.693590f, 0.693800f, 0.694010f, 0.694220f, 0.694430f, 0.694640f, 0.694850f, 0.695060f,
        0.695270f, 0.695480f, 0.695690f, 0.695900f, 0.696110f, 0.696320f, 0.696530f, 0.696740f, 0.696950f, 0.697150f,
        0.697360f, 0.697570f, 0.697780f, 0.697990f, 0.698200f, 0.698410f, 0.698620f, 0.698820f, 0.699030f, 0.699240f,
        0.699450f, 0.699660f, 0.699870f, 0.700070f, 0.700280f, 0.700490f, 0.700700f, 0.700900f, 0.701110f, 0.701320f,
        0.701530f, 0.701730f, 0.701940f});
    static constexpr Spectrum CES13({0.043374f, 0.043365f, 0.043356f, 0.043347f, 0.043338f, 0.043329f, 0.043319f,
        0.043310f, 0.043301f, 0.043292f, 0.043283f, 0.043274f, 0.043265f, 0.043256f, 0.043247f, 0.043238f, 0.043229f,
        0.043220f, 0.043211f, 0.043202f, 0.043269f, 0.043231f, 0.043196f, 0.043164f, 0.043136f, 0.043112f, 0.043092f,
        0.043078f, 0.043069f, 0.043066f, 0.043070f, 0.043081f, 0.043099f, 0.043125f, 0.043159f, 0.043202f, 0.043255f,
        0.043317f, 0.043389f, 0.043472f, 0.043566f, 0.043672f, 0.043788f, 0.043913f, 0.044047f, 0.044187f, 0.044333f,
        0.044485f, 0.044640f, 0.044797f, 0.044956f, 0.045115f, 0.045272f, 0.045426f, 0.045574f, 0.045715f, 0.045847f,
        0.045968f, 0.046076f, 0.046169f, 0.046246f, 0.046305f, 0.046348f, 0.046377f, 0.046397f, 0.046408f, 0.046415f,
        0.046419f, 0.046424f, 0.046431f, 0.046444f, 0.046465f, 0.046490f, 0.046515f, 0.046536f, 0.046549f, 0.046549f,
        0.046533f, 0.046497f, 0.046435f, 0.046345f, 0.046224f, 0.046078f, 0.045914f, 0.045741f, 0.045567f, 0.045400f,
        0.045247f, 0.045117f, 0.045017f, 0.044956f, 0.044941f, 0.044983f, 0.045093f, 0.045280f, 0.045555f, 0.045930f,
        0.046414f, 0.047018f, 0.047752f, 0.048628f, 0.049650f, 0.050810f, 0.052093f, 0.053485f, 0.054970f, 0.056536f,
        0.058168f, 0.059851f, 0.061572f, 0.063315f, 0.065067f, 0.066815f, 0.068546f, 0.070247f, 0.071905f, 0.073506f,
        0.075039f, 0.076490f, 0.077846f, 0.079094f, 0.080226f, 0.081246f, 0.082164f, 0.082990f, 0.083734f, 0.084403f,
        0.085009f, 0.085561f, 0.086067f, 0.086537f, 0.086979f, 0.087391f, 0.087769f, 0.088108f, 0.088406f, 0.088657f,
        0.088857f, 0.089004f, 0.089092f, 0.089118f, 0.089079f, 0.088978f, 0.088821f, 0.088612f, 0.088356f, 0.088057f,
        0.087721f, 0.087353f, 0.086957f, 0.086537f, 0.086099f, 0.085641f, 0.085164f, 0.084665f, 0.084144f, 0.083601f,
        0.083033f, 0.082441f, 0.081823f, 0.081178f, 0.080507f, 0.079814f, 0.079106f, 0.078389f, 0.077668f, 0.076951f,
        0.076242f, 0.075549f, 0.074876f, 0.074232f, 0.073619f, 0.073039f, 0.072491f, 0.071974f, 0.071486f, 0.071028f,
        0.070598f, 0.070195f, 0.069819f, 0.069468f, 0.069143f, 0.068847f, 0.068587f, 0.068367f, 0.068192f, 0.068068f,
        0.067999f, 0.067991f, 0.068049f, 0.068178f, 0.068388f, 0.068710f, 0.069178f, 0.069827f, 0.070693f, 0.071810f,
        0.073213f, 0.074938f, 0.077019f, 0.079491f, 0.082379f, 0.085660f, 0.089302f, 0.093272f, 0.097538f, 0.102070f,
        0.106830f, 0.111790f, 0.116910f, 0.122160f, 0.127520f, 0.132950f, 0.138410f, 0.143880f, 0.149330f, 0.154730f,
        0.160040f, 0.165230f, 0.170280f, 0.175160f, 0.179830f, 0.184310f, 0.188610f, 0.192750f, 0.196720f, 0.200560f,
        0.204270f, 0.207860f, 0.211350f, 0.214760f, 0.218080f, 0.221310f, 0.224450f, 0.227470f, 0.230370f, 0.233120f,
        0.235730f, 0.238180f, 0.240460f, 0.242540f, 0.244440f, 0.246140f, 0.247680f, 0.249060f, 0.250300f, 0.251410f,
        0.252400f, 0.253300f, 0.254110f, 0.254850f, 0.255530f, 0.256160f, 0.256740f, 0.257290f, 0.257810f, 0.258310f,
        0.258790f, 0.259260f, 0.259730f, 0.260210f, 0.260690f, 0.261180f, 0.261680f, 0.262170f, 0.262660f, 0.263140f,
        0.263610f, 0.264050f, 0.264480f, 0.264870f, 0.265240f, 0.265570f, 0.265880f, 0.266160f, 0.266410f, 0.266640f,
        0.266850f, 0.267040f, 0.267200f, 0.267350f, 0.267490f, 0.267600f, 0.267700f, 0.267780f, 0.267850f, 0.267900f,
        0.267940f, 0.267960f, 0.267960f, 0.267950f, 0.267920f, 0.267880f, 0.267830f, 0.267770f, 0.267700f, 0.267630f,
        0.267560f, 0.267480f, 0.267420f, 0.267350f, 0.267300f, 0.267250f, 0.267200f, 0.267170f, 0.267130f, 0.267100f,
        0.267060f, 0.267030f, 0.266990f, 0.266960f, 0.266910f, 0.266870f, 0.266810f, 0.266750f, 0.266680f, 0.266600f,
        0.266510f, 0.266410f, 0.266290f, 0.266160f, 0.266230f, 0.266170f, 0.266110f, 0.266040f, 0.265980f, 0.265920f,
        0.265850f, 0.265790f, 0.265730f, 0.265660f, 0.265600f, 0.265540f, 0.265470f, 0.265410f, 0.265350f, 0.265280f,
        0.265220f, 0.265160f, 0.265090f, 0.265030f, 0.264970f, 0.264900f, 0.264840f, 0.264780f, 0.264710f, 0.264650f,
        0.264590f, 0.264520f, 0.264460f, 0.264400f, 0.264330f, 0.264270f, 0.264210f, 0.264140f, 0.264080f, 0.264020f,
        0.263950f, 0.263890f, 0.263830f, 0.263760f, 0.263700f, 0.263640f, 0.263570f, 0.263510f, 0.263450f, 0.263380f,
        0.263320f, 0.263260f, 0.263190f, 0.263130f, 0.263070f, 0.263010f, 0.262940f, 0.262880f, 0.262820f, 0.262750f,
        0.262690f, 0.262630f, 0.262560f, 0.262500f, 0.262440f, 0.262380f, 0.262310f, 0.262250f, 0.262190f, 0.262120f,
        0.262060f, 0.262000f, 0.261940f, 0.261870f, 0.261810f, 0.261750f, 0.261680f, 0.261620f, 0.261560f, 0.261500f,
        0.261430f, 0.261370f, 0.261310f});
    static constexpr Spectrum CES14({0.335320f, 0.335620f, 0.335920f, 0.336220f, 0.336530f, 0.336830f, 0.337130f,
        0.337430f, 0.337730f, 0.338030f, 0.338330f, 0.338630f, 0.338930f, 0.339240f, 0.339540f, 0.339840f, 0.340140f,
        0.340440f, 0.340750f, 0.341050f, 0.341400f, 0.341630f, 0.341910f, 0.342240f, 0.342610f, 0.343010f, 0.343430f,
        0.343850f, 0.344280f, 0.344700f, 0.345100f, 0.345470f, 0.345820f, 0.346150f, 0.346460f, 0.346760f, 0.347050f,
        0.347340f, 0.347620f, 0.347910f, 0.348200f, 0.348500f, 0.348810f, 0.349120f, 0.349440f, 0.349760f, 0.350080f,
        0.350400f, 0.350710f, 0.351010f, 0.351300f, 0.351580f, 0.351850f, 0.352110f, 0.352370f, 0.352610f, 0.352860f,
        0.353090f, 0.353330f, 0.353560f, 0.353800f, 0.354040f, 0.354270f, 0.354510f, 0.354750f, 0.354990f, 0.355240f,
        0.355480f, 0.355720f, 0.355960f, 0.356200f, 0.356440f, 0.356680f, 0.356910f, 0.357150f, 0.357380f, 0.357610f,
        0.357840f, 0.358060f, 0.358280f, 0.358500f, 0.358710f, 0.358920f, 0.359130f, 0.359340f, 0.359550f, 0.359770f,
        0.359990f, 0.360220f, 0.360450f, 0.360700f, 0.360960f, 0.361220f, 0.361490f, 0.361770f, 0.362050f, 0.362330f,
        0.362610f, 0.362880f, 0.363140f, 0.363400f, 0.363650f, 0.363880f, 0.364110f, 0.364330f, 0.364540f, 0.364760f,
        0.364970f, 0.365180f, 0.365390f, 0.365600f, 0.365820f, 0.366040f, 0.366270f, 0.366510f, 0.366760f, 0.367020f,
        0.367300f, 0.367580f, 0.367880f, 0.368200f, 0.368540f, 0.368890f, 0.369270f, 0.369670f, 0.370090f, 0.370550f,
        0.371040f, 0.371550f, 0.372110f, 0.372700f, 0.373330f, 0.374000f, 0.374710f, 0.375460f, 0.376250f, 0.377080f,
        0.377950f, 0.378860f, 0.379810f, 0.380800f, 0.381830f, 0.382890f, 0.383970f, 0.385060f, 0.386140f, 0.387210f,
        0.388250f, 0.389250f, 0.390210f, 0.391100f, 0.391930f, 0.392680f, 0.393370f, 0.394000f, 0.394570f, 0.395080f,
        0.395530f, 0.395940f, 0.396290f, 0.396600f, 0.396860f, 0.397090f, 0.397270f, 0.397420f, 0.397530f, 0.397600f,
        0.397640f, 0.397660f, 0.397640f, 0.397600f, 0.397540f, 0.397460f, 0.397380f, 0.397320f, 0.397280f, 0.397280f,
        0.397330f, 0.397440f, 0.397630f, 0.397900f, 0.398270f, 0.398760f, 0.399370f, 0.400120f, 0.401030f, 0.402100f,
        0.403340f, 0.404790f, 0.406430f, 0.408300f, 0.410390f, 0.412690f, 0.415170f, 0.417810f, 0.420580f, 0.423450f,
        0.426410f, 0.429420f, 0.432460f, 0.435500f, 0.438520f, 0.441510f, 0.444440f, 0.447310f, 0.450090f, 0.452770f,
        0.455340f, 0.457780f, 0.460070f, 0.462200f, 0.464160f, 0.465950f, 0.467580f, 0.469070f, 0.470430f, 0.471670f,
        0.472790f, 0.473810f, 0.474750f, 0.475600f, 0.476390f, 0.477110f, 0.477770f, 0.478390f, 0.478950f, 0.479470f,
        0.479950f, 0.480390f, 0.480810f, 0.481200f, 0.481570f, 0.481920f, 0.482250f, 0.482570f, 0.482860f, 0.483140f,
        0.483410f, 0.483650f, 0.483880f, 0.484100f, 0.484300f, 0.484490f, 0.484670f, 0.484830f, 0.484980f, 0.485130f,
        0.485260f, 0.485380f, 0.485490f, 0.485600f, 0.485700f, 0.485790f, 0.485890f, 0.485980f, 0.486070f, 0.486180f,
        0.486290f, 0.486410f, 0.486550f, 0.486700f, 0.486870f, 0.487060f, 0.487260f, 0.487470f, 0.487680f, 0.487890f,
        0.488110f, 0.488310f, 0.488510f, 0.488700f, 0.488870f, 0.489030f, 0.489170f, 0.489300f, 0.489410f, 0.489510f,
        0.489600f, 0.489680f, 0.489740f, 0.489800f, 0.489850f, 0.489890f, 0.489920f, 0.489960f, 0.489990f, 0.490030f,
        0.490080f, 0.490140f, 0.490210f, 0.490300f, 0.490410f, 0.490530f, 0.490670f, 0.490820f, 0.490980f, 0.491160f,
        0.491340f, 0.491520f, 0.491710f, 0.491900f, 0.492090f, 0.492280f, 0.492460f, 0.492650f, 0.492830f, 0.493010f,
        0.493190f, 0.493360f, 0.493530f, 0.493700f, 0.493860f, 0.494020f, 0.494170f, 0.494310f, 0.494440f, 0.494560f,
        0.494670f, 0.494760f, 0.494840f, 0.494900f, 0.495000f, 0.495090f, 0.495170f, 0.495260f, 0.495340f, 0.495430f,
        0.495510f, 0.495600f, 0.495680f, 0.495770f, 0.495850f, 0.495940f, 0.496020f, 0.496110f, 0.496190f, 0.496280f,
        0.496360f, 0.496450f, 0.496530f, 0.496620f, 0.496700f, 0.496790f, 0.496870f, 0.496960f, 0.497040f, 0.497130f,
        0.497210f, 0.497300f, 0.497380f, 0.497470f, 0.497550f, 0.497640f, 0.497720f, 0.497810f, 0.497890f, 0.497980f,
        0.498060f, 0.498150f, 0.498230f, 0.498320f, 0.498400f, 0.498490f, 0.498570f, 0.498660f, 0.498740f, 0.498830f,
        0.498910f, 0.499000f, 0.499080f, 0.499170f, 0.499250f, 0.499340f, 0.499420f, 0.499510f, 0.499590f, 0.499680f,
        0.499760f, 0.499850f, 0.499930f, 0.500020f, 0.500100f, 0.500190f, 0.500270f, 0.500360f, 0.500440f, 0.500530f,
        0.500610f, 0.500700f, 0.500780f, 0.500870f, 0.500950f, 0.501040f, 0.501120f, 0.501210f, 0.501290f, 0.501380f,
        0.501460f, 0.501550f, 0.501630f});
    static constexpr Spectrum CES15({0.169920f, 0.170520f, 0.171130f, 0.171730f, 0.172340f, 0.172950f, 0.173560f,
        0.174170f, 0.174790f, 0.175400f, 0.176020f, 0.176640f, 0.177260f, 0.177880f, 0.178510f, 0.179130f, 0.179760f,
        0.180390f, 0.181020f, 0.181660f, 0.181410f, 0.182510f, 0.183490f, 0.184390f, 0.185200f, 0.185930f, 0.186600f,
        0.187210f, 0.187780f, 0.188320f, 0.188830f, 0.189330f, 0.189820f, 0.190320f, 0.190830f, 0.191370f, 0.191940f,
        0.192560f, 0.193240f, 0.193970f, 0.194790f, 0.195690f, 0.196680f, 0.197770f, 0.198970f, 0.200290f, 0.201730f,
        0.203310f, 0.205020f, 0.206880f, 0.208900f, 0.211070f, 0.213390f, 0.215830f, 0.218390f, 0.221030f, 0.223740f,
        0.226500f, 0.229300f, 0.232110f, 0.234920f, 0.237700f, 0.240460f, 0.243190f, 0.245870f, 0.248510f, 0.251100f,
        0.253630f, 0.256100f, 0.258510f, 0.260830f, 0.263080f, 0.265250f, 0.267350f, 0.269380f, 0.271350f, 0.273240f,
        0.275080f, 0.276860f, 0.278590f, 0.280270f, 0.281900f, 0.283490f, 0.285030f, 0.286530f, 0.287990f, 0.289420f,
        0.290810f, 0.292170f, 0.293500f, 0.294800f, 0.296070f, 0.297320f, 0.298540f, 0.299750f, 0.300930f, 0.302100f,
        0.303250f, 0.304380f, 0.305500f, 0.306600f, 0.307700f, 0.308790f, 0.309880f, 0.310980f, 0.312080f, 0.313210f,
        0.314360f, 0.315530f, 0.316740f, 0.318000f, 0.319290f, 0.320620f, 0.321980f, 0.323370f, 0.324760f, 0.326170f,
        0.327570f, 0.328960f, 0.330340f, 0.331680f, 0.333000f, 0.334270f, 0.335480f, 0.336620f, 0.337670f, 0.338620f,
        0.339460f, 0.340180f, 0.340760f, 0.341190f, 0.341460f, 0.341570f, 0.341530f, 0.341330f, 0.340980f, 0.340490f,
        0.339850f, 0.339070f, 0.338160f, 0.337120f, 0.335950f, 0.334670f, 0.333310f, 0.331880f, 0.330410f, 0.328910f,
        0.327410f, 0.325930f, 0.324490f, 0.323120f, 0.321820f, 0.320610f, 0.319490f, 0.318470f, 0.317560f, 0.316760f,
        0.316070f, 0.315500f, 0.315060f, 0.314760f, 0.314580f, 0.314540f, 0.314620f, 0.314820f, 0.315130f, 0.315540f,
        0.316050f, 0.316650f, 0.317340f, 0.318100f, 0.318930f, 0.319800f, 0.320690f, 0.321570f, 0.322410f, 0.323180f,
        0.323870f, 0.324430f, 0.324850f, 0.325100f, 0.325160f, 0.325050f, 0.324810f, 0.324480f, 0.324070f, 0.323640f,
        0.323210f, 0.322820f, 0.322500f, 0.322280f, 0.322200f, 0.322330f, 0.322710f, 0.323410f, 0.324490f, 0.326020f,
        0.328050f, 0.330630f, 0.333850f, 0.337750f, 0.342370f, 0.347660f, 0.353540f, 0.359960f, 0.366820f, 0.374050f,
        0.381580f, 0.389350f, 0.397260f, 0.405250f, 0.413260f, 0.421230f, 0.429130f, 0.436920f, 0.444570f, 0.452030f,
        0.459280f, 0.466280f, 0.472970f, 0.479340f, 0.485360f, 0.491020f, 0.496360f, 0.501390f, 0.506140f, 0.510630f,
        0.514880f, 0.518900f, 0.522720f, 0.526370f, 0.529850f, 0.533180f, 0.536360f, 0.539390f, 0.542270f, 0.545020f,
        0.547630f, 0.550110f, 0.552460f, 0.554690f, 0.556800f, 0.558790f, 0.560680f, 0.562480f, 0.564190f, 0.565820f,
        0.567380f, 0.568880f, 0.570320f, 0.571720f, 0.573080f, 0.574410f, 0.575710f, 0.576970f, 0.578200f, 0.579400f,
        0.580580f, 0.581730f, 0.582850f, 0.583950f, 0.585030f, 0.586080f, 0.587120f, 0.588130f, 0.589110f, 0.590070f,
        0.591010f, 0.591920f, 0.592810f, 0.593670f, 0.594500f, 0.595310f, 0.596110f, 0.596880f, 0.597650f, 0.598420f,
        0.599180f, 0.599950f, 0.600720f, 0.601500f, 0.602310f, 0.603120f, 0.603950f, 0.604790f, 0.605630f, 0.606490f,
        0.607350f, 0.608220f, 0.609100f, 0.609970f, 0.610840f, 0.611720f, 0.612590f, 0.613450f, 0.614310f, 0.615160f,
        0.616000f, 0.616830f, 0.617640f, 0.618430f, 0.619210f, 0.619970f, 0.620710f, 0.621440f, 0.622140f, 0.622820f,
        0.623480f, 0.624120f, 0.624740f, 0.625330f, 0.625900f, 0.626450f, 0.626980f, 0.627470f, 0.627950f, 0.628400f,
        0.628820f, 0.629220f, 0.629590f, 0.629930f, 0.630810f, 0.631320f, 0.631830f, 0.632330f, 0.632840f, 0.633340f,
        0.633850f, 0.634350f, 0.634860f, 0.635360f, 0.635870f, 0.636370f, 0.636870f, 0.637380f, 0.637880f, 0.638380f,
        0.638880f, 0.639380f, 0.639890f, 0.640390f, 0.640890f, 0.641390f, 0.641890f, 0.642390f, 0.642890f, 0.643390f,
        0.643890f, 0.644390f, 0.644890f, 0.645380f, 0.645880f, 0.646380f, 0.646880f, 0.647370f, 0.647870f, 0.648370f,
        0.648860f, 0.649360f, 0.649850f, 0.650350f, 0.650840f, 0.651340f, 0.651830f, 0.652320f, 0.652820f, 0.653310f,
        0.653800f, 0.654300f, 0.654790f, 0.655280f, 0.655770f, 0.656260f, 0.656750f, 0.657240f, 0.657730f, 0.658220f,
        0.658710f, 0.659200f, 0.659690f, 0.660180f, 0.660670f, 0.661150f, 0.661640f, 0.662130f, 0.662620f, 0.663100f,
        0.663590f, 0.664070f, 0.664560f, 0.665040f, 0.665530f, 0.666010f, 0.666500f, 0.666980f, 0.667460f, 0.667950f,
        0.668430f, 0.668910f, 0.669390f});
    static constexpr Spectrum CES16({0.034659f, 0.034799f, 0.034941f, 0.035082f, 0.035224f, 0.035367f, 0.035510f,
        0.035654f, 0.035799f, 0.035943f, 0.036089f, 0.036235f, 0.036382f, 0.036529f, 0.036677f, 0.036825f, 0.036974f,
        0.037124f, 0.037274f, 0.037424f, 0.037670f, 0.037778f, 0.037895f, 0.038020f, 0.038153f, 0.038294f, 0.038442f,
        0.038597f, 0.038759f, 0.038926f, 0.039099f, 0.039277f, 0.039460f, 0.039647f, 0.039837f, 0.040031f, 0.040229f,
        0.040428f, 0.040630f, 0.040834f, 0.041039f, 0.041245f, 0.041454f, 0.041668f, 0.041890f, 0.042121f, 0.042363f,
        0.042619f, 0.042890f, 0.043180f, 0.043489f, 0.043819f, 0.044169f, 0.044536f, 0.044917f, 0.045310f, 0.045713f,
        0.046123f, 0.046537f, 0.046953f, 0.047368f, 0.047781f, 0.048193f, 0.048605f, 0.049020f, 0.049438f, 0.049863f,
        0.050295f, 0.050737f, 0.051190f, 0.051656f, 0.052136f, 0.052631f, 0.053138f, 0.053657f, 0.054187f, 0.054727f,
        0.055275f, 0.055832f, 0.056395f, 0.056964f, 0.057538f, 0.058118f, 0.058702f, 0.059293f, 0.059891f, 0.060495f,
        0.061106f, 0.061725f, 0.062352f, 0.062987f, 0.063631f, 0.064284f, 0.064944f, 0.065610f, 0.066284f, 0.066963f,
        0.067647f, 0.068336f, 0.069029f, 0.069725f, 0.070424f, 0.071124f, 0.071825f, 0.072526f, 0.073226f, 0.073923f,
        0.074616f, 0.075306f, 0.075990f, 0.076667f, 0.077337f, 0.078001f, 0.078659f, 0.079313f, 0.079963f, 0.080611f,
        0.081258f, 0.081904f, 0.082551f, 0.083200f, 0.083852f, 0.084505f, 0.085161f, 0.085816f, 0.086473f, 0.087128f,
        0.087783f, 0.088436f, 0.089086f, 0.089734f, 0.090377f, 0.091016f, 0.091648f, 0.092271f, 0.092885f, 0.093488f,
        0.094078f, 0.094654f, 0.095214f, 0.095757f, 0.096282f, 0.096791f, 0.097285f, 0.097768f, 0.098242f, 0.098709f,
        0.099172f, 0.099632f, 0.100090f, 0.100560f, 0.101020f, 0.101490f, 0.101970f, 0.102450f, 0.102930f, 0.103410f,
        0.103900f, 0.104380f, 0.104870f, 0.105350f, 0.105840f, 0.106320f, 0.106810f, 0.107290f, 0.107770f, 0.108250f,
        0.108730f, 0.109200f, 0.109680f, 0.110150f, 0.110620f, 0.111090f, 0.111560f, 0.112030f, 0.112500f, 0.112970f,
        0.113440f, 0.113910f, 0.114380f, 0.114850f, 0.115320f, 0.115800f, 0.116290f, 0.116790f, 0.117300f, 0.117840f,
        0.118400f, 0.118990f, 0.119600f, 0.120260f, 0.120950f, 0.121680f, 0.122460f, 0.123280f, 0.124150f, 0.125070f,
        0.126050f, 0.127090f, 0.128180f, 0.129340f, 0.130570f, 0.131860f, 0.133230f, 0.134680f, 0.136220f, 0.137830f,
        0.139550f, 0.141350f, 0.143260f, 0.145270f, 0.147390f, 0.149610f, 0.151950f, 0.154400f, 0.156960f, 0.159620f,
        0.162400f, 0.165290f, 0.168290f, 0.171400f, 0.174630f, 0.177960f, 0.181390f, 0.184910f, 0.188520f, 0.192210f,
        0.195970f, 0.199800f, 0.203700f, 0.207640f, 0.211640f, 0.215690f, 0.219780f, 0.223920f, 0.228100f, 0.232330f,
        0.236610f, 0.240930f, 0.245290f, 0.249700f, 0.254150f, 0.258640f, 0.263150f, 0.267680f, 0.272230f, 0.276770f,
        0.281310f, 0.285840f, 0.290350f, 0.294830f, 0.299270f, 0.303670f, 0.308050f, 0.312380f, 0.316690f, 0.320970f,
        0.325220f, 0.329440f, 0.333630f, 0.337800f, 0.341950f, 0.346080f, 0.350200f, 0.354300f, 0.358410f, 0.362510f,
        0.366620f, 0.370750f, 0.374890f, 0.379050f, 0.383230f, 0.387430f, 0.391640f, 0.395840f, 0.400030f, 0.404190f,
        0.408320f, 0.412400f, 0.416430f, 0.420390f, 0.424280f, 0.428090f, 0.431840f, 0.435530f, 0.439170f, 0.442760f,
        0.446300f, 0.449810f, 0.453290f, 0.456730f, 0.460160f, 0.463570f, 0.466960f, 0.470340f, 0.473710f, 0.477060f,
        0.480410f, 0.483750f, 0.487090f, 0.490420f, 0.493760f, 0.497090f, 0.500410f, 0.503730f, 0.507030f, 0.510310f,
        0.513570f, 0.516810f, 0.520020f, 0.523190f, 0.526330f, 0.529430f, 0.532490f, 0.535500f, 0.538460f, 0.541370f,
        0.544220f, 0.547010f, 0.549730f, 0.552390f, 0.556110f, 0.559100f, 0.562080f, 0.565060f, 0.568040f, 0.571010f,
        0.573970f, 0.576930f, 0.579890f, 0.582830f, 0.585780f, 0.588710f, 0.591640f, 0.594570f, 0.597480f, 0.600390f,
        0.603290f, 0.606190f, 0.609080f, 0.611960f, 0.614830f, 0.617700f, 0.620550f, 0.623400f, 0.626240f, 0.629070f,
        0.631890f, 0.634710f, 0.637510f, 0.640300f, 0.643090f, 0.645870f, 0.648630f, 0.651390f, 0.654130f, 0.656870f,
        0.659590f, 0.662310f, 0.665010f, 0.667700f, 0.670390f, 0.673060f, 0.675720f, 0.678370f, 0.681000f, 0.683630f,
        0.686240f, 0.688850f, 0.691440f, 0.694010f, 0.696580f, 0.699140f, 0.701680f, 0.704210f, 0.706720f, 0.709230f,
        0.711720f, 0.714200f, 0.716660f, 0.719120f, 0.721560f, 0.723990f, 0.726400f, 0.728800f, 0.731190f, 0.733560f,
        0.735920f, 0.738270f, 0.740600f, 0.742920f, 0.745230f, 0.747520f, 0.749800f, 0.752070f, 0.754320f, 0.756560f,
        0.758780f, 0.760990f, 0.763190f});
    static constexpr Spectrum CES17({0.030798f, 0.031245f, 0.031698f, 0.032158f, 0.032624f, 0.033096f, 0.033576f,
        0.034062f, 0.034554f, 0.035054f, 0.035560f, 0.036074f, 0.036595f, 0.037123f, 0.037658f, 0.038201f, 0.038751f,
        0.039308f, 0.039874f, 0.040447f, 0.040310f, 0.041225f, 0.042094f, 0.042917f, 0.043696f, 0.044432f, 0.045126f,
        0.045779f, 0.046393f, 0.046967f, 0.047505f, 0.048005f, 0.048470f, 0.048902f, 0.049299f, 0.049665f, 0.050000f,
        0.050305f, 0.050582f, 0.050830f, 0.051053f, 0.051250f, 0.051426f, 0.051584f, 0.051729f, 0.051864f, 0.051993f,
        0.052120f, 0.052250f, 0.052385f, 0.052531f, 0.052690f, 0.052862f, 0.053048f, 0.053247f, 0.053458f, 0.053682f,
        0.053919f, 0.054167f, 0.054428f, 0.054699f, 0.054982f, 0.055276f, 0.055582f, 0.055899f, 0.056228f, 0.056569f,
        0.056921f, 0.057285f, 0.057662f, 0.058050f, 0.058451f, 0.058863f, 0.059286f, 0.059719f, 0.060161f, 0.060612f,
        0.061070f, 0.061536f, 0.062008f, 0.062485f, 0.062967f, 0.063453f, 0.063941f, 0.064430f, 0.064919f, 0.065407f,
        0.065891f, 0.066371f, 0.066846f, 0.067315f, 0.067775f, 0.068227f, 0.068669f, 0.069102f, 0.069524f, 0.069935f,
        0.070334f, 0.070721f, 0.071094f, 0.071454f, 0.071799f, 0.072131f, 0.072451f, 0.072758f, 0.073055f, 0.073342f,
        0.073620f, 0.073891f, 0.074154f, 0.074411f, 0.074662f, 0.074907f, 0.075144f, 0.075372f, 0.075590f, 0.075796f,
        0.075990f, 0.076169f, 0.076333f, 0.076480f, 0.076610f, 0.076726f, 0.076829f, 0.076923f, 0.077011f, 0.077096f,
        0.077181f, 0.077269f, 0.077363f, 0.077466f, 0.077580f, 0.077707f, 0.077847f, 0.078000f, 0.078168f, 0.078350f,
        0.078548f, 0.078762f, 0.078992f, 0.079240f, 0.079506f, 0.079796f, 0.080113f, 0.080465f, 0.080855f, 0.081290f,
        0.081774f, 0.082313f, 0.082912f, 0.083576f, 0.084312f, 0.085126f, 0.086024f, 0.087015f, 0.088106f, 0.089303f,
        0.090615f, 0.092047f, 0.093608f, 0.095305f, 0.097142f, 0.099116f, 0.101220f, 0.103450f, 0.105800f, 0.108270f,
        0.110840f, 0.113520f, 0.116290f, 0.119160f, 0.122100f, 0.125120f, 0.128210f, 0.131330f, 0.134500f, 0.137680f,
        0.140880f, 0.144070f, 0.147250f, 0.150400f, 0.153510f, 0.156570f, 0.159580f, 0.162530f, 0.165410f, 0.168220f,
        0.170950f, 0.173590f, 0.176140f, 0.178590f, 0.180930f, 0.183170f, 0.185300f, 0.187340f, 0.189280f, 0.191130f,
        0.192880f, 0.194540f, 0.196120f, 0.197610f, 0.199010f, 0.200350f, 0.201610f, 0.202810f, 0.203960f, 0.205060f,
        0.206120f, 0.207150f, 0.208150f, 0.209140f, 0.210110f, 0.211070f, 0.212000f, 0.212920f, 0.213810f, 0.214670f,
        0.215490f, 0.216270f, 0.217020f, 0.217710f, 0.218360f, 0.218950f, 0.219500f, 0.219980f, 0.220410f, 0.220790f,
        0.221100f, 0.221350f, 0.221530f, 0.221660f, 0.221710f, 0.221710f, 0.221660f, 0.221580f, 0.221460f, 0.221320f,
        0.221180f, 0.221030f, 0.220890f, 0.220770f, 0.220670f, 0.220620f, 0.220620f, 0.220680f, 0.220810f, 0.221030f,
        0.221350f, 0.221790f, 0.222340f, 0.223030f, 0.223870f, 0.224860f, 0.226010f, 0.227330f, 0.228810f, 0.230480f,
        0.232320f, 0.234350f, 0.236580f, 0.239000f, 0.241630f, 0.244480f, 0.247570f, 0.250910f, 0.254510f, 0.258390f,
        0.262570f, 0.267070f, 0.271880f, 0.277040f, 0.282560f, 0.288420f, 0.294620f, 0.301150f, 0.308000f, 0.315160f,
        0.322620f, 0.330370f, 0.338410f, 0.346720f, 0.355300f, 0.364090f, 0.373060f, 0.382170f, 0.391380f, 0.400630f,
        0.409890f, 0.419120f, 0.428270f, 0.437300f, 0.446170f, 0.454890f, 0.463450f, 0.471880f, 0.480160f, 0.488320f,
        0.496360f, 0.504270f, 0.512080f, 0.519790f, 0.527400f, 0.534900f, 0.542280f, 0.549540f, 0.556660f, 0.563620f,
        0.570430f, 0.577060f, 0.583510f, 0.589770f, 0.595820f, 0.601650f, 0.607260f, 0.612630f, 0.617750f, 0.622610f,
        0.627200f, 0.631510f, 0.635530f, 0.639240f, 0.647780f, 0.652950f, 0.658080f, 0.663180f, 0.668230f, 0.673250f,
        0.678230f, 0.683170f, 0.688070f, 0.692920f, 0.697740f, 0.702510f, 0.707240f, 0.711920f, 0.716560f, 0.721150f,
        0.725700f, 0.730200f, 0.734650f, 0.739050f, 0.743410f, 0.747720f, 0.751990f, 0.756200f, 0.760370f, 0.764480f,
        0.768550f, 0.772570f, 0.776530f, 0.780450f, 0.784320f, 0.788140f, 0.791910f, 0.795630f, 0.799300f, 0.802920f,
        0.806490f, 0.810010f, 0.813490f, 0.816910f, 0.820280f, 0.823610f, 0.826890f, 0.830110f, 0.833300f, 0.836430f,
        0.839510f, 0.842550f, 0.845540f, 0.848490f, 0.851390f, 0.854240f, 0.857040f, 0.859810f, 0.862520f, 0.865200f,
        0.867820f, 0.870410f, 0.872950f, 0.875450f, 0.877910f, 0.880320f, 0.882690f, 0.885030f, 0.887320f, 0.889570f,
        0.891780f, 0.893960f, 0.896090f, 0.898190f, 0.900250f, 0.902270f, 0.904260f, 0.906210f, 0.908120f, 0.910000f,
        0.911840f, 0.913650f, 0.915430f});
    static constexpr Spectrum CES18({0.064984f, 0.065366f, 0.065750f, 0.066136f, 0.066524f, 0.066915f, 0.067307f,
        0.067702f, 0.068099f, 0.068497f, 0.068899f, 0.069302f, 0.069707f, 0.070115f, 0.070525f, 0.070937f, 0.071351f,
        0.071767f, 0.072186f, 0.072607f, 0.073100f, 0.073417f, 0.073816f, 0.074285f, 0.074812f, 0.075384f, 0.075989f,
        0.076616f, 0.077251f, 0.077883f, 0.078500f, 0.079092f, 0.079662f, 0.080215f, 0.080758f, 0.081296f, 0.081833f,
        0.082377f, 0.082932f, 0.083505f, 0.084100f, 0.084723f, 0.085375f, 0.086058f, 0.086773f, 0.087521f, 0.088304f,
        0.089121f, 0.089976f, 0.090868f, 0.091800f, 0.092771f, 0.093777f, 0.094813f, 0.095875f, 0.096957f, 0.098055f,
        0.099163f, 0.100280f, 0.101390f, 0.102500f, 0.103600f, 0.104690f, 0.105770f, 0.106840f, 0.107890f, 0.108920f,
        0.109940f, 0.110950f, 0.111930f, 0.112900f, 0.113850f, 0.114770f, 0.115690f, 0.116580f, 0.117470f, 0.118350f,
        0.119220f, 0.120080f, 0.120940f, 0.121800f, 0.122660f, 0.123520f, 0.124380f, 0.125250f, 0.126110f, 0.126970f,
        0.127830f, 0.128690f, 0.129550f, 0.130400f, 0.131250f, 0.132100f, 0.132950f, 0.133810f, 0.134660f, 0.135520f,
        0.136380f, 0.137250f, 0.138120f, 0.139000f, 0.139890f, 0.140790f, 0.141700f, 0.142620f, 0.143550f, 0.144490f,
        0.145450f, 0.146420f, 0.147400f, 0.148400f, 0.149410f, 0.150440f, 0.151480f, 0.152530f, 0.153590f, 0.154660f,
        0.155740f, 0.156820f, 0.157910f, 0.159000f, 0.160090f, 0.161180f, 0.162250f, 0.163300f, 0.164330f, 0.165320f,
        0.166270f, 0.167170f, 0.168020f, 0.168800f, 0.169510f, 0.170150f, 0.170720f, 0.171230f, 0.171660f, 0.172020f,
        0.172320f, 0.172550f, 0.172710f, 0.172800f, 0.172830f, 0.172800f, 0.172720f, 0.172600f, 0.172460f, 0.172290f,
        0.172110f, 0.171930f, 0.171760f, 0.171600f, 0.171470f, 0.171370f, 0.171310f, 0.171290f, 0.171320f, 0.171420f,
        0.171570f, 0.171800f, 0.172110f, 0.172500f, 0.172980f, 0.173540f, 0.174170f, 0.174850f, 0.175590f, 0.176360f,
        0.177160f, 0.177970f, 0.178790f, 0.179600f, 0.180400f, 0.181170f, 0.181910f, 0.182620f, 0.183290f, 0.183910f,
        0.184480f, 0.184990f, 0.185430f, 0.185800f, 0.186100f, 0.186330f, 0.186530f, 0.186690f, 0.186850f, 0.187020f,
        0.187210f, 0.187440f, 0.187730f, 0.188100f, 0.188560f, 0.189130f, 0.189840f, 0.190690f, 0.191700f, 0.192900f,
        0.194300f, 0.195920f, 0.197780f, 0.199900f, 0.202280f, 0.204920f, 0.207780f, 0.210850f, 0.214090f, 0.217500f,
        0.221030f, 0.224680f, 0.228410f, 0.232200f, 0.236030f, 0.239890f, 0.243740f, 0.247590f, 0.251400f, 0.255170f,
        0.258870f, 0.262480f, 0.266000f, 0.269400f, 0.272670f, 0.275810f, 0.278820f, 0.281720f, 0.284500f, 0.287170f,
        0.289750f, 0.292220f, 0.294600f, 0.296900f, 0.299110f, 0.301250f, 0.303310f, 0.305290f, 0.307200f, 0.309030f,
        0.310800f, 0.312500f, 0.314130f, 0.315700f, 0.317210f, 0.318650f, 0.320050f, 0.321410f, 0.322720f, 0.324000f,
        0.325250f, 0.326480f, 0.327700f, 0.328900f, 0.330100f, 0.331290f, 0.332480f, 0.333670f, 0.334860f, 0.336050f,
        0.337230f, 0.338420f, 0.339610f, 0.340800f, 0.341990f, 0.343190f, 0.344380f, 0.345560f, 0.346740f, 0.347900f,
        0.349050f, 0.350190f, 0.351310f, 0.352400f, 0.353470f, 0.354520f, 0.355550f, 0.356570f, 0.357580f, 0.358580f,
        0.359580f, 0.360580f, 0.361590f, 0.362600f, 0.363630f, 0.364660f, 0.365710f, 0.366770f, 0.367830f, 0.368900f,
        0.369970f, 0.371050f, 0.372120f, 0.373200f, 0.374280f, 0.375350f, 0.376420f, 0.377480f, 0.378530f, 0.379570f,
        0.380600f, 0.381620f, 0.382620f, 0.383600f, 0.384560f, 0.385510f, 0.386440f, 0.387370f, 0.388290f, 0.389220f,
        0.390140f, 0.391080f, 0.392030f, 0.393000f, 0.393990f, 0.394980f, 0.395960f, 0.396910f, 0.397830f, 0.398680f,
        0.399460f, 0.400150f, 0.400740f, 0.401200f, 0.401940f, 0.402580f, 0.403210f, 0.403840f, 0.404480f, 0.405110f,
        0.405750f, 0.406380f, 0.407020f, 0.407650f, 0.408290f, 0.408930f, 0.409560f, 0.410200f, 0.410840f, 0.411470f,
        0.412110f, 0.412750f, 0.413390f, 0.414030f, 0.414670f, 0.415310f, 0.415950f, 0.416590f, 0.417230f, 0.417870f,
        0.418510f, 0.419150f, 0.419790f, 0.420430f, 0.421070f, 0.421710f, 0.422360f, 0.423000f, 0.423640f, 0.424290f,
        0.424930f, 0.425570f, 0.426220f, 0.426860f, 0.427500f, 0.428150f, 0.428790f, 0.429440f, 0.430080f, 0.430730f,
        0.431380f, 0.432020f, 0.432670f, 0.433320f, 0.433960f, 0.434610f, 0.435260f, 0.435900f, 0.436550f, 0.437200f,
        0.437850f, 0.438500f, 0.439140f, 0.439790f, 0.440440f, 0.441090f, 0.441740f, 0.442390f, 0.443040f, 0.443690f,
        0.444340f, 0.444990f, 0.445640f, 0.446290f, 0.446940f, 0.447590f, 0.448240f, 0.448890f, 0.449550f, 0.450200f,
        0.450850f, 0.451500f, 0.452150f});
    static constexpr Spectrum CES19({0.109270f, 0.109790f, 0.110320f, 0.110840f, 0.111370f, 0.111900f, 0.112430f,
        0.112970f, 0.113500f, 0.114040f, 0.114580f, 0.115130f, 0.115670f, 0.116220f, 0.116770f, 0.117330f, 0.117880f,
        0.118440f, 0.119000f, 0.119560f, 0.120200f, 0.120640f, 0.121190f, 0.121800f, 0.122470f, 0.123160f, 0.123840f,
        0.124480f, 0.125050f, 0.125540f, 0.125900f, 0.126120f, 0.126210f, 0.126190f, 0.126080f, 0.125890f, 0.125650f,
        0.125370f, 0.125070f, 0.124780f, 0.124500f, 0.124260f, 0.124050f, 0.123870f, 0.123720f, 0.123590f, 0.123490f,
        0.123400f, 0.123320f, 0.123260f, 0.123200f, 0.123150f, 0.123090f, 0.123050f, 0.123010f, 0.122980f, 0.122960f,
        0.122950f, 0.122950f, 0.122970f, 0.123000f, 0.123050f, 0.123100f, 0.123170f, 0.123250f, 0.123340f, 0.123430f,
        0.123520f, 0.123620f, 0.123710f, 0.123800f, 0.123880f, 0.123960f, 0.124040f, 0.124110f, 0.124170f, 0.124240f,
        0.124300f, 0.124370f, 0.124430f, 0.124500f, 0.124570f, 0.124640f, 0.124710f, 0.124770f, 0.124840f, 0.124900f,
        0.124960f, 0.125010f, 0.125060f, 0.125100f, 0.125130f, 0.125160f, 0.125170f, 0.125180f, 0.125180f, 0.125180f,
        0.125170f, 0.125150f, 0.125130f, 0.125100f, 0.125070f, 0.125040f, 0.125010f, 0.125000f, 0.125010f, 0.125050f,
        0.125120f, 0.125230f, 0.125390f, 0.125600f, 0.125870f, 0.126220f, 0.126680f, 0.127270f, 0.128010f, 0.128920f,
        0.130030f, 0.131370f, 0.132950f, 0.134800f, 0.136930f, 0.139340f, 0.142010f, 0.144920f, 0.148050f, 0.151380f,
        0.154900f, 0.158590f, 0.162430f, 0.166400f, 0.170490f, 0.174640f, 0.178830f, 0.183000f, 0.187110f, 0.191120f,
        0.194970f, 0.198630f, 0.202060f, 0.205200f, 0.208030f, 0.210560f, 0.212810f, 0.214830f, 0.216620f, 0.218230f,
        0.219670f, 0.220980f, 0.222180f, 0.223300f, 0.224370f, 0.225430f, 0.226540f, 0.227730f, 0.229050f, 0.230560f,
        0.232290f, 0.234290f, 0.236610f, 0.239300f, 0.242400f, 0.245930f, 0.249920f, 0.254390f, 0.259360f, 0.264850f,
        0.270890f, 0.277500f, 0.284690f, 0.292500f, 0.300920f, 0.309900f, 0.319340f, 0.329160f, 0.339280f, 0.349620f,
        0.360080f, 0.370590f, 0.381050f, 0.391400f, 0.401540f, 0.411440f, 0.421050f, 0.430330f, 0.439250f, 0.447750f,
        0.455810f, 0.463380f, 0.470430f, 0.476900f, 0.482780f, 0.488090f, 0.492860f, 0.497150f, 0.500970f, 0.504380f,
        0.507410f, 0.510100f, 0.512480f, 0.514600f, 0.516490f, 0.518160f, 0.519650f, 0.520950f, 0.522100f, 0.523100f,
        0.523970f, 0.524740f, 0.525410f, 0.526000f, 0.526530f, 0.527000f, 0.527420f, 0.527790f, 0.528100f, 0.528370f,
        0.528590f, 0.528770f, 0.528900f, 0.529000f, 0.529060f, 0.529090f, 0.529090f, 0.529070f, 0.529020f, 0.528960f,
        0.528880f, 0.528790f, 0.528700f, 0.528600f, 0.528500f, 0.528400f, 0.528300f, 0.528200f, 0.528100f, 0.527990f,
        0.527870f, 0.527760f, 0.527630f, 0.527500f, 0.527360f, 0.527220f, 0.527070f, 0.526910f, 0.526760f, 0.526600f,
        0.526440f, 0.526290f, 0.526140f, 0.526000f, 0.525860f, 0.525730f, 0.525600f, 0.525470f, 0.525340f, 0.525210f,
        0.525070f, 0.524920f, 0.524770f, 0.524600f, 0.524420f, 0.524230f, 0.524020f, 0.523820f, 0.523610f, 0.523400f,
        0.523190f, 0.522980f, 0.522790f, 0.522600f, 0.522430f, 0.522260f, 0.522110f, 0.521960f, 0.521810f, 0.521670f,
        0.521530f, 0.521390f, 0.521250f, 0.521100f, 0.520940f, 0.520780f, 0.520620f, 0.520450f, 0.520290f, 0.520120f,
        0.519960f, 0.519800f, 0.519640f, 0.519500f, 0.519360f, 0.519240f, 0.519110f, 0.518990f, 0.518870f, 0.518750f,
        0.518620f, 0.518490f, 0.518350f, 0.518200f, 0.518040f, 0.517860f, 0.517680f, 0.517490f, 0.517290f, 0.517080f,
        0.516860f, 0.516650f, 0.516420f, 0.516200f, 0.515980f, 0.515750f, 0.515530f, 0.515320f, 0.515120f, 0.514940f,
        0.514770f, 0.514620f, 0.514500f, 0.514400f, 0.514240f, 0.514110f, 0.513970f, 0.513840f, 0.513700f, 0.513570f,
        0.513430f, 0.513300f, 0.513160f, 0.513030f, 0.512890f, 0.512760f, 0.512620f, 0.512490f, 0.512350f, 0.512220f,
        0.512080f, 0.511950f, 0.511810f, 0.511680f, 0.511540f, 0.511410f, 0.511270f, 0.511140f, 0.511000f, 0.510870f,
        0.510730f, 0.510600f, 0.510460f, 0.510330f, 0.510190f, 0.510050f, 0.509920f, 0.509780f, 0.509650f, 0.509510f,
        0.509380f, 0.509240f, 0.509110f, 0.508970f, 0.508840f, 0.508700f, 0.508570f, 0.508430f, 0.508300f, 0.508160f,
        0.508030f, 0.507890f, 0.507760f, 0.507620f, 0.507490f, 0.507350f, 0.507220f, 0.507080f, 0.506950f, 0.506810f,
        0.506680f, 0.506540f, 0.506410f, 0.506270f, 0.506140f, 0.506000f, 0.505870f, 0.505730f, 0.505600f, 0.505460f,
        0.505330f, 0.505190f, 0.505060f, 0.504920f, 0.504790f, 0.504650f, 0.504520f, 0.504380f, 0.504250f, 0.504110f,
        0.503980f, 0.503840f, 0.503710f});
    static constexpr Spectrum CES20({0.025086f, 0.025460f, 0.025839f, 0.026223f, 0.026613f, 0.027009f, 0.027411f,
        0.027818f, 0.028231f, 0.028650f, 0.029075f, 0.029506f, 0.029944f, 0.030387f, 0.030837f, 0.031294f, 0.031757f,
        0.032227f, 0.032703f, 0.033186f, 0.031971f, 0.033232f, 0.034378f, 0.035412f, 0.036332f, 0.037141f, 0.037837f,
        0.038423f, 0.038898f, 0.039262f, 0.039517f, 0.039663f, 0.039700f, 0.039629f, 0.039450f, 0.039164f, 0.038771f,
        0.038271f, 0.037667f, 0.036956f, 0.036142f, 0.035226f, 0.034227f, 0.033165f, 0.032061f, 0.030935f, 0.029810f,
        0.028704f, 0.027639f, 0.026637f, 0.025716f, 0.024894f, 0.024164f, 0.023515f, 0.022937f, 0.022420f, 0.021951f,
        0.021521f, 0.021119f, 0.020734f, 0.020354f, 0.019974f, 0.019598f, 0.019235f, 0.018897f, 0.018590f, 0.018326f,
        0.018112f, 0.017959f, 0.017876f, 0.017872f, 0.017957f, 0.018143f, 0.018441f, 0.018865f, 0.019426f, 0.020137f,
        0.021010f, 0.022057f, 0.023291f, 0.024723f, 0.026371f, 0.028264f, 0.030439f, 0.032930f, 0.035774f, 0.039004f,
        0.042656f, 0.046765f, 0.051367f, 0.056496f, 0.062177f, 0.068393f, 0.075115f, 0.082313f, 0.089961f, 0.098030f,
        0.106490f, 0.115320f, 0.124480f, 0.133940f, 0.143670f, 0.153550f, 0.163470f, 0.173280f, 0.182880f, 0.192130f,
        0.200920f, 0.209110f, 0.216580f, 0.223200f, 0.228890f, 0.233680f, 0.237620f, 0.240770f, 0.243200f, 0.244970f,
        0.246140f, 0.246760f, 0.246910f, 0.246640f, 0.246000f, 0.245050f, 0.243830f, 0.242360f, 0.240700f, 0.238890f,
        0.236950f, 0.234940f, 0.232890f, 0.230850f, 0.228840f, 0.226880f, 0.224950f, 0.223060f, 0.221220f, 0.219400f,
        0.217620f, 0.215880f, 0.214160f, 0.212480f, 0.210820f, 0.209200f, 0.207610f, 0.206060f, 0.204550f, 0.203080f,
        0.201650f, 0.200280f, 0.198950f, 0.197690f, 0.196480f, 0.195330f, 0.194260f, 0.193260f, 0.192340f, 0.191520f,
        0.190790f, 0.190170f, 0.189650f, 0.189250f, 0.188970f, 0.188810f, 0.188780f, 0.188890f, 0.189130f, 0.189520f,
        0.190040f, 0.190720f, 0.191540f, 0.192520f, 0.193660f, 0.194950f, 0.196400f, 0.198000f, 0.199740f, 0.201630f,
        0.203670f, 0.205840f, 0.208150f, 0.210590f, 0.213170f, 0.215880f, 0.218700f, 0.221640f, 0.224700f, 0.227860f,
        0.231120f, 0.234480f, 0.237940f, 0.241470f, 0.245090f, 0.248810f, 0.252630f, 0.256560f, 0.260620f, 0.264810f,
        0.269150f, 0.273650f, 0.278320f, 0.283170f, 0.288220f, 0.293480f, 0.298970f, 0.304720f, 0.310750f, 0.317080f,
        0.323730f, 0.330730f, 0.338080f, 0.345830f, 0.353970f, 0.362520f, 0.371450f, 0.380770f, 0.390450f, 0.400500f,
        0.410900f, 0.421640f, 0.432720f, 0.444120f, 0.455830f, 0.467800f, 0.479950f, 0.492230f, 0.504580f, 0.516930f,
        0.529210f, 0.541370f, 0.553340f, 0.565060f, 0.576470f, 0.587540f, 0.598260f, 0.608590f, 0.618520f, 0.628010f,
        0.637060f, 0.645630f, 0.653710f, 0.661270f, 0.668300f, 0.674820f, 0.680860f, 0.686460f, 0.691660f, 0.696480f,
        0.700960f, 0.705140f, 0.709040f, 0.712700f, 0.716160f, 0.719420f, 0.722500f, 0.725420f, 0.728180f, 0.730810f,
        0.733310f, 0.735700f, 0.738000f, 0.740210f, 0.742340f, 0.744420f, 0.746440f, 0.748410f, 0.750350f, 0.752260f,
        0.754150f, 0.756020f, 0.757890f, 0.759770f, 0.761650f, 0.763550f, 0.765480f, 0.767420f, 0.769390f, 0.771380f,
        0.773410f, 0.775480f, 0.777580f, 0.779720f, 0.781910f, 0.784140f, 0.786390f, 0.788680f, 0.790980f, 0.793290f,
        0.795620f, 0.797940f, 0.800260f, 0.802560f, 0.804850f, 0.807120f, 0.809360f, 0.811590f, 0.813790f, 0.815960f,
        0.818120f, 0.820240f, 0.822340f, 0.824400f, 0.826440f, 0.828440f, 0.830410f, 0.832340f, 0.834230f, 0.836080f,
        0.837880f, 0.839630f, 0.841330f, 0.842970f, 0.844560f, 0.846090f, 0.847560f, 0.848960f, 0.850300f, 0.851560f,
        0.852760f, 0.853880f, 0.854920f, 0.855880f, 0.858060f, 0.859390f, 0.860710f, 0.862010f, 0.863310f, 0.864600f,
        0.865870f, 0.867140f, 0.868390f, 0.869640f, 0.870880f, 0.872100f, 0.873320f, 0.874520f, 0.875720f, 0.876910f,
        0.878080f, 0.879250f, 0.880410f, 0.881560f, 0.882690f, 0.883820f, 0.884940f, 0.886050f, 0.887150f, 0.888240f,
        0.889320f, 0.890400f, 0.891460f, 0.892520f, 0.893560f, 0.894600f, 0.895630f, 0.896640f, 0.897650f, 0.898660f,
        0.899650f, 0.900630f, 0.901610f, 0.902570f, 0.903530f, 0.904480f, 0.905420f, 0.906360f, 0.907280f, 0.908200f,
        0.909110f, 0.910010f, 0.910900f, 0.911780f, 0.912660f, 0.913530f, 0.914390f, 0.915240f, 0.916090f, 0.916930f,
        0.917760f, 0.918580f, 0.919390f, 0.920200f, 0.921000f, 0.921790f, 0.922580f, 0.923360f, 0.924130f, 0.924890f,
        0.925650f, 0.926400f, 0.927140f, 0.927880f, 0.928610f, 0.929330f, 0.930040f, 0.930750f, 0.931460f, 0.932150f,
        0.932840f, 0.933520f, 0.934200f});
    static constexpr Spectrum CES21({0.123560f, 0.124960f, 0.126380f, 0.127810f, 0.129260f, 0.130710f, 0.132190f,
        0.133670f, 0.135170f, 0.136690f, 0.138220f, 0.139760f, 0.141320f, 0.142890f, 0.144480f, 0.146080f, 0.147690f,
        0.149320f, 0.150970f, 0.152630f, 0.148740f, 0.153090f, 0.156940f, 0.160320f, 0.163260f, 0.165780f, 0.167910f,
        0.169690f, 0.171140f, 0.172290f, 0.173160f, 0.173790f, 0.174210f, 0.174430f, 0.174500f, 0.174440f, 0.174270f,
        0.174030f, 0.173750f, 0.173450f, 0.173160f, 0.172910f, 0.172690f, 0.172500f, 0.172340f, 0.172200f, 0.172080f,
        0.171970f, 0.171870f, 0.171770f, 0.171670f, 0.171570f, 0.171470f, 0.171360f, 0.171250f, 0.171140f, 0.171040f,
        0.170940f, 0.170850f, 0.170760f, 0.170680f, 0.170610f, 0.170550f, 0.170510f, 0.170480f, 0.170470f, 0.170470f,
        0.170490f, 0.170530f, 0.170590f, 0.170680f, 0.170790f, 0.170920f, 0.171090f, 0.171280f, 0.171500f, 0.171760f,
        0.172050f, 0.172380f, 0.172750f, 0.173160f, 0.173620f, 0.174110f, 0.174650f, 0.175230f, 0.175860f, 0.176530f,
        0.177230f, 0.177990f, 0.178780f, 0.179620f, 0.180490f, 0.181400f, 0.182350f, 0.183320f, 0.184300f, 0.185310f,
        0.186320f, 0.187330f, 0.188340f, 0.189350f, 0.190340f, 0.191320f, 0.192300f, 0.193270f, 0.194250f, 0.195230f,
        0.196220f, 0.197220f, 0.198240f, 0.199270f, 0.200340f, 0.201430f, 0.202560f, 0.203750f, 0.204990f, 0.206300f,
        0.207690f, 0.209150f, 0.210720f, 0.212380f, 0.214160f, 0.216070f, 0.218130f, 0.220360f, 0.222780f, 0.225420f,
        0.228280f, 0.231410f, 0.234800f, 0.238490f, 0.242500f, 0.246820f, 0.251460f, 0.256440f, 0.261740f, 0.267380f,
        0.273370f, 0.279690f, 0.286370f, 0.293400f, 0.300790f, 0.308510f, 0.316560f, 0.324900f, 0.333530f, 0.342430f,
        0.351570f, 0.360940f, 0.370510f, 0.380280f, 0.390220f, 0.400300f, 0.410520f, 0.420850f, 0.431260f, 0.441740f,
        0.452270f, 0.462830f, 0.473390f, 0.483940f, 0.494450f, 0.504910f, 0.515300f, 0.525590f, 0.535770f, 0.545820f,
        0.555720f, 0.565450f, 0.574990f, 0.584320f, 0.593430f, 0.602330f, 0.611040f, 0.619560f, 0.627910f, 0.636100f,
        0.644160f, 0.652090f, 0.659900f, 0.667630f, 0.675260f, 0.682800f, 0.690240f, 0.697570f, 0.704760f, 0.711830f,
        0.718750f, 0.725510f, 0.732100f, 0.738520f, 0.744750f, 0.750800f, 0.756670f, 0.762370f, 0.767920f, 0.773300f,
        0.778540f, 0.783640f, 0.788600f, 0.793430f, 0.798130f, 0.802700f, 0.807130f, 0.811410f, 0.815530f, 0.819470f,
        0.823240f, 0.826810f, 0.830180f, 0.833340f, 0.836280f, 0.839030f, 0.841600f, 0.844010f, 0.846280f, 0.848440f,
        0.850510f, 0.852510f, 0.854460f, 0.856380f, 0.858280f, 0.860170f, 0.862030f, 0.863860f, 0.865650f, 0.867370f,
        0.869040f, 0.870630f, 0.872140f, 0.873550f, 0.874870f, 0.876100f, 0.877240f, 0.878310f, 0.879320f, 0.880280f,
        0.881190f, 0.882070f, 0.882930f, 0.883780f, 0.884620f, 0.885460f, 0.886280f, 0.887090f, 0.887890f, 0.888670f,
        0.889430f, 0.890150f, 0.890860f, 0.891520f, 0.892160f, 0.892760f, 0.893340f, 0.893890f, 0.894420f, 0.894930f,
        0.895430f, 0.895920f, 0.896400f, 0.896890f, 0.897370f, 0.897850f, 0.898330f, 0.898810f, 0.899280f, 0.899740f,
        0.900190f, 0.900630f, 0.901050f, 0.901450f, 0.901840f, 0.902210f, 0.902550f, 0.902880f, 0.903190f, 0.903480f,
        0.903750f, 0.903990f, 0.904220f, 0.904430f, 0.904620f, 0.904790f, 0.904950f, 0.905090f, 0.905210f, 0.905330f,
        0.905440f, 0.905540f, 0.905630f, 0.905720f, 0.905810f, 0.905900f, 0.905990f, 0.906080f, 0.906170f, 0.906270f,
        0.906370f, 0.906480f, 0.906590f, 0.906720f, 0.906850f, 0.906990f, 0.907140f, 0.907290f, 0.907450f, 0.907620f,
        0.907780f, 0.907960f, 0.908130f, 0.908300f, 0.908480f, 0.908650f, 0.908820f, 0.908990f, 0.909160f, 0.909320f,
        0.909470f, 0.909620f, 0.909760f, 0.909890f, 0.910110f, 0.910270f, 0.910440f, 0.910600f, 0.910760f, 0.910920f,
        0.911080f, 0.911250f, 0.911410f, 0.911570f, 0.911730f, 0.911890f, 0.912050f, 0.912210f, 0.912370f, 0.912520f,
        0.912680f, 0.912840f, 0.913000f, 0.913160f, 0.913320f, 0.913470f, 0.913630f, 0.913790f, 0.913940f, 0.914100f,
        0.914260f, 0.914410f, 0.914570f, 0.914720f, 0.914880f, 0.915030f, 0.915190f, 0.915340f, 0.915500f, 0.915650f,
        0.915810f, 0.915960f, 0.916110f, 0.916260f, 0.916420f, 0.916570f, 0.916720f, 0.916870f, 0.917030f, 0.917180f,
        0.917330f, 0.917480f, 0.917630f, 0.917780f, 0.917930f, 0.918080f, 0.918230f, 0.918380f, 0.918530f, 0.918680f,
        0.918830f, 0.918970f, 0.919120f, 0.919270f, 0.919420f, 0.919570f, 0.919710f, 0.919860f, 0.920010f, 0.920150f,
        0.920300f, 0.920440f, 0.920590f, 0.920740f, 0.920880f, 0.921030f, 0.921170f, 0.921310f, 0.921460f, 0.921600f,
        0.921750f, 0.921890f, 0.922030f});
    static constexpr Spectrum CES22({0.000123f, 0.002255f, 0.004368f, 0.006464f, 0.008546f, 0.010615f, 0.012673f,
        0.014723f, 0.016767f, 0.018807f, 0.020844f, 0.022881f, 0.024921f, 0.026965f, 0.029015f, 0.031073f, 0.033142f,
        0.035224f, 0.037320f, 0.039433f, 0.041565f, 0.043716f, 0.045877f, 0.048039f, 0.050190f, 0.052320f, 0.054418f,
        0.056474f, 0.058478f, 0.060419f, 0.062286f, 0.064071f, 0.065769f, 0.067381f, 0.068903f, 0.070335f, 0.071674f,
        0.072919f, 0.074068f, 0.075120f, 0.076073f, 0.076928f, 0.077695f, 0.078389f, 0.079023f, 0.079611f, 0.080166f,
        0.080702f, 0.081233f, 0.081772f, 0.082333f, 0.082927f, 0.083556f, 0.084219f, 0.084915f, 0.085644f, 0.086405f,
        0.087197f, 0.088020f, 0.088873f, 0.089755f, 0.090665f, 0.091602f, 0.092565f, 0.093551f, 0.094558f, 0.095586f,
        0.096632f, 0.097695f, 0.098773f, 0.099864f, 0.100970f, 0.102090f, 0.103220f, 0.104370f, 0.105540f, 0.106730f,
        0.107950f, 0.109190f, 0.110470f, 0.111770f, 0.113100f, 0.114470f, 0.115870f, 0.117300f, 0.118770f, 0.120270f,
        0.121800f, 0.123370f, 0.124980f, 0.126620f, 0.128300f, 0.130000f, 0.131740f, 0.133500f, 0.135290f, 0.137090f,
        0.138910f, 0.140730f, 0.142560f, 0.144400f, 0.146240f, 0.148080f, 0.149920f, 0.151780f, 0.153650f, 0.155550f,
        0.157460f, 0.159400f, 0.161370f, 0.163370f, 0.165410f, 0.167490f, 0.169610f, 0.171770f, 0.173970f, 0.176220f,
        0.178510f, 0.180860f, 0.183250f, 0.185690f, 0.188190f, 0.190740f, 0.193350f, 0.196020f, 0.198760f, 0.201570f,
        0.204450f, 0.207400f, 0.210440f, 0.213560f, 0.216760f, 0.220050f, 0.223420f, 0.226870f, 0.230390f, 0.233990f,
        0.237660f, 0.241390f, 0.245200f, 0.249060f, 0.252990f, 0.256970f, 0.261030f, 0.265160f, 0.269350f, 0.273630f,
        0.277980f, 0.282410f, 0.286930f, 0.291540f, 0.296240f, 0.301030f, 0.305910f, 0.310880f, 0.315940f, 0.321090f,
        0.326330f, 0.331660f, 0.337070f, 0.342580f, 0.348170f, 0.353850f, 0.359600f, 0.365420f, 0.371310f, 0.377240f,
        0.383230f, 0.389260f, 0.395330f, 0.401420f, 0.407540f, 0.413690f, 0.419870f, 0.426080f, 0.432350f, 0.438660f,
        0.445020f, 0.451450f, 0.457940f, 0.464500f, 0.471130f, 0.477820f, 0.484560f, 0.491310f, 0.498080f, 0.504830f,
        0.511560f, 0.518240f, 0.524860f, 0.531410f, 0.537860f, 0.544230f, 0.550500f, 0.556700f, 0.562810f, 0.568860f,
        0.574830f, 0.580730f, 0.586570f, 0.592360f, 0.598090f, 0.603760f, 0.609360f, 0.614890f, 0.620340f, 0.625710f,
        0.630990f, 0.636160f, 0.641240f, 0.646200f, 0.651050f, 0.655780f, 0.660400f, 0.664900f, 0.669280f, 0.673540f,
        0.677700f, 0.681730f, 0.685650f, 0.689450f, 0.693140f, 0.696710f, 0.700180f, 0.703550f, 0.706820f, 0.710010f,
        0.713110f, 0.716120f, 0.719070f, 0.721940f, 0.724750f, 0.727490f, 0.730160f, 0.732770f, 0.735310f, 0.737780f,
        0.740180f, 0.742510f, 0.744770f, 0.746950f, 0.749060f, 0.751100f, 0.753080f, 0.755010f, 0.756890f, 0.758720f,
        0.760520f, 0.762280f, 0.764030f, 0.765750f, 0.767460f, 0.769150f, 0.770820f, 0.772460f, 0.774070f, 0.775640f,
        0.777150f, 0.778620f, 0.780030f, 0.781370f, 0.782650f, 0.783860f, 0.785020f, 0.786130f, 0.787200f, 0.788240f,
        0.789250f, 0.790240f, 0.791230f, 0.792210f, 0.793200f, 0.794190f, 0.795190f, 0.796200f, 0.797230f, 0.798270f,
        0.799320f, 0.800400f, 0.801500f, 0.802620f, 0.803760f, 0.804930f, 0.806090f, 0.807260f, 0.808420f, 0.809560f,
        0.810670f, 0.811750f, 0.812770f, 0.813750f, 0.814670f, 0.815520f, 0.816340f, 0.817110f, 0.817850f, 0.818560f,
        0.819260f, 0.819950f, 0.820640f, 0.821330f, 0.822030f, 0.822750f, 0.823480f, 0.824210f, 0.824950f, 0.825680f,
        0.826420f, 0.827160f, 0.827890f, 0.828620f, 0.829340f, 0.830050f, 0.830760f, 0.831470f, 0.832180f, 0.832900f,
        0.833620f, 0.834360f, 0.835100f, 0.835870f, 0.836650f, 0.837450f, 0.838240f, 0.839020f, 0.839790f, 0.840530f,
        0.841240f, 0.841900f, 0.842510f, 0.843050f, 0.843520f, 0.843930f, 0.844290f, 0.844590f, 0.844860f, 0.845100f,
        0.845310f, 0.845500f, 0.845690f, 0.845880f, 0.846080f, 0.846290f, 0.846530f, 0.846810f, 0.847120f, 0.847490f,
        0.847910f, 0.848400f, 0.848960f, 0.849600f, 0.849240f, 0.849540f, 0.849830f, 0.850130f, 0.850430f, 0.850720f,
        0.851020f, 0.851310f, 0.851600f, 0.851900f, 0.852190f, 0.852480f, 0.852770f, 0.853060f, 0.853350f, 0.853640f,
        0.853930f, 0.854220f, 0.854510f, 0.854800f, 0.855090f, 0.855380f, 0.855660f, 0.855950f, 0.856240f, 0.856520f,
        0.856810f, 0.857090f, 0.857370f, 0.857660f, 0.857940f, 0.858220f, 0.858510f, 0.858790f, 0.859070f, 0.859350f,
        0.859630f, 0.859910f, 0.860190f, 0.860470f, 0.860750f, 0.861020f, 0.861300f, 0.861580f, 0.861860f, 0.862130f,
        0.862410f, 0.862680f, 0.862960f});
    static constexpr Spectrum CES23({0.378530f, 0.383530f, 0.388220f, 0.392610f, 0.396730f, 0.400610f, 0.404260f,
        0.407720f, 0.411010f, 0.414150f, 0.417160f, 0.420070f, 0.422880f, 0.425610f, 0.428250f, 0.430820f, 0.433310f,
        0.435750f, 0.438130f, 0.440460f, 0.442740f, 0.444990f, 0.447200f, 0.449380f, 0.451520f, 0.453640f, 0.455730f,
        0.457800f, 0.459850f, 0.461870f, 0.463880f, 0.465870f, 0.467850f, 0.469800f, 0.471740f, 0.473660f, 0.475560f,
        0.477440f, 0.479300f, 0.481130f, 0.482940f, 0.484730f, 0.486490f, 0.488240f, 0.489980f, 0.491710f, 0.493440f,
        0.495170f, 0.496910f, 0.498650f, 0.500410f, 0.502190f, 0.503970f, 0.505770f, 0.507570f, 0.509370f, 0.511160f,
        0.512930f, 0.514690f, 0.516430f, 0.518140f, 0.519820f, 0.521470f, 0.523090f, 0.524700f, 0.526300f, 0.527890f,
        0.529480f, 0.531060f, 0.532660f, 0.534270f, 0.535900f, 0.537540f, 0.539200f, 0.540890f, 0.542600f, 0.544330f,
        0.546080f, 0.547870f, 0.549670f, 0.551510f, 0.553370f, 0.555260f, 0.557160f, 0.559080f, 0.560990f, 0.562910f,
        0.564810f, 0.566700f, 0.568570f, 0.570410f, 0.572220f, 0.573990f, 0.575740f, 0.577460f, 0.579150f, 0.580820f,
        0.582470f, 0.584100f, 0.585720f, 0.587320f, 0.588910f, 0.590500f, 0.592080f, 0.593670f, 0.595260f, 0.596880f,
        0.598510f, 0.600160f, 0.601850f, 0.603570f, 0.605330f, 0.607130f, 0.608970f, 0.610860f, 0.612780f, 0.614750f,
        0.616750f, 0.618800f, 0.620900f, 0.623030f, 0.625210f, 0.627420f, 0.629670f, 0.631950f, 0.634250f, 0.636570f,
        0.638910f, 0.641270f, 0.643630f, 0.645990f, 0.648350f, 0.650720f, 0.653100f, 0.655490f, 0.657890f, 0.660320f,
        0.662780f, 0.665270f, 0.667790f, 0.670360f, 0.672970f, 0.675640f, 0.678370f, 0.681170f, 0.684060f, 0.687030f,
        0.690100f, 0.693270f, 0.696560f, 0.699970f, 0.703510f, 0.707160f, 0.710900f, 0.714700f, 0.718560f, 0.722440f,
        0.726330f, 0.730210f, 0.734050f, 0.737840f, 0.741550f, 0.745180f, 0.748720f, 0.752150f, 0.755480f, 0.758690f,
        0.761770f, 0.764720f, 0.767530f, 0.770190f, 0.772700f, 0.775050f, 0.777270f, 0.779370f, 0.781350f, 0.783220f,
        0.785000f, 0.786690f, 0.788310f, 0.789860f, 0.791360f, 0.792810f, 0.794210f, 0.795570f, 0.796890f, 0.798180f,
        0.799450f, 0.800680f, 0.801900f, 0.803100f, 0.804290f, 0.805460f, 0.806600f, 0.807710f, 0.808790f, 0.809820f,
        0.810810f, 0.811740f, 0.812610f, 0.813410f, 0.814140f, 0.814800f, 0.815390f, 0.815920f, 0.816400f, 0.816810f,
        0.817180f, 0.817500f, 0.817770f, 0.818010f, 0.818210f, 0.818380f, 0.818520f, 0.818650f, 0.818760f, 0.818870f,
        0.818970f, 0.819080f, 0.819200f, 0.819340f, 0.819500f, 0.819670f, 0.819870f, 0.820080f, 0.820310f, 0.820560f,
        0.820810f, 0.821080f, 0.821350f, 0.821630f, 0.821920f, 0.822210f, 0.822500f, 0.822800f, 0.823110f, 0.823410f,
        0.823720f, 0.824040f, 0.824350f, 0.824670f, 0.824990f, 0.825310f, 0.825650f, 0.826000f, 0.826360f, 0.826740f,
        0.827150f, 0.827580f, 0.828050f, 0.828550f, 0.829090f, 0.829660f, 0.830260f, 0.830880f, 0.831520f, 0.832180f,
        0.832840f, 0.833500f, 0.834160f, 0.834810f, 0.835450f, 0.836080f, 0.836700f, 0.837320f, 0.837940f, 0.838560f,
        0.839190f, 0.839820f, 0.840470f, 0.841140f, 0.841820f, 0.842530f, 0.843240f, 0.843970f, 0.844710f, 0.845450f,
        0.846200f, 0.846950f, 0.847700f, 0.848440f, 0.849180f, 0.849910f, 0.850640f, 0.851370f, 0.852090f, 0.852820f,
        0.853540f, 0.854260f, 0.854990f, 0.855720f, 0.856450f, 0.857190f, 0.857930f, 0.858680f, 0.859430f, 0.860190f,
        0.860960f, 0.861730f, 0.862510f, 0.863290f, 0.864080f, 0.864880f, 0.865670f, 0.866460f, 0.867230f, 0.868000f,
        0.868740f, 0.869460f, 0.870160f, 0.870820f, 0.871450f, 0.872050f, 0.872620f, 0.873180f, 0.873730f, 0.874280f,
        0.874830f, 0.875390f, 0.875970f, 0.876570f, 0.877200f, 0.877850f, 0.878520f, 0.879180f, 0.879840f, 0.880480f,
        0.881100f, 0.881690f, 0.882220f, 0.882710f, 0.883140f, 0.883510f, 0.883840f, 0.884140f, 0.884410f, 0.884670f,
        0.884930f, 0.885200f, 0.885480f, 0.885790f, 0.886130f, 0.886510f, 0.886900f, 0.887300f, 0.887700f, 0.888100f,
        0.888490f, 0.888850f, 0.889170f, 0.889460f, 0.889700f, 0.889900f, 0.890060f, 0.890190f, 0.890300f, 0.890390f,
        0.890470f, 0.890550f, 0.890640f, 0.890740f, 0.890850f, 0.890980f, 0.891120f, 0.891270f, 0.891410f, 0.891550f,
        0.891690f, 0.891810f, 0.891920f, 0.892000f, 0.892060f, 0.892100f, 0.892120f, 0.892130f, 0.892120f, 0.892100f,
        0.892080f, 0.892050f, 0.892020f, 0.892000f, 0.891980f, 0.891970f, 0.891960f, 0.891960f, 0.891960f, 0.891960f,
        0.891970f, 0.891980f, 0.891990f, 0.892000f, 0.892010f, 0.892020f, 0.892030f, 0.892040f, 0.892040f, 0.892040f,
        0.892040f, 0.892030f, 0.892020f});
    static constexpr Spectrum CES24({0.200260f, 0.202600f, 0.204970f, 0.207360f, 0.209770f, 0.212200f, 0.214650f,
        0.217130f, 0.219620f, 0.222130f, 0.224660f, 0.227220f, 0.229790f, 0.232380f, 0.235000f, 0.237630f, 0.240290f,
        0.242960f, 0.245660f, 0.248370f, 0.251500f, 0.253630f, 0.256270f, 0.259300f, 0.262600f, 0.266050f, 0.269520f,
        0.272890f, 0.276050f, 0.278850f, 0.281200f, 0.282990f, 0.284260f, 0.285080f, 0.285510f, 0.285630f, 0.285500f,
        0.285180f, 0.284750f, 0.284260f, 0.283800f, 0.283410f, 0.283090f, 0.282840f, 0.282650f, 0.282490f, 0.282370f,
        0.282270f, 0.282180f, 0.282100f, 0.282000f, 0.281890f, 0.281760f, 0.281620f, 0.281480f, 0.281350f, 0.281230f,
        0.281130f, 0.281050f, 0.281010f, 0.281000f, 0.281040f, 0.281120f, 0.281240f, 0.281400f, 0.281610f, 0.281850f,
        0.282130f, 0.282450f, 0.282810f, 0.283200f, 0.283630f, 0.284090f, 0.284580f, 0.285110f, 0.285670f, 0.286270f,
        0.286900f, 0.287570f, 0.288270f, 0.289000f, 0.289770f, 0.290580f, 0.291430f, 0.292340f, 0.293300f, 0.294320f,
        0.295400f, 0.296560f, 0.297790f, 0.299100f, 0.300500f, 0.301980f, 0.303550f, 0.305210f, 0.306950f, 0.308780f,
        0.310710f, 0.312710f, 0.314810f, 0.317000f, 0.319280f, 0.321650f, 0.324120f, 0.326700f, 0.329390f, 0.332200f,
        0.335130f, 0.338180f, 0.341370f, 0.344700f, 0.348170f, 0.351780f, 0.355530f, 0.359410f, 0.363430f, 0.367570f,
        0.371840f, 0.376240f, 0.380760f, 0.385400f, 0.390160f, 0.395030f, 0.400030f, 0.405140f, 0.410390f, 0.415770f,
        0.421270f, 0.426910f, 0.432690f, 0.438600f, 0.444660f, 0.450860f, 0.457210f, 0.463710f, 0.470370f, 0.477200f,
        0.484190f, 0.491350f, 0.498690f, 0.506200f, 0.513890f, 0.521760f, 0.529780f, 0.537960f, 0.546260f, 0.554690f,
        0.563230f, 0.571870f, 0.580600f, 0.589400f, 0.598260f, 0.607180f, 0.616130f, 0.625110f, 0.634100f, 0.643090f,
        0.652070f, 0.661020f, 0.669940f, 0.678800f, 0.687600f, 0.696320f, 0.704930f, 0.713420f, 0.721760f, 0.729940f,
        0.737930f, 0.745720f, 0.753280f, 0.760600f, 0.767650f, 0.774430f, 0.780940f, 0.787180f, 0.793140f, 0.798820f,
        0.804220f, 0.809330f, 0.814160f, 0.818700f, 0.822950f, 0.826920f, 0.830620f, 0.834070f, 0.837270f, 0.840240f,
        0.842980f, 0.845520f, 0.847850f, 0.850000f, 0.851970f, 0.853780f, 0.855430f, 0.856930f, 0.858300f, 0.859540f,
        0.860660f, 0.861670f, 0.862580f, 0.863400f, 0.864140f, 0.864800f, 0.865400f, 0.865920f, 0.866370f, 0.866770f,
        0.867100f, 0.867390f, 0.867620f, 0.867800f, 0.867940f, 0.868040f, 0.868100f, 0.868130f, 0.868120f, 0.868080f,
        0.868020f, 0.867930f, 0.867830f, 0.867700f, 0.867560f, 0.867400f, 0.867240f, 0.867060f, 0.866880f, 0.866690f,
        0.866490f, 0.866290f, 0.866100f, 0.865900f, 0.865710f, 0.865510f, 0.865320f, 0.865130f, 0.864940f, 0.864740f,
        0.864540f, 0.864330f, 0.864120f, 0.863900f, 0.863670f, 0.863430f, 0.863180f, 0.862910f, 0.862630f, 0.862340f,
        0.862030f, 0.861700f, 0.861360f, 0.861000f, 0.860620f, 0.860220f, 0.859820f, 0.859400f, 0.858980f, 0.858570f,
        0.858160f, 0.857760f, 0.857370f, 0.857000f, 0.856650f, 0.856320f, 0.856010f, 0.855710f, 0.855430f, 0.855150f,
        0.854890f, 0.854620f, 0.854360f, 0.854100f, 0.853840f, 0.853560f, 0.853280f, 0.852970f, 0.852640f, 0.852270f,
        0.851870f, 0.851430f, 0.850940f, 0.850400f, 0.849800f, 0.849160f, 0.848480f, 0.847780f, 0.847080f, 0.846380f,
        0.845700f, 0.845050f, 0.844450f, 0.843900f, 0.843420f, 0.843010f, 0.842660f, 0.842370f, 0.842120f, 0.841930f,
        0.841770f, 0.841650f, 0.841560f, 0.841500f, 0.841460f, 0.841430f, 0.841410f, 0.841380f, 0.841340f, 0.841290f,
        0.841200f, 0.841080f, 0.840920f, 0.840700f, 0.840430f, 0.840110f, 0.839750f, 0.839380f, 0.838990f, 0.838610f,
        0.838250f, 0.837920f, 0.837630f, 0.837400f, 0.837050f, 0.836740f, 0.836440f, 0.836130f, 0.835820f, 0.835520f,
        0.835210f, 0.834900f, 0.834590f, 0.834280f, 0.833970f, 0.833660f, 0.833350f, 0.833040f, 0.832730f, 0.832410f,
        0.832100f, 0.831790f, 0.831480f, 0.831160f, 0.830850f, 0.830530f, 0.830220f, 0.829900f, 0.829580f, 0.829270f,
        0.828950f, 0.828630f, 0.828310f, 0.828000f, 0.827680f, 0.827360f, 0.827040f, 0.826720f, 0.826400f, 0.826070f,
        0.825750f, 0.825430f, 0.825110f, 0.824780f, 0.824460f, 0.824140f, 0.823810f, 0.823490f, 0.823160f, 0.822830f,
        0.822510f, 0.822180f, 0.821850f, 0.821520f, 0.821200f, 0.820870f, 0.820540f, 0.820210f, 0.819880f, 0.819550f,
        0.819210f, 0.818880f, 0.818550f, 0.818220f, 0.817880f, 0.817550f, 0.817220f, 0.816880f, 0.816550f, 0.816210f,
        0.815870f, 0.815540f, 0.815200f, 0.814860f, 0.814520f, 0.814190f, 0.813850f, 0.813510f, 0.813170f, 0.812830f,
        0.812490f, 0.812140f, 0.811800f});
    static constexpr Spectrum CES25({0.102810f, 0.100680f, 0.098595f, 0.096582f, 0.094694f, 0.092924f, 0.091152f,
        0.089347f, 0.087542f, 0.085757f, 0.084042f, 0.082434f, 0.080878f, 0.079330f, 0.077787f, 0.076213f, 0.074576f,
        0.072926f, 0.071304f, 0.069711f, 0.068141f, 0.066611f, 0.065182f, 0.063954f, 0.063010f, 0.062404f, 0.062137f,
        0.062086f, 0.062064f, 0.061918f, 0.061578f, 0.061078f, 0.060493f, 0.059873f, 0.059228f, 0.058536f, 0.057750f,
        0.056884f, 0.056009f, 0.055235f, 0.054642f, 0.054270f, 0.054119f, 0.054210f, 0.054513f, 0.055016f, 0.055679f,
        0.056420f, 0.057135f, 0.057791f, 0.058324f, 0.058714f, 0.058950f, 0.059003f, 0.058809f, 0.058404f, 0.057798f,
        0.057046f, 0.056238f, 0.055429f, 0.054616f, 0.053889f, 0.053318f, 0.052943f, 0.052809f, 0.052897f, 0.053070f,
        0.053216f, 0.053276f, 0.053201f, 0.053036f, 0.052882f, 0.052749f, 0.052584f, 0.052381f, 0.052101f, 0.051730f,
        0.051343f, 0.051033f, 0.050820f, 0.050747f, 0.050788f, 0.050871f, 0.050965f, 0.051104f, 0.051275f, 0.051515f,
        0.051854f, 0.052262f, 0.052658f, 0.052996f, 0.053223f, 0.053334f, 0.053340f, 0.053275f, 0.053163f, 0.053005f,
        0.052776f, 0.052489f, 0.052176f, 0.051878f, 0.051655f, 0.051578f, 0.051681f, 0.051963f, 0.052396f, 0.052925f,
        0.053494f, 0.054048f, 0.054577f, 0.055098f, 0.055654f, 0.056272f, 0.056967f, 0.057722f, 0.058514f, 0.059328f,
        0.060194f, 0.061150f, 0.062227f, 0.063430f, 0.064750f, 0.066140f, 0.067574f, 0.069053f, 0.070617f, 0.072327f,
        0.074253f, 0.076470f, 0.079084f, 0.082176f, 0.085793f, 0.089978f, 0.094738f, 0.100040f, 0.105890f, 0.112310f,
        0.119300f, 0.126910f, 0.135120f, 0.143880f, 0.153160f, 0.162940f, 0.173200f, 0.183910f, 0.195060f, 0.206540f,
        0.218210f, 0.229970f, 0.241710f, 0.253320f, 0.264740f, 0.275930f, 0.286820f, 0.297350f, 0.307500f, 0.317280f,
        0.326670f, 0.335710f, 0.344400f, 0.352710f, 0.360630f, 0.368150f, 0.375270f, 0.382020f, 0.388450f, 0.394600f,
        0.400480f, 0.406110f, 0.411490f, 0.416590f, 0.421390f, 0.425880f, 0.430050f, 0.433920f, 0.437530f, 0.440870f,
        0.443980f, 0.446870f, 0.449550f, 0.452060f, 0.454410f, 0.456680f, 0.458890f, 0.461060f, 0.463150f, 0.465150f,
        0.467030f, 0.468770f, 0.470390f, 0.471930f, 0.473410f, 0.474840f, 0.476250f, 0.477630f, 0.479000f, 0.480400f,
        0.481830f, 0.483330f, 0.484880f, 0.486440f, 0.487930f, 0.489330f, 0.490630f, 0.491790f, 0.492820f, 0.493750f,
        0.494570f, 0.495240f, 0.495800f, 0.496290f, 0.496810f, 0.497450f, 0.498250f, 0.499180f, 0.500200f, 0.501240f,
        0.502210f, 0.503080f, 0.503860f, 0.504560f, 0.505110f, 0.505540f, 0.505860f, 0.506090f, 0.506260f, 0.506410f,
        0.506610f, 0.506890f, 0.507270f, 0.507790f, 0.508400f, 0.509050f, 0.509660f, 0.510140f, 0.510420f, 0.510590f,
        0.510670f, 0.510670f, 0.510630f, 0.510550f, 0.510400f, 0.510210f, 0.510040f, 0.509910f, 0.509840f, 0.509820f,
        0.509800f, 0.509750f, 0.509670f, 0.509530f, 0.509390f, 0.509270f, 0.509190f, 0.509200f, 0.509350f, 0.509700f,
        0.510240f, 0.510960f, 0.511800f, 0.512700f, 0.513550f, 0.514260f, 0.514750f, 0.514990f, 0.514960f, 0.514650f,
        0.514110f, 0.513430f, 0.512690f, 0.511940f, 0.511320f, 0.510910f, 0.510740f, 0.510760f, 0.510940f, 0.511160f,
        0.511320f, 0.511340f, 0.511250f, 0.511070f, 0.510850f, 0.510620f, 0.510380f, 0.510120f, 0.509850f, 0.509550f,
        0.509190f, 0.508780f, 0.508260f, 0.507580f, 0.506690f, 0.505620f, 0.504390f, 0.503150f, 0.502010f, 0.501050f,
        0.500290f, 0.499770f, 0.499470f, 0.499370f, 0.499490f, 0.499820f, 0.500300f, 0.500830f, 0.501300f, 0.501700f,
        0.502100f, 0.502610f, 0.503360f, 0.504440f, 0.505820f, 0.507380f, 0.508940f, 0.510470f, 0.513900f, 0.516520f,
        0.518210f, 0.518210f, 0.518210f, 0.518210f, 0.521570f, 0.523020f, 0.524460f, 0.525900f, 0.527350f, 0.528790f,
        0.530230f, 0.531670f, 0.533110f, 0.534550f, 0.535990f, 0.537430f, 0.538870f, 0.540310f, 0.541750f, 0.543190f,
        0.544620f, 0.546060f, 0.547490f, 0.548930f, 0.550360f, 0.551790f, 0.553220f, 0.554650f, 0.556080f, 0.557510f,
        0.558940f, 0.560370f, 0.561790f, 0.563220f, 0.564640f, 0.566060f, 0.567480f, 0.568900f, 0.570320f, 0.571740f,
        0.573160f, 0.574570f, 0.575990f, 0.577400f, 0.578810f, 0.580220f, 0.581630f, 0.583040f, 0.584450f, 0.585850f,
        0.587260f, 0.588660f, 0.590060f, 0.591460f, 0.592860f, 0.594260f, 0.595650f, 0.597040f, 0.598440f, 0.599830f,
        0.601220f, 0.602600f, 0.603990f, 0.605370f, 0.606750f, 0.608130f, 0.609510f, 0.610890f, 0.612270f, 0.613640f,
        0.615010f, 0.616380f, 0.617750f, 0.619110f, 0.620480f, 0.621840f, 0.623200f, 0.624560f, 0.625920f, 0.627270f,
        0.628620f, 0.629970f, 0.631320f});
    static constexpr Spectrum CES26({0.181510f, 0.180860f, 0.180220f, 0.179570f, 0.178930f, 0.178290f, 0.177650f,
        0.177010f, 0.176380f, 0.175750f, 0.175110f, 0.174480f, 0.173860f, 0.173230f, 0.172610f, 0.171980f, 0.171360f,
        0.170740f, 0.170130f, 0.169510f, 0.168800f, 0.168340f, 0.167770f, 0.167110f, 0.166370f, 0.165580f, 0.164770f,
        0.163940f, 0.163120f, 0.162330f, 0.161600f, 0.160940f, 0.160350f, 0.159830f, 0.159380f, 0.159000f, 0.158690f,
        0.158450f, 0.158270f, 0.158150f, 0.158100f, 0.158110f, 0.158180f, 0.158320f, 0.158530f, 0.158810f, 0.159150f,
        0.159580f, 0.160070f, 0.160650f, 0.161300f, 0.162030f, 0.162840f, 0.163730f, 0.164690f, 0.165710f, 0.166790f,
        0.167940f, 0.169140f, 0.170400f, 0.171700f, 0.173050f, 0.174450f, 0.175900f, 0.177400f, 0.178950f, 0.180560f,
        0.182230f, 0.183960f, 0.185750f, 0.187600f, 0.189520f, 0.191500f, 0.193540f, 0.195630f, 0.197770f, 0.199950f,
        0.202170f, 0.204420f, 0.206700f, 0.209000f, 0.211320f, 0.213680f, 0.216090f, 0.218560f, 0.221120f, 0.223760f,
        0.226520f, 0.229400f, 0.232420f, 0.235600f, 0.238940f, 0.242440f, 0.246060f, 0.249790f, 0.253620f, 0.257510f,
        0.261460f, 0.265440f, 0.269420f, 0.273400f, 0.277350f, 0.281300f, 0.285260f, 0.289250f, 0.293290f, 0.297400f,
        0.301610f, 0.305930f, 0.310390f, 0.315000f, 0.319780f, 0.324730f, 0.329830f, 0.335090f, 0.340480f, 0.346000f,
        0.351640f, 0.357390f, 0.363250f, 0.369200f, 0.375240f, 0.381350f, 0.387550f, 0.393820f, 0.400170f, 0.406580f,
        0.413070f, 0.419620f, 0.426230f, 0.432900f, 0.439630f, 0.446420f, 0.453270f, 0.460190f, 0.467170f, 0.474230f,
        0.481360f, 0.488560f, 0.495840f, 0.503200f, 0.510640f, 0.518150f, 0.525700f, 0.533280f, 0.540880f, 0.548460f,
        0.556030f, 0.563550f, 0.571010f, 0.578400f, 0.585690f, 0.592900f, 0.600010f, 0.607040f, 0.614000f, 0.620870f,
        0.627680f, 0.634420f, 0.641090f, 0.647700f, 0.654250f, 0.660750f, 0.667170f, 0.673520f, 0.679790f, 0.685970f,
        0.692050f, 0.698040f, 0.703930f, 0.709700f, 0.715360f, 0.720890f, 0.726310f, 0.731600f, 0.736770f, 0.741820f,
        0.746730f, 0.751520f, 0.756180f, 0.760700f, 0.765090f, 0.769350f, 0.773490f, 0.777510f, 0.781410f, 0.785210f,
        0.788900f, 0.792490f, 0.795990f, 0.799400f, 0.802730f, 0.805970f, 0.809130f, 0.812220f, 0.815230f, 0.818160f,
        0.821020f, 0.823820f, 0.826540f, 0.829200f, 0.831790f, 0.834320f, 0.836790f, 0.839200f, 0.841550f, 0.843830f,
        0.846060f, 0.848230f, 0.850340f, 0.852400f, 0.854400f, 0.856350f, 0.858240f, 0.860070f, 0.861850f, 0.863580f,
        0.865240f, 0.866850f, 0.868400f, 0.869900f, 0.871340f, 0.872720f, 0.874060f, 0.875350f, 0.876590f, 0.877800f,
        0.878970f, 0.880110f, 0.881220f, 0.882300f, 0.883360f, 0.884400f, 0.885410f, 0.886400f, 0.887350f, 0.888280f,
        0.889160f, 0.890020f, 0.890830f, 0.891600f, 0.892330f, 0.893020f, 0.893670f, 0.894290f, 0.894870f, 0.895430f,
        0.895950f, 0.896460f, 0.896940f, 0.897400f, 0.897850f, 0.898280f, 0.898690f, 0.899080f, 0.899460f, 0.899830f,
        0.900170f, 0.900500f, 0.900810f, 0.901100f, 0.901370f, 0.901630f, 0.901880f, 0.902110f, 0.902320f, 0.902530f,
        0.902730f, 0.902930f, 0.903110f, 0.903300f, 0.903480f, 0.903660f, 0.903840f, 0.904020f, 0.904190f, 0.904360f,
        0.904530f, 0.904690f, 0.904850f, 0.905000f, 0.905150f, 0.905290f, 0.905430f, 0.905570f, 0.905700f, 0.905830f,
        0.905950f, 0.906070f, 0.906190f, 0.906300f, 0.906410f, 0.906510f, 0.906610f, 0.906710f, 0.906800f, 0.906890f,
        0.906970f, 0.907050f, 0.907130f, 0.907200f, 0.907270f, 0.907330f, 0.907400f, 0.907460f, 0.907520f, 0.907590f,
        0.907660f, 0.907730f, 0.907810f, 0.907900f, 0.907990f, 0.908090f, 0.908200f, 0.908300f, 0.908400f, 0.908500f,
        0.908590f, 0.908670f, 0.908740f, 0.908800f, 0.908880f, 0.908960f, 0.909030f, 0.909110f, 0.909180f, 0.909260f,
        0.909330f, 0.909410f, 0.909480f, 0.909560f, 0.909630f, 0.909700f, 0.909780f, 0.909850f, 0.909930f, 0.910000f,
        0.910070f, 0.910150f, 0.910220f, 0.910300f, 0.910370f, 0.910440f, 0.910520f, 0.910590f, 0.910660f, 0.910740f,
        0.910810f, 0.910890f, 0.910960f, 0.911030f, 0.911100f, 0.911180f, 0.911250f, 0.911320f, 0.911400f, 0.911470f,
        0.911540f, 0.911620f, 0.911690f, 0.911760f, 0.911830f, 0.911910f, 0.911980f, 0.912050f, 0.912120f, 0.912200f,
        0.912270f, 0.912340f, 0.912410f, 0.912490f, 0.912560f, 0.912630f, 0.912700f, 0.912770f, 0.912850f, 0.912920f,
        0.912990f, 0.913060f, 0.913130f, 0.913200f, 0.913280f, 0.913350f, 0.913420f, 0.913490f, 0.913560f, 0.913630f,
        0.913700f, 0.913780f, 0.913850f, 0.913920f, 0.913990f, 0.914060f, 0.914130f, 0.914200f, 0.914270f, 0.914340f,
        0.914410f, 0.914490f, 0.914560f});
    static constexpr Spectrum CES27({0.064624f, 0.064350f, 0.064299f, 0.064431f, 0.064685f, 0.064995f, 0.065259f,
        0.065449f, 0.065645f, 0.065884f, 0.066204f, 0.066605f, 0.066987f, 0.067268f, 0.067461f, 0.067575f, 0.067627f,
        0.067711f, 0.067868f, 0.068095f, 0.068390f, 0.068748f, 0.069193f, 0.069808f, 0.070678f, 0.071852f, 0.073348f,
        0.075113f, 0.077052f, 0.079087f, 0.081167f, 0.083330f, 0.085659f, 0.088183f, 0.090913f, 0.093890f, 0.097141f,
        0.100710f, 0.104720f, 0.109250f, 0.114340f, 0.119930f, 0.125940f, 0.132240f, 0.138770f, 0.145480f, 0.152330f,
        0.159270f, 0.166230f, 0.173120f, 0.179800f, 0.186180f, 0.192190f, 0.197760f, 0.202880f, 0.207580f, 0.211880f,
        0.215830f, 0.219460f, 0.222800f, 0.225870f, 0.228780f, 0.231650f, 0.234550f, 0.237560f, 0.240720f, 0.243910f,
        0.246980f, 0.249870f, 0.252500f, 0.254860f, 0.256940f, 0.258770f, 0.260320f, 0.261630f, 0.262720f, 0.263640f,
        0.264530f, 0.265510f, 0.266650f, 0.268050f, 0.269740f, 0.271680f, 0.273790f, 0.275990f, 0.278120f, 0.280150f,
        0.282060f, 0.283830f, 0.285420f, 0.286810f, 0.287910f, 0.288780f, 0.289460f, 0.290000f, 0.290480f, 0.290960f,
        0.291340f, 0.291610f, 0.291810f, 0.292000f, 0.292270f, 0.292700f, 0.293320f, 0.294150f, 0.295160f, 0.296330f,
        0.297690f, 0.299240f, 0.300970f, 0.302880f, 0.304970f, 0.307230f, 0.309620f, 0.312110f, 0.314680f, 0.317330f,
        0.320070f, 0.322920f, 0.325920f, 0.329070f, 0.332290f, 0.335510f, 0.338670f, 0.341730f, 0.344630f, 0.347360f,
        0.349910f, 0.352310f, 0.354620f, 0.356910f, 0.359250f, 0.361700f, 0.364270f, 0.366960f, 0.369740f, 0.372600f,
        0.375510f, 0.378460f, 0.381380f, 0.384220f, 0.386960f, 0.389590f, 0.392100f, 0.394520f, 0.396860f, 0.399050f,
        0.401080f, 0.402910f, 0.404520f, 0.405920f, 0.407160f, 0.408220f, 0.409130f, 0.409870f, 0.410450f, 0.410870f,
        0.411210f, 0.411520f, 0.411890f, 0.412370f, 0.413000f, 0.413740f, 0.414590f, 0.415490f, 0.416460f, 0.417480f,
        0.418580f, 0.419740f, 0.420930f, 0.422070f, 0.423130f, 0.424090f, 0.424960f, 0.425780f, 0.426570f, 0.427330f,
        0.428020f, 0.428600f, 0.429060f, 0.429410f, 0.429690f, 0.429960f, 0.430260f, 0.430590f, 0.430900f, 0.431180f,
        0.431400f, 0.431560f, 0.431690f, 0.431860f, 0.432110f, 0.432470f, 0.432940f, 0.433510f, 0.434170f, 0.434880f,
        0.435620f, 0.436400f, 0.437180f, 0.437990f, 0.438780f, 0.439540f, 0.440240f, 0.440850f, 0.441340f, 0.441730f,
        0.442050f, 0.442320f, 0.442540f, 0.442760f, 0.442960f, 0.443180f, 0.443420f, 0.443730f, 0.444100f, 0.444500f,
        0.444880f, 0.445220f, 0.445500f, 0.445740f, 0.445970f, 0.446210f, 0.446410f, 0.446550f, 0.446560f, 0.446430f,
        0.446240f, 0.446090f, 0.446040f, 0.446170f, 0.446470f, 0.446850f, 0.447230f, 0.447620f, 0.447940f, 0.448210f,
        0.448480f, 0.448750f, 0.448950f, 0.449090f, 0.449150f, 0.449150f, 0.449110f, 0.449100f, 0.449190f, 0.449450f,
        0.449830f, 0.450330f, 0.450910f, 0.451550f, 0.452180f, 0.452800f, 0.453400f, 0.453990f, 0.454560f, 0.455110f,
        0.455620f, 0.456100f, 0.456510f, 0.456810f, 0.456940f, 0.456870f, 0.456580f, 0.456040f, 0.455280f, 0.454350f,
        0.453320f, 0.452220f, 0.451110f, 0.450030f, 0.449030f, 0.448160f, 0.447470f, 0.446970f, 0.446660f, 0.446510f,
        0.446400f, 0.446280f, 0.446080f, 0.445790f, 0.445450f, 0.445120f, 0.444800f, 0.444570f, 0.444450f, 0.444400f,
        0.444370f, 0.444400f, 0.444400f, 0.444300f, 0.444100f, 0.443800f, 0.443420f, 0.443040f, 0.442720f, 0.442450f,
        0.442270f, 0.442200f, 0.442250f, 0.442420f, 0.442860f, 0.443590f, 0.444580f, 0.445750f, 0.447040f, 0.448340f,
        0.449620f, 0.450890f, 0.452190f, 0.453580f, 0.455040f, 0.456520f, 0.457880f, 0.459090f, 0.461210f, 0.462480f,
        0.463570f, 0.463570f, 0.463570f, 0.463570f, 0.467160f, 0.468380f, 0.469600f, 0.470820f, 0.472050f, 0.473270f,
        0.474490f, 0.475720f, 0.476940f, 0.478160f, 0.479390f, 0.480610f, 0.481840f, 0.483060f, 0.484290f, 0.485510f,
        0.486740f, 0.487960f, 0.489190f, 0.490420f, 0.491640f, 0.492870f, 0.494090f, 0.495320f, 0.496550f, 0.497770f,
        0.499000f, 0.500230f, 0.501450f, 0.502680f, 0.503910f, 0.505130f, 0.506360f, 0.507580f, 0.508810f, 0.510040f,
        0.511260f, 0.512490f, 0.513710f, 0.514940f, 0.516160f, 0.517390f, 0.518610f, 0.519840f, 0.521060f, 0.522290f,
        0.523510f, 0.524730f, 0.525960f, 0.527180f, 0.528400f, 0.529630f, 0.530850f, 0.532070f, 0.533290f, 0.534510f,
        0.535730f, 0.536950f, 0.538170f, 0.539390f, 0.540610f, 0.541830f, 0.543040f, 0.544260f, 0.545480f, 0.546690f,
        0.547910f, 0.549120f, 0.550340f, 0.551550f, 0.552770f, 0.553980f, 0.555190f, 0.556400f, 0.557610f, 0.558820f,
        0.560030f, 0.561240f, 0.562450f});
    static constexpr Spectrum CES28({0.045800f, 0.043900f, 0.047500f, 0.045014f, 0.044729f, 0.044886f, 0.044414f,
        0.044271f, 0.044186f, 0.043786f, 0.043943f, 0.043800f, 0.043157f, 0.043100f, 0.043271f, 0.043171f, 0.043429f,
        0.043629f, 0.043800f, 0.044329f, 0.044671f, 0.044686f, 0.045329f, 0.045571f, 0.045714f, 0.046500f, 0.047071f,
        0.047200f, 0.047971f, 0.048043f, 0.047929f, 0.048114f, 0.047771f, 0.047300f, 0.047486f, 0.046986f, 0.046786f,
        0.047143f, 0.046686f, 0.046114f, 0.046214f, 0.045714f, 0.045357f, 0.045014f, 0.044086f, 0.044000f, 0.044443f,
        0.044086f, 0.044100f, 0.044429f, 0.044557f, 0.045171f, 0.045486f, 0.045314f, 0.045971f, 0.046529f, 0.046443f,
        0.046771f, 0.046714f, 0.046686f, 0.047114f, 0.047114f, 0.046871f, 0.047500f, 0.047743f, 0.047943f, 0.048671f,
        0.049000f, 0.048843f, 0.048814f, 0.048014f, 0.047271f, 0.047129f, 0.046300f, 0.045586f, 0.045371f, 0.045000f,
        0.045000f, 0.045657f, 0.045514f, 0.045557f, 0.045914f, 0.046000f, 0.045986f, 0.046200f, 0.045743f, 0.046000f,
        0.046429f, 0.046343f, 0.046514f, 0.047043f, 0.047257f, 0.048129f, 0.048586f, 0.048414f, 0.048886f, 0.049400f,
        0.049586f, 0.050071f, 0.050114f, 0.049986f, 0.050843f, 0.051114f, 0.050786f, 0.050943f, 0.051129f, 0.051171f,
        0.051757f, 0.051286f, 0.051129f, 0.051843f, 0.051786f, 0.051500f, 0.051971f, 0.052000f, 0.052657f, 0.053614f,
        0.054100f, 0.055286f, 0.056814f, 0.057800f, 0.059114f, 0.060529f, 0.061714f, 0.063843f, 0.065829f, 0.067700f,
        0.070186f, 0.073129f, 0.075829f, 0.079229f, 0.082057f, 0.084914f, 0.088086f, 0.091014f, 0.093157f, 0.096086f,
        0.098414f, 0.100470f, 0.102660f, 0.104430f, 0.106000f, 0.108240f, 0.109560f, 0.110690f, 0.112370f, 0.113840f,
        0.115030f, 0.116310f, 0.116760f, 0.117170f, 0.118160f, 0.118200f, 0.118030f, 0.118170f, 0.117990f, 0.118300f,
        0.118770f, 0.118130f, 0.118110f, 0.118430f, 0.118490f, 0.118770f, 0.118840f, 0.118260f, 0.118690f, 0.119090f,
        0.118830f, 0.118940f, 0.119030f, 0.118770f, 0.119210f, 0.118990f, 0.118390f, 0.118260f, 0.117400f, 0.115770f,
        0.114610f, 0.113730f, 0.114290f, 0.115030f, 0.115840f, 0.116540f, 0.117970f, 0.119310f, 0.120310f, 0.120240f,
        0.119960f, 0.119700f, 0.120340f, 0.120270f, 0.120290f, 0.120730f, 0.120700f, 0.120640f, 0.120870f, 0.120310f,
        0.120100f, 0.120170f, 0.119910f, 0.120090f, 0.120630f, 0.120630f, 0.121160f, 0.121990f, 0.122140f, 0.122340f,
        0.122390f, 0.122130f, 0.122360f, 0.122370f, 0.121810f, 0.121810f, 0.121830f, 0.121330f, 0.121160f, 0.120610f,
        0.119890f, 0.120090f, 0.120000f, 0.119510f, 0.119600f, 0.119130f, 0.118790f, 0.119190f, 0.118960f, 0.118410f,
        0.118490f, 0.118190f, 0.118260f, 0.118640f, 0.118440f, 0.118190f, 0.118390f, 0.118200f, 0.118260f, 0.118560f,
        0.118230f, 0.117990f, 0.118130f, 0.118130f, 0.118290f, 0.118630f, 0.118410f, 0.118790f, 0.119370f, 0.119370f,
        0.119340f, 0.119340f, 0.118900f, 0.118810f, 0.118490f, 0.117590f, 0.117230f, 0.116910f, 0.116340f, 0.116460f,
        0.116360f, 0.115810f, 0.116090f, 0.116210f, 0.116060f, 0.116310f, 0.116190f, 0.115890f, 0.116370f, 0.116560f,
        0.116560f, 0.117010f, 0.116910f, 0.116560f, 0.116840f, 0.116540f, 0.116100f, 0.116140f, 0.115710f, 0.115490f,
        0.115810f, 0.115640f, 0.115410f, 0.115370f, 0.114800f, 0.114500f, 0.114310f, 0.113330f, 0.112670f, 0.112860f,
        0.112970f, 0.113590f, 0.114100f, 0.114200f, 0.115100f, 0.116030f, 0.116060f, 0.116190f, 0.116240f, 0.116170f,
        0.116610f, 0.116630f, 0.116160f, 0.116360f, 0.116430f, 0.116010f, 0.115710f, 0.114590f, 0.113200f, 0.113040f,
        0.113490f, 0.113570f, 0.114140f, 0.114560f, 0.115740f, 0.117510f, 0.118090f, 0.117690f, 0.117810f, 0.117290f,
        0.117230f, 0.117140f, 0.116470f, 0.116270f, 0.116270f, 0.115630f, 0.115340f, 0.115060f, 0.114430f, 0.114540f,
        0.114760f, 0.114630f, 0.114990f, 0.115590f, 0.115860f, 0.116740f, 0.117010f, 0.116800f, 0.116970f, 0.117060f,
        0.116870f, 0.117230f, 0.117160f, 0.117110f, 0.117640f, 0.117800f, 0.117970f, 0.118560f, 0.118710f, 0.118890f,
        0.119310f, 0.119260f, 0.119410f, 0.119790f, 0.119470f, 0.119240f, 0.119390f, 0.119100f, 0.118840f, 0.118730f,
        0.118240f, 0.118310f, 0.118440f, 0.117800f, 0.117660f, 0.118030f, 0.117940f, 0.118190f, 0.118330f, 0.117990f,
        0.118440f, 0.118940f, 0.118660f, 0.118800f, 0.118910f, 0.118800f, 0.119210f, 0.119230f, 0.119060f, 0.119660f,
        0.119930f, 0.119930f, 0.120330f, 0.120300f, 0.120260f, 0.120640f, 0.120500f, 0.120530f, 0.121060f, 0.121010f,
        0.121070f, 0.121400f, 0.121010f, 0.120990f, 0.121190f, 0.120810f, 0.120730f, 0.121140f, 0.121200f, 0.121310f,
        0.121440f, 0.121140f, 0.121090f});
    static constexpr Spectrum CES29({0.039813f, 0.040155f, 0.040499f, 0.040847f, 0.041197f, 0.041550f, 0.041907f,
        0.042266f, 0.042628f, 0.042993f, 0.043360f, 0.043731f, 0.044105f, 0.044482f, 0.044862f, 0.045246f, 0.045632f,
        0.046021f, 0.046414f, 0.046810f, 0.047981f, 0.048044f, 0.048164f, 0.048341f, 0.048575f, 0.048865f, 0.049211f,
        0.049612f, 0.050068f, 0.050580f, 0.051145f, 0.051765f, 0.052438f, 0.053165f, 0.053945f, 0.054778f, 0.055663f,
        0.056600f, 0.057589f, 0.058629f, 0.059720f, 0.060864f, 0.062066f, 0.063336f, 0.064683f, 0.066113f, 0.067637f,
        0.069262f, 0.070996f, 0.072849f, 0.074829f, 0.076943f, 0.079197f, 0.081593f, 0.084136f, 0.086830f, 0.089679f,
        0.092686f, 0.095856f, 0.099192f, 0.102700f, 0.106380f, 0.110230f, 0.114260f, 0.118450f, 0.122810f, 0.127340f,
        0.132040f, 0.136900f, 0.141920f, 0.147110f, 0.152450f, 0.157940f, 0.163580f, 0.169340f, 0.175230f, 0.181220f,
        0.187320f, 0.193510f, 0.199770f, 0.206110f, 0.212510f, 0.218950f, 0.225430f, 0.231920f, 0.238420f, 0.244900f,
        0.251360f, 0.257770f, 0.264130f, 0.270430f, 0.276640f, 0.282770f, 0.288810f, 0.294780f, 0.300660f, 0.306470f,
        0.312190f, 0.317850f, 0.323420f, 0.328920f, 0.334350f, 0.339700f, 0.344980f, 0.350190f, 0.355330f, 0.360400f,
        0.365390f, 0.370320f, 0.375180f, 0.379970f, 0.384690f, 0.389370f, 0.394010f, 0.398640f, 0.403280f, 0.407950f,
        0.412660f, 0.417430f, 0.422280f, 0.427230f, 0.432300f, 0.437480f, 0.442780f, 0.448180f, 0.453690f, 0.459310f,
        0.465030f, 0.470850f, 0.476760f, 0.482770f, 0.488860f, 0.495030f, 0.501250f, 0.507510f, 0.513780f, 0.520050f,
        0.526300f, 0.532510f, 0.538660f, 0.544730f, 0.550710f, 0.556600f, 0.562390f, 0.568080f, 0.573670f, 0.579150f,
        0.584540f, 0.589820f, 0.594990f, 0.600060f, 0.605020f, 0.609870f, 0.614590f, 0.619190f, 0.623650f, 0.627960f,
        0.632130f, 0.636140f, 0.639980f, 0.643650f, 0.647150f, 0.650480f, 0.653650f, 0.656670f, 0.659560f, 0.662310f,
        0.664950f, 0.667470f, 0.669900f, 0.672240f, 0.674500f, 0.676680f, 0.678790f, 0.680830f, 0.682800f, 0.684720f,
        0.686580f, 0.688380f, 0.690130f, 0.691840f, 0.693500f, 0.695120f, 0.696700f, 0.698240f, 0.699750f, 0.701210f,
        0.702640f, 0.704040f, 0.705410f, 0.706740f, 0.708050f, 0.709330f, 0.710590f, 0.711830f, 0.713050f, 0.714270f,
        0.715470f, 0.716680f, 0.717890f, 0.719100f, 0.720310f, 0.721540f, 0.722760f, 0.723990f, 0.725210f, 0.726420f,
        0.727630f, 0.728820f, 0.729990f, 0.731140f, 0.732270f, 0.733380f, 0.734480f, 0.735560f, 0.736630f, 0.737700f,
        0.738760f, 0.739820f, 0.740890f, 0.741960f, 0.743050f, 0.744140f, 0.745250f, 0.746370f, 0.747510f, 0.748670f,
        0.749850f, 0.751040f, 0.752260f, 0.753500f, 0.754760f, 0.756050f, 0.757350f, 0.758670f, 0.759990f, 0.761330f,
        0.762670f, 0.764010f, 0.765340f, 0.766670f, 0.767990f, 0.769290f, 0.770580f, 0.771840f, 0.773080f, 0.774300f,
        0.775480f, 0.776630f, 0.777750f, 0.778820f, 0.779840f, 0.780830f, 0.781760f, 0.782660f, 0.783510f, 0.784310f,
        0.785070f, 0.785790f, 0.786460f, 0.787090f, 0.787670f, 0.788210f, 0.788710f, 0.789180f, 0.789620f, 0.790030f,
        0.790420f, 0.790790f, 0.791140f, 0.791480f, 0.791800f, 0.792130f, 0.792470f, 0.792830f, 0.793210f, 0.793640f,
        0.794120f, 0.794660f, 0.795270f, 0.795970f, 0.796750f, 0.797610f, 0.798540f, 0.799540f, 0.800600f, 0.801700f,
        0.802840f, 0.804000f, 0.805190f, 0.806380f, 0.807580f, 0.808770f, 0.809970f, 0.811170f, 0.812360f, 0.813560f,
        0.814750f, 0.815950f, 0.817140f, 0.818320f, 0.819510f, 0.820690f, 0.821860f, 0.823030f, 0.824190f, 0.825330f,
        0.826460f, 0.827580f, 0.828680f, 0.829760f, 0.830820f, 0.831860f, 0.832870f, 0.833860f, 0.834820f, 0.835750f,
        0.836650f, 0.837520f, 0.838350f, 0.839150f, 0.840450f, 0.841400f, 0.842350f, 0.843290f, 0.844230f, 0.845160f,
        0.846090f, 0.847010f, 0.847930f, 0.848850f, 0.849760f, 0.850660f, 0.851570f, 0.852460f, 0.853350f, 0.854240f,
        0.855130f, 0.856000f, 0.856880f, 0.857750f, 0.858610f, 0.859480f, 0.860330f, 0.861180f, 0.862030f, 0.862880f,
        0.863720f, 0.864550f, 0.865380f, 0.866210f, 0.867030f, 0.867850f, 0.868660f, 0.869470f, 0.870280f, 0.871080f,
        0.871870f, 0.872670f, 0.873450f, 0.874240f, 0.875020f, 0.875790f, 0.876560f, 0.877330f, 0.878100f, 0.878850f,
        0.879610f, 0.880360f, 0.881110f, 0.881850f, 0.882590f, 0.883330f, 0.884060f, 0.884780f, 0.885510f, 0.886230f,
        0.886940f, 0.887650f, 0.888360f, 0.889060f, 0.889760f, 0.890460f, 0.891150f, 0.891840f, 0.892520f, 0.893200f,
        0.893880f, 0.894550f, 0.895220f, 0.895890f, 0.896550f, 0.897200f, 0.897860f, 0.898510f, 0.899160f, 0.899800f,
        0.900440f, 0.901070f, 0.901710f});
    static constexpr Spectrum CES30({0.030430f, 0.030767f, 0.031096f, 0.031410f, 0.031715f, 0.031976f, 0.032081f,
        0.032014f, 0.031841f, 0.031608f, 0.031346f, 0.031072f, 0.030760f, 0.030398f, 0.029997f, 0.029564f, 0.029106f,
        0.028686f, 0.028354f, 0.028124f, 0.027991f, 0.027991f, 0.028144f, 0.028445f, 0.028871f, 0.029427f, 0.030051f,
        0.030652f, 0.031146f, 0.031473f, 0.031587f, 0.031529f, 0.031337f, 0.031059f, 0.030771f, 0.030541f, 0.030379f,
        0.030318f, 0.030356f, 0.030475f, 0.030682f, 0.031018f, 0.031522f, 0.032335f, 0.033535f, 0.035104f, 0.036995f,
        0.039153f, 0.041434f, 0.043778f, 0.046178f, 0.048632f, 0.051116f, 0.053631f, 0.056110f, 0.058577f, 0.061110f,
        0.063807f, 0.066759f, 0.070079f, 0.073778f, 0.077809f, 0.082134f, 0.086771f, 0.091747f, 0.097106f, 0.102870f,
        0.109000f, 0.115380f, 0.121890f, 0.128430f, 0.134930f, 0.141400f, 0.147860f, 0.154310f, 0.160710f, 0.166990f,
        0.173060f, 0.178820f, 0.184250f, 0.189450f, 0.194430f, 0.199240f, 0.203960f, 0.208600f, 0.213010f, 0.217130f,
        0.220950f, 0.224410f, 0.227490f, 0.230250f, 0.232660f, 0.234700f, 0.236360f, 0.237650f, 0.238600f, 0.239400f,
        0.240180f, 0.241030f, 0.242010f, 0.243130f, 0.244340f, 0.245590f, 0.246900f, 0.248270f, 0.249700f, 0.251170f,
        0.252690f, 0.254190f, 0.255680f, 0.257150f, 0.258630f, 0.260080f, 0.261530f, 0.262940f, 0.264320f, 0.265650f,
        0.266930f, 0.268150f, 0.269230f, 0.270110f, 0.270740f, 0.271130f, 0.271350f, 0.271460f, 0.271550f, 0.271670f,
        0.271790f, 0.271890f, 0.271950f, 0.272000f, 0.272080f, 0.272270f, 0.272580f, 0.273010f, 0.273530f, 0.274090f,
        0.274640f, 0.275200f, 0.275750f, 0.276300f, 0.276830f, 0.277330f, 0.277710f, 0.278000f, 0.278190f, 0.278310f,
        0.278410f, 0.278500f, 0.278570f, 0.278620f, 0.278650f, 0.278660f, 0.278680f, 0.278760f, 0.278890f, 0.279040f,
        0.279190f, 0.279310f, 0.279410f, 0.279550f, 0.279770f, 0.280110f, 0.280580f, 0.281130f, 0.281730f, 0.282360f,
        0.283020f, 0.283720f, 0.284470f, 0.285240f, 0.285980f, 0.286670f, 0.287270f, 0.287810f, 0.288330f, 0.288880f,
        0.289490f, 0.290180f, 0.290950f, 0.291770f, 0.292640f, 0.293550f, 0.294530f, 0.295570f, 0.296640f, 0.297720f,
        0.298780f, 0.299770f, 0.300700f, 0.301540f, 0.302320f, 0.303050f, 0.303740f, 0.304400f, 0.305070f, 0.305750f,
        0.306440f, 0.307150f, 0.307920f, 0.308710f, 0.309540f, 0.310430f, 0.311370f, 0.312350f, 0.313400f, 0.314470f,
        0.315580f, 0.316720f, 0.317880f, 0.319030f, 0.320240f, 0.321470f, 0.322700f, 0.323910f, 0.325090f, 0.326180f,
        0.327180f, 0.328090f, 0.328940f, 0.329740f, 0.330540f, 0.331320f, 0.332120f, 0.332960f, 0.333880f, 0.334860f,
        0.335900f, 0.336970f, 0.338000f, 0.338940f, 0.339740f, 0.340430f, 0.341040f, 0.341630f, 0.342200f, 0.342770f,
        0.343320f, 0.343800f, 0.344180f, 0.344510f, 0.344840f, 0.345210f, 0.345630f, 0.346120f, 0.346610f, 0.347020f,
        0.347290f, 0.347480f, 0.347630f, 0.347790f, 0.348080f, 0.348500f, 0.349030f, 0.349680f, 0.350400f, 0.351160f,
        0.351980f, 0.352860f, 0.353670f, 0.354380f, 0.354930f, 0.355200f, 0.355220f, 0.355090f, 0.354850f, 0.354550f,
        0.354250f, 0.353950f, 0.353640f, 0.353380f, 0.353190f, 0.353150f, 0.353310f, 0.353630f, 0.354030f, 0.354480f,
        0.354880f, 0.355150f, 0.355370f, 0.355610f, 0.355920f, 0.356420f, 0.357060f, 0.357790f, 0.358520f, 0.359170f,
        0.359580f, 0.359790f, 0.359840f, 0.359710f, 0.359440f, 0.359100f, 0.358650f, 0.358070f, 0.357390f, 0.356640f,
        0.355890f, 0.355250f, 0.354770f, 0.354490f, 0.354370f, 0.354270f, 0.354110f, 0.353910f, 0.353730f, 0.353660f,
        0.353830f, 0.354220f, 0.354760f, 0.355400f, 0.356090f, 0.356820f, 0.357700f, 0.358790f, 0.360750f, 0.361580f,
        0.362290f, 0.363410f, 0.365140f, 0.365140f, 0.365920f, 0.366830f, 0.367750f, 0.368670f, 0.369590f, 0.370510f,
        0.371430f, 0.372350f, 0.373270f, 0.374190f, 0.375120f, 0.376040f, 0.376970f, 0.377900f, 0.378820f, 0.379750f,
        0.380680f, 0.381610f, 0.382540f, 0.383480f, 0.384410f, 0.385340f, 0.386280f, 0.387210f, 0.388150f, 0.389090f,
        0.390030f, 0.390960f, 0.391900f, 0.392840f, 0.393790f, 0.394730f, 0.395670f, 0.396610f, 0.397560f, 0.398500f,
        0.399450f, 0.400400f, 0.401340f, 0.402290f, 0.403240f, 0.404190f, 0.405140f, 0.406090f, 0.407040f, 0.408000f,
        0.408950f, 0.409900f, 0.410860f, 0.411810f, 0.412770f, 0.413720f, 0.414680f, 0.415640f, 0.416600f, 0.417560f,
        0.418520f, 0.419480f, 0.420440f, 0.421400f, 0.422360f, 0.423320f, 0.424290f, 0.425250f, 0.426220f, 0.427180f,
        0.428150f, 0.429110f, 0.430080f, 0.431050f, 0.432010f, 0.432980f, 0.433950f, 0.434920f, 0.435890f, 0.436860f,
        0.437830f, 0.438800f, 0.439770f});
    static constexpr Spectrum CES31({0.044010f, 0.043824f, 0.043705f, 0.043651f, 0.043661f, 0.043734f, 0.043868f,
        0.044062f, 0.044314f, 0.044624f, 0.044990f, 0.045410f, 0.045884f, 0.046409f, 0.046984f, 0.047609f, 0.048281f,
        0.049000f, 0.049763f, 0.050571f, 0.051420f, 0.052310f, 0.053239f, 0.054205f, 0.055205f, 0.056238f, 0.057302f,
        0.058394f, 0.059512f, 0.060655f, 0.061820f, 0.063007f, 0.064220f, 0.065467f, 0.066754f, 0.068088f, 0.069475f,
        0.070922f, 0.072436f, 0.074023f, 0.075690f, 0.077442f, 0.079281f, 0.081207f, 0.083219f, 0.085318f, 0.087504f,
        0.089777f, 0.092138f, 0.094585f, 0.097120f, 0.099742f, 0.102450f, 0.105240f, 0.108120f, 0.111080f, 0.114130f,
        0.117250f, 0.120460f, 0.123740f, 0.127100f, 0.130540f, 0.134060f, 0.137660f, 0.141350f, 0.145130f, 0.149000f,
        0.152980f, 0.157050f, 0.161220f, 0.165510f, 0.169900f, 0.174410f, 0.179010f, 0.183710f, 0.188500f, 0.193370f,
        0.198330f, 0.203360f, 0.208470f, 0.213640f, 0.218880f, 0.224190f, 0.229610f, 0.235140f, 0.240810f, 0.246640f,
        0.252640f, 0.258840f, 0.265250f, 0.271890f, 0.278780f, 0.285880f, 0.293160f, 0.300590f, 0.308130f, 0.315740f,
        0.323390f, 0.331040f, 0.338660f, 0.346210f, 0.353670f, 0.361040f, 0.368370f, 0.375680f, 0.382990f, 0.390330f,
        0.397730f, 0.405210f, 0.412810f, 0.420550f, 0.428450f, 0.436470f, 0.444590f, 0.452770f, 0.460970f, 0.469160f,
        0.477290f, 0.485330f, 0.493250f, 0.501010f, 0.508580f, 0.515950f, 0.523140f, 0.530160f, 0.537020f, 0.543710f,
        0.550260f, 0.556660f, 0.562930f, 0.569080f, 0.575110f, 0.581010f, 0.586770f, 0.592390f, 0.597860f, 0.603170f,
        0.608310f, 0.613260f, 0.618030f, 0.622600f, 0.626970f, 0.631140f, 0.635140f, 0.638970f, 0.642640f, 0.646180f,
        0.649590f, 0.652890f, 0.656090f, 0.659210f, 0.662250f, 0.665220f, 0.668130f, 0.670960f, 0.673730f, 0.676420f,
        0.679060f, 0.681620f, 0.684130f, 0.686570f, 0.688950f, 0.691270f, 0.693530f, 0.695740f, 0.697900f, 0.700010f,
        0.702080f, 0.704110f, 0.706090f, 0.708040f, 0.709960f, 0.711840f, 0.713700f, 0.715530f, 0.717340f, 0.719130f,
        0.720900f, 0.722660f, 0.724410f, 0.726150f, 0.727890f, 0.729610f, 0.731320f, 0.733010f, 0.734670f, 0.736300f,
        0.737900f, 0.739450f, 0.740960f, 0.742410f, 0.743810f, 0.745150f, 0.746460f, 0.747720f, 0.748950f, 0.750150f,
        0.751340f, 0.752510f, 0.753680f, 0.754840f, 0.756010f, 0.757180f, 0.758360f, 0.759530f, 0.760700f, 0.761870f,
        0.763040f, 0.764190f, 0.765340f, 0.766480f, 0.767600f, 0.768710f, 0.769800f, 0.770860f, 0.771900f, 0.772910f,
        0.773880f, 0.774810f, 0.775700f, 0.776550f, 0.777350f, 0.778110f, 0.778820f, 0.779510f, 0.780170f, 0.780810f,
        0.781430f, 0.782040f, 0.782650f, 0.783260f, 0.783870f, 0.784500f, 0.785120f, 0.785760f, 0.786410f, 0.787060f,
        0.787720f, 0.788400f, 0.789080f, 0.789780f, 0.790490f, 0.791200f, 0.791930f, 0.792650f, 0.793380f, 0.794110f,
        0.794830f, 0.795540f, 0.796250f, 0.796940f, 0.797620f, 0.798270f, 0.798910f, 0.799530f, 0.800120f, 0.800690f,
        0.801220f, 0.801730f, 0.802200f, 0.802640f, 0.803040f, 0.803420f, 0.803770f, 0.804110f, 0.804450f, 0.804790f,
        0.805150f, 0.805530f, 0.805950f, 0.806400f, 0.806900f, 0.807440f, 0.808000f, 0.808580f, 0.809160f, 0.809730f,
        0.810280f, 0.810800f, 0.811280f, 0.811700f, 0.812060f, 0.812360f, 0.812620f, 0.812840f, 0.813020f, 0.813200f,
        0.813360f, 0.813520f, 0.813690f, 0.813880f, 0.814100f, 0.814340f, 0.814620f, 0.814910f, 0.815240f, 0.815590f,
        0.815960f, 0.816360f, 0.816780f, 0.817220f, 0.817680f, 0.818160f, 0.818650f, 0.819160f, 0.819660f, 0.820170f,
        0.820680f, 0.821180f, 0.821670f, 0.822150f, 0.822610f, 0.823050f, 0.823480f, 0.823890f, 0.824290f, 0.824680f,
        0.825050f, 0.825410f, 0.825770f, 0.826110f, 0.826450f, 0.826770f, 0.827090f, 0.827400f, 0.827690f, 0.827980f,
        0.828250f, 0.828510f, 0.828760f, 0.828990f, 0.829210f, 0.829410f, 0.829600f, 0.829770f, 0.829930f, 0.830070f,
        0.830200f, 0.830310f, 0.830400f, 0.830480f, 0.830540f, 0.830580f, 0.830610f, 0.830610f, 0.830600f, 0.830570f,
        0.830520f, 0.830450f, 0.830370f, 0.830260f, 0.830550f, 0.830560f, 0.830580f, 0.830590f, 0.830600f, 0.830620f,
        0.830630f, 0.830650f, 0.830660f, 0.830670f, 0.830690f, 0.830700f, 0.830720f, 0.830730f, 0.830740f, 0.830760f,
        0.830770f, 0.830790f, 0.830800f, 0.830810f, 0.830830f, 0.830840f, 0.830850f, 0.830870f, 0.830880f, 0.830900f,
        0.830910f, 0.830920f, 0.830940f, 0.830950f, 0.830970f, 0.830980f, 0.830990f, 0.831010f, 0.831020f, 0.831030f,
        0.831050f, 0.831060f, 0.831080f, 0.831090f, 0.831100f, 0.831120f, 0.831130f, 0.831150f, 0.831160f, 0.831170f,
        0.831190f, 0.831200f, 0.831210f});
    static constexpr Spectrum CES32({0.046277f, 0.046104f, 0.045931f, 0.045759f, 0.045587f, 0.045416f, 0.045246f,
        0.045077f, 0.044908f, 0.044739f, 0.044571f, 0.044404f, 0.044237f, 0.044071f, 0.043906f, 0.043741f, 0.043577f,
        0.043413f, 0.043250f, 0.043088f, 0.042900f, 0.042779f, 0.042628f, 0.042453f, 0.042261f, 0.042058f, 0.041851f,
        0.041645f, 0.041447f, 0.041263f, 0.041100f, 0.040963f, 0.040852f, 0.040770f, 0.040715f, 0.040689f, 0.040692f,
        0.040724f, 0.040786f, 0.040878f, 0.041000f, 0.041153f, 0.041338f, 0.041555f, 0.041804f, 0.042086f, 0.042400f,
        0.042749f, 0.043131f, 0.043548f, 0.044000f, 0.044487f, 0.045010f, 0.045568f, 0.046163f, 0.046793f, 0.047461f,
        0.048165f, 0.048906f, 0.049684f, 0.050500f, 0.051354f, 0.052249f, 0.053189f, 0.054177f, 0.055216f, 0.056310f,
        0.057461f, 0.058675f, 0.059953f, 0.061300f, 0.062718f, 0.064209f, 0.065774f, 0.067414f, 0.069131f, 0.070925f,
        0.072798f, 0.074751f, 0.076784f, 0.078900f, 0.081099f, 0.083384f, 0.085756f, 0.088219f, 0.090773f, 0.093421f,
        0.096166f, 0.099009f, 0.101950f, 0.105000f, 0.108150f, 0.111410f, 0.114770f, 0.118230f, 0.121800f, 0.125480f,
        0.129250f, 0.133130f, 0.137110f, 0.141200f, 0.145390f, 0.149670f, 0.154050f, 0.158520f, 0.163070f, 0.167700f,
        0.172400f, 0.177170f, 0.182000f, 0.186900f, 0.191850f, 0.196860f, 0.201910f, 0.207020f, 0.212180f, 0.217390f,
        0.222650f, 0.227950f, 0.233300f, 0.238700f, 0.244140f, 0.249630f, 0.255160f, 0.260740f, 0.266360f, 0.272030f,
        0.277750f, 0.283520f, 0.289340f, 0.295200f, 0.301110f, 0.307070f, 0.313070f, 0.319110f, 0.325170f, 0.331260f,
        0.337380f, 0.343510f, 0.349650f, 0.355800f, 0.361950f, 0.368100f, 0.374230f, 0.380350f, 0.386440f, 0.392490f,
        0.398500f, 0.404460f, 0.410360f, 0.416200f, 0.421960f, 0.427640f, 0.433220f, 0.438680f, 0.444020f, 0.449210f,
        0.454250f, 0.459120f, 0.463810f, 0.468300f, 0.472580f, 0.476650f, 0.480510f, 0.484140f, 0.487560f, 0.490760f,
        0.493730f, 0.496480f, 0.499010f, 0.501300f, 0.503360f, 0.505200f, 0.506820f, 0.508220f, 0.509410f, 0.510400f,
        0.511180f, 0.511780f, 0.512180f, 0.512400f, 0.512440f, 0.512320f, 0.512060f, 0.511650f, 0.511130f, 0.510500f,
        0.509770f, 0.508970f, 0.508110f, 0.507200f, 0.506250f, 0.505290f, 0.504300f, 0.503320f, 0.502350f, 0.501390f,
        0.500470f, 0.499590f, 0.498760f, 0.498000f, 0.497310f, 0.496680f, 0.496110f, 0.495600f, 0.495120f, 0.494680f,
        0.494270f, 0.493870f, 0.493480f, 0.493100f, 0.492710f, 0.492300f, 0.491880f, 0.491420f, 0.490930f, 0.490380f,
        0.489790f, 0.489130f, 0.488400f, 0.487600f, 0.486710f, 0.485750f, 0.484710f, 0.483620f, 0.482480f, 0.481290f,
        0.480070f, 0.478830f, 0.477570f, 0.476300f, 0.475040f, 0.473790f, 0.472570f, 0.471410f, 0.470300f, 0.469280f,
        0.468350f, 0.467540f, 0.466850f, 0.466300f, 0.465910f, 0.465690f, 0.465660f, 0.465840f, 0.466230f, 0.466860f,
        0.467730f, 0.468870f, 0.470290f, 0.472000f, 0.474020f, 0.476340f, 0.478980f, 0.481930f, 0.485180f, 0.488750f,
        0.492620f, 0.496800f, 0.501300f, 0.506100f, 0.511210f, 0.516620f, 0.522310f, 0.528260f, 0.534470f, 0.540920f,
        0.547590f, 0.554470f, 0.561540f, 0.568800f, 0.576220f, 0.583780f, 0.591450f, 0.599200f, 0.606990f, 0.614810f,
        0.622620f, 0.630390f, 0.638090f, 0.645700f, 0.653180f, 0.660540f, 0.667750f, 0.674820f, 0.681730f, 0.688490f,
        0.695080f, 0.701500f, 0.707740f, 0.713800f, 0.719670f, 0.725360f, 0.730890f, 0.736270f, 0.741510f, 0.746640f,
        0.751650f, 0.756580f, 0.761420f, 0.766200f, 0.770930f, 0.775620f, 0.780280f, 0.784920f, 0.789560f, 0.794210f,
        0.798880f, 0.803570f, 0.808310f, 0.813100f, 0.817940f, 0.822780f, 0.827530f, 0.832140f, 0.836520f, 0.840610f,
        0.844340f, 0.847620f, 0.850400f, 0.852600f, 0.855940f, 0.858780f, 0.861560f, 0.864300f, 0.866990f, 0.869640f,
        0.872240f, 0.874800f, 0.877320f, 0.879790f, 0.882220f, 0.884600f, 0.886950f, 0.889250f, 0.891510f, 0.893730f,
        0.895910f, 0.898050f, 0.900150f, 0.902210f, 0.904230f, 0.906220f, 0.908170f, 0.910080f, 0.911960f, 0.913800f,
        0.915610f, 0.917380f, 0.919120f, 0.920830f, 0.922500f, 0.924140f, 0.925750f, 0.927320f, 0.928870f, 0.930380f,
        0.931870f, 0.933320f, 0.934750f, 0.936150f, 0.937520f, 0.938860f, 0.940180f, 0.941470f, 0.942730f, 0.943970f,
        0.945180f, 0.946370f, 0.947530f, 0.948670f, 0.949790f, 0.950880f, 0.951950f, 0.953000f, 0.954030f, 0.955030f,
        0.956010f, 0.956980f, 0.957920f, 0.958850f, 0.959750f, 0.960640f, 0.961500f, 0.962350f, 0.963180f, 0.963990f,
        0.964790f, 0.965570f, 0.966330f, 0.967070f, 0.967800f, 0.968520f, 0.969220f, 0.969900f, 0.970570f, 0.971220f,
        0.971860f, 0.972490f, 0.973100f});
    static constexpr Spectrum CES33({0.173150f, 0.174480f, 0.175790f, 0.177060f, 0.178310f, 0.179550f, 0.180760f,
        0.181960f, 0.183150f, 0.184330f, 0.185510f, 0.186690f, 0.187870f, 0.189060f, 0.190260f, 0.191470f, 0.192710f,
        0.193960f, 0.195230f, 0.196540f, 0.197870f, 0.199240f, 0.200620f, 0.202000f, 0.203370f, 0.204700f, 0.205980f,
        0.207190f, 0.208310f, 0.209330f, 0.210230f, 0.211000f, 0.211660f, 0.212230f, 0.212740f, 0.213210f, 0.213670f,
        0.214150f, 0.214680f, 0.215270f, 0.215950f, 0.216750f, 0.217660f, 0.218700f, 0.219850f, 0.221120f, 0.222500f,
        0.224010f, 0.225630f, 0.227360f, 0.229210f, 0.231180f, 0.233280f, 0.235550f, 0.237990f, 0.240630f, 0.243500f,
        0.246620f, 0.249990f, 0.253660f, 0.257640f, 0.261950f, 0.266600f, 0.271630f, 0.277030f, 0.282830f, 0.289050f,
        0.295700f, 0.302800f, 0.310360f, 0.318410f, 0.326940f, 0.335920f, 0.345270f, 0.354940f, 0.364870f, 0.375000f,
        0.385270f, 0.395620f, 0.405990f, 0.416310f, 0.426540f, 0.436670f, 0.446670f, 0.456560f, 0.466320f, 0.475930f,
        0.485410f, 0.494730f, 0.503900f, 0.512900f, 0.521720f, 0.530340f, 0.538730f, 0.546840f, 0.554660f, 0.562140f,
        0.569250f, 0.575970f, 0.582260f, 0.588090f, 0.593440f, 0.598320f, 0.602760f, 0.606810f, 0.610480f, 0.613800f,
        0.616810f, 0.619530f, 0.622000f, 0.624240f, 0.626290f, 0.628160f, 0.629910f, 0.631540f, 0.633100f, 0.634610f,
        0.636110f, 0.637620f, 0.639170f, 0.640790f, 0.642510f, 0.644320f, 0.646210f, 0.648160f, 0.650180f, 0.652230f,
        0.654320f, 0.656440f, 0.658560f, 0.660680f, 0.662790f, 0.664880f, 0.666950f, 0.668990f, 0.671010f, 0.672990f,
        0.674930f, 0.676830f, 0.678680f, 0.680480f, 0.682230f, 0.683930f, 0.685580f, 0.687210f, 0.688810f, 0.690400f,
        0.691970f, 0.693550f, 0.695130f, 0.696720f, 0.698330f, 0.699960f, 0.701600f, 0.703240f, 0.704890f, 0.706520f,
        0.708140f, 0.709740f, 0.711320f, 0.712870f, 0.714380f, 0.715850f, 0.717290f, 0.718670f, 0.720010f, 0.721300f,
        0.722530f, 0.723710f, 0.724820f, 0.725880f, 0.726870f, 0.727810f, 0.728720f, 0.729600f, 0.730480f, 0.731360f,
        0.732260f, 0.733190f, 0.734170f, 0.735220f, 0.736340f, 0.737520f, 0.738740f, 0.739980f, 0.741230f, 0.742480f,
        0.743690f, 0.744860f, 0.745960f, 0.746990f, 0.747920f, 0.748770f, 0.749540f, 0.750250f, 0.750910f, 0.751530f,
        0.752120f, 0.752690f, 0.753260f, 0.753830f, 0.754420f, 0.755020f, 0.755640f, 0.756280f, 0.756920f, 0.757570f,
        0.758220f, 0.758880f, 0.759550f, 0.760210f, 0.760870f, 0.761540f, 0.762210f, 0.762900f, 0.763610f, 0.764350f,
        0.765110f, 0.765910f, 0.766750f, 0.767630f, 0.768560f, 0.769540f, 0.770550f, 0.771590f, 0.772650f, 0.773720f,
        0.774790f, 0.775870f, 0.776930f, 0.777970f, 0.778990f, 0.779970f, 0.780930f, 0.781840f, 0.782710f, 0.783530f,
        0.784310f, 0.785020f, 0.785680f, 0.786270f, 0.786800f, 0.787260f, 0.787690f, 0.788070f, 0.788440f, 0.788790f,
        0.789130f, 0.789490f, 0.789870f, 0.790280f, 0.790730f, 0.791210f, 0.791720f, 0.792250f, 0.792780f, 0.793310f,
        0.793840f, 0.794350f, 0.794830f, 0.795280f, 0.795690f, 0.796060f, 0.796410f, 0.796740f, 0.797060f, 0.797370f,
        0.797690f, 0.798030f, 0.798380f, 0.798770f, 0.799190f, 0.799650f, 0.800130f, 0.800640f, 0.801180f, 0.801740f,
        0.802320f, 0.802910f, 0.803510f, 0.804120f, 0.804740f, 0.805350f, 0.805970f, 0.806580f, 0.807180f, 0.807770f,
        0.808350f, 0.808900f, 0.809440f, 0.809950f, 0.810430f, 0.810890f, 0.811330f, 0.811750f, 0.812150f, 0.812540f,
        0.812920f, 0.813290f, 0.813660f, 0.814030f, 0.814400f, 0.814770f, 0.815150f, 0.815540f, 0.815930f, 0.816340f,
        0.816750f, 0.817180f, 0.817620f, 0.818070f, 0.818540f, 0.819030f, 0.819520f, 0.820030f, 0.820550f, 0.821080f,
        0.821620f, 0.822160f, 0.822710f, 0.823260f, 0.823810f, 0.824360f, 0.824900f, 0.825430f, 0.825930f, 0.826420f,
        0.826880f, 0.827300f, 0.827690f, 0.828040f, 0.828340f, 0.828610f, 0.828840f, 0.829050f, 0.829240f, 0.829430f,
        0.829620f, 0.829820f, 0.830030f, 0.830270f, 0.830540f, 0.830860f, 0.831220f, 0.831640f, 0.832130f, 0.832680f,
        0.833320f, 0.834050f, 0.834880f, 0.835810f, 0.835220f, 0.835650f, 0.836080f, 0.836500f, 0.836930f, 0.837350f,
        0.837780f, 0.838200f, 0.838620f, 0.839040f, 0.839460f, 0.839880f, 0.840300f, 0.840720f, 0.841140f, 0.841550f,
        0.841970f, 0.842380f, 0.842790f, 0.843210f, 0.843620f, 0.844030f, 0.844440f, 0.844850f, 0.845250f, 0.845660f,
        0.846070f, 0.846470f, 0.846880f, 0.847280f, 0.847680f, 0.848090f, 0.848490f, 0.848890f, 0.849290f, 0.849680f,
        0.850080f, 0.850480f, 0.850870f, 0.851270f, 0.851660f, 0.852060f, 0.852450f, 0.852840f, 0.853230f, 0.853620f,
        0.854010f, 0.854400f, 0.854780f});
    static constexpr Spectrum CES34({0.043542f, 0.043559f, 0.043576f, 0.043593f, 0.043610f, 0.043627f, 0.043644f,
        0.043661f, 0.043678f, 0.043695f, 0.043712f, 0.043729f, 0.043746f, 0.043763f, 0.043780f, 0.043797f, 0.043814f,
        0.043831f, 0.043848f, 0.043866f, 0.043885f, 0.043898f, 0.043915f, 0.043933f, 0.043953f, 0.043972f, 0.043991f,
        0.044009f, 0.044023f, 0.044034f, 0.044040f, 0.044041f, 0.044037f, 0.044031f, 0.044022f, 0.044012f, 0.044003f,
        0.043996f, 0.043992f, 0.043991f, 0.043996f, 0.044007f, 0.044024f, 0.044047f, 0.044075f, 0.044108f, 0.044145f,
        0.044187f, 0.044232f, 0.044280f, 0.044331f, 0.044384f, 0.044441f, 0.044501f, 0.044567f, 0.044638f, 0.044716f,
        0.044802f, 0.044896f, 0.044999f, 0.045112f, 0.045237f, 0.045375f, 0.045528f, 0.045699f, 0.045889f, 0.046101f,
        0.046337f, 0.046599f, 0.046890f, 0.047210f, 0.047563f, 0.047944f, 0.048352f, 0.048781f, 0.049230f, 0.049694f,
        0.050171f, 0.050656f, 0.051146f, 0.051639f, 0.052132f, 0.052632f, 0.053149f, 0.053691f, 0.054266f, 0.054885f,
        0.055555f, 0.056285f, 0.057085f, 0.057963f, 0.058926f, 0.059974f, 0.061105f, 0.062318f, 0.063611f, 0.064982f,
        0.066429f, 0.067951f, 0.069545f, 0.071210f, 0.072945f, 0.074752f, 0.076632f, 0.078589f, 0.080625f, 0.082742f,
        0.084942f, 0.087228f, 0.089603f, 0.092068f, 0.094626f, 0.097274f, 0.100010f, 0.102830f, 0.105730f, 0.108720f,
        0.111770f, 0.114910f, 0.118110f, 0.121380f, 0.124720f, 0.128140f, 0.131660f, 0.135290f, 0.139060f, 0.142980f,
        0.147060f, 0.151330f, 0.155790f, 0.160480f, 0.165390f, 0.170530f, 0.175890f, 0.181480f, 0.187270f, 0.193280f,
        0.199490f, 0.205910f, 0.212520f, 0.219320f, 0.226300f, 0.233400f, 0.240550f, 0.247690f, 0.254740f, 0.261640f,
        0.268330f, 0.274730f, 0.280790f, 0.286420f, 0.291580f, 0.296280f, 0.300530f, 0.304370f, 0.307800f, 0.310860f,
        0.313560f, 0.315930f, 0.317980f, 0.319750f, 0.321240f, 0.322490f, 0.323520f, 0.324340f, 0.324970f, 0.325450f,
        0.325800f, 0.326020f, 0.326150f, 0.326210f, 0.326220f, 0.326180f, 0.326090f, 0.325960f, 0.325790f, 0.325590f,
        0.325350f, 0.325080f, 0.324780f, 0.324460f, 0.324110f, 0.323750f, 0.323360f, 0.322950f, 0.322530f, 0.322090f,
        0.321630f, 0.321160f, 0.320690f, 0.320200f, 0.319700f, 0.319190f, 0.318670f, 0.318150f, 0.317620f, 0.317080f,
        0.316530f, 0.315980f, 0.315420f, 0.314860f, 0.314280f, 0.313700f, 0.313110f, 0.312500f, 0.311890f, 0.311250f,
        0.310590f, 0.309920f, 0.309210f, 0.308490f, 0.307730f, 0.306950f, 0.306140f, 0.305300f, 0.304430f, 0.303540f,
        0.302620f, 0.301680f, 0.300710f, 0.299710f, 0.298680f, 0.297640f, 0.296580f, 0.295510f, 0.294430f, 0.293360f,
        0.292290f, 0.291230f, 0.290190f, 0.289170f, 0.288170f, 0.287200f, 0.286260f, 0.285340f, 0.284460f, 0.283590f,
        0.282760f, 0.281960f, 0.281180f, 0.280430f, 0.279710f, 0.279020f, 0.278350f, 0.277700f, 0.277070f, 0.276460f,
        0.275860f, 0.275280f, 0.274700f, 0.274140f, 0.273580f, 0.273030f, 0.272480f, 0.271950f, 0.271410f, 0.270890f,
        0.270370f, 0.269860f, 0.269350f, 0.268850f, 0.268360f, 0.267880f, 0.267400f, 0.266940f, 0.266490f, 0.266050f,
        0.265640f, 0.265230f, 0.264850f, 0.264490f, 0.264150f, 0.263840f, 0.263540f, 0.263270f, 0.263030f, 0.262800f,
        0.262600f, 0.262420f, 0.262260f, 0.262130f, 0.262020f, 0.261940f, 0.261880f, 0.261850f, 0.261860f, 0.261900f,
        0.261990f, 0.262110f, 0.262270f, 0.262480f, 0.262740f, 0.263030f, 0.263360f, 0.263720f, 0.264100f, 0.264500f,
        0.264900f, 0.265310f, 0.265700f, 0.266090f, 0.266460f, 0.266820f, 0.267160f, 0.267510f, 0.267860f, 0.268210f,
        0.268580f, 0.268970f, 0.269380f, 0.269810f, 0.270280f, 0.270780f, 0.271280f, 0.271780f, 0.272280f, 0.272750f,
        0.273180f, 0.273570f, 0.273910f, 0.274170f, 0.274590f, 0.274950f, 0.275300f, 0.275660f, 0.276020f, 0.276380f,
        0.276740f, 0.277100f, 0.277460f, 0.277820f, 0.278180f, 0.278540f, 0.278900f, 0.279270f, 0.279630f, 0.279990f,
        0.280350f, 0.280710f, 0.281080f, 0.281440f, 0.281800f, 0.282170f, 0.282530f, 0.282900f, 0.283260f, 0.283630f,
        0.283990f, 0.284360f, 0.284720f, 0.285090f, 0.285450f, 0.285820f, 0.286190f, 0.286560f, 0.286920f, 0.287290f,
        0.287660f, 0.288030f, 0.288400f, 0.288760f, 0.289130f, 0.289500f, 0.289870f, 0.290240f, 0.290610f, 0.290980f,
        0.291350f, 0.291730f, 0.292100f, 0.292470f, 0.292840f, 0.293210f, 0.293590f, 0.293960f, 0.294330f, 0.294700f,
        0.295080f, 0.295450f, 0.295830f, 0.296200f, 0.296580f, 0.296950f, 0.297330f, 0.297700f, 0.298080f, 0.298450f,
        0.298830f, 0.299210f, 0.299580f, 0.299960f, 0.300340f, 0.300720f, 0.301090f, 0.301470f, 0.301850f, 0.302230f,
        0.302610f, 0.302990f, 0.303370f});
    static constexpr Spectrum CES35({0.052651f, 0.052641f, 0.052631f, 0.052621f, 0.052611f, 0.052601f, 0.052591f,
        0.052581f, 0.052571f, 0.052561f, 0.052552f, 0.052542f, 0.052532f, 0.052522f, 0.052512f, 0.052502f, 0.052492f,
        0.052482f, 0.052472f, 0.052462f, 0.052449f, 0.052444f, 0.052436f, 0.052425f, 0.052409f, 0.052388f, 0.052360f,
        0.052325f, 0.052282f, 0.052229f, 0.052166f, 0.052092f, 0.052008f, 0.051914f, 0.051812f, 0.051703f, 0.051587f,
        0.051466f, 0.051341f, 0.051212f, 0.051080f, 0.050947f, 0.050814f, 0.050680f, 0.050547f, 0.050415f, 0.050285f,
        0.050159f, 0.050035f, 0.049916f, 0.049802f, 0.049694f, 0.049591f, 0.049494f, 0.049403f, 0.049318f, 0.049240f,
        0.049168f, 0.049102f, 0.049042f, 0.048989f, 0.048943f, 0.048904f, 0.048872f, 0.048846f, 0.048828f, 0.048817f,
        0.048813f, 0.048816f, 0.048827f, 0.048845f, 0.048871f, 0.048903f, 0.048941f, 0.048984f, 0.049031f, 0.049081f,
        0.049133f, 0.049186f, 0.049239f, 0.049292f, 0.049343f, 0.049397f, 0.049455f, 0.049523f, 0.049602f, 0.049696f,
        0.049809f, 0.049944f, 0.050105f, 0.050294f, 0.050515f, 0.050772f, 0.051066f, 0.051400f, 0.051777f, 0.052199f,
        0.052670f, 0.053192f, 0.053767f, 0.054399f, 0.055089f, 0.055836f, 0.056640f, 0.057498f, 0.058408f, 0.059371f,
        0.060383f, 0.061444f, 0.062552f, 0.063705f, 0.064901f, 0.066134f, 0.067397f, 0.068683f, 0.069985f, 0.071295f,
        0.072607f, 0.073913f, 0.075207f, 0.076481f, 0.077730f, 0.078951f, 0.080144f, 0.081308f, 0.082441f, 0.083544f,
        0.084614f, 0.085652f, 0.086657f, 0.087627f, 0.088565f, 0.089482f, 0.090393f, 0.091313f, 0.092257f, 0.093240f,
        0.094276f, 0.095381f, 0.096568f, 0.097853f, 0.099246f, 0.100740f, 0.102320f, 0.103980f, 0.105700f, 0.107480f,
        0.109290f, 0.111130f, 0.112990f, 0.114850f, 0.116700f, 0.118530f, 0.120320f, 0.122060f, 0.123720f, 0.125300f,
        0.126780f, 0.128140f, 0.129370f, 0.130450f, 0.131380f, 0.132150f, 0.132780f, 0.133290f, 0.133680f, 0.133960f,
        0.134150f, 0.134260f, 0.134300f, 0.134280f, 0.134210f, 0.134100f, 0.133950f, 0.133770f, 0.133550f, 0.133310f,
        0.133050f, 0.132760f, 0.132460f, 0.132160f, 0.131840f, 0.131520f, 0.131200f, 0.130870f, 0.130540f, 0.130210f,
        0.129870f, 0.129530f, 0.129200f, 0.128860f, 0.128530f, 0.128200f, 0.127870f, 0.127530f, 0.127200f, 0.126870f,
        0.126540f, 0.126210f, 0.125870f, 0.125540f, 0.125200f, 0.124860f, 0.124510f, 0.124160f, 0.123810f, 0.123450f,
        0.123090f, 0.122720f, 0.122340f, 0.121960f, 0.121570f, 0.121170f, 0.120760f, 0.120350f, 0.119920f, 0.119490f,
        0.119050f, 0.118600f, 0.118140f, 0.117670f, 0.117190f, 0.116700f, 0.116210f, 0.115710f, 0.115220f, 0.114730f,
        0.114240f, 0.113760f, 0.113290f, 0.112830f, 0.112380f, 0.111950f, 0.111530f, 0.111130f, 0.110740f, 0.110360f,
        0.109990f, 0.109630f, 0.109280f, 0.108940f, 0.108610f, 0.108290f, 0.107980f, 0.107680f, 0.107380f, 0.107090f,
        0.106800f, 0.106520f, 0.106250f, 0.105980f, 0.105710f, 0.105450f, 0.105190f, 0.104930f, 0.104680f, 0.104430f,
        0.104190f, 0.103950f, 0.103710f, 0.103490f, 0.103260f, 0.103040f, 0.102830f, 0.102620f, 0.102420f, 0.102230f,
        0.102040f, 0.101860f, 0.101700f, 0.101540f, 0.101390f, 0.101250f, 0.101120f, 0.100990f, 0.100880f, 0.100770f,
        0.100670f, 0.100580f, 0.100490f, 0.100410f, 0.100340f, 0.100270f, 0.100210f, 0.100160f, 0.100120f, 0.100100f,
        0.100100f, 0.100110f, 0.100140f, 0.100190f, 0.100260f, 0.100360f, 0.100460f, 0.100580f, 0.100710f, 0.100840f,
        0.100970f, 0.101090f, 0.101210f, 0.101310f, 0.101400f, 0.101480f, 0.101550f, 0.101620f, 0.101680f, 0.101750f,
        0.101820f, 0.101910f, 0.102010f, 0.102130f, 0.102270f, 0.102420f, 0.102590f, 0.102760f, 0.102930f, 0.103100f,
        0.103270f, 0.103410f, 0.103540f, 0.103640f, 0.103800f, 0.103930f, 0.104070f, 0.104200f, 0.104340f, 0.104480f,
        0.104610f, 0.104750f, 0.104890f, 0.105020f, 0.105160f, 0.105300f, 0.105430f, 0.105570f, 0.105710f, 0.105850f,
        0.105990f, 0.106120f, 0.106260f, 0.106400f, 0.106540f, 0.106680f, 0.106820f, 0.106960f, 0.107090f, 0.107230f,
        0.107370f, 0.107510f, 0.107650f, 0.107790f, 0.107930f, 0.108070f, 0.108210f, 0.108350f, 0.108500f, 0.108640f,
        0.108780f, 0.108920f, 0.109060f, 0.109200f, 0.109340f, 0.109490f, 0.109630f, 0.109770f, 0.109910f, 0.110050f,
        0.110200f, 0.110340f, 0.110480f, 0.110630f, 0.110770f, 0.110910f, 0.111060f, 0.111200f, 0.111350f, 0.111490f,
        0.111630f, 0.111780f, 0.111920f, 0.112070f, 0.112210f, 0.112360f, 0.112500f, 0.112650f, 0.112790f, 0.112940f,
        0.113090f, 0.113230f, 0.113380f, 0.113530f, 0.113670f, 0.113820f, 0.113970f, 0.114110f, 0.114260f, 0.114410f,
        0.114560f, 0.114700f, 0.114850f});
    static constexpr Spectrum CES36({0.034314f, 0.034654f, 0.035058f, 0.035570f, 0.036236f, 0.037075f, 0.038061f,
        0.039149f, 0.040302f, 0.041518f, 0.042834f, 0.044321f, 0.046071f, 0.048169f, 0.050668f, 0.053584f, 0.056925f,
        0.060705f, 0.064953f, 0.069717f, 0.075062f, 0.079028f, 0.083284f, 0.087829f, 0.092658f, 0.097763f, 0.103136f,
        0.108773f, 0.114675f, 0.120846f, 0.127288f, 0.134004f, 0.140999f, 0.148279f, 0.155848f, 0.163707f, 0.171860f,
        0.180305f, 0.189035f, 0.198031f, 0.207270f, 0.216728f, 0.226378f, 0.236189f, 0.246131f, 0.256177f, 0.266301f,
        0.276478f, 0.286688f, 0.296910f, 0.307123f, 0.317310f, 0.327453f, 0.337533f, 0.347522f, 0.357395f, 0.367130f,
        0.376703f, 0.386085f, 0.395235f, 0.404130f, 0.412744f, 0.421050f, 0.429029f, 0.436667f, 0.443969f, 0.450934f,
        0.457566f, 0.463867f, 0.469840f, 0.475494f, 0.480833f, 0.485866f, 0.490598f, 0.495033f, 0.499170f, 0.503006f,
        0.506530f, 0.509734f, 0.512618f, 0.515188f, 0.517463f, 0.519459f, 0.521201f, 0.522710f, 0.524002f, 0.525096f,
        0.526007f, 0.526755f, 0.527349f, 0.527811f, 0.528164f, 0.528422f, 0.528594f, 0.528680f, 0.528701f, 0.528670f,
        0.528596f, 0.528495f, 0.528391f, 0.528316f, 0.528274f, 0.528273f, 0.528315f, 0.528403f, 0.528534f, 0.528696f,
        0.528888f, 0.529106f, 0.529343f, 0.529591f, 0.529834f, 0.530062f, 0.530262f, 0.530426f, 0.530552f, 0.530643f,
        0.530711f, 0.530768f, 0.530833f, 0.530922f, 0.531046f, 0.531211f, 0.531417f, 0.531659f, 0.531933f, 0.532232f,
        0.532547f, 0.532868f, 0.533179f, 0.533468f, 0.533715f, 0.533902f, 0.534016f, 0.534054f, 0.534019f, 0.533911f,
        0.533737f, 0.533498f, 0.533195f, 0.532830f, 0.532405f, 0.531930f, 0.531421f, 0.530899f, 0.530387f, 0.529907f,
        0.529479f, 0.529113f, 0.528818f, 0.528600f, 0.528457f, 0.528382f, 0.528368f, 0.528406f, 0.528482f, 0.528580f,
        0.528684f, 0.528780f, 0.528855f, 0.528897f, 0.528894f, 0.528840f, 0.528734f, 0.528576f, 0.528370f, 0.528123f,
        0.527848f, 0.527559f, 0.527270f, 0.526998f, 0.526762f, 0.526578f, 0.526462f, 0.526430f, 0.526495f, 0.526666f,
        0.526948f, 0.527340f, 0.527839f, 0.528432f, 0.529100f, 0.529822f, 0.530578f, 0.531347f, 0.532103f, 0.532828f,
        0.533506f, 0.534125f, 0.534677f, 0.535156f, 0.535565f, 0.535912f, 0.536210f, 0.536469f, 0.536704f, 0.536931f,
        0.537164f, 0.537410f, 0.537679f, 0.537978f, 0.538313f, 0.538688f, 0.539105f, 0.539566f, 0.540071f, 0.540616f,
        0.541198f, 0.541810f, 0.542446f, 0.543099f, 0.543758f, 0.544413f, 0.545054f, 0.545668f, 0.546240f, 0.546757f,
        0.547212f, 0.547599f, 0.547917f, 0.548170f, 0.548373f, 0.548540f, 0.548687f, 0.548831f, 0.548987f, 0.549171f,
        0.549395f, 0.549664f, 0.549982f, 0.550349f, 0.550760f, 0.551208f, 0.551680f, 0.552167f, 0.552663f, 0.553161f,
        0.553660f, 0.554161f, 0.554667f, 0.555178f, 0.555689f, 0.556189f, 0.556663f, 0.557094f, 0.557462f, 0.557751f,
        0.557953f, 0.558062f, 0.558079f, 0.558008f, 0.557863f, 0.557667f, 0.557441f, 0.557216f, 0.557021f, 0.556886f,
        0.556831f, 0.556862f, 0.556979f, 0.557174f, 0.557433f, 0.557733f, 0.558053f, 0.558377f, 0.558688f, 0.558970f,
        0.559212f, 0.559409f, 0.559560f, 0.559666f, 0.559731f, 0.559768f, 0.559794f, 0.559832f, 0.559894f, 0.559995f,
        0.560142f, 0.560338f, 0.560570f, 0.560821f, 0.561081f, 0.561339f, 0.561578f, 0.561778f, 0.561924f, 0.562006f,
        0.562016f, 0.561954f, 0.561827f, 0.561651f, 0.561449f, 0.561230f, 0.560997f, 0.560751f, 0.560412f, 0.560168f,
        0.559972f, 0.559766f, 0.559481f, 0.559120f, 0.558430f, 0.557644f, 0.556779f, 0.555848f, 0.554862f, 0.553832f,
        0.552765f, 0.551671f, 0.550553f, 0.549410f, 0.548247f, 0.547062f, 0.545851f, 0.544600f, 0.543296f, 0.541927f,
        0.540482f, 0.538951f, 0.537331f, 0.535626f, 0.533850f, 0.532009f, 0.530121f, 0.528208f, 0.526304f, 0.524443f,
        0.522651f, 0.520952f, 0.519354f, 0.517858f, 0.516445f, 0.515090f, 0.513772f, 0.512467f, 0.511156f, 0.509893f,
        0.508466f, 0.506897f, 0.505230f, 0.503531f, 0.501796f, 0.500284f, 0.498773f, 0.497261f, 0.495750f, 0.494239f,
        0.492728f, 0.491217f, 0.489706f, 0.488195f, 0.486684f, 0.485174f, 0.483665f, 0.482155f, 0.480646f, 0.479137f,
        0.477628f, 0.476120f, 0.474612f, 0.473105f, 0.471599f, 0.470092f, 0.468587f, 0.467081f, 0.465577f, 0.464073f,
        0.462570f, 0.461067f, 0.459565f, 0.458064f, 0.456564f, 0.455064f, 0.453565f, 0.452067f, 0.450570f, 0.449074f,
        0.447579f, 0.446084f, 0.444591f, 0.443099f, 0.441600f, 0.440103f, 0.438607f, 0.437112f, 0.435618f, 0.434126f,
        0.432634f, 0.431144f, 0.429656f, 0.428169f, 0.426683f, 0.425199f, 0.423717f, 0.422237f, 0.420758f, 0.419281f,
        0.417806f, 0.416332f, 0.414860f});
    static constexpr Spectrum CES37({0.050280f, 0.050270f, 0.050261f, 0.050251f, 0.050242f, 0.050232f, 0.050223f,
        0.050213f, 0.050204f, 0.050195f, 0.050185f, 0.050176f, 0.050166f, 0.050157f, 0.050147f, 0.050138f, 0.050128f,
        0.050119f, 0.050110f, 0.050100f, 0.050000f, 0.050041f, 0.050068f, 0.050085f, 0.050091f, 0.050089f, 0.050080f,
        0.050065f, 0.050046f, 0.050024f, 0.050000f, 0.049976f, 0.049954f, 0.049935f, 0.049920f, 0.049911f, 0.049909f,
        0.049915f, 0.049932f, 0.049959f, 0.050000f, 0.050054f, 0.050122f, 0.050202f, 0.050292f, 0.050392f, 0.050501f,
        0.050617f, 0.050740f, 0.050868f, 0.051000f, 0.051136f, 0.051278f, 0.051429f, 0.051592f, 0.051770f, 0.051967f,
        0.052186f, 0.052428f, 0.052699f, 0.053000f, 0.053333f, 0.053696f, 0.054081f, 0.054485f, 0.054901f, 0.055326f,
        0.055753f, 0.056178f, 0.056596f, 0.057000f, 0.057387f, 0.057756f, 0.058106f, 0.058437f, 0.058749f, 0.059041f,
        0.059312f, 0.059562f, 0.059792f, 0.060000f, 0.060186f, 0.060352f, 0.060496f, 0.060621f, 0.060728f, 0.060816f,
        0.060886f, 0.060940f, 0.060978f, 0.061000f, 0.061009f, 0.061014f, 0.061022f, 0.061045f, 0.061091f, 0.061169f,
        0.061290f, 0.061462f, 0.061696f, 0.062000f, 0.062381f, 0.062834f, 0.063349f, 0.063919f, 0.064535f, 0.065188f,
        0.065869f, 0.066571f, 0.067284f, 0.068000f, 0.068709f, 0.069395f, 0.070042f, 0.070631f, 0.071145f, 0.071568f,
        0.071882f, 0.072071f, 0.072116f, 0.072000f, 0.071723f, 0.071345f, 0.070944f, 0.070598f, 0.070384f, 0.070379f,
        0.070661f, 0.071307f, 0.072394f, 0.074000f, 0.076182f, 0.078913f, 0.082148f, 0.085840f, 0.089943f, 0.094411f,
        0.099196f, 0.104250f, 0.109540f, 0.115000f, 0.120590f, 0.126260f, 0.131940f, 0.137570f, 0.143090f, 0.148450f,
        0.153580f, 0.158420f, 0.162910f, 0.167000f, 0.170640f, 0.173870f, 0.176740f, 0.179330f, 0.181680f, 0.183860f,
        0.185920f, 0.187930f, 0.189930f, 0.192000f, 0.194170f, 0.196400f, 0.198640f, 0.200830f, 0.202930f, 0.204870f,
        0.206600f, 0.208070f, 0.209220f, 0.210000f, 0.210360f, 0.210330f, 0.209920f, 0.209170f, 0.208110f, 0.206750f,
        0.205130f, 0.203290f, 0.201230f, 0.199000f, 0.196620f, 0.194120f, 0.191520f, 0.188850f, 0.186150f, 0.183440f,
        0.180730f, 0.178080f, 0.175490f, 0.173000f, 0.170630f, 0.168390f, 0.166270f, 0.164280f, 0.162420f, 0.160680f,
        0.159070f, 0.157580f, 0.156230f, 0.155000f, 0.153900f, 0.152950f, 0.152150f, 0.151510f, 0.151060f, 0.150810f,
        0.150760f, 0.150930f, 0.151340f, 0.152000f, 0.152910f, 0.154000f, 0.155220f, 0.156480f, 0.157710f, 0.158850f,
        0.159820f, 0.160550f, 0.160970f, 0.161000f, 0.160600f, 0.159800f, 0.158660f, 0.157240f, 0.155590f, 0.153780f,
        0.151850f, 0.149880f, 0.147910f, 0.146000f, 0.144200f, 0.142510f, 0.140910f, 0.139380f, 0.137920f, 0.136500f,
        0.135110f, 0.133750f, 0.132380f, 0.131000f, 0.129600f, 0.128170f, 0.126740f, 0.125300f, 0.123860f, 0.122440f,
        0.121030f, 0.119650f, 0.118300f, 0.117000f, 0.115740f, 0.114520f, 0.113340f, 0.112170f, 0.111010f, 0.109850f,
        0.108680f, 0.107490f, 0.106260f, 0.105000f, 0.103690f, 0.102360f, 0.101030f, 0.099727f, 0.098480f, 0.097315f,
        0.096258f, 0.095335f, 0.094574f, 0.094000f, 0.093639f, 0.093510f, 0.093633f, 0.094026f, 0.094708f, 0.095697f,
        0.097012f, 0.098672f, 0.100700f, 0.103100f, 0.105900f, 0.109090f, 0.112640f, 0.116550f, 0.120800f, 0.125370f,
        0.130250f, 0.135430f, 0.140880f, 0.146600f, 0.152570f, 0.158770f, 0.165210f, 0.171880f, 0.178770f, 0.185870f,
        0.193180f, 0.200700f, 0.208400f, 0.216300f, 0.224380f, 0.232630f, 0.241040f, 0.249610f, 0.258330f, 0.267180f,
        0.276160f, 0.285270f, 0.294480f, 0.303800f, 0.313210f, 0.322710f, 0.332280f, 0.341920f, 0.351620f, 0.361370f,
        0.371160f, 0.380990f, 0.390840f, 0.400700f, 0.412320f, 0.422850f, 0.433460f, 0.444120f, 0.454840f, 0.465600f,
        0.476390f, 0.487210f, 0.498030f, 0.508860f, 0.519680f, 0.530480f, 0.541250f, 0.551990f, 0.562670f, 0.573300f,
        0.583860f, 0.594340f, 0.604740f, 0.615050f, 0.625250f, 0.635340f, 0.645320f, 0.655170f, 0.664890f, 0.674470f,
        0.683900f, 0.693190f, 0.702330f, 0.711300f, 0.720110f, 0.728760f, 0.737240f, 0.745540f, 0.753670f, 0.761620f,
        0.769400f, 0.776990f, 0.784410f, 0.791640f, 0.798700f, 0.805570f, 0.812260f, 0.818780f, 0.825120f, 0.831280f,
        0.837270f, 0.843080f, 0.848730f, 0.854210f, 0.859520f, 0.864670f, 0.869660f, 0.874490f, 0.879170f, 0.883690f,
        0.888070f, 0.892300f, 0.896400f, 0.900350f, 0.904170f, 0.907860f, 0.911420f, 0.914850f, 0.918170f, 0.921360f,
        0.924440f, 0.927410f, 0.930280f, 0.933030f, 0.935690f, 0.938250f, 0.940710f, 0.943080f, 0.945360f, 0.947560f,
        0.949670f, 0.951700f, 0.953650f});
    static constexpr Spectrum CES38({0.113230f, 0.112997f, 0.113218f, 0.113931f, 0.114694f, 0.115251f, 0.115559f,
        0.115345f, 0.114266f, 0.112351f, 0.109791f, 0.106768f, 0.103524f, 0.100469f, 0.097758f, 0.095408f, 0.093312f,
        0.091402f, 0.089653f, 0.088110f, 0.086886f, 0.085770f, 0.084746f, 0.083812f, 0.082972f, 0.082241f, 0.081656f,
        0.081253f, 0.081070f, 0.081123f, 0.081421f, 0.081956f, 0.082728f, 0.083753f, 0.085097f, 0.086814f, 0.088936f,
        0.091465f, 0.094410f, 0.097757f, 0.101512f, 0.105705f, 0.110375f, 0.115547f, 0.121235f, 0.127414f, 0.134031f,
        0.141047f, 0.148443f, 0.156198f, 0.164312f, 0.172801f, 0.181655f, 0.190844f, 0.200349f, 0.210142f, 0.220191f,
        0.230490f, 0.241021f, 0.251752f, 0.262648f, 0.273683f, 0.284826f, 0.296077f, 0.307436f, 0.318892f, 0.330418f,
        0.341963f, 0.353471f, 0.364903f, 0.376225f, 0.387415f, 0.398481f, 0.409400f, 0.420120f, 0.430598f, 0.440799f,
        0.450684f, 0.460254f, 0.469511f, 0.478444f, 0.487061f, 0.495372f, 0.503345f, 0.510981f, 0.518293f, 0.525276f,
        0.531938f, 0.538319f, 0.544430f, 0.550291f, 0.555929f, 0.561330f, 0.566472f, 0.571363f, 0.575990f, 0.580352f,
        0.584485f, 0.588400f, 0.592100f, 0.595617f, 0.598954f, 0.602103f, 0.605079f, 0.607880f, 0.610494f, 0.612926f,
        0.615168f, 0.617214f, 0.619084f, 0.620773f, 0.622269f, 0.623602f, 0.624807f, 0.625919f, 0.626983f, 0.628050f,
        0.629135f, 0.630269f, 0.631458f, 0.632704f, 0.634029f, 0.635463f, 0.636974f, 0.638535f, 0.640113f, 0.641670f,
        0.643159f, 0.644582f, 0.645926f, 0.647191f, 0.648359f, 0.649408f, 0.650311f, 0.651073f, 0.651709f, 0.652249f,
        0.652711f, 0.653110f, 0.653449f, 0.653734f, 0.653950f, 0.654111f, 0.654242f, 0.654358f, 0.654470f, 0.654599f,
        0.654745f, 0.654915f, 0.655127f, 0.655376f, 0.655664f, 0.656026f, 0.656450f, 0.656907f, 0.657392f, 0.657892f,
        0.658369f, 0.658809f, 0.659203f, 0.659541f, 0.659831f, 0.660072f, 0.660223f, 0.660286f, 0.660268f, 0.660163f,
        0.659977f, 0.659758f, 0.659509f, 0.659234f, 0.658941f, 0.658629f, 0.658301f, 0.658009f, 0.657789f, 0.657654f,
        0.657615f, 0.657675f, 0.657804f, 0.657988f, 0.658223f, 0.658493f, 0.658781f, 0.659090f, 0.659395f, 0.659677f,
        0.659928f, 0.660158f, 0.660346f, 0.660493f, 0.660605f, 0.660691f, 0.660769f, 0.660857f, 0.660942f, 0.661013f,
        0.661075f, 0.661128f, 0.661169f, 0.661224f, 0.661312f, 0.661429f, 0.661566f, 0.661725f, 0.661901f, 0.662100f,
        0.662331f, 0.662588f, 0.662862f, 0.663148f, 0.663438f, 0.663725f, 0.664015f, 0.664306f, 0.664581f, 0.664827f,
        0.665038f, 0.665209f, 0.665339f, 0.665444f, 0.665534f, 0.665604f, 0.665671f, 0.665756f, 0.665853f, 0.665972f,
        0.666133f, 0.666329f, 0.666560f, 0.666844f, 0.667176f, 0.667552f, 0.667973f, 0.668402f, 0.668811f, 0.669212f,
        0.669602f, 0.669959f, 0.670290f, 0.670597f, 0.670864f, 0.671071f, 0.671221f, 0.671319f, 0.671382f, 0.671408f,
        0.671404f, 0.671371f, 0.671304f, 0.671187f, 0.671040f, 0.670880f, 0.670721f, 0.670584f, 0.670487f, 0.670419f,
        0.670378f, 0.670359f, 0.670362f, 0.670409f, 0.670511f, 0.670645f, 0.670796f, 0.670950f, 0.671076f, 0.671146f,
        0.671174f, 0.671169f, 0.671133f, 0.671061f, 0.670952f, 0.670808f, 0.670637f, 0.670441f, 0.670241f, 0.670053f,
        0.669883f, 0.669726f, 0.669587f, 0.669443f, 0.669272f, 0.669063f, 0.668844f, 0.668632f, 0.668458f, 0.668346f,
        0.668307f, 0.668326f, 0.668392f, 0.668510f, 0.668710f, 0.669008f, 0.669394f, 0.669837f, 0.670283f, 0.670612f,
        0.670871f, 0.671101f, 0.671307f, 0.671493f, 0.671745f, 0.671963f, 0.672163f, 0.672366f, 0.672582f, 0.672836f,
        0.673150f, 0.673508f, 0.673888f, 0.674283f, 0.674689f, 0.675103f, 0.675525f, 0.675951f, 0.676380f, 0.676810f,
        0.677238f, 0.677666f, 0.678109f, 0.678573f, 0.679083f, 0.679671f, 0.680346f, 0.681093f, 0.681900f, 0.682734f,
        0.683556f, 0.684358f, 0.685145f, 0.685924f, 0.686703f, 0.687470f, 0.688197f, 0.688876f, 0.689522f, 0.690160f,
        0.690888f, 0.691667f, 0.692464f, 0.693278f, 0.694111f, 0.694861f, 0.695610f, 0.696359f, 0.697106f, 0.697852f,
        0.698597f, 0.699340f, 0.700083f, 0.700825f, 0.701566f, 0.702306f, 0.703045f, 0.703783f, 0.704519f, 0.705255f,
        0.705989f, 0.706723f, 0.707455f, 0.708186f, 0.708916f, 0.709646f, 0.710374f, 0.711101f, 0.711827f, 0.712551f,
        0.713275f, 0.713998f, 0.714719f, 0.715440f, 0.716159f, 0.716877f, 0.717594f, 0.718310f, 0.719025f, 0.719739f,
        0.720452f, 0.721163f, 0.721874f, 0.722583f, 0.723299f, 0.724013f, 0.724726f, 0.725437f, 0.726147f, 0.726855f,
        0.727562f, 0.728268f, 0.728971f, 0.729674f, 0.730375f, 0.731074f, 0.731772f, 0.732468f, 0.733163f, 0.733857f,
        0.734549f, 0.735238f, 0.735927f});
    static constexpr Spectrum CES39({0.045658f, 0.045673f, 0.045687f, 0.045702f, 0.045716f, 0.045731f, 0.045745f,
        0.045760f, 0.045775f, 0.045789f, 0.045804f, 0.045818f, 0.045833f, 0.045848f, 0.045862f, 0.045877f, 0.045892f,
        0.045906f, 0.045921f, 0.045935f, 0.046062f, 0.046028f, 0.046002f, 0.045984f, 0.045973f, 0.045969f, 0.045973f,
        0.045984f, 0.046003f, 0.046029f, 0.046062f, 0.046102f, 0.046150f, 0.046204f, 0.046266f, 0.046335f, 0.046411f,
        0.046494f, 0.046583f, 0.046680f, 0.046783f, 0.046894f, 0.047011f, 0.047135f, 0.047267f, 0.047407f, 0.047554f,
        0.047710f, 0.047873f, 0.048045f, 0.048226f, 0.048415f, 0.048613f, 0.048817f, 0.049029f, 0.049246f, 0.049468f,
        0.049694f, 0.049923f, 0.050156f, 0.050390f, 0.050626f, 0.050863f, 0.051104f, 0.051349f, 0.051600f, 0.051857f,
        0.052120f, 0.052393f, 0.052674f, 0.052966f, 0.053269f, 0.053579f, 0.053894f, 0.054211f, 0.054526f, 0.054835f,
        0.055136f, 0.055425f, 0.055699f, 0.055954f, 0.056189f, 0.056405f, 0.056607f, 0.056798f, 0.056981f, 0.057161f,
        0.057340f, 0.057523f, 0.057712f, 0.057912f, 0.058128f, 0.058366f, 0.058638f, 0.058952f, 0.059317f, 0.059743f,
        0.060239f, 0.060814f, 0.061478f, 0.062240f, 0.063106f, 0.064068f, 0.065115f, 0.066236f, 0.067419f, 0.068654f,
        0.069929f, 0.071234f, 0.072556f, 0.073885f, 0.075210f, 0.076522f, 0.077815f, 0.079079f, 0.080307f, 0.081491f,
        0.082624f, 0.083697f, 0.084702f, 0.085632f, 0.086481f, 0.087252f, 0.087950f, 0.088583f, 0.089155f, 0.089671f,
        0.090138f, 0.090562f, 0.090947f, 0.091299f, 0.091623f, 0.091916f, 0.092172f, 0.092386f, 0.092555f, 0.092672f,
        0.092734f, 0.092735f, 0.092671f, 0.092536f, 0.092328f, 0.092051f, 0.091710f, 0.091311f, 0.090859f, 0.090359f,
        0.089818f, 0.089242f, 0.088634f, 0.088002f, 0.087350f, 0.086682f, 0.086000f, 0.085308f, 0.084608f, 0.083903f,
        0.083197f, 0.082492f, 0.081791f, 0.081098f, 0.080414f, 0.079739f, 0.079073f, 0.078415f, 0.077765f, 0.077122f,
        0.076485f, 0.075854f, 0.075227f, 0.074606f, 0.073989f, 0.073381f, 0.072785f, 0.072206f, 0.071650f, 0.071120f,
        0.070620f, 0.070156f, 0.069731f, 0.069350f, 0.069017f, 0.068730f, 0.068485f, 0.068279f, 0.068110f, 0.067974f,
        0.067868f, 0.067790f, 0.067735f, 0.067702f, 0.067687f, 0.067692f, 0.067718f, 0.067766f, 0.067839f, 0.067939f,
        0.068066f, 0.068222f, 0.068409f, 0.068629f, 0.068883f, 0.069167f, 0.069481f, 0.069820f, 0.070182f, 0.070564f,
        0.070963f, 0.071377f, 0.071802f, 0.072236f, 0.072675f, 0.073118f, 0.073562f, 0.074003f, 0.074440f, 0.074869f,
        0.075289f, 0.075695f, 0.076087f, 0.076461f, 0.076815f, 0.077149f, 0.077465f, 0.077763f, 0.078045f, 0.078312f,
        0.078564f, 0.078802f, 0.079028f, 0.079243f, 0.079447f, 0.079641f, 0.079826f, 0.080001f, 0.080168f, 0.080327f,
        0.080479f, 0.080623f, 0.080760f, 0.080892f, 0.081017f, 0.081138f, 0.081252f, 0.081362f, 0.081467f, 0.081567f,
        0.081662f, 0.081753f, 0.081840f, 0.081922f, 0.082001f, 0.082079f, 0.082158f, 0.082240f, 0.082329f, 0.082426f,
        0.082535f, 0.082657f, 0.082795f, 0.082953f, 0.083131f, 0.083327f, 0.083541f, 0.083770f, 0.084011f, 0.084262f,
        0.084520f, 0.084785f, 0.085053f, 0.085323f, 0.085591f, 0.085856f, 0.086116f, 0.086368f, 0.086610f, 0.086840f,
        0.087055f, 0.087253f, 0.087432f, 0.087590f, 0.087725f, 0.087838f, 0.087931f, 0.088007f, 0.088068f, 0.088115f,
        0.088150f, 0.088176f, 0.088195f, 0.088208f, 0.088218f, 0.088224f, 0.088228f, 0.088229f, 0.088229f, 0.088227f,
        0.088223f, 0.088219f, 0.088213f, 0.088208f, 0.088203f, 0.088197f, 0.088191f, 0.088184f, 0.088175f, 0.088166f,
        0.088154f, 0.088140f, 0.088124f, 0.088105f, 0.088083f, 0.088057f, 0.088028f, 0.087995f, 0.087957f, 0.087915f,
        0.087868f, 0.087815f, 0.087757f, 0.087693f, 0.087727f, 0.087694f, 0.087661f, 0.087629f, 0.087596f, 0.087563f,
        0.087530f, 0.087497f, 0.087464f, 0.087432f, 0.087399f, 0.087366f, 0.087333f, 0.087301f, 0.087268f, 0.087235f,
        0.087203f, 0.087170f, 0.087137f, 0.087105f, 0.087072f, 0.087039f, 0.087007f, 0.086974f, 0.086941f, 0.086909f,
        0.086876f, 0.086844f, 0.086811f, 0.086779f, 0.086746f, 0.086713f, 0.086681f, 0.086648f, 0.086616f, 0.086583f,
        0.086551f, 0.086519f, 0.086486f, 0.086454f, 0.086421f, 0.086389f, 0.086356f, 0.086324f, 0.086292f, 0.086259f,
        0.086227f, 0.086195f, 0.086162f, 0.086130f, 0.086098f, 0.086065f, 0.086033f, 0.086001f, 0.085968f, 0.085936f,
        0.085904f, 0.085872f, 0.085839f, 0.085807f, 0.085775f, 0.085743f, 0.085711f, 0.085678f, 0.085646f, 0.085614f,
        0.085582f, 0.085550f, 0.085518f, 0.085485f, 0.085453f, 0.085421f, 0.085389f, 0.085357f, 0.085325f, 0.085293f,
        0.085261f, 0.085229f, 0.085197f});
    static constexpr Spectrum CES40({0.075004f, 0.075007f, 0.075009f, 0.075011f, 0.075014f, 0.075016f, 0.075018f,
        0.075020f, 0.075023f, 0.075025f, 0.075027f, 0.075030f, 0.075032f, 0.075034f, 0.075036f, 0.075039f, 0.075041f,
        0.075043f, 0.075046f, 0.075048f, 0.075224f, 0.075152f, 0.075091f, 0.075041f, 0.075004f, 0.074977f, 0.074963f,
        0.074959f, 0.074968f, 0.074987f, 0.075018f, 0.075060f, 0.075114f, 0.075179f, 0.075255f, 0.075342f, 0.075441f,
        0.075550f, 0.075671f, 0.075803f, 0.075945f, 0.076099f, 0.076264f, 0.076441f, 0.076628f, 0.076828f, 0.077039f,
        0.077263f, 0.077498f, 0.077746f, 0.078006f, 0.078279f, 0.078565f, 0.078863f, 0.079174f, 0.079498f, 0.079834f,
        0.080183f, 0.080544f, 0.080918f, 0.081304f, 0.081703f, 0.082117f, 0.082548f, 0.082999f, 0.083472f, 0.083969f,
        0.084493f, 0.085047f, 0.085631f, 0.086250f, 0.086903f, 0.087580f, 0.088269f, 0.088959f, 0.089637f, 0.090291f,
        0.090910f, 0.091481f, 0.091993f, 0.092433f, 0.092794f, 0.093087f, 0.093326f, 0.093528f, 0.093707f, 0.093880f,
        0.094060f, 0.094264f, 0.094507f, 0.094803f, 0.095171f, 0.095639f, 0.096236f, 0.096992f, 0.097937f, 0.099101f,
        0.100510f, 0.102210f, 0.104210f, 0.106550f, 0.109250f, 0.112280f, 0.115610f, 0.119200f, 0.123020f, 0.127030f,
        0.131200f, 0.135480f, 0.139850f, 0.144270f, 0.148700f, 0.153120f, 0.157520f, 0.161870f, 0.166170f, 0.170370f,
        0.174480f, 0.178470f, 0.182310f, 0.186000f, 0.189510f, 0.192860f, 0.196050f, 0.199080f, 0.201980f, 0.204740f,
        0.207380f, 0.209900f, 0.212320f, 0.214650f, 0.216880f, 0.218990f, 0.220960f, 0.222760f, 0.224350f, 0.225710f,
        0.226810f, 0.227620f, 0.228110f, 0.228250f, 0.228020f, 0.227440f, 0.226550f, 0.225370f, 0.223930f, 0.222260f,
        0.220400f, 0.218370f, 0.216200f, 0.213930f, 0.211570f, 0.209150f, 0.206670f, 0.204140f, 0.201580f, 0.198990f,
        0.196390f, 0.193770f, 0.191170f, 0.188580f, 0.186010f, 0.183460f, 0.180940f, 0.178440f, 0.175960f, 0.173500f,
        0.171060f, 0.168630f, 0.166230f, 0.163840f, 0.161480f, 0.159140f, 0.156840f, 0.154610f, 0.152440f, 0.150360f,
        0.148370f, 0.146490f, 0.144740f, 0.143130f, 0.141670f, 0.140350f, 0.139160f, 0.138100f, 0.137160f, 0.136330f,
        0.135610f, 0.134970f, 0.134430f, 0.133960f, 0.133570f, 0.133250f, 0.133020f, 0.132880f, 0.132860f, 0.132940f,
        0.133140f, 0.133480f, 0.133960f, 0.134580f, 0.135360f, 0.136280f, 0.137330f, 0.138500f, 0.139780f, 0.141150f,
        0.142610f, 0.144140f, 0.145720f, 0.147360f, 0.149020f, 0.150710f, 0.152410f, 0.154110f, 0.155790f, 0.157460f,
        0.159080f, 0.160660f, 0.162190f, 0.163640f, 0.165010f, 0.166310f, 0.167540f, 0.168700f, 0.169800f, 0.170840f,
        0.171830f, 0.172780f, 0.173690f, 0.174560f, 0.175400f, 0.176210f, 0.176990f, 0.177740f, 0.178460f, 0.179140f,
        0.179800f, 0.180420f, 0.181010f, 0.181570f, 0.182090f, 0.182590f, 0.183060f, 0.183500f, 0.183910f, 0.184300f,
        0.184680f, 0.185030f, 0.185370f, 0.185690f, 0.186000f, 0.186310f, 0.186610f, 0.186930f, 0.187250f, 0.187600f,
        0.187970f, 0.188370f, 0.188810f, 0.189300f, 0.189830f, 0.190410f, 0.191040f, 0.191690f, 0.192380f, 0.193100f,
        0.193830f, 0.194580f, 0.195340f, 0.196100f, 0.196860f, 0.197610f, 0.198350f, 0.199070f, 0.199760f, 0.200420f,
        0.201040f, 0.201610f, 0.202130f, 0.202590f, 0.202980f, 0.203320f, 0.203590f, 0.203820f, 0.204000f, 0.204130f,
        0.204230f, 0.204290f, 0.204330f, 0.204340f, 0.204340f, 0.204310f, 0.204270f, 0.204210f, 0.204130f, 0.204040f,
        0.203930f, 0.203810f, 0.203670f, 0.203520f, 0.203350f, 0.203170f, 0.202980f, 0.202770f, 0.202540f, 0.202300f,
        0.202040f, 0.201760f, 0.201460f, 0.201150f, 0.200810f, 0.200460f, 0.200080f, 0.199680f, 0.199270f, 0.198820f,
        0.198360f, 0.197870f, 0.197360f, 0.196820f, 0.196780f, 0.196390f, 0.196010f, 0.195620f, 0.195240f, 0.194860f,
        0.194470f, 0.194090f, 0.193710f, 0.193330f, 0.192950f, 0.192570f, 0.192190f, 0.191810f, 0.191440f, 0.191060f,
        0.190680f, 0.190310f, 0.189930f, 0.189560f, 0.189180f, 0.188810f, 0.188430f, 0.188060f, 0.187690f, 0.187320f,
        0.186950f, 0.186580f, 0.186210f, 0.185840f, 0.185470f, 0.185100f, 0.184730f, 0.184370f, 0.184000f, 0.183640f,
        0.183270f, 0.182910f, 0.182540f, 0.182180f, 0.181810f, 0.181450f, 0.181090f, 0.180730f, 0.180370f, 0.180010f,
        0.179650f, 0.179290f, 0.178930f, 0.178570f, 0.178220f, 0.177860f, 0.177500f, 0.177150f, 0.176790f, 0.176440f,
        0.176080f, 0.175730f, 0.175380f, 0.175020f, 0.174670f, 0.174320f, 0.173970f, 0.173620f, 0.173270f, 0.172920f,
        0.172570f, 0.172230f, 0.171880f, 0.171530f, 0.171190f, 0.170840f, 0.170490f, 0.170150f, 0.169810f, 0.169460f,
        0.169120f, 0.168780f, 0.168440f});
    static constexpr Spectrum CES41({0.379730f, 0.381200f, 0.382680f, 0.384160f, 0.385630f, 0.387120f, 0.388600f,
        0.390080f, 0.391570f, 0.393060f, 0.394550f, 0.396040f, 0.397540f, 0.399040f, 0.400530f, 0.402040f, 0.403540f,
        0.405040f, 0.406550f, 0.408060f, 0.409800f, 0.410940f, 0.412370f, 0.414020f, 0.415830f, 0.417750f, 0.419720f,
        0.421680f, 0.423560f, 0.425330f, 0.426900f, 0.428250f, 0.429390f, 0.430390f, 0.431270f, 0.432090f, 0.432890f,
        0.433720f, 0.434610f, 0.435630f, 0.436800f, 0.438170f, 0.439730f, 0.441480f, 0.443390f, 0.445460f, 0.447670f,
        0.450010f, 0.452470f, 0.455040f, 0.457700f, 0.460440f, 0.463270f, 0.466190f, 0.469200f, 0.472310f, 0.475520f,
        0.478840f, 0.482280f, 0.485830f, 0.489500f, 0.493300f, 0.497220f, 0.501280f, 0.505460f, 0.509770f, 0.514220f,
        0.518790f, 0.523490f, 0.528330f, 0.533300f, 0.538400f, 0.543630f, 0.548970f, 0.554420f, 0.559960f, 0.565590f,
        0.571300f, 0.577080f, 0.582910f, 0.588800f, 0.594720f, 0.600680f, 0.606640f, 0.612600f, 0.618540f, 0.624450f,
        0.630320f, 0.636120f, 0.641860f, 0.647500f, 0.653040f, 0.658470f, 0.663770f, 0.668940f, 0.673960f, 0.678830f,
        0.683530f, 0.688040f, 0.692370f, 0.696500f, 0.700420f, 0.704130f, 0.707620f, 0.710910f, 0.713980f, 0.716850f,
        0.719500f, 0.721940f, 0.724180f, 0.726200f, 0.728020f, 0.729630f, 0.731050f, 0.732290f, 0.733360f, 0.734260f,
        0.735010f, 0.735610f, 0.736070f, 0.736400f, 0.736610f, 0.736710f, 0.736700f, 0.736590f, 0.736380f, 0.736090f,
        0.735720f, 0.735280f, 0.734770f, 0.734200f, 0.733580f, 0.732900f, 0.732180f, 0.731420f, 0.730610f, 0.729770f,
        0.728890f, 0.727990f, 0.727060f, 0.726100f, 0.725120f, 0.724130f, 0.723110f, 0.722070f, 0.721020f, 0.719950f,
        0.718860f, 0.717760f, 0.716640f, 0.715500f, 0.714350f, 0.713190f, 0.712020f, 0.710840f, 0.709670f, 0.708500f,
        0.707330f, 0.706170f, 0.705030f, 0.703900f, 0.702790f, 0.701700f, 0.700640f, 0.699610f, 0.698600f, 0.697630f,
        0.696690f, 0.695790f, 0.694920f, 0.694100f, 0.693320f, 0.692580f, 0.691890f, 0.691240f, 0.690620f, 0.690040f,
        0.689500f, 0.689000f, 0.688530f, 0.688100f, 0.687700f, 0.687320f, 0.686970f, 0.686630f, 0.686300f, 0.685980f,
        0.685650f, 0.685320f, 0.684970f, 0.684600f, 0.684210f, 0.683790f, 0.683360f, 0.682900f, 0.682430f, 0.681940f,
        0.681450f, 0.680940f, 0.680420f, 0.679900f, 0.679380f, 0.678850f, 0.678330f, 0.677810f, 0.677290f, 0.676790f,
        0.676290f, 0.675810f, 0.675350f, 0.674900f, 0.674470f, 0.674070f, 0.673690f, 0.673330f, 0.673000f, 0.672700f,
        0.672430f, 0.672190f, 0.671980f, 0.671800f, 0.671660f, 0.671550f, 0.671460f, 0.671410f, 0.671370f, 0.671350f,
        0.671350f, 0.671360f, 0.671380f, 0.671400f, 0.671430f, 0.671460f, 0.671500f, 0.671560f, 0.671640f, 0.671740f,
        0.671870f, 0.672040f, 0.672250f, 0.672500f, 0.672800f, 0.673160f, 0.673580f, 0.674070f, 0.674640f, 0.675290f,
        0.676040f, 0.676880f, 0.677830f, 0.678900f, 0.680090f, 0.681400f, 0.682850f, 0.684430f, 0.686160f, 0.688030f,
        0.690050f, 0.692240f, 0.694590f, 0.697100f, 0.699790f, 0.702640f, 0.705660f, 0.708840f, 0.712180f, 0.715660f,
        0.719290f, 0.723060f, 0.726960f, 0.731000f, 0.735160f, 0.739440f, 0.743820f, 0.748290f, 0.752840f, 0.757460f,
        0.762140f, 0.766860f, 0.771620f, 0.776400f, 0.781190f, 0.785980f, 0.790760f, 0.795520f, 0.800250f, 0.804930f,
        0.809560f, 0.814120f, 0.818600f, 0.823000f, 0.827300f, 0.831500f, 0.835590f, 0.839570f, 0.843430f, 0.847180f,
        0.850810f, 0.854300f, 0.857670f, 0.860900f, 0.863990f, 0.866960f, 0.869800f, 0.872540f, 0.875180f, 0.877720f,
        0.880190f, 0.882580f, 0.884920f, 0.887200f, 0.889440f, 0.891610f, 0.893700f, 0.895680f, 0.897540f, 0.899240f,
        0.900770f, 0.902100f, 0.903220f, 0.904100f, 0.905490f, 0.906650f, 0.907810f, 0.908950f, 0.910080f, 0.911200f,
        0.912310f, 0.913400f, 0.914480f, 0.915550f, 0.916610f, 0.917650f, 0.918680f, 0.919700f, 0.920710f, 0.921710f,
        0.922690f, 0.923670f, 0.924630f, 0.925580f, 0.926520f, 0.927450f, 0.928370f, 0.929280f, 0.930180f, 0.931060f,
        0.931940f, 0.932810f, 0.933660f, 0.934510f, 0.935340f, 0.936170f, 0.936990f, 0.937790f, 0.938590f, 0.939370f,
        0.940150f, 0.940920f, 0.941680f, 0.942430f, 0.943170f, 0.943900f, 0.944630f, 0.945340f, 0.946040f, 0.946740f,
        0.947430f, 0.948110f, 0.948780f, 0.949450f, 0.950100f, 0.950750f, 0.951390f, 0.952020f, 0.952640f, 0.953260f,
        0.953870f, 0.954470f, 0.955060f, 0.955650f, 0.956220f, 0.956800f, 0.957360f, 0.957920f, 0.958470f, 0.959010f,
        0.959550f, 0.960080f, 0.960600f, 0.961120f, 0.961630f, 0.962130f, 0.962630f, 0.963120f, 0.963600f, 0.964080f,
        0.964560f, 0.965020f, 0.965480f});
    static constexpr Spectrum CES42({0.065979f, 0.065857f, 0.065736f, 0.065615f, 0.065494f, 0.065373f, 0.065252f,
        0.065132f, 0.065012f, 0.064892f, 0.064772f, 0.064652f, 0.064533f, 0.064414f, 0.064295f, 0.064176f, 0.064057f,
        0.063939f, 0.063821f, 0.063703f, 0.063744f, 0.063565f, 0.063393f, 0.063228f, 0.063071f, 0.062924f, 0.062787f,
        0.062661f, 0.062546f, 0.062443f, 0.062354f, 0.062279f, 0.062218f, 0.062172f, 0.062143f, 0.062131f, 0.062137f,
        0.062161f, 0.062205f, 0.062269f, 0.062354f, 0.062460f, 0.062587f, 0.062732f, 0.062894f, 0.063072f, 0.063264f,
        0.063468f, 0.063684f, 0.063908f, 0.064141f, 0.064381f, 0.064627f, 0.064879f, 0.065137f, 0.065402f, 0.065674f,
        0.065951f, 0.066235f, 0.066526f, 0.066822f, 0.067125f, 0.067437f, 0.067760f, 0.068095f, 0.068445f, 0.068812f,
        0.069199f, 0.069607f, 0.070038f, 0.070496f, 0.070981f, 0.071499f, 0.072052f, 0.072645f, 0.073281f, 0.073964f,
        0.074699f, 0.075488f, 0.076336f, 0.077248f, 0.078224f, 0.079265f, 0.080369f, 0.081532f, 0.082754f, 0.084032f,
        0.085363f, 0.086746f, 0.088179f, 0.089659f, 0.091185f, 0.092763f, 0.094398f, 0.096095f, 0.097861f, 0.099702f,
        0.101620f, 0.103630f, 0.105730f, 0.107930f, 0.110230f, 0.112670f, 0.115260f, 0.118020f, 0.120990f, 0.124200f,
        0.127650f, 0.131390f, 0.135430f, 0.139800f, 0.144530f, 0.149640f, 0.155170f, 0.161150f, 0.167610f, 0.174590f,
        0.182100f, 0.190190f, 0.198880f, 0.208210f, 0.218200f, 0.228810f, 0.240000f, 0.251740f, 0.263990f, 0.276700f,
        0.289820f, 0.303330f, 0.317180f, 0.331330f, 0.345720f, 0.360260f, 0.374810f, 0.389270f, 0.403520f, 0.417440f,
        0.430910f, 0.443810f, 0.456030f, 0.467460f, 0.477980f, 0.487610f, 0.496350f, 0.504220f, 0.511240f, 0.517410f,
        0.522770f, 0.527330f, 0.531090f, 0.534080f, 0.536320f, 0.537850f, 0.538710f, 0.538950f, 0.538610f, 0.537750f,
        0.536400f, 0.534620f, 0.532440f, 0.529910f, 0.527080f, 0.523960f, 0.520570f, 0.516920f, 0.513040f, 0.508930f,
        0.504620f, 0.500110f, 0.495430f, 0.490590f, 0.485610f, 0.480500f, 0.475290f, 0.470000f, 0.464660f, 0.459270f,
        0.453870f, 0.448470f, 0.443100f, 0.437770f, 0.432510f, 0.427310f, 0.422180f, 0.417110f, 0.412110f, 0.407160f,
        0.402270f, 0.397430f, 0.392650f, 0.387930f, 0.383250f, 0.378640f, 0.374100f, 0.369640f, 0.365270f, 0.361010f,
        0.356870f, 0.352850f, 0.348970f, 0.345230f, 0.341650f, 0.338220f, 0.334940f, 0.331790f, 0.328770f, 0.325880f,
        0.323090f, 0.320410f, 0.317830f, 0.315340f, 0.312940f, 0.310610f, 0.308360f, 0.306170f, 0.304040f, 0.301970f,
        0.299940f, 0.297960f, 0.296010f, 0.294100f, 0.292200f, 0.290330f, 0.288480f, 0.286640f, 0.284800f, 0.282970f,
        0.281150f, 0.279310f, 0.277480f, 0.275630f, 0.273770f, 0.271890f, 0.270020f, 0.268130f, 0.266250f, 0.264370f,
        0.262500f, 0.260640f, 0.258790f, 0.256960f, 0.255160f, 0.253390f, 0.251660f, 0.249990f, 0.248390f, 0.246860f,
        0.245430f, 0.244090f, 0.242870f, 0.241770f, 0.240810f, 0.239980f, 0.239300f, 0.238770f, 0.238400f, 0.238190f,
        0.238150f, 0.238280f, 0.238590f, 0.239090f, 0.239780f, 0.240650f, 0.241700f, 0.242940f, 0.244360f, 0.245950f,
        0.247710f, 0.249640f, 0.251730f, 0.253980f, 0.256390f, 0.258940f, 0.261580f, 0.264310f, 0.267090f, 0.269900f,
        0.272710f, 0.275500f, 0.278230f, 0.280890f, 0.283450f, 0.285910f, 0.288260f, 0.290510f, 0.292670f, 0.294720f,
        0.296670f, 0.298530f, 0.300280f, 0.301940f, 0.303510f, 0.305010f, 0.306470f, 0.307920f, 0.309390f, 0.310900f,
        0.312500f, 0.314200f, 0.316030f, 0.318030f, 0.320210f, 0.322600f, 0.325190f, 0.328000f, 0.331030f, 0.334290f,
        0.337790f, 0.341540f, 0.345540f, 0.349800f, 0.354330f, 0.359140f, 0.364230f, 0.369610f, 0.375290f, 0.381290f,
        0.387590f, 0.394220f, 0.401180f, 0.408480f, 0.409770f, 0.415280f, 0.420800f, 0.426350f, 0.431910f, 0.437490f,
        0.443090f, 0.448700f, 0.454320f, 0.459960f, 0.465600f, 0.471260f, 0.476920f, 0.482590f, 0.488260f, 0.493940f,
        0.499610f, 0.505290f, 0.510970f, 0.516640f, 0.522310f, 0.527970f, 0.533630f, 0.539270f, 0.544910f, 0.550540f,
        0.556150f, 0.561750f, 0.567330f, 0.572900f, 0.578450f, 0.583970f, 0.589480f, 0.594960f, 0.600420f, 0.605860f,
        0.611270f, 0.616650f, 0.622010f, 0.627330f, 0.632620f, 0.637890f, 0.643120f, 0.648310f, 0.653470f, 0.658600f,
        0.663680f, 0.668730f, 0.673740f, 0.678720f, 0.683650f, 0.688540f, 0.693390f, 0.698200f, 0.702960f, 0.707680f,
        0.712360f, 0.716990f, 0.721570f, 0.726110f, 0.730600f, 0.735050f, 0.739450f, 0.743800f, 0.748100f, 0.752360f,
        0.756570f, 0.760730f, 0.764830f, 0.768890f, 0.772910f, 0.776870f, 0.780780f, 0.784640f, 0.788450f, 0.792220f,
        0.795930f, 0.799590f, 0.803210f});
    static constexpr Spectrum CES43({0.088623f, 0.087423f, 0.086238f, 0.085067f, 0.083911f, 0.082769f, 0.081641f,
        0.080528f, 0.079428f, 0.078342f, 0.077269f, 0.076210f, 0.075164f, 0.074132f, 0.073112f, 0.072106f, 0.071112f,
        0.070131f, 0.069162f, 0.068206f, 0.067100f, 0.066418f, 0.065564f, 0.064569f, 0.063462f, 0.062275f, 0.061038f,
        0.059781f, 0.058536f, 0.057332f, 0.056200f, 0.055165f, 0.054228f, 0.053385f, 0.052630f, 0.051959f, 0.051368f,
        0.050851f, 0.050404f, 0.050022f, 0.049700f, 0.049434f, 0.049221f, 0.049055f, 0.048933f, 0.048850f, 0.048804f,
        0.048790f, 0.048804f, 0.048842f, 0.048900f, 0.048975f, 0.049070f, 0.049187f, 0.049330f, 0.049502f, 0.049705f,
        0.049943f, 0.050220f, 0.050538f, 0.050900f, 0.051310f, 0.051775f, 0.052299f, 0.052891f, 0.053555f, 0.054300f,
        0.055130f, 0.056052f, 0.057074f, 0.058200f, 0.059436f, 0.060782f, 0.062232f, 0.063786f, 0.065439f, 0.067189f,
        0.069033f, 0.070969f, 0.072992f, 0.075100f, 0.077297f, 0.079611f, 0.082078f, 0.084734f, 0.087613f, 0.090751f,
        0.094183f, 0.097945f, 0.102070f, 0.106600f, 0.111550f, 0.116910f, 0.122650f, 0.128740f, 0.135150f, 0.141850f,
        0.148810f, 0.156010f, 0.163420f, 0.171000f, 0.178740f, 0.186680f, 0.194840f, 0.203280f, 0.212030f, 0.221130f,
        0.230630f, 0.240560f, 0.250970f, 0.261900f, 0.273370f, 0.285320f, 0.297680f, 0.310360f, 0.323300f, 0.336410f,
        0.349630f, 0.362870f, 0.376050f, 0.389100f, 0.401950f, 0.414540f, 0.426820f, 0.438730f, 0.450240f, 0.461270f,
        0.471790f, 0.481730f, 0.491050f, 0.499700f, 0.507630f, 0.514870f, 0.521430f, 0.527340f, 0.532630f, 0.537330f,
        0.541460f, 0.545040f, 0.548120f, 0.550700f, 0.552820f, 0.554510f, 0.555780f, 0.556680f, 0.557220f, 0.557430f,
        0.557330f, 0.556960f, 0.556340f, 0.555500f, 0.554460f, 0.553230f, 0.551830f, 0.550260f, 0.548540f, 0.546670f,
        0.544670f, 0.542560f, 0.540330f, 0.538000f, 0.535580f, 0.533080f, 0.530490f, 0.527820f, 0.525060f, 0.522240f,
        0.519330f, 0.516360f, 0.513310f, 0.510200f, 0.507020f, 0.503790f, 0.500500f, 0.497180f, 0.493820f, 0.490430f,
        0.487030f, 0.483620f, 0.480210f, 0.476800f, 0.473410f, 0.470030f, 0.466670f, 0.463340f, 0.460020f, 0.456740f,
        0.453480f, 0.450250f, 0.447060f, 0.443900f, 0.440780f, 0.437690f, 0.434640f, 0.431610f, 0.428610f, 0.425640f,
        0.422680f, 0.419740f, 0.416810f, 0.413900f, 0.410990f, 0.408070f, 0.405120f, 0.402120f, 0.399050f, 0.395890f,
        0.392630f, 0.389240f, 0.385700f, 0.382000f, 0.378130f, 0.374110f, 0.370000f, 0.365820f, 0.361630f, 0.357470f,
        0.353370f, 0.349380f, 0.345550f, 0.341900f, 0.338480f, 0.335300f, 0.332350f, 0.329640f, 0.327170f, 0.324930f,
        0.322940f, 0.321180f, 0.319670f, 0.318400f, 0.317370f, 0.316550f, 0.315930f, 0.315460f, 0.315120f, 0.314880f,
        0.314720f, 0.314610f, 0.314510f, 0.314400f, 0.314260f, 0.314080f, 0.313870f, 0.313630f, 0.313360f, 0.313070f,
        0.312750f, 0.312420f, 0.312070f, 0.311700f, 0.311320f, 0.310930f, 0.310530f, 0.310110f, 0.309680f, 0.309230f,
        0.308770f, 0.308300f, 0.307810f, 0.307300f, 0.306770f, 0.306170f, 0.305460f, 0.304580f, 0.303490f, 0.302130f,
        0.300470f, 0.298440f, 0.296000f, 0.293100f, 0.289710f, 0.285860f, 0.281600f, 0.276990f, 0.272080f, 0.266920f,
        0.261560f, 0.256050f, 0.250450f, 0.244800f, 0.239160f, 0.233590f, 0.228130f, 0.222840f, 0.217770f, 0.212980f,
        0.208510f, 0.204420f, 0.200770f, 0.197600f, 0.194960f, 0.192810f, 0.191120f, 0.189850f, 0.188940f, 0.188360f,
        0.188070f, 0.188030f, 0.188180f, 0.188500f, 0.188960f, 0.189660f, 0.190710f, 0.192250f, 0.194370f, 0.197210f,
        0.200880f, 0.205510f, 0.211210f, 0.218100f, 0.226230f, 0.235360f, 0.245180f, 0.255390f, 0.265670f, 0.275710f,
        0.285200f, 0.293830f, 0.301310f, 0.307300f, 0.316930f, 0.325300f, 0.333780f, 0.342370f, 0.351070f, 0.359860f,
        0.368750f, 0.377730f, 0.386800f, 0.395940f, 0.405160f, 0.414440f, 0.423790f, 0.433190f, 0.442640f, 0.452130f,
        0.461650f, 0.471210f, 0.480780f, 0.490370f, 0.499970f, 0.509560f, 0.519150f, 0.528730f, 0.538280f, 0.547800f,
        0.557300f, 0.566740f, 0.576150f, 0.585490f, 0.594780f, 0.603990f, 0.613140f, 0.622200f, 0.631180f, 0.640070f,
        0.648870f, 0.657570f, 0.666160f, 0.674640f, 0.683010f, 0.691260f, 0.699390f, 0.707400f, 0.715280f, 0.723040f,
        0.730660f, 0.738150f, 0.745500f, 0.752710f, 0.759790f, 0.766720f, 0.773520f, 0.780170f, 0.786690f, 0.793060f,
        0.799290f, 0.805380f, 0.811320f, 0.817130f, 0.822790f, 0.828320f, 0.833710f, 0.838970f, 0.844090f, 0.849070f,
        0.853930f, 0.858650f, 0.863240f, 0.867710f, 0.872060f, 0.876280f, 0.880380f, 0.884370f, 0.888230f, 0.891990f,
        0.895630f, 0.899170f, 0.902590f});
    static constexpr Spectrum CES44({0.041634f, 0.041743f, 0.041853f, 0.041964f, 0.042074f, 0.042185f, 0.042296f,
        0.042408f, 0.042520f, 0.042632f, 0.042744f, 0.042857f, 0.042969f, 0.043083f, 0.043196f, 0.043310f, 0.043424f,
        0.043538f, 0.043653f, 0.043767f, 0.043900f, 0.043988f, 0.044097f, 0.044223f, 0.044362f, 0.044509f, 0.044658f,
        0.044807f, 0.044950f, 0.045082f, 0.045200f, 0.045299f, 0.045381f, 0.045449f, 0.045503f, 0.045548f, 0.045585f,
        0.045616f, 0.045644f, 0.045671f, 0.045700f, 0.045732f, 0.045767f, 0.045805f, 0.045845f, 0.045887f, 0.045930f,
        0.045973f, 0.046016f, 0.046059f, 0.046100f, 0.046140f, 0.046178f, 0.046216f, 0.046253f, 0.046291f, 0.046330f,
        0.046369f, 0.046410f, 0.046454f, 0.046500f, 0.046549f, 0.046600f, 0.046653f, 0.046707f, 0.046761f, 0.046814f,
        0.046865f, 0.046913f, 0.046959f, 0.047000f, 0.047037f, 0.047069f, 0.047097f, 0.047120f, 0.047141f, 0.047158f,
        0.047172f, 0.047184f, 0.047193f, 0.047200f, 0.047206f, 0.047211f, 0.047218f, 0.047226f, 0.047239f, 0.047256f,
        0.047279f, 0.047310f, 0.047350f, 0.047400f, 0.047461f, 0.047532f, 0.047612f, 0.047699f, 0.047792f, 0.047890f,
        0.047991f, 0.048094f, 0.048197f, 0.048300f, 0.048401f, 0.048500f, 0.048598f, 0.048695f, 0.048793f, 0.048891f,
        0.048990f, 0.049091f, 0.049194f, 0.049300f, 0.049409f, 0.049521f, 0.049634f, 0.049748f, 0.049862f, 0.049975f,
        0.050086f, 0.050194f, 0.050299f, 0.050400f, 0.050496f, 0.050586f, 0.050670f, 0.050750f, 0.050823f, 0.050891f,
        0.050952f, 0.051008f, 0.051057f, 0.051100f, 0.051137f, 0.051167f, 0.051191f, 0.051209f, 0.051221f, 0.051228f,
        0.051229f, 0.051224f, 0.051215f, 0.051200f, 0.051180f, 0.051156f, 0.051129f, 0.051099f, 0.051067f, 0.051033f,
        0.050999f, 0.050965f, 0.050932f, 0.050900f, 0.050870f, 0.050843f, 0.050816f, 0.050789f, 0.050762f, 0.050734f,
        0.050705f, 0.050673f, 0.050638f, 0.050600f, 0.050558f, 0.050512f, 0.050463f, 0.050412f, 0.050360f, 0.050307f,
        0.050253f, 0.050201f, 0.050149f, 0.050100f, 0.050053f, 0.050008f, 0.049963f, 0.049919f, 0.049873f, 0.049826f,
        0.049776f, 0.049722f, 0.049664f, 0.049600f, 0.049530f, 0.049456f, 0.049378f, 0.049299f, 0.049221f, 0.049145f,
        0.049074f, 0.049008f, 0.048949f, 0.048900f, 0.048861f, 0.048833f, 0.048812f, 0.048799f, 0.048792f, 0.048789f,
        0.048790f, 0.048793f, 0.048797f, 0.048800f, 0.048802f, 0.048802f, 0.048802f, 0.048801f, 0.048799f, 0.048798f,
        0.048797f, 0.048797f, 0.048798f, 0.048800f, 0.048804f, 0.048809f, 0.048814f, 0.048819f, 0.048823f, 0.048825f,
        0.048825f, 0.048821f, 0.048813f, 0.048800f, 0.048782f, 0.048760f, 0.048735f, 0.048709f, 0.048683f, 0.048658f,
        0.048636f, 0.048618f, 0.048605f, 0.048600f, 0.048603f, 0.048613f, 0.048630f, 0.048654f, 0.048684f, 0.048719f,
        0.048759f, 0.048803f, 0.048850f, 0.048900f, 0.048952f, 0.049007f, 0.049063f, 0.049121f, 0.049181f, 0.049242f,
        0.049305f, 0.049369f, 0.049434f, 0.049500f, 0.049567f, 0.049635f, 0.049703f, 0.049773f, 0.049842f, 0.049913f,
        0.049984f, 0.050056f, 0.050128f, 0.050200f, 0.050273f, 0.050346f, 0.050421f, 0.050497f, 0.050574f, 0.050654f,
        0.050736f, 0.050821f, 0.050909f, 0.051000f, 0.051095f, 0.051193f, 0.051293f, 0.051394f, 0.051497f, 0.051600f,
        0.051703f, 0.051804f, 0.051903f, 0.052000f, 0.052094f, 0.052183f, 0.052269f, 0.052349f, 0.052424f, 0.052494f,
        0.052556f, 0.052612f, 0.052660f, 0.052700f, 0.052732f, 0.052755f, 0.052770f, 0.052779f, 0.052780f, 0.052775f,
        0.052764f, 0.052748f, 0.052726f, 0.052700f, 0.052670f, 0.052636f, 0.052598f, 0.052559f, 0.052517f, 0.052474f,
        0.052431f, 0.052387f, 0.052343f, 0.052300f, 0.052259f, 0.052219f, 0.052181f, 0.052146f, 0.052113f, 0.052084f,
        0.052057f, 0.052034f, 0.052015f, 0.052000f, 0.051975f, 0.051954f, 0.051933f, 0.051912f, 0.051891f, 0.051870f,
        0.051849f, 0.051828f, 0.051807f, 0.051787f, 0.051766f, 0.051745f, 0.051724f, 0.051703f, 0.051682f, 0.051661f,
        0.051640f, 0.051620f, 0.051599f, 0.051578f, 0.051557f, 0.051536f, 0.051516f, 0.051495f, 0.051474f, 0.051453f,
        0.051432f, 0.051412f, 0.051391f, 0.051370f, 0.051349f, 0.051329f, 0.051308f, 0.051287f, 0.051266f, 0.051246f,
        0.051225f, 0.051204f, 0.051184f, 0.051163f, 0.051142f, 0.051122f, 0.051101f, 0.051080f, 0.051060f, 0.051039f,
        0.051019f, 0.050998f, 0.050977f, 0.050957f, 0.050936f, 0.050916f, 0.050895f, 0.050875f, 0.050854f, 0.050833f,
        0.050813f, 0.050792f, 0.050772f, 0.050751f, 0.050731f, 0.050710f, 0.050690f, 0.050669f, 0.050649f, 0.050628f,
        0.050608f, 0.050588f, 0.050567f, 0.050547f, 0.050526f, 0.050506f, 0.050485f, 0.050465f, 0.050445f, 0.050424f,
        0.050404f, 0.050383f, 0.050363f});
    static constexpr Spectrum CES45({0.055400f, 0.054600f, 0.055000f, 0.055186f, 0.055200f, 0.055186f, 0.055257f,
        0.055314f, 0.055643f, 0.055871f, 0.056086f, 0.056300f, 0.056586f, 0.057043f, 0.057243f, 0.057229f, 0.057486f,
        0.057700f, 0.057643f, 0.057657f, 0.057600f, 0.057557f, 0.057600f, 0.057500f, 0.057500f, 0.057729f, 0.057843f,
        0.057771f, 0.057986f, 0.058171f, 0.058414f, 0.058586f, 0.058786f, 0.058814f, 0.058914f, 0.059000f, 0.059400f,
        0.059771f, 0.060100f, 0.060300f, 0.060757f, 0.061257f, 0.061543f, 0.061614f, 0.061714f, 0.061814f, 0.062186f,
        0.062600f, 0.062800f, 0.063186f, 0.063643f, 0.063914f, 0.064200f, 0.064529f, 0.064857f, 0.065357f, 0.065814f,
        0.066214f, 0.066914f, 0.067686f, 0.068414f, 0.069000f, 0.069500f, 0.069929f, 0.070371f, 0.070629f, 0.071086f,
        0.071429f, 0.072000f, 0.072571f, 0.073329f, 0.074300f, 0.075400f, 0.076157f, 0.077257f, 0.078314f, 0.079414f,
        0.080586f, 0.081957f, 0.083143f, 0.084700f, 0.086186f, 0.087786f, 0.089543f, 0.091457f, 0.093114f, 0.095071f,
        0.097029f, 0.099043f, 0.101030f, 0.103130f, 0.105110f, 0.107410f, 0.109890f, 0.112530f, 0.115230f, 0.118230f,
        0.121060f, 0.124030f, 0.127270f, 0.130640f, 0.133810f, 0.137600f, 0.141240f, 0.145270f, 0.149610f, 0.154060f,
        0.158330f, 0.163340f, 0.168230f, 0.173590f, 0.179030f, 0.184700f, 0.190410f, 0.197200f, 0.203690f, 0.210590f,
        0.217630f, 0.225070f, 0.232260f, 0.239960f, 0.246900f, 0.254090f, 0.260840f, 0.267790f, 0.274140f, 0.280800f,
        0.286870f, 0.293000f, 0.298660f, 0.304430f, 0.309260f, 0.314010f, 0.318200f, 0.322060f, 0.325370f, 0.328490f,
        0.330900f, 0.333590f, 0.335910f, 0.337560f, 0.338990f, 0.340490f, 0.341370f, 0.342070f, 0.342530f, 0.342740f,
        0.343060f, 0.343040f, 0.342700f, 0.342390f, 0.342090f, 0.341230f, 0.340210f, 0.339040f, 0.337970f, 0.336630f,
        0.335330f, 0.333590f, 0.332010f, 0.330440f, 0.328910f, 0.326910f, 0.325290f, 0.323240f, 0.321470f, 0.319940f,
        0.318460f, 0.316270f, 0.314500f, 0.312290f, 0.309900f, 0.307460f, 0.305000f, 0.301890f, 0.299140f, 0.295790f,
        0.292490f, 0.289700f, 0.287400f, 0.284790f, 0.283070f, 0.281410f, 0.279760f, 0.277810f, 0.276390f, 0.274310f,
        0.271860f, 0.269010f, 0.266140f, 0.263070f, 0.260570f, 0.257590f, 0.254260f, 0.251290f, 0.248610f, 0.245660f,
        0.243090f, 0.240190f, 0.237160f, 0.234260f, 0.231210f, 0.227860f, 0.225040f, 0.222190f, 0.219430f, 0.216690f,
        0.213960f, 0.211170f, 0.208630f, 0.205760f, 0.202740f, 0.199610f, 0.196730f, 0.193800f, 0.191140f, 0.188630f,
        0.186200f, 0.183970f, 0.182090f, 0.179900f, 0.177960f, 0.176170f, 0.174460f, 0.172870f, 0.171590f, 0.170070f,
        0.169060f, 0.168160f, 0.167030f, 0.165710f, 0.164740f, 0.163670f, 0.162600f, 0.161490f, 0.160470f, 0.159630f,
        0.159070f, 0.158430f, 0.157790f, 0.157390f, 0.156940f, 0.156270f, 0.155600f, 0.155030f, 0.154540f, 0.154310f,
        0.154100f, 0.153660f, 0.153610f, 0.153640f, 0.153410f, 0.153200f, 0.152810f, 0.152260f, 0.152040f, 0.151700f,
        0.151170f, 0.150730f, 0.150460f, 0.150240f, 0.150300f, 0.150290f, 0.150230f, 0.150330f, 0.150390f, 0.150390f,
        0.150470f, 0.150440f, 0.150390f, 0.150370f, 0.150330f, 0.150210f, 0.150100f, 0.150240f, 0.150410f, 0.150590f,
        0.150900f, 0.151530f, 0.151960f, 0.152340f, 0.152490f, 0.152800f, 0.153200f, 0.153590f, 0.153670f, 0.154730f,
        0.155910f, 0.157100f, 0.158130f, 0.159240f, 0.160410f, 0.161730f, 0.162530f, 0.163230f, 0.163940f, 0.164700f,
        0.165360f, 0.166100f, 0.166840f, 0.167370f, 0.167860f, 0.168470f, 0.168960f, 0.169030f, 0.169140f, 0.169600f,
        0.170600f, 0.171470f, 0.172300f, 0.173300f, 0.174660f, 0.175840f, 0.176660f, 0.177040f, 0.177730f, 0.178400f,
        0.178940f, 0.179500f, 0.179830f, 0.179870f, 0.180130f, 0.180200f, 0.179940f, 0.180000f, 0.180200f, 0.180260f,
        0.180540f, 0.180500f, 0.180300f, 0.180000f, 0.179730f, 0.179160f, 0.178590f, 0.178060f, 0.177730f, 0.177370f,
        0.177100f, 0.176500f, 0.175960f, 0.175860f, 0.175700f, 0.175460f, 0.175540f, 0.175990f, 0.176440f, 0.177060f,
        0.177460f, 0.177810f, 0.178170f, 0.178470f, 0.178610f, 0.179040f, 0.179690f, 0.180370f, 0.181160f, 0.182230f,
        0.183390f, 0.184770f, 0.186160f, 0.187470f, 0.188760f, 0.190410f, 0.192030f, 0.193640f, 0.195260f, 0.196870f,
        0.198470f, 0.200300f, 0.201930f, 0.203190f, 0.204490f, 0.205740f, 0.206890f, 0.208140f, 0.209300f, 0.210270f,
        0.211340f, 0.212340f, 0.213130f, 0.214070f, 0.214970f, 0.215810f, 0.216630f, 0.217670f, 0.218640f, 0.219730f,
        0.220490f, 0.220460f, 0.220200f, 0.219970f, 0.219500f, 0.218790f, 0.218070f, 0.217660f, 0.217840f, 0.218110f,
        0.218400f, 0.218340f, 0.218670f});
    static constexpr Spectrum CES46({0.099822f, 0.100420f, 0.101030f, 0.101640f, 0.102250f, 0.102870f, 0.103490f,
        0.104110f, 0.104730f, 0.105360f, 0.105990f, 0.106630f, 0.107270f, 0.107910f, 0.108560f, 0.109200f, 0.109860f,
        0.110510f, 0.111170f, 0.111830f, 0.112600f, 0.113110f, 0.113750f, 0.114480f, 0.115290f, 0.116150f, 0.117030f,
        0.117890f, 0.118730f, 0.119510f, 0.120200f, 0.120790f, 0.121280f, 0.121680f, 0.122020f, 0.122310f, 0.122560f,
        0.122790f, 0.123010f, 0.123250f, 0.123500f, 0.123790f, 0.124120f, 0.124480f, 0.124880f, 0.125320f, 0.125790f,
        0.126300f, 0.126830f, 0.127400f, 0.128000f, 0.128630f, 0.129290f, 0.129990f, 0.130740f, 0.131520f, 0.132350f,
        0.133230f, 0.134160f, 0.135150f, 0.136200f, 0.137310f, 0.138490f, 0.139740f, 0.141060f, 0.142450f, 0.143930f,
        0.145490f, 0.147130f, 0.148870f, 0.150700f, 0.152630f, 0.154660f, 0.156790f, 0.159040f, 0.161410f, 0.163900f,
        0.166520f, 0.169270f, 0.172160f, 0.175200f, 0.178380f, 0.181720f, 0.185200f, 0.188830f, 0.192610f, 0.196530f,
        0.200600f, 0.204820f, 0.209190f, 0.213700f, 0.218350f, 0.223140f, 0.228060f, 0.233090f, 0.238220f, 0.243450f,
        0.248760f, 0.254150f, 0.259600f, 0.265100f, 0.270640f, 0.276210f, 0.281790f, 0.287350f, 0.292890f, 0.298380f,
        0.303800f, 0.309140f, 0.314380f, 0.319500f, 0.324490f, 0.329330f, 0.334030f, 0.338580f, 0.342970f, 0.347200f,
        0.351260f, 0.355150f, 0.358870f, 0.362400f, 0.365750f, 0.368920f, 0.371900f, 0.374710f, 0.377350f, 0.379810f,
        0.382100f, 0.384230f, 0.386200f, 0.388000f, 0.389650f, 0.391140f, 0.392490f, 0.393700f, 0.394770f, 0.395710f,
        0.396530f, 0.397230f, 0.397820f, 0.398300f, 0.398680f, 0.398960f, 0.399160f, 0.399270f, 0.399310f, 0.399280f,
        0.399190f, 0.399040f, 0.398840f, 0.398600f, 0.398320f, 0.398000f, 0.397650f, 0.397260f, 0.396820f, 0.396340f,
        0.395820f, 0.395260f, 0.394650f, 0.394000f, 0.393300f, 0.392550f, 0.391760f, 0.390920f, 0.390030f, 0.389100f,
        0.388120f, 0.387090f, 0.386020f, 0.384900f, 0.383740f, 0.382530f, 0.381280f, 0.380010f, 0.378700f, 0.377370f,
        0.376020f, 0.374660f, 0.373280f, 0.371900f, 0.370510f, 0.369120f, 0.367730f, 0.366320f, 0.364900f, 0.363460f,
        0.362010f, 0.360530f, 0.359030f, 0.357500f, 0.355940f, 0.354350f, 0.352730f, 0.351070f, 0.349380f, 0.347650f,
        0.345880f, 0.344060f, 0.342200f, 0.340300f, 0.338350f, 0.336350f, 0.334290f, 0.332180f, 0.330000f, 0.327760f,
        0.325460f, 0.323080f, 0.320630f, 0.318100f, 0.315490f, 0.312810f, 0.310080f, 0.307290f, 0.304470f, 0.301630f,
        0.298760f, 0.295900f, 0.293040f, 0.290200f, 0.287390f, 0.284610f, 0.281880f, 0.279190f, 0.276560f, 0.273980f,
        0.271470f, 0.269040f, 0.266680f, 0.264400f, 0.262210f, 0.260110f, 0.258100f, 0.256170f, 0.254330f, 0.252560f,
        0.250880f, 0.249280f, 0.247750f, 0.246300f, 0.244920f, 0.243620f, 0.242370f, 0.241180f, 0.240040f, 0.238950f,
        0.237890f, 0.236870f, 0.235880f, 0.234900f, 0.233940f, 0.233000f, 0.232070f, 0.231150f, 0.230250f, 0.229370f,
        0.228500f, 0.227650f, 0.226820f, 0.226000f, 0.225200f, 0.224410f, 0.223650f, 0.222900f, 0.222180f, 0.221470f,
        0.220790f, 0.220140f, 0.219510f, 0.218900f, 0.218320f, 0.217770f, 0.217260f, 0.216770f, 0.216330f, 0.215920f,
        0.215550f, 0.215220f, 0.214940f, 0.214700f, 0.214510f, 0.214370f, 0.214270f, 0.214220f, 0.214220f, 0.214250f,
        0.214330f, 0.214450f, 0.214610f, 0.214800f, 0.215030f, 0.215300f, 0.215600f, 0.215920f, 0.216280f, 0.216660f,
        0.217070f, 0.217490f, 0.217940f, 0.218400f, 0.218880f, 0.219370f, 0.219880f, 0.220420f, 0.220980f, 0.221560f,
        0.222170f, 0.222810f, 0.223490f, 0.224200f, 0.224940f, 0.225710f, 0.226490f, 0.227250f, 0.227990f, 0.228690f,
        0.229340f, 0.229920f, 0.230410f, 0.230800f, 0.231420f, 0.231950f, 0.232490f, 0.233020f, 0.233560f, 0.234090f,
        0.234630f, 0.235160f, 0.235700f, 0.236240f, 0.236780f, 0.237320f, 0.237860f, 0.238410f, 0.238950f, 0.239490f,
        0.240040f, 0.240580f, 0.241130f, 0.241680f, 0.242230f, 0.242770f, 0.243320f, 0.243880f, 0.244430f, 0.244980f,
        0.245530f, 0.246090f, 0.246640f, 0.247200f, 0.247750f, 0.248310f, 0.248870f, 0.249430f, 0.249990f, 0.250550f,
        0.251110f, 0.251670f, 0.252240f, 0.252800f, 0.253370f, 0.253930f, 0.254500f, 0.255070f, 0.255640f, 0.256210f,
        0.256780f, 0.257350f, 0.257920f, 0.258490f, 0.259060f, 0.259640f, 0.260210f, 0.260790f, 0.261370f, 0.261940f,
        0.262520f, 0.263100f, 0.263680f, 0.264260f, 0.264840f, 0.265420f, 0.266010f, 0.266590f, 0.267180f, 0.267760f,
        0.268350f, 0.268940f, 0.269520f, 0.270110f, 0.270700f, 0.271290f, 0.271880f, 0.272480f, 0.273070f, 0.273660f,
        0.274260f, 0.274850f, 0.275450f});
    static constexpr Spectrum CES47({0.142540f, 0.142640f, 0.142740f, 0.142840f, 0.142950f, 0.143050f, 0.143150f,
        0.143250f, 0.143350f, 0.143460f, 0.143560f, 0.143660f, 0.143760f, 0.143870f, 0.143970f, 0.144070f, 0.144170f,
        0.144280f, 0.144380f, 0.144480f, 0.144600f, 0.144680f, 0.144770f, 0.144890f, 0.145010f, 0.145140f, 0.145280f,
        0.145420f, 0.145550f, 0.145680f, 0.145800f, 0.145910f, 0.146010f, 0.146130f, 0.146280f, 0.146460f, 0.146690f,
        0.146990f, 0.147360f, 0.147830f, 0.148400f, 0.149080f, 0.149880f, 0.150790f, 0.151810f, 0.152940f, 0.154170f,
        0.155500f, 0.156940f, 0.158470f, 0.160100f, 0.161820f, 0.163640f, 0.165560f, 0.167580f, 0.169700f, 0.171930f,
        0.174280f, 0.176730f, 0.179310f, 0.182000f, 0.184820f, 0.187760f, 0.190830f, 0.194030f, 0.197370f, 0.200850f,
        0.204460f, 0.208230f, 0.212140f, 0.216200f, 0.220420f, 0.224790f, 0.229320f, 0.234010f, 0.238860f, 0.243870f,
        0.249050f, 0.254400f, 0.259910f, 0.265600f, 0.271460f, 0.277480f, 0.283670f, 0.290010f, 0.296500f, 0.303130f,
        0.309900f, 0.316810f, 0.323840f, 0.331000f, 0.338270f, 0.345630f, 0.353050f, 0.360520f, 0.368000f, 0.375480f,
        0.382920f, 0.390310f, 0.397610f, 0.404800f, 0.411860f, 0.418750f, 0.425450f, 0.431910f, 0.438110f, 0.444010f,
        0.449590f, 0.454800f, 0.459610f, 0.464000f, 0.467930f, 0.471420f, 0.474470f, 0.477100f, 0.479330f, 0.481160f,
        0.482610f, 0.483690f, 0.484410f, 0.484800f, 0.484860f, 0.484620f, 0.484090f, 0.483310f, 0.482290f, 0.481060f,
        0.479650f, 0.478060f, 0.476340f, 0.474500f, 0.472560f, 0.470540f, 0.468430f, 0.466260f, 0.464020f, 0.461730f,
        0.459400f, 0.457020f, 0.454620f, 0.452200f, 0.449760f, 0.447320f, 0.444870f, 0.442410f, 0.439950f, 0.437490f,
        0.435030f, 0.432580f, 0.430130f, 0.427700f, 0.425280f, 0.422870f, 0.420480f, 0.418120f, 0.415770f, 0.413460f,
        0.411170f, 0.408910f, 0.406690f, 0.404500f, 0.402350f, 0.400240f, 0.398170f, 0.396130f, 0.394120f, 0.392150f,
        0.390200f, 0.388270f, 0.386380f, 0.384500f, 0.382650f, 0.380810f, 0.379000f, 0.377210f, 0.375430f, 0.373680f,
        0.371960f, 0.370250f, 0.368560f, 0.366900f, 0.365260f, 0.363640f, 0.362050f, 0.360480f, 0.358950f, 0.357450f,
        0.355980f, 0.354550f, 0.353150f, 0.351800f, 0.350490f, 0.349220f, 0.348010f, 0.346850f, 0.345750f, 0.344710f,
        0.343750f, 0.342850f, 0.342040f, 0.341300f, 0.340650f, 0.340070f, 0.339560f, 0.339110f, 0.338700f, 0.338330f,
        0.337980f, 0.337650f, 0.337330f, 0.337000f, 0.336660f, 0.336300f, 0.335910f, 0.335500f, 0.335060f, 0.334580f,
        0.334060f, 0.333490f, 0.332870f, 0.332200f, 0.331470f, 0.330690f, 0.329860f, 0.328990f, 0.328080f, 0.327140f,
        0.326180f, 0.325200f, 0.324200f, 0.323200f, 0.322200f, 0.321210f, 0.320250f, 0.319330f, 0.318480f, 0.317690f,
        0.317000f, 0.316410f, 0.315940f, 0.315600f, 0.315410f, 0.315380f, 0.315530f, 0.315850f, 0.316360f, 0.317070f,
        0.317990f, 0.319120f, 0.320490f, 0.322100f, 0.323960f, 0.326060f, 0.328430f, 0.331060f, 0.333950f, 0.337100f,
        0.340540f, 0.344240f, 0.348230f, 0.352500f, 0.357060f, 0.361900f, 0.367020f, 0.372420f, 0.378090f, 0.384020f,
        0.390230f, 0.396700f, 0.403420f, 0.410400f, 0.417630f, 0.425100f, 0.432800f, 0.440710f, 0.448840f, 0.457160f,
        0.465670f, 0.474350f, 0.483200f, 0.492200f, 0.501340f, 0.510610f, 0.519980f, 0.529440f, 0.538970f, 0.548540f,
        0.558130f, 0.567740f, 0.577340f, 0.586900f, 0.596420f, 0.605870f, 0.615250f, 0.624550f, 0.633750f, 0.642840f,
        0.651810f, 0.660660f, 0.669360f, 0.677900f, 0.686280f, 0.694510f, 0.702590f, 0.710530f, 0.718340f, 0.726030f,
        0.733620f, 0.741100f, 0.748490f, 0.755800f, 0.763020f, 0.770090f, 0.776930f, 0.783450f, 0.789580f, 0.795230f,
        0.800330f, 0.804800f, 0.808550f, 0.811500f, 0.816070f, 0.819910f, 0.823700f, 0.827420f, 0.831070f, 0.834670f,
        0.838200f, 0.841680f, 0.845090f, 0.848440f, 0.851730f, 0.854970f, 0.858140f, 0.861260f, 0.864310f, 0.867320f,
        0.870260f, 0.873150f, 0.875980f, 0.878760f, 0.881490f, 0.884160f, 0.886780f, 0.889350f, 0.891870f, 0.894330f,
        0.896750f, 0.899120f, 0.901440f, 0.903710f, 0.905930f, 0.908110f, 0.910240f, 0.912330f, 0.914380f, 0.916380f,
        0.918340f, 0.920250f, 0.922130f, 0.923960f, 0.925760f, 0.927510f, 0.929230f, 0.930910f, 0.932560f, 0.934160f,
        0.935730f, 0.937270f, 0.938770f, 0.940240f, 0.941680f, 0.943080f, 0.944450f, 0.945790f, 0.947100f, 0.948380f,
        0.949630f, 0.950850f, 0.952040f, 0.953210f, 0.954350f, 0.955460f, 0.956550f, 0.957610f, 0.958640f, 0.959660f,
        0.960640f, 0.961610f, 0.962550f, 0.963470f, 0.964370f, 0.965250f, 0.966100f, 0.966940f, 0.967760f, 0.968550f,
        0.969330f, 0.970090f, 0.970830f});
    static constexpr Spectrum CES48({0.087762f, 0.087867f, 0.087989f, 0.088125f, 0.088273f, 0.088432f, 0.088600f,
        0.088775f, 0.088955f, 0.089139f, 0.089325f, 0.089511f, 0.089695f, 0.089875f, 0.090050f, 0.090218f, 0.090377f,
        0.090525f, 0.090661f, 0.090783f, 0.090888f, 0.090977f, 0.091058f, 0.091138f, 0.091228f, 0.091335f, 0.091469f,
        0.091639f, 0.091853f, 0.092121f, 0.092451f, 0.092856f, 0.093362f, 0.093998f, 0.094796f, 0.095783f, 0.096990f,
        0.098446f, 0.100180f, 0.102230f, 0.104610f, 0.107350f, 0.110440f, 0.113850f, 0.117560f, 0.121540f, 0.125770f,
        0.130230f, 0.134890f, 0.139740f, 0.144740f, 0.149880f, 0.155150f, 0.160560f, 0.166110f, 0.171790f, 0.177610f,
        0.183560f, 0.189660f, 0.195890f, 0.202270f, 0.208790f, 0.215440f, 0.222220f, 0.229130f, 0.236160f, 0.243300f,
        0.250550f, 0.257910f, 0.265370f, 0.272920f, 0.280560f, 0.288260f, 0.296010f, 0.303770f, 0.311530f, 0.319260f,
        0.326930f, 0.334530f, 0.342030f, 0.349400f, 0.356630f, 0.363710f, 0.370650f, 0.377440f, 0.384080f, 0.390580f,
        0.396940f, 0.403150f, 0.409230f, 0.415160f, 0.420950f, 0.426610f, 0.432140f, 0.437540f, 0.442820f, 0.447980f,
        0.453020f, 0.457960f, 0.462800f, 0.467530f, 0.472170f, 0.476710f, 0.481150f, 0.485480f, 0.489710f, 0.493830f,
        0.497840f, 0.501740f, 0.505520f, 0.509180f, 0.512720f, 0.516150f, 0.519470f, 0.522690f, 0.525830f, 0.528880f,
        0.531860f, 0.534770f, 0.537620f, 0.540420f, 0.543180f, 0.545890f, 0.548570f, 0.551220f, 0.553830f, 0.556410f,
        0.558960f, 0.561480f, 0.563980f, 0.566460f, 0.568920f, 0.571350f, 0.573740f, 0.576090f, 0.578380f, 0.580600f,
        0.582750f, 0.584820f, 0.586790f, 0.588660f, 0.590420f, 0.592080f, 0.593630f, 0.595100f, 0.596480f, 0.597790f,
        0.599020f, 0.600190f, 0.601300f, 0.602360f, 0.603370f, 0.604340f, 0.605240f, 0.606070f, 0.606830f, 0.607510f,
        0.608090f, 0.608570f, 0.608950f, 0.609210f, 0.609350f, 0.609350f, 0.609230f, 0.608970f, 0.608570f, 0.608020f,
        0.607320f, 0.606460f, 0.605450f, 0.604270f, 0.602930f, 0.601420f, 0.599770f, 0.597980f, 0.596070f, 0.594030f,
        0.591890f, 0.589640f, 0.587310f, 0.584900f, 0.582420f, 0.579860f, 0.577200f, 0.574450f, 0.571590f, 0.568600f,
        0.565480f, 0.562210f, 0.558800f, 0.555210f, 0.551450f, 0.547530f, 0.543470f, 0.539270f, 0.534950f, 0.530530f,
        0.526020f, 0.521440f, 0.516800f, 0.512110f, 0.507390f, 0.502640f, 0.497860f, 0.493050f, 0.488210f, 0.483330f,
        0.478430f, 0.473500f, 0.468530f, 0.463540f, 0.458520f, 0.453470f, 0.448400f, 0.443320f, 0.438240f, 0.433160f,
        0.428080f, 0.423030f, 0.417990f, 0.412980f, 0.408000f, 0.403060f, 0.398160f, 0.393290f, 0.388460f, 0.383660f,
        0.378900f, 0.374180f, 0.369500f, 0.364860f, 0.360260f, 0.355700f, 0.351190f, 0.346720f, 0.342310f, 0.337960f,
        0.333670f, 0.329430f, 0.325270f, 0.321170f, 0.317140f, 0.313190f, 0.309310f, 0.305500f, 0.301770f, 0.298110f,
        0.294520f, 0.291010f, 0.287570f, 0.284200f, 0.280910f, 0.277690f, 0.274540f, 0.271450f, 0.268440f, 0.265490f,
        0.262600f, 0.259770f, 0.257000f, 0.254290f, 0.251630f, 0.249030f, 0.246490f, 0.244010f, 0.241580f, 0.239220f,
        0.236920f, 0.234680f, 0.232510f, 0.230400f, 0.228360f, 0.226380f, 0.224460f, 0.222600f, 0.220800f, 0.219050f,
        0.217340f, 0.215680f, 0.214060f, 0.212480f, 0.210940f, 0.209430f, 0.207950f, 0.206520f, 0.205110f, 0.203740f,
        0.202410f, 0.201110f, 0.199850f, 0.198620f, 0.197430f, 0.196270f, 0.195140f, 0.194050f, 0.193000f, 0.191980f,
        0.191000f, 0.190050f, 0.189140f, 0.188270f, 0.187430f, 0.186640f, 0.185870f, 0.185140f, 0.184440f, 0.183770f,
        0.183130f, 0.182520f, 0.181930f, 0.181370f, 0.180830f, 0.180320f, 0.179830f, 0.179370f, 0.178930f, 0.178520f,
        0.178130f, 0.177770f, 0.177440f, 0.177130f, 0.176850f, 0.176600f, 0.176370f, 0.176170f, 0.175990f, 0.175850f,
        0.175720f, 0.175620f, 0.175550f, 0.175500f, 0.175470f, 0.175470f, 0.175500f, 0.175550f, 0.175640f, 0.175750f,
        0.175890f, 0.176070f, 0.176280f, 0.176520f, 0.176800f, 0.177120f, 0.177470f, 0.177860f, 0.178300f, 0.178770f,
        0.179290f, 0.179850f, 0.180460f, 0.181110f, 0.180880f, 0.181270f, 0.181650f, 0.182040f, 0.182430f, 0.182810f,
        0.183200f, 0.183590f, 0.183980f, 0.184370f, 0.184760f, 0.185150f, 0.185540f, 0.185940f, 0.186330f, 0.186720f,
        0.187120f, 0.187510f, 0.187910f, 0.188310f, 0.188700f, 0.189100f, 0.189500f, 0.189900f, 0.190300f, 0.190700f,
        0.191100f, 0.191500f, 0.191900f, 0.192300f, 0.192710f, 0.193110f, 0.193520f, 0.193920f, 0.194330f, 0.194730f,
        0.195140f, 0.195550f, 0.195960f, 0.196370f, 0.196780f, 0.197190f, 0.197600f, 0.198010f, 0.198420f, 0.198840f,
        0.199250f, 0.199660f, 0.200080f});
    static constexpr Spectrum CES49({0.000572f, 0.001977f, 0.003359f, 0.004720f, 0.006063f, 0.007391f, 0.008705f,
        0.010009f, 0.011306f, 0.012597f, 0.013886f, 0.015175f, 0.016466f, 0.017763f, 0.019067f, 0.020382f, 0.021709f,
        0.023052f, 0.024413f, 0.025795f, 0.027200f, 0.028628f, 0.030069f, 0.031511f, 0.032939f, 0.034341f, 0.035705f,
        0.037017f, 0.038265f, 0.039434f, 0.040514f, 0.041493f, 0.042373f, 0.043157f, 0.043851f, 0.044457f, 0.044981f,
        0.045426f, 0.045796f, 0.046095f, 0.046328f, 0.046499f, 0.046618f, 0.046693f, 0.046734f, 0.046751f, 0.046754f,
        0.046751f, 0.046753f, 0.046769f, 0.046809f, 0.046880f, 0.046984f, 0.047119f, 0.047286f, 0.047482f, 0.047709f,
        0.047964f, 0.048247f, 0.048557f, 0.048894f, 0.049257f, 0.049650f, 0.050074f, 0.050533f, 0.051029f, 0.051566f,
        0.052147f, 0.052774f, 0.053451f, 0.054181f, 0.054967f, 0.055814f, 0.056730f, 0.057720f, 0.058791f, 0.059950f,
        0.061203f, 0.062555f, 0.064014f, 0.065586f, 0.067276f, 0.069086f, 0.071017f, 0.073068f, 0.075241f, 0.077536f,
        0.079953f, 0.082494f, 0.085159f, 0.087948f, 0.090861f, 0.093892f, 0.097036f, 0.100280f, 0.103630f, 0.107070f,
        0.110600f, 0.114210f, 0.117880f, 0.121630f, 0.125430f, 0.129280f, 0.133140f, 0.137010f, 0.140870f, 0.144680f,
        0.148440f, 0.152130f, 0.155730f, 0.159210f, 0.162560f, 0.165790f, 0.168880f, 0.171850f, 0.174690f, 0.177400f,
        0.179990f, 0.182450f, 0.184790f, 0.187010f, 0.189110f, 0.191080f, 0.192930f, 0.194660f, 0.196270f, 0.197750f,
        0.199100f, 0.200330f, 0.201430f, 0.202410f, 0.203260f, 0.203980f, 0.204580f, 0.205060f, 0.205420f, 0.205660f,
        0.205790f, 0.205810f, 0.205730f, 0.205530f, 0.205230f, 0.204840f, 0.204350f, 0.203780f, 0.203130f, 0.202410f,
        0.201620f, 0.200780f, 0.199890f, 0.198950f, 0.197970f, 0.196960f, 0.195900f, 0.194820f, 0.193690f, 0.192530f,
        0.191330f, 0.190100f, 0.188840f, 0.187540f, 0.186210f, 0.184840f, 0.183450f, 0.182020f, 0.180560f, 0.179070f,
        0.177540f, 0.175990f, 0.174410f, 0.172790f, 0.171150f, 0.169470f, 0.167780f, 0.166080f, 0.164360f, 0.162630f,
        0.160900f, 0.159170f, 0.157450f, 0.155740f, 0.154040f, 0.152360f, 0.150700f, 0.149050f, 0.147410f, 0.145800f,
        0.144190f, 0.142600f, 0.141030f, 0.139470f, 0.137930f, 0.136400f, 0.134900f, 0.133420f, 0.131970f, 0.130540f,
        0.129160f, 0.127810f, 0.126500f, 0.125230f, 0.124010f, 0.122830f, 0.121690f, 0.120590f, 0.119510f, 0.118470f,
        0.117450f, 0.116450f, 0.115460f, 0.114490f, 0.113530f, 0.112570f, 0.111630f, 0.110690f, 0.109770f, 0.108860f,
        0.107960f, 0.107060f, 0.106190f, 0.105320f, 0.104470f, 0.103630f, 0.102800f, 0.101990f, 0.101200f, 0.100420f,
        0.099663f, 0.098923f, 0.098203f, 0.097504f, 0.096826f, 0.096171f, 0.095540f, 0.094935f, 0.094355f, 0.093804f,
        0.093282f, 0.092790f, 0.092329f, 0.091902f, 0.091508f, 0.091149f, 0.090824f, 0.090532f, 0.090275f, 0.090052f,
        0.089862f, 0.089707f, 0.089585f, 0.089497f, 0.089442f, 0.089419f, 0.089424f, 0.089455f, 0.089510f, 0.089585f,
        0.089679f, 0.089789f, 0.089911f, 0.090044f, 0.090185f, 0.090335f, 0.090493f, 0.090661f, 0.090838f, 0.091025f,
        0.091222f, 0.091430f, 0.091649f, 0.091879f, 0.092121f, 0.092374f, 0.092636f, 0.092907f, 0.093185f, 0.093470f,
        0.093759f, 0.094053f, 0.094350f, 0.094649f, 0.094949f, 0.095249f, 0.095548f, 0.095846f, 0.096143f, 0.096438f,
        0.096729f, 0.097017f, 0.097302f, 0.097581f, 0.097855f, 0.098125f, 0.098392f, 0.098656f, 0.098917f, 0.099178f,
        0.099439f, 0.099700f, 0.099963f, 0.100230f, 0.100500f, 0.100770f, 0.101050f, 0.101340f, 0.101640f, 0.101950f,
        0.102280f, 0.102620f, 0.102980f, 0.103360f, 0.103760f, 0.104180f, 0.104620f, 0.105070f, 0.105530f, 0.105990f,
        0.106470f, 0.106940f, 0.107420f, 0.107890f, 0.108360f, 0.108820f, 0.109280f, 0.109730f, 0.110190f, 0.110640f,
        0.111090f, 0.111550f, 0.112000f, 0.112460f, 0.112920f, 0.113390f, 0.113860f, 0.114350f, 0.114840f, 0.115350f,
        0.115880f, 0.116420f, 0.116980f, 0.117570f, 0.118180f, 0.118820f, 0.119480f, 0.120170f, 0.120900f, 0.121660f,
        0.122460f, 0.123290f, 0.124170f, 0.125090f, 0.125350f, 0.126070f, 0.126790f, 0.127520f, 0.128250f, 0.128980f,
        0.129720f, 0.130460f, 0.131200f, 0.131950f, 0.132700f, 0.133450f, 0.134210f, 0.134970f, 0.135740f, 0.136510f,
        0.137280f, 0.138060f, 0.138840f, 0.139620f, 0.140410f, 0.141200f, 0.141990f, 0.142790f, 0.143590f, 0.144400f,
        0.145210f, 0.146020f, 0.146840f, 0.147660f, 0.148490f, 0.149320f, 0.150150f, 0.150980f, 0.151820f, 0.152670f,
        0.153520f, 0.154370f, 0.155220f, 0.156080f, 0.156950f, 0.157810f, 0.158680f, 0.159560f, 0.160440f, 0.161320f,
        0.162210f, 0.163100f, 0.163990f});
    static constexpr Spectrum CES50({0.078830f, 0.079503f, 0.080181f, 0.080865f, 0.081553f, 0.082247f, 0.082947f,
        0.083652f, 0.084362f, 0.085077f, 0.085799f, 0.086525f, 0.087258f, 0.087996f, 0.088739f, 0.089488f, 0.090243f,
        0.091004f, 0.091770f, 0.092542f, 0.093385f, 0.094134f, 0.094898f, 0.095675f, 0.096464f, 0.097266f, 0.098079f,
        0.098902f, 0.099734f, 0.100570f, 0.101420f, 0.102280f, 0.103140f, 0.104010f, 0.104880f, 0.105750f, 0.106630f,
        0.107510f, 0.108390f, 0.109280f, 0.110160f, 0.111040f, 0.111920f, 0.112800f, 0.113690f, 0.114580f, 0.115480f,
        0.116390f, 0.117310f, 0.118240f, 0.119190f, 0.120150f, 0.121130f, 0.122120f, 0.123120f, 0.124130f, 0.125150f,
        0.126180f, 0.127220f, 0.128260f, 0.129310f, 0.130360f, 0.131420f, 0.132490f, 0.133570f, 0.134670f, 0.135780f,
        0.136920f, 0.138090f, 0.139290f, 0.140520f, 0.141790f, 0.143090f, 0.144420f, 0.145770f, 0.147140f, 0.148530f,
        0.149920f, 0.151320f, 0.152720f, 0.154120f, 0.155510f, 0.156900f, 0.158300f, 0.159730f, 0.161180f, 0.162680f,
        0.164230f, 0.165850f, 0.167530f, 0.169300f, 0.171170f, 0.173150f, 0.175260f, 0.177540f, 0.180000f, 0.182660f,
        0.185550f, 0.188690f, 0.192100f, 0.195800f, 0.199810f, 0.204100f, 0.208640f, 0.213390f, 0.218320f, 0.223410f,
        0.228600f, 0.233870f, 0.239190f, 0.244530f, 0.249840f, 0.255100f, 0.260260f, 0.265310f, 0.270190f, 0.274890f,
        0.279360f, 0.283570f, 0.287480f, 0.291070f, 0.294310f, 0.297200f, 0.299790f, 0.302090f, 0.304120f, 0.305920f,
        0.307500f, 0.308900f, 0.310130f, 0.311220f, 0.312190f, 0.313030f, 0.313750f, 0.314330f, 0.314770f, 0.315050f,
        0.315170f, 0.315120f, 0.314900f, 0.314490f, 0.313900f, 0.313120f, 0.312170f, 0.311060f, 0.309790f, 0.308380f,
        0.306840f, 0.305170f, 0.303380f, 0.301490f, 0.299500f, 0.297410f, 0.295210f, 0.292900f, 0.290470f, 0.287920f,
        0.285240f, 0.282430f, 0.279480f, 0.276380f, 0.273150f, 0.269770f, 0.266280f, 0.262670f, 0.258960f, 0.255170f,
        0.251300f, 0.247370f, 0.243390f, 0.239370f, 0.235320f, 0.231280f, 0.227250f, 0.223260f, 0.219330f, 0.215490f,
        0.211750f, 0.208130f, 0.204660f, 0.201360f, 0.198240f, 0.195300f, 0.192530f, 0.189900f, 0.187420f, 0.185060f,
        0.182820f, 0.180680f, 0.178620f, 0.176650f, 0.174740f, 0.172890f, 0.171120f, 0.169420f, 0.167800f, 0.166270f,
        0.164830f, 0.163480f, 0.162220f, 0.161070f, 0.160020f, 0.159070f, 0.158210f, 0.157440f, 0.156760f, 0.156160f,
        0.155620f, 0.155160f, 0.154760f, 0.154420f, 0.154130f, 0.153880f, 0.153660f, 0.153460f, 0.153280f, 0.153100f,
        0.152920f, 0.152720f, 0.152490f, 0.152230f, 0.151930f, 0.151590f, 0.151220f, 0.150820f, 0.150400f, 0.149970f,
        0.149540f, 0.149100f, 0.148680f, 0.148260f, 0.147870f, 0.147500f, 0.147170f, 0.146870f, 0.146620f, 0.146410f,
        0.146260f, 0.146160f, 0.146140f, 0.146180f, 0.146300f, 0.146490f, 0.146750f, 0.147070f, 0.147450f, 0.147880f,
        0.148360f, 0.148890f, 0.149450f, 0.150050f, 0.150680f, 0.151350f, 0.152050f, 0.152790f, 0.153580f, 0.154410f,
        0.155300f, 0.156240f, 0.157230f, 0.158290f, 0.159410f, 0.160580f, 0.161800f, 0.163070f, 0.164360f, 0.165690f,
        0.167030f, 0.168380f, 0.169740f, 0.171090f, 0.172430f, 0.173740f, 0.175020f, 0.176260f, 0.177430f, 0.178530f,
        0.179550f, 0.180480f, 0.181300f, 0.182010f, 0.182590f, 0.183040f, 0.183400f, 0.183650f, 0.183810f, 0.183890f,
        0.183890f, 0.183840f, 0.183740f, 0.183590f, 0.183420f, 0.183210f, 0.182970f, 0.182700f, 0.182400f, 0.182070f,
        0.181720f, 0.181350f, 0.180940f, 0.180520f, 0.180070f, 0.179600f, 0.179100f, 0.178570f, 0.178030f, 0.177450f,
        0.176850f, 0.176210f, 0.175550f, 0.174860f, 0.174140f, 0.173390f, 0.172600f, 0.171780f, 0.170930f, 0.170040f,
        0.169120f, 0.168160f, 0.167160f, 0.166130f, 0.165900f, 0.165110f, 0.164340f, 0.163560f, 0.162790f, 0.162010f,
        0.161250f, 0.160480f, 0.159720f, 0.158960f, 0.158210f, 0.157450f, 0.156700f, 0.155960f, 0.155210f, 0.154470f,
        0.153730f, 0.153000f, 0.152260f, 0.151540f, 0.150810f, 0.150080f, 0.149360f, 0.148650f, 0.147930f, 0.147220f,
        0.146510f, 0.145800f, 0.145100f, 0.144400f, 0.143700f, 0.143000f, 0.142310f, 0.141620f, 0.140930f, 0.140250f,
        0.139570f, 0.138890f, 0.138210f, 0.137540f, 0.136870f, 0.136200f, 0.135530f, 0.134870f, 0.134210f, 0.133560f,
        0.132900f, 0.132250f, 0.131600f, 0.130960f, 0.130310f, 0.129670f, 0.129030f, 0.128400f, 0.127770f, 0.127140f,
        0.126510f, 0.125880f, 0.125260f, 0.124640f, 0.124030f, 0.123410f, 0.122800f, 0.122190f, 0.121590f, 0.120980f,
        0.120380f, 0.119780f, 0.119190f, 0.118590f, 0.118000f, 0.117410f, 0.116830f, 0.116250f, 0.115660f, 0.115090f,
        0.114510f, 0.113940f, 0.113370f});
    static constexpr Spectrum CES51({0.099501f, 0.100400f, 0.101300f, 0.102210f, 0.103130f, 0.104060f, 0.104990f,
        0.105930f, 0.106880f, 0.107840f, 0.108800f, 0.109770f, 0.110750f, 0.111740f, 0.112730f, 0.113730f, 0.114740f,
        0.115760f, 0.116780f, 0.117810f, 0.118890f, 0.119920f, 0.120960f, 0.122010f, 0.123080f, 0.124160f, 0.125250f,
        0.126350f, 0.127460f, 0.128580f, 0.129710f, 0.130840f, 0.131980f, 0.133130f, 0.134280f, 0.135430f, 0.136580f,
        0.137740f, 0.138900f, 0.140060f, 0.141220f, 0.142380f, 0.143530f, 0.144690f, 0.145860f, 0.147030f, 0.148200f,
        0.149390f, 0.150590f, 0.151800f, 0.153030f, 0.154270f, 0.155530f, 0.156800f, 0.158080f, 0.159370f, 0.160660f,
        0.161960f, 0.163250f, 0.164540f, 0.165830f, 0.167110f, 0.168390f, 0.169670f, 0.170960f, 0.172260f, 0.173580f,
        0.174920f, 0.176290f, 0.177690f, 0.179130f, 0.180610f, 0.182120f, 0.183660f, 0.185220f, 0.186800f, 0.188390f,
        0.189980f, 0.191570f, 0.193150f, 0.194710f, 0.196250f, 0.197780f, 0.199320f, 0.200870f, 0.202440f, 0.204060f,
        0.205730f, 0.207460f, 0.209280f, 0.211180f, 0.213190f, 0.215340f, 0.217640f, 0.220130f, 0.222830f, 0.225780f,
        0.229000f, 0.232520f, 0.236360f, 0.240560f, 0.245130f, 0.250040f, 0.255250f, 0.260720f, 0.266400f, 0.272260f,
        0.278260f, 0.284340f, 0.290480f, 0.296630f, 0.302740f, 0.308780f, 0.314690f, 0.320440f, 0.325970f, 0.331240f,
        0.336210f, 0.340830f, 0.345050f, 0.348830f, 0.352140f, 0.354990f, 0.357430f, 0.359490f, 0.361190f, 0.362590f,
        0.363700f, 0.364570f, 0.365220f, 0.365700f, 0.366030f, 0.366210f, 0.366250f, 0.366130f, 0.365860f, 0.365420f,
        0.364830f, 0.364060f, 0.363130f, 0.362030f, 0.360750f, 0.359300f, 0.357690f, 0.355920f, 0.354010f, 0.351950f,
        0.349770f, 0.347450f, 0.345020f, 0.342480f, 0.339830f, 0.337060f, 0.334180f, 0.331160f, 0.328000f, 0.324700f,
        0.321240f, 0.317610f, 0.313810f, 0.309830f, 0.305660f, 0.301310f, 0.296810f, 0.292180f, 0.287420f, 0.282570f,
        0.277630f, 0.272620f, 0.267570f, 0.262490f, 0.257400f, 0.252330f, 0.247320f, 0.242390f, 0.237570f, 0.232910f,
        0.228420f, 0.224140f, 0.220100f, 0.216340f, 0.212880f, 0.209700f, 0.206790f, 0.204120f, 0.201670f, 0.199430f,
        0.197380f, 0.195490f, 0.193750f, 0.192130f, 0.190620f, 0.189220f, 0.187940f, 0.186780f, 0.185750f, 0.184860f,
        0.184120f, 0.183520f, 0.183080f, 0.182800f, 0.182690f, 0.182730f, 0.182930f, 0.183270f, 0.183730f, 0.184320f,
        0.185010f, 0.185810f, 0.186700f, 0.187660f, 0.188700f, 0.189780f, 0.190880f, 0.191970f, 0.193040f, 0.194040f,
        0.194970f, 0.195780f, 0.196470f, 0.196990f, 0.197340f, 0.197520f, 0.197550f, 0.197460f, 0.197250f, 0.196950f,
        0.196580f, 0.196150f, 0.195690f, 0.195210f, 0.194720f, 0.194250f, 0.193800f, 0.193390f, 0.193020f, 0.192720f,
        0.192480f, 0.192330f, 0.192280f, 0.192330f, 0.192490f, 0.192770f, 0.193150f, 0.193620f, 0.194180f, 0.194810f,
        0.195500f, 0.196250f, 0.197050f, 0.197880f, 0.198750f, 0.199650f, 0.200600f, 0.201590f, 0.202640f, 0.203750f,
        0.204920f, 0.206160f, 0.207490f, 0.208900f, 0.210400f, 0.211980f, 0.213630f, 0.215330f, 0.217090f, 0.218880f,
        0.220690f, 0.222520f, 0.224350f, 0.226170f, 0.227970f, 0.229730f, 0.231450f, 0.233090f, 0.234650f, 0.236120f,
        0.237470f, 0.238680f, 0.239750f, 0.240660f, 0.241390f, 0.241960f, 0.242370f, 0.242650f, 0.242810f, 0.242860f,
        0.242810f, 0.242680f, 0.242490f, 0.242240f, 0.241960f, 0.241630f, 0.241270f, 0.240870f, 0.240440f, 0.239980f,
        0.239470f, 0.238940f, 0.238380f, 0.237780f, 0.237150f, 0.236490f, 0.235800f, 0.235080f, 0.234320f, 0.233520f,
        0.232680f, 0.231810f, 0.230890f, 0.229940f, 0.228940f, 0.227900f, 0.226810f, 0.225680f, 0.224500f, 0.223270f,
        0.221990f, 0.220650f, 0.219270f, 0.217830f, 0.217520f, 0.216430f, 0.215350f, 0.214270f, 0.213200f, 0.212130f,
        0.211070f, 0.210000f, 0.208950f, 0.207890f, 0.206840f, 0.205800f, 0.204750f, 0.203720f, 0.202680f, 0.201650f,
        0.200630f, 0.199600f, 0.198580f, 0.197570f, 0.196560f, 0.195550f, 0.194550f, 0.193550f, 0.192560f, 0.191560f,
        0.190580f, 0.189590f, 0.188610f, 0.187640f, 0.186670f, 0.185700f, 0.184740f, 0.183780f, 0.182820f, 0.181870f,
        0.180920f, 0.179970f, 0.179030f, 0.178100f, 0.177160f, 0.176230f, 0.175310f, 0.174390f, 0.173470f, 0.172550f,
        0.171640f, 0.170740f, 0.169840f, 0.168940f, 0.168040f, 0.167150f, 0.166260f, 0.165380f, 0.164500f, 0.163620f,
        0.162750f, 0.161880f, 0.161020f, 0.160160f, 0.159300f, 0.158450f, 0.157600f, 0.156750f, 0.155910f, 0.155070f,
        0.154230f, 0.153400f, 0.152580f, 0.151750f, 0.150930f, 0.150110f, 0.149300f, 0.148490f, 0.147690f, 0.146880f,
        0.146080f, 0.145290f, 0.144500f});
    static constexpr Spectrum CES52({0.058595f, 0.058559f, 0.058524f, 0.058488f, 0.058453f, 0.058417f, 0.058382f,
        0.058346f, 0.058311f, 0.058275f, 0.058240f, 0.058205f, 0.058169f, 0.058134f, 0.058098f, 0.058063f, 0.058028f,
        0.057993f, 0.057957f, 0.057922f, 0.058015f, 0.057923f, 0.057840f, 0.057767f, 0.057703f, 0.057649f, 0.057603f,
        0.057565f, 0.057536f, 0.057514f, 0.057500f, 0.057493f, 0.057493f, 0.057500f, 0.057513f, 0.057532f, 0.057556f,
        0.057587f, 0.057622f, 0.057662f, 0.057706f, 0.057755f, 0.057809f, 0.057869f, 0.057936f, 0.058010f, 0.058094f,
        0.058186f, 0.058289f, 0.058404f, 0.058531f, 0.058670f, 0.058823f, 0.058991f, 0.059172f, 0.059369f, 0.059581f,
        0.059809f, 0.060053f, 0.060313f, 0.060592f, 0.060888f, 0.061203f, 0.061541f, 0.061904f, 0.062293f, 0.062710f,
        0.063159f, 0.063641f, 0.064158f, 0.064713f, 0.065306f, 0.065927f, 0.066563f, 0.067205f, 0.067839f, 0.068453f,
        0.069037f, 0.069578f, 0.070064f, 0.070484f, 0.070830f, 0.071114f, 0.071351f, 0.071559f, 0.071751f, 0.071946f,
        0.072158f, 0.072404f, 0.072699f, 0.073060f, 0.073506f, 0.074068f, 0.074780f, 0.075677f, 0.076793f, 0.078163f,
        0.079821f, 0.081802f, 0.084139f, 0.086868f, 0.090011f, 0.093537f, 0.097405f, 0.101570f, 0.106000f, 0.110650f,
        0.115470f, 0.120420f, 0.125470f, 0.130560f, 0.135670f, 0.140760f, 0.145810f, 0.150800f, 0.155710f, 0.160510f,
        0.165180f, 0.169690f, 0.174030f, 0.178170f, 0.182090f, 0.185790f, 0.189300f, 0.192610f, 0.195750f, 0.198710f,
        0.201510f, 0.204170f, 0.206690f, 0.209080f, 0.211350f, 0.213470f, 0.215390f, 0.217090f, 0.218530f, 0.219670f,
        0.220480f, 0.220910f, 0.220940f, 0.220520f, 0.219630f, 0.218300f, 0.216570f, 0.214470f, 0.212040f, 0.209320f,
        0.206350f, 0.203160f, 0.199800f, 0.196300f, 0.192710f, 0.189030f, 0.185290f, 0.181510f, 0.177700f, 0.173900f,
        0.170100f, 0.166340f, 0.162640f, 0.159000f, 0.155450f, 0.151990f, 0.148590f, 0.145270f, 0.142000f, 0.138780f,
        0.135600f, 0.132460f, 0.129340f, 0.126230f, 0.123140f, 0.120070f, 0.117040f, 0.114070f, 0.111150f, 0.108320f,
        0.105590f, 0.102960f, 0.100460f, 0.098101f, 0.095888f, 0.093821f, 0.091891f, 0.090089f, 0.088407f, 0.086837f,
        0.085371f, 0.084000f, 0.082716f, 0.081510f, 0.080376f, 0.079310f, 0.078309f, 0.077372f, 0.076497f, 0.075680f,
        0.074920f, 0.074214f, 0.073561f, 0.072957f, 0.072401f, 0.071892f, 0.071428f, 0.071008f, 0.070631f, 0.070295f,
        0.070001f, 0.069746f, 0.069530f, 0.069350f, 0.069207f, 0.069095f, 0.069011f, 0.068950f, 0.068909f, 0.068882f,
        0.068865f, 0.068855f, 0.068846f, 0.068835f, 0.068819f, 0.068797f, 0.068772f, 0.068745f, 0.068717f, 0.068691f,
        0.068667f, 0.068648f, 0.068635f, 0.068629f, 0.068632f, 0.068643f, 0.068662f, 0.068687f, 0.068719f, 0.068755f,
        0.068796f, 0.068841f, 0.068889f, 0.068938f, 0.068990f, 0.069043f, 0.069100f, 0.069160f, 0.069226f, 0.069297f,
        0.069376f, 0.069461f, 0.069556f, 0.069660f, 0.069774f, 0.069903f, 0.070048f, 0.070214f, 0.070403f, 0.070620f,
        0.070867f, 0.071147f, 0.071465f, 0.071824f, 0.072224f, 0.072665f, 0.073140f, 0.073646f, 0.074178f, 0.074732f,
        0.075303f, 0.075888f, 0.076481f, 0.077079f, 0.077677f, 0.078269f, 0.078852f, 0.079419f, 0.079967f, 0.080490f,
        0.080983f, 0.081440f, 0.081858f, 0.082231f, 0.082556f, 0.082834f, 0.083068f, 0.083263f, 0.083421f, 0.083546f,
        0.083640f, 0.083708f, 0.083753f, 0.083777f, 0.083784f, 0.083774f, 0.083747f, 0.083702f, 0.083640f, 0.083560f,
        0.083461f, 0.083345f, 0.083210f, 0.083056f, 0.082883f, 0.082691f, 0.082480f, 0.082251f, 0.082002f, 0.081735f,
        0.081449f, 0.081145f, 0.080821f, 0.080480f, 0.080119f, 0.079740f, 0.079343f, 0.078927f, 0.078493f, 0.078041f,
        0.077570f, 0.077081f, 0.076574f, 0.076049f, 0.075958f, 0.075567f, 0.075178f, 0.074791f, 0.074405f, 0.074022f,
        0.073640f, 0.073260f, 0.072882f, 0.072505f, 0.072131f, 0.071758f, 0.071387f, 0.071018f, 0.070650f, 0.070285f,
        0.069921f, 0.069558f, 0.069198f, 0.068839f, 0.068482f, 0.068127f, 0.067773f, 0.067421f, 0.067071f, 0.066722f,
        0.066376f, 0.066030f, 0.065687f, 0.065345f, 0.065005f, 0.064666f, 0.064330f, 0.063994f, 0.063661f, 0.063329f,
        0.062998f, 0.062670f, 0.062342f, 0.062017f, 0.061693f, 0.061370f, 0.061049f, 0.060730f, 0.060413f, 0.060096f,
        0.059782f, 0.059469f, 0.059157f, 0.058847f, 0.058539f, 0.058232f, 0.057926f, 0.057622f, 0.057320f, 0.057019f,
        0.056720f, 0.056422f, 0.056125f, 0.055830f, 0.055536f, 0.055244f, 0.054954f, 0.054664f, 0.054377f, 0.054090f,
        0.053805f, 0.053522f, 0.053240f, 0.052959f, 0.052679f, 0.052402f, 0.052125f, 0.051850f, 0.051576f, 0.051304f,
        0.051032f, 0.050763f, 0.050494f});
    static constexpr Spectrum CES53({0.072779f, 0.073026f, 0.073273f, 0.073521f, 0.073770f, 0.074019f, 0.074270f,
        0.074521f, 0.074773f, 0.075025f, 0.075279f, 0.075533f, 0.075788f, 0.076044f, 0.076300f, 0.076558f, 0.076816f,
        0.077075f, 0.077335f, 0.077595f, 0.077900f, 0.078095f, 0.078340f, 0.078628f, 0.078951f, 0.079302f, 0.079672f,
        0.080055f, 0.080442f, 0.080826f, 0.081200f, 0.081558f, 0.081903f, 0.082239f, 0.082573f, 0.082909f, 0.083252f,
        0.083608f, 0.083981f, 0.084377f, 0.084800f, 0.085255f, 0.085740f, 0.086255f, 0.086796f, 0.087362f, 0.087951f,
        0.088562f, 0.089191f, 0.089838f, 0.090500f, 0.091176f, 0.091866f, 0.092572f, 0.093293f, 0.094031f, 0.094786f,
        0.095560f, 0.096353f, 0.097166f, 0.098000f, 0.098857f, 0.099741f, 0.100660f, 0.101620f, 0.102630f, 0.103690f,
        0.104810f, 0.106000f, 0.107260f, 0.108600f, 0.110030f, 0.111540f, 0.113130f, 0.114820f, 0.116580f, 0.118440f,
        0.120370f, 0.122400f, 0.124510f, 0.126700f, 0.128980f, 0.131350f, 0.133800f, 0.136350f, 0.139000f, 0.141750f,
        0.144600f, 0.147560f, 0.150620f, 0.153800f, 0.157090f, 0.160500f, 0.164040f, 0.167710f, 0.171520f, 0.175470f,
        0.179570f, 0.183820f, 0.188230f, 0.192800f, 0.197540f, 0.202460f, 0.207540f, 0.212800f, 0.218220f, 0.223820f,
        0.229580f, 0.235520f, 0.241620f, 0.247900f, 0.254340f, 0.260910f, 0.267580f, 0.274300f, 0.281050f, 0.287770f,
        0.294430f, 0.301000f, 0.307440f, 0.313700f, 0.319750f, 0.325560f, 0.331080f, 0.336280f, 0.341110f, 0.345540f,
        0.349530f, 0.353050f, 0.356050f, 0.358500f, 0.360370f, 0.361670f, 0.362440f, 0.362710f, 0.362510f, 0.361860f,
        0.360790f, 0.359340f, 0.357530f, 0.355400f, 0.352970f, 0.350250f, 0.347270f, 0.344040f, 0.340580f, 0.336890f,
        0.333010f, 0.328940f, 0.324700f, 0.320300f, 0.315760f, 0.311100f, 0.306320f, 0.301440f, 0.296460f, 0.291400f,
        0.286280f, 0.281090f, 0.275860f, 0.270600f, 0.265310f, 0.260010f, 0.254710f, 0.249410f, 0.244120f, 0.238850f,
        0.233620f, 0.228430f, 0.223280f, 0.218200f, 0.213180f, 0.208240f, 0.203370f, 0.198580f, 0.193870f, 0.189250f,
        0.184710f, 0.180280f, 0.175940f, 0.171700f, 0.167570f, 0.163540f, 0.159600f, 0.155770f, 0.152040f, 0.148390f,
        0.144840f, 0.141370f, 0.138000f, 0.134700f, 0.131490f, 0.128350f, 0.125290f, 0.122310f, 0.119400f, 0.116570f,
        0.113800f, 0.111100f, 0.108470f, 0.105900f, 0.103390f, 0.100950f, 0.098578f, 0.096269f, 0.094028f, 0.091858f,
        0.089758f, 0.087731f, 0.085778f, 0.083900f, 0.082099f, 0.080376f, 0.078729f, 0.077161f, 0.075671f, 0.074259f,
        0.072926f, 0.071671f, 0.070496f, 0.069400f, 0.068383f, 0.067442f, 0.066574f, 0.065774f, 0.065038f, 0.064364f,
        0.063746f, 0.063182f, 0.062668f, 0.062200f, 0.061774f, 0.061388f, 0.061038f, 0.060723f, 0.060439f, 0.060184f,
        0.059955f, 0.059750f, 0.059566f, 0.059400f, 0.059250f, 0.059114f, 0.058992f, 0.058881f, 0.058781f, 0.058691f,
        0.058609f, 0.058534f, 0.058465f, 0.058400f, 0.058339f, 0.058281f, 0.058226f, 0.058173f, 0.058123f, 0.058075f,
        0.058029f, 0.057985f, 0.057942f, 0.057900f, 0.057860f, 0.057822f, 0.057791f, 0.057766f, 0.057752f, 0.057749f,
        0.057760f, 0.057787f, 0.057833f, 0.057900f, 0.057989f, 0.058099f, 0.058230f, 0.058379f, 0.058546f, 0.058729f,
        0.058928f, 0.059140f, 0.059364f, 0.059600f, 0.059846f, 0.060101f, 0.060367f, 0.060642f, 0.060927f, 0.061222f,
        0.061527f, 0.061841f, 0.062166f, 0.062500f, 0.062844f, 0.063195f, 0.063554f, 0.063917f, 0.064283f, 0.064650f,
        0.065017f, 0.065382f, 0.065744f, 0.066100f, 0.066450f, 0.066793f, 0.067133f, 0.067470f, 0.067804f, 0.068139f,
        0.068474f, 0.068812f, 0.069154f, 0.069500f, 0.069852f, 0.070204f, 0.070552f, 0.070890f, 0.071212f, 0.071514f,
        0.071789f, 0.072032f, 0.072237f, 0.072400f, 0.072663f, 0.072887f, 0.073112f, 0.073338f, 0.073564f, 0.073791f,
        0.074018f, 0.074246f, 0.074475f, 0.074705f, 0.074935f, 0.075165f, 0.075397f, 0.075629f, 0.075862f, 0.076095f,
        0.076329f, 0.076564f, 0.076799f, 0.077035f, 0.077272f, 0.077509f, 0.077747f, 0.077986f, 0.078225f, 0.078465f,
        0.078706f, 0.078947f, 0.079189f, 0.079432f, 0.079675f, 0.079919f, 0.080164f, 0.080409f, 0.080656f, 0.080902f,
        0.081150f, 0.081398f, 0.081647f, 0.081897f, 0.082147f, 0.082398f, 0.082649f, 0.082902f, 0.083155f, 0.083409f,
        0.083663f, 0.083918f, 0.084174f, 0.084431f, 0.084688f, 0.084946f, 0.085205f, 0.085464f, 0.085724f, 0.085985f,
        0.086247f, 0.086509f, 0.086772f, 0.087036f, 0.087300f, 0.087566f, 0.087832f, 0.088098f, 0.088366f, 0.088634f,
        0.088903f, 0.089172f, 0.089443f, 0.089714f, 0.089986f, 0.090258f, 0.090531f, 0.090805f, 0.091080f, 0.091356f,
        0.091632f, 0.091909f, 0.092187f});
    static constexpr Spectrum CES54({0.408570f, 0.411660f, 0.414750f, 0.417860f, 0.420960f, 0.424080f, 0.427200f,
        0.430320f, 0.433450f, 0.436590f, 0.439730f, 0.442880f, 0.446030f, 0.449180f, 0.452340f, 0.455510f, 0.458670f,
        0.461840f, 0.465020f, 0.468190f, 0.471900f, 0.474260f, 0.477220f, 0.480700f, 0.484600f, 0.488830f, 0.493290f,
        0.497890f, 0.502540f, 0.507140f, 0.511600f, 0.515850f, 0.519890f, 0.523770f, 0.527510f, 0.531140f, 0.534690f,
        0.538200f, 0.541700f, 0.545230f, 0.548800f, 0.552450f, 0.556170f, 0.559930f, 0.563730f, 0.567540f, 0.571350f,
        0.575130f, 0.578880f, 0.582580f, 0.586200f, 0.589740f, 0.593200f, 0.596610f, 0.599980f, 0.603330f, 0.606680f,
        0.610040f, 0.613440f, 0.616890f, 0.620400f, 0.623990f, 0.627670f, 0.631410f, 0.635210f, 0.639070f, 0.642970f,
        0.646910f, 0.650890f, 0.654890f, 0.658900f, 0.662920f, 0.666940f, 0.670960f, 0.674950f, 0.678910f, 0.682840f,
        0.686720f, 0.690550f, 0.694310f, 0.698000f, 0.701610f, 0.705140f, 0.708580f, 0.711950f, 0.715230f, 0.718440f,
        0.721570f, 0.724620f, 0.727600f, 0.730500f, 0.733320f, 0.736070f, 0.738730f, 0.741310f, 0.743790f, 0.746180f,
        0.748470f, 0.750660f, 0.752730f, 0.754700f, 0.756550f, 0.758290f, 0.759920f, 0.761440f, 0.762850f, 0.764160f,
        0.765360f, 0.766470f, 0.767480f, 0.768400f, 0.769220f, 0.769960f, 0.770610f, 0.771180f, 0.771660f, 0.772070f,
        0.772410f, 0.772670f, 0.772870f, 0.773000f, 0.773070f, 0.773090f, 0.773060f, 0.772990f, 0.772890f, 0.772770f,
        0.772630f, 0.772480f, 0.772340f, 0.772200f, 0.772080f, 0.771960f, 0.771850f, 0.771750f, 0.771640f, 0.771520f,
        0.771400f, 0.771250f, 0.771090f, 0.770900f, 0.770680f, 0.770430f, 0.770150f, 0.769840f, 0.769500f, 0.769130f,
        0.768720f, 0.768280f, 0.767810f, 0.767300f, 0.766760f, 0.766180f, 0.765560f, 0.764910f, 0.764210f, 0.763460f,
        0.762670f, 0.761830f, 0.760940f, 0.760000f, 0.759000f, 0.757950f, 0.756830f, 0.755640f, 0.754370f, 0.753030f,
        0.751610f, 0.750100f, 0.748500f, 0.746800f, 0.745000f, 0.743110f, 0.741120f, 0.739060f, 0.736910f, 0.734690f,
        0.732410f, 0.730060f, 0.727650f, 0.725200f, 0.722700f, 0.720160f, 0.717570f, 0.714930f, 0.712260f, 0.709530f,
        0.706770f, 0.703960f, 0.701100f, 0.698200f, 0.695250f, 0.692250f, 0.689180f, 0.686030f, 0.682800f, 0.679460f,
        0.676020f, 0.672450f, 0.668750f, 0.664900f, 0.660900f, 0.656750f, 0.652440f, 0.647970f, 0.643360f, 0.638590f,
        0.633670f, 0.628600f, 0.623370f, 0.618000f, 0.612480f, 0.606850f, 0.601130f, 0.595370f, 0.589590f, 0.583830f,
        0.578130f, 0.572520f, 0.567030f, 0.561700f, 0.556560f, 0.551610f, 0.546870f, 0.542330f, 0.538020f, 0.533920f,
        0.530060f, 0.526430f, 0.523040f, 0.519900f, 0.517010f, 0.514360f, 0.511950f, 0.509770f, 0.507800f, 0.506040f,
        0.504480f, 0.503110f, 0.501920f, 0.500900f, 0.500040f, 0.499330f, 0.498760f, 0.498310f, 0.497980f, 0.497750f,
        0.497610f, 0.497540f, 0.497550f, 0.497600f, 0.497700f, 0.497830f, 0.497990f, 0.498180f, 0.498390f, 0.498610f,
        0.498830f, 0.499060f, 0.499290f, 0.499500f, 0.499700f, 0.499900f, 0.500120f, 0.500360f, 0.500650f, 0.500990f,
        0.501400f, 0.501890f, 0.502490f, 0.503200f, 0.504040f, 0.505000f, 0.506110f, 0.507350f, 0.508740f, 0.510280f,
        0.511970f, 0.513820f, 0.515830f, 0.518000f, 0.520340f, 0.522830f, 0.525460f, 0.528200f, 0.531040f, 0.533960f,
        0.536950f, 0.539980f, 0.543030f, 0.546100f, 0.549160f, 0.552210f, 0.555240f, 0.558240f, 0.561220f, 0.564160f,
        0.567070f, 0.569930f, 0.572740f, 0.575500f, 0.578200f, 0.580850f, 0.583470f, 0.586050f, 0.588620f, 0.591170f,
        0.593730f, 0.596300f, 0.598880f, 0.601500f, 0.604150f, 0.606800f, 0.609410f, 0.611940f, 0.614350f, 0.616600f,
        0.618650f, 0.620460f, 0.621990f, 0.623200f, 0.625130f, 0.626770f, 0.628420f, 0.630050f, 0.631690f, 0.633320f,
        0.634950f, 0.636570f, 0.638200f, 0.639820f, 0.641430f, 0.643050f, 0.644660f, 0.646260f, 0.647870f, 0.649470f,
        0.651060f, 0.652660f, 0.654250f, 0.655830f, 0.657410f, 0.658990f, 0.660570f, 0.662140f, 0.663710f, 0.665280f,
        0.666840f, 0.668400f, 0.669950f, 0.671500f, 0.673050f, 0.674590f, 0.676130f, 0.677660f, 0.679200f, 0.680720f,
        0.682250f, 0.683770f, 0.685280f, 0.686800f, 0.688300f, 0.689810f, 0.691310f, 0.692800f, 0.694300f, 0.695780f,
        0.697270f, 0.698750f, 0.700220f, 0.701700f, 0.703160f, 0.704630f, 0.706090f, 0.707540f, 0.708990f, 0.710440f,
        0.711880f, 0.713320f, 0.714750f, 0.716180f, 0.717600f, 0.719020f, 0.720440f, 0.721850f, 0.723260f, 0.724660f,
        0.726060f, 0.727460f, 0.728850f, 0.730230f, 0.731610f, 0.732990f, 0.734360f, 0.735730f, 0.737090f, 0.738450f,
        0.739800f, 0.741150f, 0.742490f});
    static constexpr Spectrum CES55({0.300730f, 0.306210f, 0.311520f, 0.316680f, 0.321710f, 0.326610f, 0.331420f,
        0.336140f, 0.340780f, 0.345370f, 0.349920f, 0.354440f, 0.358940f, 0.363420f, 0.367870f, 0.372300f, 0.376710f,
        0.381100f, 0.385480f, 0.389830f, 0.394160f, 0.398470f, 0.402770f, 0.407050f, 0.411300f, 0.415530f, 0.419740f,
        0.423930f, 0.428080f, 0.432220f, 0.436320f, 0.440400f, 0.444450f, 0.448480f, 0.452510f, 0.456520f, 0.460530f,
        0.464550f, 0.468580f, 0.472620f, 0.476680f, 0.480770f, 0.484870f, 0.488990f, 0.493120f, 0.497250f, 0.501380f,
        0.505500f, 0.509600f, 0.513690f, 0.517750f, 0.521780f, 0.525790f, 0.529780f, 0.533760f, 0.537740f, 0.541730f,
        0.545730f, 0.549750f, 0.553800f, 0.557890f, 0.562020f, 0.566190f, 0.570410f, 0.574670f, 0.578980f, 0.583340f,
        0.587740f, 0.592200f, 0.596700f, 0.601250f, 0.605850f, 0.610490f, 0.615170f, 0.619860f, 0.624580f, 0.629290f,
        0.634010f, 0.638710f, 0.643390f, 0.648040f, 0.652650f, 0.657210f, 0.661710f, 0.666140f, 0.670490f, 0.674750f,
        0.678910f, 0.682950f, 0.686870f, 0.690660f, 0.694310f, 0.697810f, 0.701180f, 0.704420f, 0.707540f, 0.710520f,
        0.713390f, 0.716150f, 0.718790f, 0.721320f, 0.723750f, 0.726080f, 0.728300f, 0.730430f, 0.732460f, 0.734380f,
        0.736220f, 0.737960f, 0.739600f, 0.741150f, 0.742610f, 0.743980f, 0.745250f, 0.746430f, 0.747520f, 0.748510f,
        0.749400f, 0.750200f, 0.750900f, 0.751510f, 0.752020f, 0.752430f, 0.752750f, 0.752980f, 0.753120f, 0.753180f,
        0.753150f, 0.753040f, 0.752850f, 0.752590f, 0.752250f, 0.751840f, 0.751350f, 0.750770f, 0.750100f, 0.749350f,
        0.748500f, 0.747560f, 0.746520f, 0.745370f, 0.744120f, 0.742760f, 0.741300f, 0.739750f, 0.738100f, 0.736360f,
        0.734530f, 0.732610f, 0.730610f, 0.728530f, 0.726370f, 0.724130f, 0.721810f, 0.719390f, 0.716890f, 0.714280f,
        0.711580f, 0.708780f, 0.705860f, 0.702840f, 0.699700f, 0.696450f, 0.693100f, 0.689640f, 0.686080f, 0.682430f,
        0.678680f, 0.674850f, 0.670940f, 0.666940f, 0.662870f, 0.658740f, 0.654550f, 0.650330f, 0.646070f, 0.641800f,
        0.637530f, 0.633260f, 0.629020f, 0.624800f, 0.620630f, 0.616500f, 0.612420f, 0.608400f, 0.604440f, 0.600530f,
        0.596700f, 0.592930f, 0.589240f, 0.585620f, 0.582080f, 0.578620f, 0.575210f, 0.571860f, 0.568560f, 0.565300f,
        0.562060f, 0.558840f, 0.555630f, 0.552420f, 0.549200f, 0.545970f, 0.542720f, 0.539440f, 0.536110f, 0.532750f,
        0.529320f, 0.525840f, 0.522280f, 0.518650f, 0.514940f, 0.511170f, 0.507360f, 0.503550f, 0.499750f, 0.496000f,
        0.492320f, 0.488740f, 0.485280f, 0.481960f, 0.478810f, 0.475840f, 0.473030f, 0.470390f, 0.467910f, 0.465590f,
        0.463430f, 0.461430f, 0.459580f, 0.457880f, 0.456330f, 0.454920f, 0.453640f, 0.452490f, 0.451460f, 0.450550f,
        0.449750f, 0.449050f, 0.448440f, 0.447920f, 0.447480f, 0.447120f, 0.446830f, 0.446610f, 0.446440f, 0.446330f,
        0.446270f, 0.446250f, 0.446260f, 0.446310f, 0.446380f, 0.446490f, 0.446620f, 0.446780f, 0.446990f, 0.447230f,
        0.447510f, 0.447830f, 0.448200f, 0.448620f, 0.449090f, 0.449620f, 0.450210f, 0.450870f, 0.451600f, 0.452410f,
        0.453300f, 0.454290f, 0.455370f, 0.456560f, 0.457850f, 0.459240f, 0.460730f, 0.462300f, 0.463960f, 0.465680f,
        0.467480f, 0.469330f, 0.471240f, 0.473190f, 0.475180f, 0.477210f, 0.479260f, 0.481330f, 0.483410f, 0.485500f,
        0.487580f, 0.489660f, 0.491720f, 0.493760f, 0.495770f, 0.497750f, 0.499700f, 0.501610f, 0.503470f, 0.505300f,
        0.507080f, 0.508810f, 0.510490f, 0.512120f, 0.513690f, 0.515200f, 0.516640f, 0.518010f, 0.519310f, 0.520530f,
        0.521660f, 0.522700f, 0.523650f, 0.524500f, 0.525250f, 0.525900f, 0.526460f, 0.526940f, 0.527350f, 0.527690f,
        0.527970f, 0.528200f, 0.528390f, 0.528540f, 0.528660f, 0.528760f, 0.528840f, 0.528920f, 0.528990f, 0.529060f,
        0.529150f, 0.529250f, 0.529380f, 0.529540f, 0.529740f, 0.529980f, 0.530280f, 0.530650f, 0.531090f, 0.531620f,
        0.532250f, 0.532970f, 0.533820f, 0.534780f, 0.535880f, 0.537100f, 0.538470f, 0.539960f, 0.541580f, 0.543340f,
        0.545230f, 0.547250f, 0.549400f, 0.551680f, 0.554100f, 0.556660f, 0.559400f, 0.562330f, 0.565480f, 0.568860f,
        0.572490f, 0.576400f, 0.580610f, 0.585140f, 0.589980f, 0.595050f, 0.600230f, 0.605400f, 0.610450f, 0.615270f,
        0.619740f, 0.623740f, 0.627170f, 0.629910f, 0.631880f, 0.633140f, 0.633790f, 0.633940f, 0.633690f, 0.633140f,
        0.632380f, 0.631520f, 0.630670f, 0.629910f, 0.629340f, 0.628940f, 0.628710f, 0.628620f, 0.628650f, 0.628780f,
        0.628990f, 0.629260f, 0.629580f, 0.629910f, 0.630240f, 0.630560f, 0.630830f, 0.631040f, 0.631170f, 0.631200f,
        0.631110f, 0.630880f, 0.630480f});
    static constexpr Spectrum CES56({0.155210f, 0.158450f, 0.161750f, 0.165110f, 0.168520f, 0.171990f, 0.175510f,
        0.179090f, 0.182730f, 0.186430f, 0.190180f, 0.193990f, 0.197850f, 0.201780f, 0.205760f, 0.209800f, 0.213900f,
        0.218060f, 0.222270f, 0.226550f, 0.231500f, 0.234890f, 0.239110f, 0.243980f, 0.249350f, 0.255030f, 0.260840f,
        0.266620f, 0.272190f, 0.277380f, 0.282010f, 0.285950f, 0.289250f, 0.291990f, 0.294260f, 0.296160f, 0.297760f,
        0.299150f, 0.300420f, 0.301670f, 0.302970f, 0.304400f, 0.305960f, 0.307640f, 0.309420f, 0.311300f, 0.313250f,
        0.315270f, 0.317340f, 0.319450f, 0.321580f, 0.323730f, 0.325890f, 0.328080f, 0.330290f, 0.332530f, 0.334810f,
        0.337130f, 0.339490f, 0.341900f, 0.344370f, 0.346900f, 0.349510f, 0.352230f, 0.355080f, 0.358080f, 0.361250f,
        0.364610f, 0.368200f, 0.372030f, 0.376120f, 0.380490f, 0.385120f, 0.389990f, 0.395070f, 0.400320f, 0.405730f,
        0.411260f, 0.416900f, 0.422610f, 0.428360f, 0.434140f, 0.439980f, 0.445910f, 0.451950f, 0.458150f, 0.464530f,
        0.471140f, 0.478000f, 0.485150f, 0.492620f, 0.500430f, 0.508530f, 0.516870f, 0.525380f, 0.533990f, 0.542660f,
        0.551310f, 0.559890f, 0.568330f, 0.576580f, 0.584580f, 0.592300f, 0.599740f, 0.606870f, 0.613680f, 0.620160f,
        0.626290f, 0.632060f, 0.637450f, 0.642440f, 0.647030f, 0.651210f, 0.654990f, 0.658380f, 0.661370f, 0.663980f,
        0.666200f, 0.668050f, 0.669530f, 0.670630f, 0.671370f, 0.671770f, 0.671840f, 0.671590f, 0.671060f, 0.670240f,
        0.669180f, 0.667870f, 0.666340f, 0.664610f, 0.662690f, 0.660600f, 0.658330f, 0.655900f, 0.653320f, 0.650590f,
        0.647720f, 0.644720f, 0.641600f, 0.638360f, 0.635010f, 0.631550f, 0.627980f, 0.624300f, 0.620520f, 0.616620f,
        0.612610f, 0.608500f, 0.604270f, 0.599940f, 0.595500f, 0.590940f, 0.586280f, 0.581510f, 0.576640f, 0.571660f,
        0.566580f, 0.561390f, 0.556100f, 0.550710f, 0.545220f, 0.539640f, 0.533970f, 0.528230f, 0.522410f, 0.516540f,
        0.510610f, 0.504640f, 0.498630f, 0.492590f, 0.486530f, 0.480450f, 0.474360f, 0.468270f, 0.462180f, 0.456090f,
        0.450020f, 0.443960f, 0.437920f, 0.431910f, 0.425930f, 0.419980f, 0.414070f, 0.408170f, 0.402310f, 0.396460f,
        0.390640f, 0.384840f, 0.379050f, 0.373280f, 0.367520f, 0.361780f, 0.356050f, 0.350340f, 0.344650f, 0.338980f,
        0.333330f, 0.327700f, 0.322090f, 0.316510f, 0.310950f, 0.305420f, 0.299930f, 0.294470f, 0.289040f, 0.283660f,
        0.278330f, 0.273040f, 0.267810f, 0.262630f, 0.257510f, 0.252460f, 0.247500f, 0.242630f, 0.237870f, 0.233220f,
        0.228710f, 0.224340f, 0.220130f, 0.216090f, 0.212230f, 0.208540f, 0.205030f, 0.201710f, 0.198560f, 0.195590f,
        0.192790f, 0.190180f, 0.187740f, 0.185470f, 0.183380f, 0.181450f, 0.179670f, 0.178040f, 0.176530f, 0.175140f,
        0.173860f, 0.172670f, 0.171570f, 0.170540f, 0.169570f, 0.168660f, 0.167800f, 0.166990f, 0.166210f, 0.165480f,
        0.164780f, 0.164100f, 0.163450f, 0.162810f, 0.162190f, 0.161580f, 0.160980f, 0.160390f, 0.159810f, 0.159240f,
        0.158680f, 0.158130f, 0.157580f, 0.157040f, 0.156510f, 0.155990f, 0.155480f, 0.155000f, 0.154550f, 0.154140f,
        0.153760f, 0.153440f, 0.153160f, 0.152950f, 0.152800f, 0.152720f, 0.152700f, 0.152750f, 0.152860f, 0.153040f,
        0.153290f, 0.153600f, 0.153980f, 0.154430f, 0.154950f, 0.155530f, 0.156190f, 0.156920f, 0.157720f, 0.158610f,
        0.159570f, 0.160610f, 0.161740f, 0.162950f, 0.164250f, 0.165610f, 0.167030f, 0.168480f, 0.169940f, 0.171410f,
        0.172850f, 0.174250f, 0.175600f, 0.176870f, 0.178050f, 0.179160f, 0.180210f, 0.181200f, 0.182160f, 0.183110f,
        0.184050f, 0.185000f, 0.185970f, 0.186990f, 0.188060f, 0.189160f, 0.190270f, 0.191360f, 0.192430f, 0.193440f,
        0.194370f, 0.195200f, 0.195910f, 0.196470f, 0.197370f, 0.198130f, 0.198910f, 0.199680f, 0.200450f, 0.201230f,
        0.202010f, 0.202790f, 0.203570f, 0.204360f, 0.205150f, 0.205940f, 0.206730f, 0.207530f, 0.208330f, 0.209120f,
        0.209930f, 0.210730f, 0.211540f, 0.212350f, 0.213160f, 0.213970f, 0.214790f, 0.215600f, 0.216420f, 0.217250f,
        0.218070f, 0.218900f, 0.219730f, 0.220560f, 0.221390f, 0.222230f, 0.223060f, 0.223900f, 0.224750f, 0.225590f,
        0.226440f, 0.227290f, 0.228140f, 0.228990f, 0.229850f, 0.230710f, 0.231570f, 0.232430f, 0.233300f, 0.234160f,
        0.235030f, 0.235900f, 0.236780f, 0.237650f, 0.238530f, 0.239410f, 0.240290f, 0.241180f, 0.242070f, 0.242960f,
        0.243850f, 0.244740f, 0.245640f, 0.246540f, 0.247440f, 0.248340f, 0.249240f, 0.250150f, 0.251060f, 0.251970f,
        0.252890f, 0.253800f, 0.254720f, 0.255640f, 0.256560f, 0.257490f, 0.258410f, 0.259340f, 0.260270f, 0.261210f,
        0.262140f, 0.263080f, 0.264020f});
    static constexpr Spectrum CES57({0.134700f, 0.134750f, 0.134820f, 0.134880f, 0.134930f, 0.134970f, 0.134990f,
        0.134980f, 0.134930f, 0.134840f, 0.134700f, 0.134510f, 0.134280f, 0.134050f, 0.133850f, 0.133700f, 0.133630f,
        0.133670f, 0.133840f, 0.134170f, 0.134700f, 0.135430f, 0.136350f, 0.137410f, 0.138590f, 0.139840f, 0.141130f,
        0.142430f, 0.143700f, 0.144900f, 0.146000f, 0.146970f, 0.147840f, 0.148610f, 0.149320f, 0.149990f, 0.150650f,
        0.151310f, 0.152010f, 0.152760f, 0.153600f, 0.154540f, 0.155580f, 0.156730f, 0.157990f, 0.159340f, 0.160810f,
        0.162370f, 0.164050f, 0.165820f, 0.167700f, 0.169680f, 0.171760f, 0.173940f, 0.176210f, 0.178570f, 0.181020f,
        0.183550f, 0.186160f, 0.188840f, 0.191600f, 0.194430f, 0.197330f, 0.200320f, 0.203380f, 0.206540f, 0.209790f,
        0.213130f, 0.216580f, 0.220130f, 0.223800f, 0.227580f, 0.231490f, 0.235510f, 0.239660f, 0.243940f, 0.248360f,
        0.252910f, 0.257590f, 0.262420f, 0.267400f, 0.272520f, 0.277780f, 0.283170f, 0.288660f, 0.294270f, 0.299960f,
        0.305730f, 0.311570f, 0.317460f, 0.323400f, 0.329370f, 0.335360f, 0.341370f, 0.347370f, 0.353360f, 0.359320f,
        0.365250f, 0.371130f, 0.376950f, 0.382700f, 0.388370f, 0.393930f, 0.399380f, 0.404680f, 0.409830f, 0.414800f,
        0.419570f, 0.424120f, 0.428440f, 0.432500f, 0.436290f, 0.439800f, 0.443040f, 0.446010f, 0.448700f, 0.451110f,
        0.453250f, 0.455110f, 0.456690f, 0.458000f, 0.459030f, 0.459800f, 0.460300f, 0.460550f, 0.460550f, 0.460310f,
        0.459840f, 0.459140f, 0.458230f, 0.457100f, 0.455770f, 0.454250f, 0.452540f, 0.450660f, 0.448610f, 0.446400f,
        0.444050f, 0.441560f, 0.438940f, 0.436200f, 0.433350f, 0.430390f, 0.427330f, 0.424150f, 0.420880f, 0.417500f,
        0.414020f, 0.410450f, 0.406770f, 0.403000f, 0.399130f, 0.395170f, 0.391100f, 0.386920f, 0.382640f, 0.378240f,
        0.373720f, 0.369070f, 0.364300f, 0.359400f, 0.354370f, 0.349210f, 0.343940f, 0.338580f, 0.333130f, 0.327610f,
        0.322040f, 0.316420f, 0.310770f, 0.305100f, 0.299420f, 0.293760f, 0.288110f, 0.282500f, 0.276940f, 0.271430f,
        0.266000f, 0.260660f, 0.255420f, 0.250300f, 0.245300f, 0.240430f, 0.235670f, 0.231040f, 0.226510f, 0.222100f,
        0.217800f, 0.213600f, 0.209500f, 0.205500f, 0.201590f, 0.197770f, 0.194030f, 0.190350f, 0.186740f, 0.183190f,
        0.179690f, 0.176230f, 0.172800f, 0.169400f, 0.166020f, 0.162640f, 0.159270f, 0.155880f, 0.152470f, 0.149030f,
        0.145540f, 0.141990f, 0.138380f, 0.134700f, 0.130940f, 0.127120f, 0.123290f, 0.119480f, 0.115720f, 0.112060f,
        0.108520f, 0.105140f, 0.101950f, 0.099000f, 0.096307f, 0.093865f, 0.091661f, 0.089678f, 0.087903f, 0.086319f,
        0.084911f, 0.083666f, 0.082567f, 0.081600f, 0.080750f, 0.080004f, 0.079350f, 0.078775f, 0.078267f, 0.077814f,
        0.077403f, 0.077022f, 0.076659f, 0.076300f, 0.075937f, 0.075569f, 0.075202f, 0.074837f, 0.074478f, 0.074129f,
        0.073794f, 0.073475f, 0.073176f, 0.072900f, 0.072649f, 0.072416f, 0.072192f, 0.071967f, 0.071733f, 0.071480f,
        0.071199f, 0.070882f, 0.070518f, 0.070100f, 0.069618f, 0.069065f, 0.068435f, 0.067721f, 0.066916f, 0.066014f,
        0.065008f, 0.063892f, 0.062658f, 0.061300f, 0.059817f, 0.058226f, 0.056551f, 0.054814f, 0.053040f, 0.051250f,
        0.049468f, 0.047717f, 0.046020f, 0.044400f, 0.042877f, 0.041455f, 0.040136f, 0.038922f, 0.037813f, 0.036811f,
        0.035917f, 0.035133f, 0.034461f, 0.033900f, 0.033452f, 0.033109f, 0.032864f, 0.032708f, 0.032634f, 0.032633f,
        0.032697f, 0.032818f, 0.032989f, 0.033200f, 0.033455f, 0.033799f, 0.034288f, 0.034979f, 0.035927f, 0.037189f,
        0.038821f, 0.040880f, 0.043421f, 0.046500f, 0.050142f, 0.054238f, 0.058649f, 0.063236f, 0.067858f, 0.072375f,
        0.076648f, 0.080536f, 0.083900f, 0.086600f, 0.088536f, 0.089770f, 0.090404f, 0.090539f, 0.090279f, 0.089724f,
        0.088977f, 0.088139f, 0.087313f, 0.086600f, 0.086081f, 0.085751f, 0.085581f, 0.085544f, 0.085614f, 0.085763f,
        0.085963f, 0.086188f, 0.086409f, 0.086600f, 0.086739f, 0.086828f, 0.086873f, 0.086883f, 0.086864f, 0.086824f,
        0.086771f, 0.086710f, 0.086651f, 0.086600f, 0.086563f, 0.086539f, 0.086527f, 0.086524f, 0.086529f, 0.086540f,
        0.086554f, 0.086570f, 0.086586f, 0.086600f, 0.086610f, 0.086616f, 0.086620f, 0.086620f, 0.086619f, 0.086616f,
        0.086612f, 0.086608f, 0.086604f, 0.086600f, 0.086597f, 0.086596f, 0.086595f, 0.086594f, 0.086595f, 0.086596f,
        0.086597f, 0.086598f, 0.086599f, 0.086600f, 0.086601f, 0.086601f, 0.086602f, 0.086602f, 0.086602f, 0.086602f,
        0.086601f, 0.086601f, 0.086600f, 0.086600f, 0.086600f, 0.086599f, 0.086599f, 0.086598f, 0.086598f, 0.086598f,
        0.086598f, 0.086599f, 0.086599f});
    static constexpr Spectrum CES58({0.106060f, 0.106870f, 0.107700f, 0.108560f, 0.109440f, 0.110350f, 0.111280f,
        0.112220f, 0.113180f, 0.114150f, 0.115140f, 0.116140f, 0.117150f, 0.118160f, 0.119180f, 0.120200f, 0.121230f,
        0.122250f, 0.123280f, 0.124300f, 0.125310f, 0.126320f, 0.127320f, 0.128320f, 0.129320f, 0.130330f, 0.131340f,
        0.132360f, 0.133390f, 0.134440f, 0.135500f, 0.136580f, 0.137690f, 0.138810f, 0.139960f, 0.141140f, 0.142340f,
        0.143580f, 0.144840f, 0.146140f, 0.147470f, 0.148840f, 0.150240f, 0.151670f, 0.153140f, 0.154640f, 0.156170f,
        0.157740f, 0.159330f, 0.160940f, 0.162590f, 0.164260f, 0.165960f, 0.167700f, 0.169470f, 0.171280f, 0.173130f,
        0.175030f, 0.176980f, 0.178980f, 0.181030f, 0.183140f, 0.185330f, 0.187590f, 0.189950f, 0.192420f, 0.195000f,
        0.197710f, 0.200560f, 0.203560f, 0.206730f, 0.210070f, 0.213580f, 0.217270f, 0.221140f, 0.225180f, 0.229400f,
        0.233800f, 0.238370f, 0.243120f, 0.248040f, 0.253140f, 0.258390f, 0.263780f, 0.269280f, 0.274860f, 0.280520f,
        0.286220f, 0.291950f, 0.297680f, 0.303390f, 0.309060f, 0.314680f, 0.320250f, 0.325740f, 0.331160f, 0.336500f,
        0.341740f, 0.346880f, 0.351900f, 0.356810f, 0.361590f, 0.366210f, 0.370670f, 0.374940f, 0.379010f, 0.382850f,
        0.386440f, 0.389780f, 0.392820f, 0.395570f, 0.398000f, 0.400120f, 0.401930f, 0.403450f, 0.404680f, 0.405640f,
        0.406330f, 0.406770f, 0.406950f, 0.406900f, 0.406620f, 0.406120f, 0.405410f, 0.404500f, 0.403410f, 0.402150f,
        0.400730f, 0.399160f, 0.397450f, 0.395620f, 0.393670f, 0.391610f, 0.389450f, 0.387180f, 0.384830f, 0.382380f,
        0.379850f, 0.377230f, 0.374550f, 0.371790f, 0.368970f, 0.366080f, 0.363140f, 0.360130f, 0.357070f, 0.353950f,
        0.350770f, 0.347540f, 0.344270f, 0.340940f, 0.337570f, 0.334140f, 0.330670f, 0.327150f, 0.323570f, 0.319950f,
        0.316270f, 0.312530f, 0.308750f, 0.304900f, 0.301000f, 0.297040f, 0.293030f, 0.288980f, 0.284880f, 0.280730f,
        0.276550f, 0.272330f, 0.268080f, 0.263800f, 0.259490f, 0.255170f, 0.250830f, 0.246490f, 0.242150f, 0.237820f,
        0.233510f, 0.229230f, 0.224980f, 0.220770f, 0.216610f, 0.212490f, 0.208410f, 0.204350f, 0.200330f, 0.196330f,
        0.192340f, 0.188360f, 0.184390f, 0.180410f, 0.176430f, 0.172440f, 0.168460f, 0.164480f, 0.160510f, 0.156550f,
        0.152610f, 0.148680f, 0.144780f, 0.140910f, 0.137070f, 0.133260f, 0.129490f, 0.125770f, 0.122090f, 0.118480f,
        0.114910f, 0.111420f, 0.107990f, 0.104630f, 0.101350f, 0.098155f, 0.095050f, 0.092042f, 0.089137f, 0.086342f,
        0.083663f, 0.081107f, 0.078681f, 0.076390f, 0.074240f, 0.072227f, 0.070346f, 0.068593f, 0.066963f, 0.065451f,
        0.064051f, 0.062759f, 0.061571f, 0.060480f, 0.059483f, 0.058573f, 0.057747f, 0.056998f, 0.056322f, 0.055713f,
        0.055167f, 0.054678f, 0.054240f, 0.053850f, 0.053502f, 0.053192f, 0.052917f, 0.052674f, 0.052459f, 0.052270f,
        0.052103f, 0.051954f, 0.051821f, 0.051700f, 0.051588f, 0.051485f, 0.051390f, 0.051303f, 0.051224f, 0.051153f,
        0.051089f, 0.051032f, 0.050983f, 0.050940f, 0.050904f, 0.050878f, 0.050865f, 0.050868f, 0.050889f, 0.050933f,
        0.051001f, 0.051099f, 0.051227f, 0.051390f, 0.051590f, 0.051827f, 0.052100f, 0.052410f, 0.052755f, 0.053134f,
        0.053548f, 0.053996f, 0.054477f, 0.054990f, 0.055535f, 0.056112f, 0.056718f, 0.057353f, 0.058016f, 0.058706f,
        0.059421f, 0.060161f, 0.060924f, 0.061710f, 0.062517f, 0.063342f, 0.064181f, 0.065033f, 0.065892f, 0.066757f,
        0.067623f, 0.068488f, 0.069348f, 0.070200f, 0.071041f, 0.071871f, 0.072689f, 0.073495f, 0.074288f, 0.075069f,
        0.075837f, 0.076591f, 0.077332f, 0.078060f, 0.078773f, 0.079468f, 0.080143f, 0.080792f, 0.081412f, 0.082001f,
        0.082554f, 0.083067f, 0.083537f, 0.083960f, 0.084334f, 0.084657f, 0.084932f, 0.085159f, 0.085337f, 0.085469f,
        0.085555f, 0.085594f, 0.085589f, 0.085540f, 0.085448f, 0.085318f, 0.085156f, 0.084969f, 0.084760f, 0.084538f,
        0.084307f, 0.084073f, 0.083842f, 0.083620f, 0.083413f, 0.083226f, 0.083066f, 0.082937f, 0.082847f, 0.082801f,
        0.082804f, 0.082863f, 0.082983f, 0.083170f, 0.082468f, 0.082356f, 0.082245f, 0.082134f, 0.082023f, 0.081912f,
        0.081801f, 0.081691f, 0.081580f, 0.081470f, 0.081360f, 0.081250f, 0.081140f, 0.081030f, 0.080920f, 0.080811f,
        0.080701f, 0.080592f, 0.080483f, 0.080374f, 0.080265f, 0.080156f, 0.080048f, 0.079939f, 0.079831f, 0.079723f,
        0.079615f, 0.079507f, 0.079399f, 0.079291f, 0.079184f, 0.079076f, 0.078969f, 0.078862f, 0.078755f, 0.078648f,
        0.078542f, 0.078435f, 0.078329f, 0.078222f, 0.078116f, 0.078010f, 0.077904f, 0.077798f, 0.077693f, 0.077587f,
        0.077482f, 0.077376f, 0.077271f});
    static constexpr Spectrum CES59({0.184710f, 0.194010f, 0.203670f, 0.213670f, 0.224030f, 0.234740f, 0.245800f,
        0.257210f, 0.268960f, 0.281040f, 0.293450f, 0.306170f, 0.319190f, 0.332500f, 0.346090f, 0.359930f, 0.374010f,
        0.388300f, 0.402790f, 0.417450f, 0.434500f, 0.445880f, 0.460060f, 0.476460f, 0.494490f, 0.513570f, 0.533110f,
        0.552530f, 0.571240f, 0.588660f, 0.604200f, 0.617420f, 0.628410f, 0.637410f, 0.644660f, 0.650380f, 0.654830f,
        0.658230f, 0.660810f, 0.662830f, 0.664500f, 0.666030f, 0.667460f, 0.668790f, 0.670030f, 0.671170f, 0.672230f,
        0.673190f, 0.674080f, 0.674880f, 0.675600f, 0.676250f, 0.676830f, 0.677360f, 0.677850f, 0.678300f, 0.678740f,
        0.679170f, 0.679590f, 0.680040f, 0.680500f, 0.681000f, 0.681530f, 0.682100f, 0.682690f, 0.683320f, 0.683980f,
        0.684670f, 0.685390f, 0.686130f, 0.686900f, 0.687700f, 0.688520f, 0.689380f, 0.690270f, 0.691190f, 0.692160f,
        0.693170f, 0.694230f, 0.695340f, 0.696500f, 0.697720f, 0.699000f, 0.700350f, 0.701790f, 0.703320f, 0.704940f,
        0.706670f, 0.708520f, 0.710490f, 0.712600f, 0.714840f, 0.717230f, 0.719740f, 0.722370f, 0.725130f, 0.728000f,
        0.730970f, 0.734050f, 0.737230f, 0.740500f, 0.743850f, 0.747290f, 0.750790f, 0.754350f, 0.757960f, 0.761620f,
        0.765310f, 0.769020f, 0.772760f, 0.776500f, 0.780240f, 0.783960f, 0.787630f, 0.791230f, 0.794730f, 0.798120f,
        0.801360f, 0.804440f, 0.807330f, 0.810000f, 0.812440f, 0.814640f, 0.816620f, 0.818370f, 0.819900f, 0.821220f,
        0.822320f, 0.823210f, 0.823910f, 0.824400f, 0.824700f, 0.824810f, 0.824740f, 0.824490f, 0.824070f, 0.823490f,
        0.822740f, 0.821840f, 0.820790f, 0.819600f, 0.818270f, 0.816800f, 0.815200f, 0.813470f, 0.811600f, 0.809600f,
        0.807470f, 0.805210f, 0.802820f, 0.800300f, 0.797660f, 0.794890f, 0.792010f, 0.789020f, 0.785910f, 0.782700f,
        0.779390f, 0.775990f, 0.772490f, 0.768900f, 0.765230f, 0.761470f, 0.757620f, 0.753670f, 0.749630f, 0.745500f,
        0.741260f, 0.736910f, 0.732460f, 0.727900f, 0.723230f, 0.718460f, 0.713620f, 0.708730f, 0.703800f, 0.698850f,
        0.693910f, 0.688990f, 0.684110f, 0.679300f, 0.674570f, 0.669920f, 0.665370f, 0.660910f, 0.656570f, 0.652330f,
        0.648210f, 0.644210f, 0.640340f, 0.636600f, 0.633000f, 0.629540f, 0.626200f, 0.622990f, 0.619900f, 0.616930f,
        0.614070f, 0.611310f, 0.608660f, 0.606100f, 0.603630f, 0.601250f, 0.598940f, 0.596700f, 0.594510f, 0.592370f,
        0.590270f, 0.588200f, 0.586140f, 0.584100f, 0.582060f, 0.580020f, 0.577990f, 0.575960f, 0.573940f, 0.571930f,
        0.569930f, 0.567940f, 0.565960f, 0.564000f, 0.562050f, 0.560130f, 0.558240f, 0.556400f, 0.554610f, 0.552880f,
        0.551230f, 0.549660f, 0.548180f, 0.546800f, 0.545530f, 0.544370f, 0.543290f, 0.542310f, 0.541400f, 0.540560f,
        0.539770f, 0.539040f, 0.538350f, 0.537700f, 0.537070f, 0.536460f, 0.535870f, 0.535300f, 0.534730f, 0.534170f,
        0.533610f, 0.533050f, 0.532480f, 0.531900f, 0.531310f, 0.530720f, 0.530150f, 0.529590f, 0.529080f, 0.528610f,
        0.528210f, 0.527880f, 0.527640f, 0.527500f, 0.527470f, 0.527550f, 0.527740f, 0.528020f, 0.528400f, 0.528880f,
        0.529440f, 0.530080f, 0.530800f, 0.531600f, 0.532460f, 0.533370f, 0.534270f, 0.535160f, 0.535990f, 0.536730f,
        0.537360f, 0.537830f, 0.538120f, 0.538200f, 0.538050f, 0.537690f, 0.537180f, 0.536540f, 0.535830f, 0.535080f,
        0.534350f, 0.533660f, 0.533060f, 0.532600f, 0.532300f, 0.532170f, 0.532180f, 0.532330f, 0.532600f, 0.532980f,
        0.533460f, 0.534010f, 0.534630f, 0.535300f, 0.536010f, 0.536710f, 0.537370f, 0.537940f, 0.538390f, 0.538660f,
        0.538720f, 0.538530f, 0.538030f, 0.537200f, 0.536010f, 0.534510f, 0.532770f, 0.530890f, 0.528920f, 0.526940f,
        0.525040f, 0.523280f, 0.521740f, 0.520500f, 0.518640f, 0.517020f, 0.515400f, 0.513780f, 0.512160f, 0.510540f,
        0.508920f, 0.507300f, 0.505680f, 0.504060f, 0.502440f, 0.500810f, 0.499190f, 0.497570f, 0.495950f, 0.494330f,
        0.492710f, 0.491090f, 0.489460f, 0.487840f, 0.486220f, 0.484600f, 0.482980f, 0.481360f, 0.479740f, 0.478130f,
        0.476510f, 0.474890f, 0.473270f, 0.471660f, 0.470040f, 0.468420f, 0.466810f, 0.465190f, 0.463580f, 0.461970f,
        0.460360f, 0.458750f, 0.457140f, 0.455530f, 0.453920f, 0.452310f, 0.450700f, 0.449100f, 0.447490f, 0.445890f,
        0.444290f, 0.442690f, 0.441090f, 0.439490f, 0.437890f, 0.436300f, 0.434700f, 0.433110f, 0.431520f, 0.429930f,
        0.428340f, 0.426750f, 0.425160f, 0.423580f, 0.422000f, 0.420410f, 0.418840f, 0.417260f, 0.415680f, 0.414110f,
        0.412530f, 0.410960f, 0.409390f, 0.407830f, 0.406260f, 0.404700f, 0.403130f, 0.401570f, 0.400020f, 0.398460f,
        0.396910f, 0.395360f, 0.393810f});
    static constexpr Spectrum CES60({0.437100f, 0.444820f, 0.452220f, 0.459310f, 0.466140f, 0.472720f, 0.479100f,
        0.485300f, 0.491350f, 0.497280f, 0.503120f, 0.508890f, 0.514600f, 0.520240f, 0.525810f, 0.531290f, 0.536690f,
        0.542000f, 0.547220f, 0.552340f, 0.557350f, 0.562260f, 0.567060f, 0.571750f, 0.576340f, 0.580830f, 0.585230f,
        0.589520f, 0.593720f, 0.597820f, 0.601830f, 0.605750f, 0.609590f, 0.613340f, 0.617020f, 0.620630f, 0.624170f,
        0.627650f, 0.631070f, 0.634440f, 0.637770f, 0.641050f, 0.644300f, 0.647500f, 0.650670f, 0.653810f, 0.656920f,
        0.660000f, 0.663050f, 0.666080f, 0.669080f, 0.672060f, 0.675020f, 0.677950f, 0.680840f, 0.683700f, 0.686510f,
        0.689270f, 0.691970f, 0.694620f, 0.697210f, 0.699730f, 0.702180f, 0.704580f, 0.706910f, 0.709180f, 0.711390f,
        0.713550f, 0.715660f, 0.717720f, 0.719730f, 0.721700f, 0.723620f, 0.725490f, 0.727330f, 0.729120f, 0.730880f,
        0.732590f, 0.734270f, 0.735910f, 0.737510f, 0.739080f, 0.740610f, 0.742110f, 0.743570f, 0.745000f, 0.746390f,
        0.747740f, 0.749060f, 0.750340f, 0.751590f, 0.752800f, 0.753970f, 0.755110f, 0.756220f, 0.757290f, 0.758330f,
        0.759350f, 0.760330f, 0.761280f, 0.762200f, 0.763100f, 0.763960f, 0.764790f, 0.765590f, 0.766360f, 0.767080f,
        0.767760f, 0.768400f, 0.768990f, 0.769530f, 0.770020f, 0.770460f, 0.770840f, 0.771170f, 0.771450f, 0.771660f,
        0.771820f, 0.771920f, 0.771960f, 0.771930f, 0.771840f, 0.771690f, 0.771470f, 0.771180f, 0.770820f, 0.770400f,
        0.769900f, 0.769330f, 0.768680f, 0.767960f, 0.767160f, 0.766290f, 0.765340f, 0.764330f, 0.763240f, 0.762080f,
        0.760860f, 0.759570f, 0.758210f, 0.756800f, 0.755330f, 0.753800f, 0.752220f, 0.750600f, 0.748950f, 0.747270f,
        0.745560f, 0.743840f, 0.742100f, 0.740360f, 0.738620f, 0.736870f, 0.735120f, 0.733360f, 0.731590f, 0.729790f,
        0.727980f, 0.726140f, 0.724270f, 0.722370f, 0.720430f, 0.718470f, 0.716470f, 0.714450f, 0.712410f, 0.710350f,
        0.708290f, 0.706210f, 0.704140f, 0.702060f, 0.699990f, 0.697930f, 0.695890f, 0.693870f, 0.691890f, 0.689950f,
        0.688050f, 0.686200f, 0.684410f, 0.682680f, 0.681020f, 0.679440f, 0.677930f, 0.676500f, 0.675150f, 0.673870f,
        0.672670f, 0.671550f, 0.670510f, 0.669560f, 0.668690f, 0.667900f, 0.667170f, 0.666510f, 0.665890f, 0.665330f,
        0.664800f, 0.664300f, 0.663820f, 0.663360f, 0.662900f, 0.662450f, 0.661990f, 0.661540f, 0.661070f, 0.660600f,
        0.660110f, 0.659610f, 0.659090f, 0.658540f, 0.657970f, 0.657380f, 0.656770f, 0.656160f, 0.655530f, 0.654910f,
        0.654290f, 0.653670f, 0.653070f, 0.652490f, 0.651930f, 0.651400f, 0.650900f, 0.650440f, 0.650020f, 0.649650f,
        0.649330f, 0.649070f, 0.648870f, 0.648740f, 0.648680f, 0.648690f, 0.648750f, 0.648870f, 0.649030f, 0.649240f,
        0.649480f, 0.649740f, 0.650020f, 0.650320f, 0.650630f, 0.650940f, 0.651270f, 0.651620f, 0.651990f, 0.652380f,
        0.652810f, 0.653270f, 0.653760f, 0.654300f, 0.654880f, 0.655510f, 0.656170f, 0.656890f, 0.657640f, 0.658440f,
        0.659270f, 0.660150f, 0.661060f, 0.662020f, 0.663010f, 0.664050f, 0.665120f, 0.666230f, 0.667380f, 0.668560f,
        0.669790f, 0.671060f, 0.672370f, 0.673720f, 0.675110f, 0.676530f, 0.677980f, 0.679440f, 0.680900f, 0.682370f,
        0.683820f, 0.685250f, 0.686660f, 0.688030f, 0.689360f, 0.690640f, 0.691870f, 0.693050f, 0.694190f, 0.695280f,
        0.696310f, 0.697300f, 0.698230f, 0.699110f, 0.699940f, 0.700710f, 0.701450f, 0.702140f, 0.702800f, 0.703440f,
        0.704050f, 0.704640f, 0.705220f, 0.705790f, 0.706360f, 0.706910f, 0.707450f, 0.707950f, 0.708420f, 0.708850f,
        0.709210f, 0.709520f, 0.709760f, 0.709910f, 0.709980f, 0.709970f, 0.709900f, 0.709770f, 0.709600f, 0.709400f,
        0.709180f, 0.708960f, 0.708730f, 0.708530f, 0.708350f, 0.708210f, 0.708110f, 0.708060f, 0.708070f, 0.708140f,
        0.708270f, 0.708490f, 0.708790f, 0.709180f, 0.709670f, 0.710250f, 0.710940f, 0.711730f, 0.712610f, 0.713600f,
        0.714690f, 0.715890f, 0.717190f, 0.718590f, 0.720100f, 0.721710f, 0.723410f, 0.725210f, 0.727100f, 0.729060f,
        0.731110f, 0.733220f, 0.735410f, 0.737660f, 0.739970f, 0.742380f, 0.744900f, 0.747580f, 0.750430f, 0.753500f,
        0.756810f, 0.760390f, 0.764280f, 0.768500f, 0.773060f, 0.777870f, 0.782810f, 0.787780f, 0.792650f, 0.797310f,
        0.801640f, 0.805530f, 0.808860f, 0.811530f, 0.813450f, 0.814670f, 0.815310f, 0.815460f, 0.815210f, 0.814670f,
        0.813940f, 0.813100f, 0.812270f, 0.811530f, 0.810970f, 0.810590f, 0.810360f, 0.810270f, 0.810300f, 0.810430f,
        0.810640f, 0.810900f, 0.811210f, 0.811530f, 0.811850f, 0.812160f, 0.812420f, 0.812630f, 0.812760f, 0.812790f,
        0.812700f, 0.812470f, 0.812090f});
    static constexpr Spectrum CES61({0.288730f, 0.292590f, 0.296480f, 0.300390f, 0.304340f, 0.308310f, 0.312310f,
        0.316340f, 0.320400f, 0.324490f, 0.328600f, 0.332740f, 0.336900f, 0.341100f, 0.345310f, 0.349550f, 0.353810f,
        0.358100f, 0.362410f, 0.366740f, 0.369770f, 0.374720f, 0.379600f, 0.384410f, 0.389160f, 0.393830f, 0.398430f,
        0.402960f, 0.407420f, 0.411800f, 0.416110f, 0.420350f, 0.424510f, 0.428600f, 0.432600f, 0.436540f, 0.440390f,
        0.444170f, 0.447860f, 0.451480f, 0.455020f, 0.458470f, 0.461850f, 0.465150f, 0.468370f, 0.471510f, 0.474590f,
        0.477590f, 0.480520f, 0.483380f, 0.486180f, 0.488910f, 0.491570f, 0.494160f, 0.496690f, 0.499150f, 0.501540f,
        0.503860f, 0.506100f, 0.508280f, 0.510390f, 0.512430f, 0.514400f, 0.516300f, 0.518140f, 0.519920f, 0.521640f,
        0.523320f, 0.524940f, 0.526520f, 0.528060f, 0.529550f, 0.531000f, 0.532390f, 0.533720f, 0.534970f, 0.536140f,
        0.537230f, 0.538220f, 0.539100f, 0.539870f, 0.540520f, 0.541060f, 0.541510f, 0.541890f, 0.542200f, 0.542480f,
        0.542720f, 0.542950f, 0.543190f, 0.543440f, 0.543730f, 0.544070f, 0.544490f, 0.545000f, 0.545630f, 0.546400f,
        0.547330f, 0.548440f, 0.549750f, 0.551280f, 0.553040f, 0.555030f, 0.557200f, 0.559530f, 0.562000f, 0.564590f,
        0.567250f, 0.569980f, 0.572730f, 0.575490f, 0.578230f, 0.580950f, 0.583620f, 0.586260f, 0.588840f, 0.591360f,
        0.593810f, 0.596190f, 0.598490f, 0.600700f, 0.602810f, 0.604830f, 0.606750f, 0.608580f, 0.610310f, 0.611950f,
        0.613490f, 0.614950f, 0.616300f, 0.617570f, 0.618740f, 0.619790f, 0.620670f, 0.621370f, 0.621840f, 0.622050f,
        0.621970f, 0.621570f, 0.620810f, 0.619650f, 0.618090f, 0.616130f, 0.613810f, 0.611170f, 0.608240f, 0.605050f,
        0.601640f, 0.598030f, 0.594270f, 0.590380f, 0.586390f, 0.582320f, 0.578170f, 0.573950f, 0.569650f, 0.565290f,
        0.560870f, 0.556400f, 0.551880f, 0.547310f, 0.542700f, 0.538070f, 0.533400f, 0.528700f, 0.523990f, 0.519260f,
        0.514520f, 0.509770f, 0.505020f, 0.500270f, 0.495530f, 0.490810f, 0.486150f, 0.481550f, 0.477040f, 0.472650f,
        0.468380f, 0.464270f, 0.460340f, 0.456600f, 0.453080f, 0.449760f, 0.446630f, 0.443670f, 0.440880f, 0.438240f,
        0.435740f, 0.433360f, 0.431090f, 0.428920f, 0.426830f, 0.424820f, 0.422900f, 0.421060f, 0.419300f, 0.417620f,
        0.416030f, 0.414520f, 0.413090f, 0.411750f, 0.410480f, 0.409300f, 0.408200f, 0.407170f, 0.406210f, 0.405320f,
        0.404500f, 0.403750f, 0.403050f, 0.402420f, 0.401840f, 0.401320f, 0.400830f, 0.400390f, 0.399980f, 0.399590f,
        0.399220f, 0.398860f, 0.398510f, 0.398150f, 0.397790f, 0.397420f, 0.397060f, 0.396700f, 0.396340f, 0.395990f,
        0.395660f, 0.395340f, 0.395050f, 0.394780f, 0.394530f, 0.394320f, 0.394140f, 0.393990f, 0.393870f, 0.393790f,
        0.393750f, 0.393760f, 0.393800f, 0.393880f, 0.394010f, 0.394190f, 0.394400f, 0.394660f, 0.394950f, 0.395270f,
        0.395630f, 0.396010f, 0.396420f, 0.396860f, 0.397320f, 0.397810f, 0.398320f, 0.398860f, 0.399430f, 0.400040f,
        0.400670f, 0.401350f, 0.402060f, 0.402820f, 0.403610f, 0.404440f, 0.405300f, 0.406190f, 0.407100f, 0.408020f,
        0.408950f, 0.409890f, 0.410820f, 0.411750f, 0.412660f, 0.413550f, 0.414420f, 0.415250f, 0.416040f, 0.416780f,
        0.417460f, 0.418080f, 0.418620f, 0.419090f, 0.419480f, 0.419780f, 0.420010f, 0.420170f, 0.420280f, 0.420320f,
        0.420320f, 0.420280f, 0.420200f, 0.420080f, 0.419950f, 0.419780f, 0.419600f, 0.419400f, 0.419170f, 0.418930f,
        0.418670f, 0.418400f, 0.418100f, 0.417800f, 0.417480f, 0.417150f, 0.416810f, 0.416450f, 0.416070f, 0.415680f,
        0.415270f, 0.414840f, 0.414400f, 0.413930f, 0.413440f, 0.412930f, 0.412400f, 0.411850f, 0.411270f, 0.410660f,
        0.410030f, 0.409370f, 0.408690f, 0.407980f, 0.407800f, 0.407250f, 0.406710f, 0.406160f, 0.405620f, 0.405070f,
        0.404530f, 0.403980f, 0.403440f, 0.402890f, 0.402350f, 0.401800f, 0.401260f, 0.400720f, 0.400170f, 0.399630f,
        0.399090f, 0.398540f, 0.398000f, 0.397460f, 0.396920f, 0.396380f, 0.395830f, 0.395290f, 0.394750f, 0.394210f,
        0.393670f, 0.393130f, 0.392590f, 0.392050f, 0.391510f, 0.390970f, 0.390440f, 0.389900f, 0.389360f, 0.388820f,
        0.388280f, 0.387750f, 0.387210f, 0.386670f, 0.386140f, 0.385600f, 0.385060f, 0.384530f, 0.383990f, 0.383460f,
        0.382920f, 0.382390f, 0.381850f, 0.381320f, 0.380790f, 0.380250f, 0.379720f, 0.379190f, 0.378650f, 0.378120f,
        0.377590f, 0.377060f, 0.376530f, 0.375990f, 0.375460f, 0.374930f, 0.374400f, 0.373870f, 0.373340f, 0.372810f,
        0.372290f, 0.371760f, 0.371230f, 0.370700f, 0.370170f, 0.369640f, 0.369120f, 0.368590f, 0.368060f, 0.367540f,
        0.367010f, 0.366490f, 0.365960f});
    static constexpr Spectrum CES62({0.170090f, 0.171900f, 0.173730f, 0.175570f, 0.177430f, 0.179300f, 0.181190f,
        0.183090f, 0.185010f, 0.186950f, 0.188900f, 0.190860f, 0.192850f, 0.194840f, 0.196860f, 0.198880f, 0.200930f,
        0.202990f, 0.205060f, 0.207160f, 0.208490f, 0.211000f, 0.213440f, 0.215820f, 0.218150f, 0.220430f, 0.222660f,
        0.224870f, 0.227040f, 0.229190f, 0.231330f, 0.233460f, 0.235580f, 0.237700f, 0.239840f, 0.241990f, 0.244160f,
        0.246350f, 0.248580f, 0.250860f, 0.253170f, 0.255540f, 0.257960f, 0.260430f, 0.262950f, 0.265510f, 0.268120f,
        0.270770f, 0.273470f, 0.276210f, 0.278990f, 0.281800f, 0.284670f, 0.287580f, 0.290540f, 0.293570f, 0.296660f,
        0.299810f, 0.303050f, 0.306360f, 0.309760f, 0.313250f, 0.316820f, 0.320470f, 0.324190f, 0.327960f, 0.331790f,
        0.335670f, 0.339580f, 0.343530f, 0.347490f, 0.351470f, 0.355460f, 0.359440f, 0.363400f, 0.367330f, 0.371230f,
        0.375080f, 0.378860f, 0.382580f, 0.386210f, 0.389750f, 0.393190f, 0.396520f, 0.399730f, 0.402810f, 0.405760f,
        0.408560f, 0.411200f, 0.413690f, 0.416000f, 0.418130f, 0.420080f, 0.421870f, 0.423480f, 0.424930f, 0.426230f,
        0.427370f, 0.428350f, 0.429200f, 0.429900f, 0.430460f, 0.430880f, 0.431160f, 0.431290f, 0.431280f, 0.431120f,
        0.430800f, 0.430330f, 0.429700f, 0.428900f, 0.427950f, 0.426840f, 0.425580f, 0.424170f, 0.422620f, 0.420940f,
        0.419140f, 0.417210f, 0.415170f, 0.413020f, 0.410770f, 0.408430f, 0.406010f, 0.403530f, 0.401010f, 0.398460f,
        0.395880f, 0.393310f, 0.390740f, 0.388200f, 0.385690f, 0.383230f, 0.380830f, 0.378500f, 0.376230f, 0.374060f,
        0.371970f, 0.369980f, 0.368110f, 0.366360f, 0.364720f, 0.363200f, 0.361770f, 0.360400f, 0.359080f, 0.357780f,
        0.356490f, 0.355190f, 0.353850f, 0.352460f, 0.350990f, 0.349460f, 0.347860f, 0.346210f, 0.344510f, 0.342780f,
        0.341010f, 0.339210f, 0.337400f, 0.335580f, 0.333750f, 0.331930f, 0.330130f, 0.328350f, 0.326600f, 0.324900f,
        0.323250f, 0.321660f, 0.320140f, 0.318700f, 0.317350f, 0.316100f, 0.314960f, 0.313930f, 0.313030f, 0.312270f,
        0.311650f, 0.311180f, 0.310880f, 0.310760f, 0.310810f, 0.311020f, 0.311360f, 0.311830f, 0.312380f, 0.313000f,
        0.313670f, 0.314360f, 0.315050f, 0.315720f, 0.316350f, 0.316920f, 0.317440f, 0.317880f, 0.318250f, 0.318540f,
        0.318740f, 0.318830f, 0.318820f, 0.318700f, 0.318460f, 0.318100f, 0.317650f, 0.317120f, 0.316510f, 0.315830f,
        0.315110f, 0.314340f, 0.313550f, 0.312740f, 0.311930f, 0.311120f, 0.310330f, 0.309550f, 0.308810f, 0.308100f,
        0.307440f, 0.306820f, 0.306270f, 0.305790f, 0.305390f, 0.305060f, 0.304800f, 0.304610f, 0.304490f, 0.304430f,
        0.304430f, 0.304500f, 0.304620f, 0.304800f, 0.305030f, 0.305300f, 0.305600f, 0.305930f, 0.306270f, 0.306600f,
        0.306940f, 0.307250f, 0.307530f, 0.307780f, 0.307980f, 0.308140f, 0.308260f, 0.308360f, 0.308430f, 0.308500f,
        0.308560f, 0.308620f, 0.308690f, 0.308770f, 0.308880f, 0.309010f, 0.309180f, 0.309390f, 0.309640f, 0.309950f,
        0.310300f, 0.310720f, 0.311200f, 0.311750f, 0.312380f, 0.313090f, 0.313890f, 0.314790f, 0.315790f, 0.316910f,
        0.318150f, 0.319520f, 0.321020f, 0.322670f, 0.324470f, 0.326430f, 0.328550f, 0.330840f, 0.333300f, 0.335950f,
        0.338790f, 0.341820f, 0.345050f, 0.348480f, 0.352130f, 0.356000f, 0.360080f, 0.364390f, 0.368930f, 0.373700f,
        0.378700f, 0.383940f, 0.389420f, 0.395150f, 0.401120f, 0.407350f, 0.413840f, 0.420580f, 0.427580f, 0.434850f,
        0.442390f, 0.450190f, 0.458270f, 0.466630f, 0.475270f, 0.484160f, 0.493290f, 0.502630f, 0.512160f, 0.521870f,
        0.531720f, 0.541700f, 0.551780f, 0.561940f, 0.572170f, 0.582430f, 0.592700f, 0.602970f, 0.613210f, 0.623400f,
        0.633520f, 0.643540f, 0.653450f, 0.663210f, 0.671790f, 0.681040f, 0.690160f, 0.699130f, 0.707940f, 0.716610f,
        0.725120f, 0.733470f, 0.741650f, 0.749670f, 0.757520f, 0.765200f, 0.772720f, 0.780060f, 0.787230f, 0.794220f,
        0.801050f, 0.807700f, 0.814180f, 0.820500f, 0.826640f, 0.832620f, 0.838430f, 0.844070f, 0.849560f, 0.854880f,
        0.860050f, 0.865060f, 0.869920f, 0.874630f, 0.879200f, 0.883620f, 0.887900f, 0.892040f, 0.896040f, 0.899910f,
        0.903660f, 0.907280f, 0.910770f, 0.914150f, 0.917410f, 0.920560f, 0.923600f, 0.926530f, 0.929360f, 0.932080f,
        0.934710f, 0.937240f, 0.939690f, 0.942040f, 0.944310f, 0.946490f, 0.948590f, 0.950610f, 0.952560f, 0.954440f,
        0.956240f, 0.957980f, 0.959640f, 0.961250f, 0.962800f, 0.964280f, 0.965710f, 0.967080f, 0.968400f, 0.969670f,
        0.970890f, 0.972060f, 0.973190f, 0.974270f, 0.975310f, 0.976310f, 0.977260f, 0.978190f, 0.979070f, 0.979920f,
        0.980730f, 0.981520f, 0.982270f});
    static constexpr Spectrum CES63({0.023525f, 0.024208f, 0.024909f, 0.025631f, 0.026372f, 0.027135f, 0.027919f,
        0.028725f, 0.029553f, 0.030405f, 0.031280f, 0.032180f, 0.033105f, 0.034055f, 0.035032f, 0.036036f, 0.037067f,
        0.038127f, 0.039216f, 0.040334f, 0.039557f, 0.041640f, 0.043603f, 0.045448f, 0.047180f, 0.048801f, 0.050316f,
        0.051728f, 0.053040f, 0.054256f, 0.055379f, 0.056414f, 0.057362f, 0.058229f, 0.059016f, 0.059729f, 0.060370f,
        0.060943f, 0.061452f, 0.061899f, 0.062289f, 0.062626f, 0.062914f, 0.063162f, 0.063375f, 0.063560f, 0.063724f,
        0.063872f, 0.064012f, 0.064150f, 0.064292f, 0.064444f, 0.064608f, 0.064786f, 0.064978f, 0.065185f, 0.065409f,
        0.065651f, 0.065912f, 0.066193f, 0.066495f, 0.066820f, 0.067168f, 0.067538f, 0.067932f, 0.068350f, 0.068791f,
        0.069257f, 0.069747f, 0.070262f, 0.070801f, 0.071367f, 0.071960f, 0.072584f, 0.073240f, 0.073932f, 0.074662f,
        0.075432f, 0.076246f, 0.077105f, 0.078012f, 0.078968f, 0.079969f, 0.081008f, 0.082079f, 0.083176f, 0.084292f,
        0.085421f, 0.086558f, 0.087695f, 0.088827f, 0.089948f, 0.091055f, 0.092146f, 0.093219f, 0.094271f, 0.095300f,
        0.096304f, 0.097280f, 0.098227f, 0.099142f, 0.100020f, 0.100860f, 0.101650f, 0.102390f, 0.103060f, 0.103660f,
        0.104190f, 0.104640f, 0.104990f, 0.105250f, 0.105410f, 0.105460f, 0.105420f, 0.105290f, 0.105080f, 0.104780f,
        0.104400f, 0.103950f, 0.103430f, 0.102850f, 0.102200f, 0.101500f, 0.100750f, 0.099938f, 0.099080f, 0.098175f,
        0.097225f, 0.096233f, 0.095202f, 0.094135f, 0.093034f, 0.091903f, 0.090746f, 0.089569f, 0.088376f, 0.087170f,
        0.085958f, 0.084742f, 0.083527f, 0.082318f, 0.081119f, 0.079930f, 0.078753f, 0.077589f, 0.076437f, 0.075300f,
        0.074177f, 0.073069f, 0.071977f, 0.070902f, 0.069844f, 0.068808f, 0.067796f, 0.066811f, 0.065859f, 0.064941f,
        0.064061f, 0.063224f, 0.062432f, 0.061688f, 0.060997f, 0.060359f, 0.059776f, 0.059247f, 0.058775f, 0.058359f,
        0.058001f, 0.057702f, 0.057462f, 0.057282f, 0.057164f, 0.057108f, 0.057114f, 0.057184f, 0.057319f, 0.057519f,
        0.057784f, 0.058117f, 0.058516f, 0.058985f, 0.059521f, 0.060122f, 0.060782f, 0.061498f, 0.062265f, 0.063077f,
        0.063931f, 0.064822f, 0.065745f, 0.066696f, 0.067669f, 0.068660f, 0.069662f, 0.070670f, 0.071677f, 0.072679f,
        0.073669f, 0.074641f, 0.075590f, 0.076510f, 0.077395f, 0.078240f, 0.079041f, 0.079793f, 0.080492f, 0.081132f,
        0.081710f, 0.082220f, 0.082658f, 0.083019f, 0.083300f, 0.083502f, 0.083628f, 0.083681f, 0.083664f, 0.083579f,
        0.083429f, 0.083217f, 0.082946f, 0.082618f, 0.082237f, 0.081808f, 0.081334f, 0.080822f, 0.080276f, 0.079700f,
        0.079101f, 0.078484f, 0.077852f, 0.077211f, 0.076565f, 0.075917f, 0.075267f, 0.074618f, 0.073970f, 0.073325f,
        0.072684f, 0.072049f, 0.071421f, 0.070801f, 0.070192f, 0.069593f, 0.069007f, 0.068436f, 0.067880f, 0.067342f,
        0.066823f, 0.066324f, 0.065847f, 0.065394f, 0.064965f, 0.064558f, 0.064173f, 0.063806f, 0.063457f, 0.063122f,
        0.062799f, 0.062488f, 0.062185f, 0.061889f, 0.061598f, 0.061315f, 0.061045f, 0.060790f, 0.060553f, 0.060340f,
        0.060153f, 0.059996f, 0.059872f, 0.059786f, 0.059739f, 0.059734f, 0.059768f, 0.059841f, 0.059954f, 0.060105f,
        0.060295f, 0.060522f, 0.060786f, 0.061088f, 0.061425f, 0.061798f, 0.062205f, 0.062644f, 0.063115f, 0.063615f,
        0.064145f, 0.064702f, 0.065286f, 0.065894f, 0.066528f, 0.067188f, 0.067878f, 0.068601f, 0.069360f, 0.070158f,
        0.070999f, 0.071885f, 0.072820f, 0.073806f, 0.074845f, 0.075935f, 0.077070f, 0.078247f, 0.079460f, 0.080706f,
        0.081980f, 0.083277f, 0.084593f, 0.085923f, 0.087264f, 0.088610f, 0.089957f, 0.091300f, 0.092636f, 0.093960f,
        0.095267f, 0.096552f, 0.097812f, 0.099042f, 0.100920f, 0.102390f, 0.103890f, 0.105410f, 0.106940f, 0.108500f,
        0.110080f, 0.111670f, 0.113290f, 0.114920f, 0.116580f, 0.118260f, 0.119950f, 0.121670f, 0.123410f, 0.125170f,
        0.126960f, 0.128760f, 0.130590f, 0.132440f, 0.134310f, 0.136200f, 0.138120f, 0.140050f, 0.142020f, 0.144000f,
        0.146010f, 0.148030f, 0.150090f, 0.152160f, 0.154260f, 0.156390f, 0.158530f, 0.160700f, 0.162900f, 0.165120f,
        0.167360f, 0.169630f, 0.171920f, 0.174240f, 0.176580f, 0.178940f, 0.181330f, 0.183750f, 0.186190f, 0.188650f,
        0.191140f, 0.193660f, 0.196200f, 0.198760f, 0.201350f, 0.203970f, 0.206610f, 0.209270f, 0.211970f, 0.214680f,
        0.217420f, 0.220190f, 0.222980f, 0.225800f, 0.228640f, 0.231510f, 0.234400f, 0.237310f, 0.240260f, 0.243220f,
        0.246210f, 0.249230f, 0.252270f, 0.255340f, 0.258430f, 0.261540f, 0.264680f, 0.267840f, 0.271020f, 0.274230f,
        0.277470f, 0.280720f, 0.284000f});
    static constexpr Spectrum CES64({0.169890f, 0.172280f, 0.174700f, 0.177140f, 0.179610f, 0.182110f, 0.184630f,
        0.187190f, 0.189760f, 0.192370f, 0.195000f, 0.197660f, 0.200350f, 0.203060f, 0.205800f, 0.208570f, 0.211370f,
        0.214190f, 0.217040f, 0.219920f, 0.223260f, 0.225500f, 0.228300f, 0.231540f, 0.235120f, 0.238930f, 0.242840f,
        0.246760f, 0.250570f, 0.254170f, 0.257430f, 0.260280f, 0.262750f, 0.264890f, 0.266760f, 0.268420f, 0.269920f,
        0.271320f, 0.272680f, 0.274040f, 0.275480f, 0.277030f, 0.278690f, 0.280460f, 0.282320f, 0.284270f, 0.286290f,
        0.288370f, 0.290500f, 0.292670f, 0.294880f, 0.297110f, 0.299370f, 0.301660f, 0.303980f, 0.306350f, 0.308770f,
        0.311250f, 0.313780f, 0.316390f, 0.319060f, 0.321810f, 0.324650f, 0.327580f, 0.330620f, 0.333770f, 0.337040f,
        0.340440f, 0.343970f, 0.347660f, 0.351500f, 0.355500f, 0.359670f, 0.363990f, 0.368470f, 0.373100f, 0.377880f,
        0.382810f, 0.387870f, 0.393080f, 0.398420f, 0.403890f, 0.409470f, 0.415150f, 0.420890f, 0.426680f, 0.432500f,
        0.438330f, 0.444140f, 0.449920f, 0.455650f, 0.461300f, 0.466850f, 0.472290f, 0.477580f, 0.482710f, 0.487650f,
        0.492390f, 0.496900f, 0.501150f, 0.505140f, 0.508830f, 0.512230f, 0.515330f, 0.518140f, 0.520640f, 0.522840f,
        0.524730f, 0.526320f, 0.527610f, 0.528580f, 0.529250f, 0.529610f, 0.529680f, 0.529460f, 0.528950f, 0.528180f,
        0.527140f, 0.525840f, 0.524290f, 0.522500f, 0.520470f, 0.518210f, 0.515740f, 0.513050f, 0.510160f, 0.507070f,
        0.503800f, 0.500350f, 0.496730f, 0.492950f, 0.489020f, 0.484930f, 0.480710f, 0.476350f, 0.471850f, 0.467230f,
        0.462500f, 0.457640f, 0.452680f, 0.447610f, 0.442440f, 0.437190f, 0.431850f, 0.426430f, 0.420950f, 0.415410f,
        0.409820f, 0.404180f, 0.398510f, 0.392810f, 0.387090f, 0.381350f, 0.375590f, 0.369820f, 0.364030f, 0.358220f,
        0.352400f, 0.346570f, 0.340720f, 0.334870f, 0.329010f, 0.323140f, 0.317270f, 0.311410f, 0.305560f, 0.299730f,
        0.293920f, 0.288140f, 0.282380f, 0.276670f, 0.271000f, 0.265380f, 0.259830f, 0.254340f, 0.248930f, 0.243620f,
        0.238400f, 0.233290f, 0.228300f, 0.223440f, 0.218710f, 0.214120f, 0.209660f, 0.205340f, 0.201150f, 0.197100f,
        0.193180f, 0.189390f, 0.185730f, 0.182210f, 0.178820f, 0.175550f, 0.172400f, 0.169360f, 0.166430f, 0.163600f,
        0.160850f, 0.158190f, 0.155610f, 0.153100f, 0.150650f, 0.148270f, 0.145940f, 0.143670f, 0.141450f, 0.139270f,
        0.137140f, 0.135050f, 0.133000f, 0.130990f, 0.129010f, 0.127060f, 0.125150f, 0.123290f, 0.121470f, 0.119700f,
        0.117990f, 0.116330f, 0.114730f, 0.113200f, 0.111740f, 0.110340f, 0.109010f, 0.107750f, 0.106570f, 0.105450f,
        0.104410f, 0.103440f, 0.102540f, 0.101710f, 0.100960f, 0.100270f, 0.099654f, 0.099094f, 0.098591f, 0.098139f,
        0.097733f, 0.097369f, 0.097043f, 0.096750f, 0.096486f, 0.096247f, 0.096033f, 0.095839f, 0.095665f, 0.095508f,
        0.095365f, 0.095234f, 0.095113f, 0.095000f, 0.094892f, 0.094791f, 0.094698f, 0.094615f, 0.094543f, 0.094483f,
        0.094438f, 0.094408f, 0.094395f, 0.094400f, 0.094425f, 0.094472f, 0.094541f, 0.094635f, 0.094753f, 0.094898f,
        0.095071f, 0.095273f, 0.095506f, 0.095770f, 0.096067f, 0.096394f, 0.096750f, 0.097131f, 0.097535f, 0.097959f,
        0.098401f, 0.098859f, 0.099329f, 0.099810f, 0.100300f, 0.100790f, 0.101290f, 0.101800f, 0.102310f, 0.102810f,
        0.103320f, 0.103830f, 0.104340f, 0.104840f, 0.105340f, 0.105830f, 0.106310f, 0.106790f, 0.107250f, 0.107700f,
        0.108130f, 0.108540f, 0.108940f, 0.109320f, 0.109680f, 0.110010f, 0.110310f, 0.110590f, 0.110840f, 0.111060f,
        0.111240f, 0.111400f, 0.111510f, 0.111590f, 0.111630f, 0.111640f, 0.111610f, 0.111570f, 0.111510f, 0.111430f,
        0.111360f, 0.111280f, 0.111210f, 0.111150f, 0.111070f, 0.111000f, 0.110930f, 0.110860f, 0.110790f, 0.110720f,
        0.110650f, 0.110580f, 0.110510f, 0.110440f, 0.110370f, 0.110300f, 0.110230f, 0.110160f, 0.110080f, 0.110010f,
        0.109940f, 0.109870f, 0.109800f, 0.109730f, 0.109660f, 0.109590f, 0.109520f, 0.109450f, 0.109380f, 0.109310f,
        0.109240f, 0.109170f, 0.109100f, 0.109030f, 0.108960f, 0.108890f, 0.108830f, 0.108760f, 0.108690f, 0.108620f,
        0.108550f, 0.108480f, 0.108410f, 0.108340f, 0.108270f, 0.108200f, 0.108130f, 0.108060f, 0.107990f, 0.107920f,
        0.107850f, 0.107790f, 0.107720f, 0.107650f, 0.107580f, 0.107510f, 0.107440f, 0.107370f, 0.107300f, 0.107230f,
        0.107170f, 0.107100f, 0.107030f, 0.106960f, 0.106890f, 0.106820f, 0.106750f, 0.106690f, 0.106620f, 0.106550f,
        0.106480f, 0.106410f, 0.106340f, 0.106280f, 0.106210f, 0.106140f, 0.106070f, 0.106000f, 0.105940f, 0.105870f,
        0.105800f, 0.105730f, 0.105660f});
    static constexpr Spectrum CES65({0.006016f, 0.006529f, 0.007086f, 0.007689f, 0.008344f, 0.009053f, 0.009823f,
        0.010657f, 0.011561f, 0.012540f, 0.013602f, 0.014752f, 0.015998f, 0.017348f, 0.018809f, 0.020390f, 0.022101f,
        0.023953f, 0.025955f, 0.028120f, 0.029489f, 0.032119f, 0.035136f, 0.038499f, 0.042163f, 0.046088f, 0.050229f,
        0.054545f, 0.058993f, 0.063529f, 0.068113f, 0.072700f, 0.077249f, 0.081717f, 0.086060f, 0.090237f, 0.094205f,
        0.097922f, 0.101340f, 0.104430f, 0.107130f, 0.109430f, 0.111350f, 0.112950f, 0.114250f, 0.115330f, 0.116210f,
        0.116940f, 0.117580f, 0.118170f, 0.118750f, 0.119360f, 0.120010f, 0.120680f, 0.121370f, 0.122070f, 0.122780f,
        0.123480f, 0.124180f, 0.124850f, 0.125500f, 0.126120f, 0.126720f, 0.127310f, 0.127900f, 0.128500f, 0.129120f,
        0.129780f, 0.130480f, 0.131230f, 0.132060f, 0.132950f, 0.133930f, 0.134980f, 0.136110f, 0.137330f, 0.138630f,
        0.140010f, 0.141480f, 0.143030f, 0.144670f, 0.146390f, 0.148200f, 0.150100f, 0.152090f, 0.154170f, 0.156320f,
        0.158570f, 0.160900f, 0.163310f, 0.165810f, 0.168390f, 0.171030f, 0.173700f, 0.176370f, 0.179020f, 0.181630f,
        0.184160f, 0.186580f, 0.188890f, 0.191030f, 0.193000f, 0.194780f, 0.196350f, 0.197690f, 0.198800f, 0.199650f,
        0.200240f, 0.200550f, 0.200560f, 0.200270f, 0.199650f, 0.198740f, 0.197550f, 0.196100f, 0.194430f, 0.192540f,
        0.190470f, 0.188240f, 0.185870f, 0.183390f, 0.180810f, 0.178160f, 0.175440f, 0.172670f, 0.169870f, 0.167050f,
        0.164230f, 0.161420f, 0.158630f, 0.155890f, 0.153190f, 0.150560f, 0.148000f, 0.145520f, 0.143130f, 0.140830f,
        0.138640f, 0.136550f, 0.134590f, 0.132750f, 0.131050f, 0.129470f, 0.128000f, 0.126640f, 0.125370f, 0.124180f,
        0.123070f, 0.122010f, 0.121010f, 0.120040f, 0.119100f, 0.118190f, 0.117290f, 0.116400f, 0.115520f, 0.114640f,
        0.113760f, 0.112860f, 0.111940f, 0.111010f, 0.110040f, 0.109060f, 0.108070f, 0.107070f, 0.106090f, 0.105120f,
        0.104170f, 0.103260f, 0.102390f, 0.101570f, 0.100820f, 0.100130f, 0.099521f, 0.098997f, 0.098566f, 0.098236f,
        0.098017f, 0.097915f, 0.097940f, 0.098098f, 0.098395f, 0.098818f, 0.099351f, 0.099976f, 0.100680f, 0.101440f,
        0.102240f, 0.103070f, 0.103920f, 0.104750f, 0.105560f, 0.106340f, 0.107060f, 0.107700f, 0.108270f, 0.108730f,
        0.109070f, 0.109270f, 0.109330f, 0.109220f, 0.108930f, 0.108470f, 0.107840f, 0.107050f, 0.106120f, 0.105040f,
        0.103830f, 0.102500f, 0.101050f, 0.099489f, 0.097830f, 0.096085f, 0.094268f, 0.092393f, 0.090474f, 0.088527f,
        0.086563f, 0.084599f, 0.082647f, 0.080723f, 0.078837f, 0.076992f, 0.075185f, 0.073417f, 0.071686f, 0.069990f,
        0.068329f, 0.066702f, 0.065108f, 0.063546f, 0.062013f, 0.060510f, 0.059034f, 0.057584f, 0.056159f, 0.054757f,
        0.053376f, 0.052015f, 0.050673f, 0.049347f, 0.048037f, 0.046743f, 0.045464f, 0.044200f, 0.042951f, 0.041718f,
        0.040500f, 0.039297f, 0.038109f, 0.036936f, 0.035778f, 0.034637f, 0.033514f, 0.032411f, 0.031331f, 0.030273f,
        0.029241f, 0.028235f, 0.027258f, 0.026312f, 0.025397f, 0.024515f, 0.023665f, 0.022848f, 0.022065f, 0.021316f,
        0.020601f, 0.019921f, 0.019276f, 0.018666f, 0.018093f, 0.017553f, 0.017046f, 0.016571f, 0.016126f, 0.015709f,
        0.015319f, 0.014955f, 0.014615f, 0.014298f, 0.014002f, 0.013725f, 0.013468f, 0.013229f, 0.013007f, 0.012801f,
        0.012609f, 0.012432f, 0.012267f, 0.012113f, 0.011971f, 0.011839f, 0.011717f, 0.011605f, 0.011502f, 0.011408f,
        0.011324f, 0.011248f, 0.011180f, 0.011120f, 0.011069f, 0.011024f, 0.010988f, 0.010958f, 0.010935f, 0.010920f,
        0.010911f, 0.010908f, 0.010912f, 0.010922f, 0.010938f, 0.010959f, 0.010987f, 0.011019f, 0.011057f, 0.011100f,
        0.011148f, 0.011201f, 0.011258f, 0.011319f, 0.011270f, 0.011300f, 0.011329f, 0.011359f, 0.011389f, 0.011418f,
        0.011448f, 0.011478f, 0.011508f, 0.011539f, 0.011569f, 0.011599f, 0.011630f, 0.011660f, 0.011691f, 0.011721f,
        0.011752f, 0.011783f, 0.011813f, 0.011844f, 0.011875f, 0.011907f, 0.011938f, 0.011969f, 0.012000f, 0.012032f,
        0.012063f, 0.012095f, 0.012127f, 0.012158f, 0.012190f, 0.012222f, 0.012254f, 0.012286f, 0.012318f, 0.012351f,
        0.012383f, 0.012415f, 0.012448f, 0.012480f, 0.012513f, 0.012546f, 0.012579f, 0.012612f, 0.012645f, 0.012678f,
        0.012711f, 0.012744f, 0.012777f, 0.012811f, 0.012844f, 0.012878f, 0.012912f, 0.012945f, 0.012979f, 0.013013f,
        0.013047f, 0.013081f, 0.013116f, 0.013150f, 0.013184f, 0.013219f, 0.013253f, 0.013288f, 0.013323f, 0.013358f,
        0.013393f, 0.013428f, 0.013463f, 0.013498f, 0.013533f, 0.013569f, 0.013604f, 0.013640f, 0.013675f, 0.013711f,
        0.013747f, 0.013783f, 0.013819f});
    static constexpr Spectrum CES66({0.099631f, 0.100770f, 0.101920f, 0.103080f, 0.104260f, 0.105440f, 0.106640f,
        0.107850f, 0.109070f, 0.110300f, 0.111550f, 0.112800f, 0.114070f, 0.115360f, 0.116650f, 0.117960f, 0.119280f,
        0.120610f, 0.121960f, 0.123320f, 0.124900f, 0.125950f, 0.127270f, 0.128810f, 0.130500f, 0.132310f, 0.134180f,
        0.136070f, 0.137910f, 0.139680f, 0.141300f, 0.142750f, 0.144040f, 0.145200f, 0.146240f, 0.147210f, 0.148130f,
        0.149020f, 0.149910f, 0.150830f, 0.151800f, 0.152850f, 0.153970f, 0.155180f, 0.156450f, 0.157800f, 0.159210f,
        0.160690f, 0.162240f, 0.163840f, 0.165500f, 0.167220f, 0.168990f, 0.170820f, 0.172710f, 0.174670f, 0.176680f,
        0.178760f, 0.180910f, 0.183120f, 0.185400f, 0.187750f, 0.190180f, 0.192690f, 0.195290f, 0.197980f, 0.200760f,
        0.203650f, 0.206650f, 0.209770f, 0.213000f, 0.216360f, 0.219840f, 0.223460f, 0.227200f, 0.231070f, 0.235070f,
        0.239200f, 0.243470f, 0.247870f, 0.252400f, 0.257070f, 0.261840f, 0.266710f, 0.271630f, 0.276600f, 0.281570f,
        0.286530f, 0.291460f, 0.296320f, 0.301100f, 0.305770f, 0.310300f, 0.314690f, 0.318920f, 0.322950f, 0.326790f,
        0.330390f, 0.333760f, 0.336870f, 0.339700f, 0.342240f, 0.344470f, 0.346420f, 0.348070f, 0.349430f, 0.350500f,
        0.351280f, 0.351770f, 0.351980f, 0.351900f, 0.351540f, 0.350900f, 0.349990f, 0.348810f, 0.347370f, 0.345680f,
        0.343740f, 0.341560f, 0.339150f, 0.336500f, 0.333630f, 0.330550f, 0.327270f, 0.323800f, 0.320160f, 0.316360f,
        0.312410f, 0.308330f, 0.304120f, 0.299800f, 0.295380f, 0.290880f, 0.286290f, 0.281640f, 0.276920f, 0.272150f,
        0.267330f, 0.262480f, 0.257600f, 0.252700f, 0.247790f, 0.242880f, 0.237980f, 0.233080f, 0.228200f, 0.223350f,
        0.218520f, 0.213730f, 0.208990f, 0.204300f, 0.199660f, 0.195090f, 0.190570f, 0.186120f, 0.181740f, 0.177430f,
        0.173180f, 0.169010f, 0.164920f, 0.160900f, 0.156960f, 0.153110f, 0.149330f, 0.145630f, 0.142020f, 0.138490f,
        0.135040f, 0.131680f, 0.128400f, 0.125200f, 0.122090f, 0.119060f, 0.116130f, 0.113290f, 0.110540f, 0.107890f,
        0.105340f, 0.102880f, 0.100540f, 0.098300f, 0.096169f, 0.094144f, 0.092220f, 0.090394f, 0.088664f, 0.087024f,
        0.085471f, 0.084002f, 0.082613f, 0.081300f, 0.080060f, 0.078888f, 0.077781f, 0.076735f, 0.075745f, 0.074807f,
        0.073918f, 0.073073f, 0.072268f, 0.071500f, 0.070764f, 0.070059f, 0.069382f, 0.068732f, 0.068107f, 0.067506f,
        0.066926f, 0.066367f, 0.065825f, 0.065300f, 0.064790f, 0.064295f, 0.063815f, 0.063350f, 0.062901f, 0.062468f,
        0.062051f, 0.061651f, 0.061267f, 0.060900f, 0.060550f, 0.060217f, 0.059902f, 0.059605f, 0.059325f, 0.059063f,
        0.058820f, 0.058594f, 0.058388f, 0.058200f, 0.058031f, 0.057879f, 0.057743f, 0.057621f, 0.057512f, 0.057414f,
        0.057325f, 0.057244f, 0.057170f, 0.057100f, 0.057034f, 0.056971f, 0.056911f, 0.056855f, 0.056802f, 0.056754f,
        0.056709f, 0.056668f, 0.056632f, 0.056600f, 0.056572f, 0.056549f, 0.056530f, 0.056515f, 0.056504f, 0.056497f,
        0.056493f, 0.056492f, 0.056495f, 0.056500f, 0.056508f, 0.056520f, 0.056536f, 0.056556f, 0.056581f, 0.056612f,
        0.056648f, 0.056692f, 0.056742f, 0.056800f, 0.056866f, 0.056939f, 0.057020f, 0.057105f, 0.057196f, 0.057292f,
        0.057390f, 0.057492f, 0.057595f, 0.057700f, 0.057805f, 0.057910f, 0.058015f, 0.058118f, 0.058221f, 0.058322f,
        0.058420f, 0.058517f, 0.058610f, 0.058700f, 0.058786f, 0.058869f, 0.058948f, 0.059024f, 0.059096f, 0.059164f,
        0.059228f, 0.059289f, 0.059346f, 0.059400f, 0.059450f, 0.059496f, 0.059538f, 0.059576f, 0.059609f, 0.059638f,
        0.059661f, 0.059680f, 0.059693f, 0.059700f, 0.059702f, 0.059698f, 0.059691f, 0.059680f, 0.059667f, 0.059653f,
        0.059638f, 0.059624f, 0.059611f, 0.059600f, 0.059585f, 0.059572f, 0.059559f, 0.059545f, 0.059532f, 0.059519f,
        0.059506f, 0.059492f, 0.059479f, 0.059466f, 0.059453f, 0.059439f, 0.059426f, 0.059413f, 0.059400f, 0.059386f,
        0.059373f, 0.059360f, 0.059347f, 0.059333f, 0.059320f, 0.059307f, 0.059294f, 0.059280f, 0.059267f, 0.059254f,
        0.059241f, 0.059227f, 0.059214f, 0.059201f, 0.059188f, 0.059175f, 0.059161f, 0.059148f, 0.059135f, 0.059122f,
        0.059109f, 0.059095f, 0.059082f, 0.059069f, 0.059056f, 0.059043f, 0.059030f, 0.059016f, 0.059003f, 0.058990f,
        0.058977f, 0.058964f, 0.058951f, 0.058937f, 0.058924f, 0.058911f, 0.058898f, 0.058885f, 0.058872f, 0.058858f,
        0.058845f, 0.058832f, 0.058819f, 0.058806f, 0.058793f, 0.058780f, 0.058767f, 0.058753f, 0.058740f, 0.058727f,
        0.058714f, 0.058701f, 0.058688f, 0.058675f, 0.058662f, 0.058649f, 0.058635f, 0.058622f, 0.058609f, 0.058596f,
        0.058583f, 0.058570f, 0.058557f});
    static constexpr Spectrum CES67({0.109270f, 0.110920f, 0.112600f, 0.114290f, 0.116010f, 0.117750f, 0.119520f,
        0.121300f, 0.123110f, 0.124950f, 0.126800f, 0.128680f, 0.130590f, 0.132510f, 0.134460f, 0.136440f, 0.138440f,
        0.140470f, 0.142510f, 0.144590f, 0.147000f, 0.148630f, 0.150660f, 0.153020f, 0.155630f, 0.158420f, 0.161300f,
        0.164190f, 0.167030f, 0.169720f, 0.172200f, 0.174400f, 0.176350f, 0.178070f, 0.179620f, 0.181020f, 0.182320f,
        0.183540f, 0.184740f, 0.185950f, 0.187200f, 0.188530f, 0.189940f, 0.191440f, 0.193010f, 0.194670f, 0.196410f,
        0.198240f, 0.200140f, 0.202130f, 0.204200f, 0.206350f, 0.208580f, 0.210880f, 0.213270f, 0.215730f, 0.218260f,
        0.220860f, 0.223540f, 0.226290f, 0.229100f, 0.231980f, 0.234940f, 0.237980f, 0.241100f, 0.244310f, 0.247630f,
        0.251050f, 0.254580f, 0.258230f, 0.262000f, 0.265900f, 0.269940f, 0.274100f, 0.278390f, 0.282810f, 0.287360f,
        0.292030f, 0.296830f, 0.301750f, 0.306800f, 0.311970f, 0.317230f, 0.322550f, 0.327920f, 0.333300f, 0.338670f,
        0.343980f, 0.349230f, 0.354380f, 0.359400f, 0.364270f, 0.368970f, 0.373480f, 0.377800f, 0.381920f, 0.385810f,
        0.389470f, 0.392880f, 0.396020f, 0.398900f, 0.401490f, 0.403780f, 0.405770f, 0.407440f, 0.408780f, 0.409790f,
        0.410460f, 0.410770f, 0.410720f, 0.410300f, 0.409500f, 0.408340f, 0.406820f, 0.404980f, 0.402820f, 0.400360f,
        0.397630f, 0.394630f, 0.391380f, 0.387900f, 0.384210f, 0.380320f, 0.376230f, 0.371970f, 0.367540f, 0.362960f,
        0.358230f, 0.353370f, 0.348390f, 0.343300f, 0.338110f, 0.332840f, 0.327490f, 0.322080f, 0.316610f, 0.311110f,
        0.305570f, 0.300020f, 0.294460f, 0.288900f, 0.283360f, 0.277840f, 0.272340f, 0.266880f, 0.261450f, 0.256070f,
        0.250740f, 0.245470f, 0.240250f, 0.235100f, 0.230020f, 0.225010f, 0.220080f, 0.215220f, 0.210430f, 0.205720f,
        0.201080f, 0.196510f, 0.192020f, 0.187600f, 0.183260f, 0.178990f, 0.174800f, 0.170690f, 0.166660f, 0.162700f,
        0.158830f, 0.155040f, 0.151330f, 0.147700f, 0.144160f, 0.140700f, 0.137330f, 0.134070f, 0.130900f, 0.127830f,
        0.124880f, 0.122030f, 0.119310f, 0.116700f, 0.114220f, 0.111860f, 0.109610f, 0.107480f, 0.105450f, 0.103530f,
        0.101710f, 0.099983f, 0.098348f, 0.096800f, 0.095334f, 0.093947f, 0.092633f, 0.091387f, 0.090206f, 0.089084f,
        0.088017f, 0.087001f, 0.086030f, 0.085100f, 0.084207f, 0.083348f, 0.082520f, 0.081721f, 0.080949f, 0.080201f,
        0.079474f, 0.078767f, 0.078076f, 0.077400f, 0.076736f, 0.076084f, 0.075446f, 0.074821f, 0.074211f, 0.073616f,
        0.073037f, 0.072474f, 0.071928f, 0.071400f, 0.070891f, 0.070401f, 0.069934f, 0.069489f, 0.069070f, 0.068676f,
        0.068311f, 0.067976f, 0.067672f, 0.067400f, 0.067162f, 0.066956f, 0.066778f, 0.066627f, 0.066498f, 0.066389f,
        0.066298f, 0.066221f, 0.066156f, 0.066100f, 0.066050f, 0.066005f, 0.065963f, 0.065925f, 0.065889f, 0.065853f,
        0.065818f, 0.065781f, 0.065742f, 0.065700f, 0.065654f, 0.065607f, 0.065559f, 0.065514f, 0.065472f, 0.065437f,
        0.065410f, 0.065394f, 0.065390f, 0.065400f, 0.065427f, 0.065469f, 0.065527f, 0.065599f, 0.065685f, 0.065784f,
        0.065896f, 0.066020f, 0.066155f, 0.066300f, 0.066455f, 0.066618f, 0.066788f, 0.066962f, 0.067138f, 0.067316f,
        0.067493f, 0.067667f, 0.067836f, 0.068000f, 0.068156f, 0.068303f, 0.068441f, 0.068569f, 0.068687f, 0.068794f,
        0.068889f, 0.068973f, 0.069043f, 0.069100f, 0.069143f, 0.069174f, 0.069193f, 0.069201f, 0.069200f, 0.069192f,
        0.069176f, 0.069155f, 0.069129f, 0.069100f, 0.069069f, 0.069035f, 0.068997f, 0.068956f, 0.068911f, 0.068861f,
        0.068805f, 0.068744f, 0.068675f, 0.068600f, 0.068517f, 0.068429f, 0.068338f, 0.068246f, 0.068155f, 0.068068f,
        0.067986f, 0.067913f, 0.067850f, 0.067800f, 0.067722f, 0.067655f, 0.067588f, 0.067521f, 0.067454f, 0.067388f,
        0.067321f, 0.067254f, 0.067188f, 0.067121f, 0.067055f, 0.066989f, 0.066922f, 0.066856f, 0.066790f, 0.066724f,
        0.066658f, 0.066592f, 0.066526f, 0.066460f, 0.066394f, 0.066328f, 0.066263f, 0.066197f, 0.066131f, 0.066066f,
        0.066000f, 0.065935f, 0.065870f, 0.065804f, 0.065739f, 0.065674f, 0.065609f, 0.065544f, 0.065479f, 0.065414f,
        0.065349f, 0.065284f, 0.065220f, 0.065155f, 0.065090f, 0.065026f, 0.064961f, 0.064897f, 0.064833f, 0.064768f,
        0.064704f, 0.064640f, 0.064576f, 0.064512f, 0.064447f, 0.064384f, 0.064320f, 0.064256f, 0.064192f, 0.064128f,
        0.064065f, 0.064001f, 0.063937f, 0.063874f, 0.063811f, 0.063747f, 0.063684f, 0.063621f, 0.063557f, 0.063494f,
        0.063431f, 0.063368f, 0.063305f, 0.063242f, 0.063179f, 0.063117f, 0.063054f, 0.062991f, 0.062929f, 0.062866f,
        0.062804f, 0.062741f, 0.062679f});
    static constexpr Spectrum CES68({0.221490f, 0.224330f, 0.227190f, 0.230080f, 0.233000f, 0.235940f, 0.238910f,
        0.241910f, 0.244930f, 0.247970f, 0.251040f, 0.254140f, 0.257260f, 0.260400f, 0.263570f, 0.266760f, 0.269980f,
        0.273230f, 0.276490f, 0.279790f, 0.282450f, 0.286090f, 0.289700f, 0.293280f, 0.296830f, 0.300360f, 0.303860f,
        0.307350f, 0.310810f, 0.314260f, 0.317690f, 0.321110f, 0.324510f, 0.327910f, 0.331290f, 0.334670f, 0.338040f,
        0.341400f, 0.344770f, 0.348130f, 0.351490f, 0.354860f, 0.358220f, 0.361580f, 0.364920f, 0.368260f, 0.371570f,
        0.374860f, 0.378130f, 0.381370f, 0.384570f, 0.387740f, 0.390860f, 0.393940f, 0.396980f, 0.399970f, 0.402910f,
        0.405800f, 0.408640f, 0.411420f, 0.414150f, 0.416810f, 0.419420f, 0.421970f, 0.424460f, 0.426900f, 0.429290f,
        0.431630f, 0.433920f, 0.436160f, 0.438360f, 0.440510f, 0.442600f, 0.444620f, 0.446550f, 0.448370f, 0.450070f,
        0.451640f, 0.453050f, 0.454300f, 0.455360f, 0.456240f, 0.456930f, 0.457460f, 0.457840f, 0.458100f, 0.458240f,
        0.458290f, 0.458260f, 0.458170f, 0.458040f, 0.457890f, 0.457730f, 0.457590f, 0.457490f, 0.457450f, 0.457500f,
        0.457650f, 0.457940f, 0.458370f, 0.458970f, 0.459760f, 0.460720f, 0.461820f, 0.463040f, 0.464350f, 0.465730f,
        0.467160f, 0.468600f, 0.470040f, 0.471440f, 0.472790f, 0.474060f, 0.475260f, 0.476350f, 0.477320f, 0.478170f,
        0.478880f, 0.479420f, 0.479800f, 0.479990f, 0.479990f, 0.479790f, 0.479420f, 0.478870f, 0.478160f, 0.477310f,
        0.476310f, 0.475180f, 0.473930f, 0.472570f, 0.471110f, 0.469520f, 0.467780f, 0.465880f, 0.463780f, 0.461480f,
        0.458940f, 0.456140f, 0.453070f, 0.449700f, 0.446010f, 0.442020f, 0.437770f, 0.433280f, 0.428560f, 0.423660f,
        0.418600f, 0.413400f, 0.408100f, 0.402710f, 0.397260f, 0.391770f, 0.386260f, 0.380730f, 0.375200f, 0.369690f,
        0.364210f, 0.358770f, 0.353400f, 0.348090f, 0.342870f, 0.337740f, 0.332680f, 0.327700f, 0.322790f, 0.317940f,
        0.313160f, 0.308430f, 0.303760f, 0.299140f, 0.294570f, 0.290060f, 0.285630f, 0.281280f, 0.277040f, 0.272920f,
        0.268940f, 0.265100f, 0.261420f, 0.257930f, 0.254620f, 0.251490f, 0.248540f, 0.245760f, 0.243130f, 0.240640f,
        0.238290f, 0.236070f, 0.233960f, 0.231960f, 0.230060f, 0.228260f, 0.226560f, 0.224970f, 0.223490f, 0.222120f,
        0.220870f, 0.219740f, 0.218730f, 0.217840f, 0.217080f, 0.216450f, 0.215920f, 0.215510f, 0.215190f, 0.214970f,
        0.214830f, 0.214770f, 0.214780f, 0.214850f, 0.214980f, 0.215160f, 0.215380f, 0.215640f, 0.215930f, 0.216250f,
        0.216580f, 0.216930f, 0.217280f, 0.217630f, 0.217980f, 0.218330f, 0.218660f, 0.219000f, 0.219340f, 0.219670f,
        0.220010f, 0.220350f, 0.220690f, 0.221040f, 0.221390f, 0.221740f, 0.222100f, 0.222460f, 0.222810f, 0.223160f,
        0.223490f, 0.223820f, 0.224140f, 0.224440f, 0.224720f, 0.224980f, 0.225240f, 0.225480f, 0.225710f, 0.225940f,
        0.226160f, 0.226370f, 0.226590f, 0.226810f, 0.227030f, 0.227260f, 0.227500f, 0.227770f, 0.228070f, 0.228400f,
        0.228780f, 0.229200f, 0.229670f, 0.230210f, 0.230810f, 0.231470f, 0.232180f, 0.232930f, 0.233730f, 0.234560f,
        0.235410f, 0.236280f, 0.237160f, 0.238040f, 0.238920f, 0.239780f, 0.240630f, 0.241450f, 0.242240f, 0.242980f,
        0.243670f, 0.244300f, 0.244860f, 0.245350f, 0.245760f, 0.246100f, 0.246360f, 0.246560f, 0.246690f, 0.246780f,
        0.246810f, 0.246810f, 0.246770f, 0.246690f, 0.246600f, 0.246480f, 0.246330f, 0.246170f, 0.245990f, 0.245790f,
        0.245570f, 0.245340f, 0.245100f, 0.244840f, 0.244570f, 0.244290f, 0.244000f, 0.243700f, 0.243380f, 0.243050f,
        0.242700f, 0.242330f, 0.241940f, 0.241540f, 0.241120f, 0.240680f, 0.240210f, 0.239720f, 0.239210f, 0.238680f,
        0.238120f, 0.237540f, 0.236930f, 0.236290f, 0.236170f, 0.235700f, 0.235220f, 0.234750f, 0.234280f, 0.233800f,
        0.233330f, 0.232860f, 0.232390f, 0.231920f, 0.231450f, 0.230980f, 0.230520f, 0.230050f, 0.229580f, 0.229120f,
        0.228650f, 0.228190f, 0.227720f, 0.227260f, 0.226800f, 0.226340f, 0.225870f, 0.225410f, 0.224950f, 0.224500f,
        0.224040f, 0.223580f, 0.223120f, 0.222670f, 0.222210f, 0.221750f, 0.221300f, 0.220850f, 0.220390f, 0.219940f,
        0.219490f, 0.219040f, 0.218590f, 0.218140f, 0.217690f, 0.217240f, 0.216790f, 0.216340f, 0.215900f, 0.215450f,
        0.215010f, 0.214560f, 0.214120f, 0.213680f, 0.213230f, 0.212790f, 0.212350f, 0.211910f, 0.211470f, 0.211030f,
        0.210590f, 0.210160f, 0.209720f, 0.209280f, 0.208850f, 0.208410f, 0.207980f, 0.207540f, 0.207110f, 0.206680f,
        0.206250f, 0.205810f, 0.205380f, 0.204950f, 0.204520f, 0.204100f, 0.203670f, 0.203240f, 0.202820f, 0.202390f,
        0.201960f, 0.201540f, 0.201120f});
    static constexpr Spectrum CES69({0.132250f, 0.133640f, 0.135050f, 0.136460f, 0.137890f, 0.139330f, 0.140780f,
        0.142250f, 0.143720f, 0.145210f, 0.146720f, 0.148230f, 0.149760f, 0.151310f, 0.152860f, 0.154430f, 0.156010f,
        0.157610f, 0.159210f, 0.160840f, 0.162200f, 0.163980f, 0.165750f, 0.167500f, 0.169250f, 0.170980f, 0.172710f,
        0.174440f, 0.176160f, 0.177890f, 0.179610f, 0.181340f, 0.183070f, 0.184810f, 0.186560f, 0.188310f, 0.190090f,
        0.191870f, 0.193670f, 0.195490f, 0.197330f, 0.199200f, 0.201080f, 0.202970f, 0.204880f, 0.206790f, 0.208700f,
        0.210610f, 0.212520f, 0.214410f, 0.216300f, 0.218160f, 0.220000f, 0.221830f, 0.223630f, 0.225400f, 0.227150f,
        0.228860f, 0.230550f, 0.232200f, 0.233810f, 0.235390f, 0.236930f, 0.238440f, 0.239920f, 0.241360f, 0.242770f,
        0.244150f, 0.245510f, 0.246830f, 0.248140f, 0.249410f, 0.250650f, 0.251840f, 0.252950f, 0.253990f, 0.254920f,
        0.255750f, 0.256450f, 0.257010f, 0.257410f, 0.257650f, 0.257730f, 0.257680f, 0.257510f, 0.257240f, 0.256880f,
        0.256460f, 0.255980f, 0.255470f, 0.254940f, 0.254410f, 0.253900f, 0.253430f, 0.253010f, 0.252670f, 0.252430f,
        0.252300f, 0.252300f, 0.252450f, 0.252770f, 0.253280f, 0.253950f, 0.254780f, 0.255730f, 0.256790f, 0.257930f,
        0.259150f, 0.260400f, 0.261690f, 0.262980f, 0.264250f, 0.265500f, 0.266720f, 0.267890f, 0.269020f, 0.270080f,
        0.271080f, 0.272010f, 0.272840f, 0.273590f, 0.274230f, 0.274780f, 0.275230f, 0.275590f, 0.275870f, 0.276060f,
        0.276170f, 0.276210f, 0.276170f, 0.276060f, 0.275890f, 0.275630f, 0.275260f, 0.274770f, 0.274130f, 0.273330f,
        0.272330f, 0.271130f, 0.269700f, 0.268020f, 0.266080f, 0.263900f, 0.261490f, 0.258880f, 0.256090f, 0.253150f,
        0.250080f, 0.246900f, 0.243630f, 0.240310f, 0.236940f, 0.233540f, 0.230120f, 0.226700f, 0.223270f, 0.219840f,
        0.216430f, 0.213050f, 0.209710f, 0.206400f, 0.203150f, 0.199950f, 0.196800f, 0.193690f, 0.190640f, 0.187630f,
        0.184660f, 0.181730f, 0.178850f, 0.176000f, 0.173200f, 0.170440f, 0.167750f, 0.165130f, 0.162600f, 0.160170f,
        0.157850f, 0.155650f, 0.153590f, 0.151680f, 0.149930f, 0.148340f, 0.146880f, 0.145570f, 0.144380f, 0.143310f,
        0.142350f, 0.141490f, 0.140720f, 0.140040f, 0.139440f, 0.138920f, 0.138510f, 0.138200f, 0.138020f, 0.137970f,
        0.138060f, 0.138310f, 0.138720f, 0.139320f, 0.140100f, 0.141050f, 0.142170f, 0.143420f, 0.144810f, 0.146310f,
        0.147910f, 0.149580f, 0.151330f, 0.153130f, 0.154960f, 0.156820f, 0.158690f, 0.160560f, 0.162420f, 0.164260f,
        0.166070f, 0.167830f, 0.169530f, 0.171160f, 0.172710f, 0.174190f, 0.175600f, 0.176930f, 0.178200f, 0.179410f,
        0.180570f, 0.181670f, 0.182720f, 0.183730f, 0.184700f, 0.185630f, 0.186510f, 0.187360f, 0.188160f, 0.188930f,
        0.189650f, 0.190330f, 0.190970f, 0.191560f, 0.192120f, 0.192630f, 0.193110f, 0.193550f, 0.193960f, 0.194340f,
        0.194690f, 0.195010f, 0.195310f, 0.195580f, 0.195840f, 0.196080f, 0.196320f, 0.196570f, 0.196820f, 0.197080f,
        0.197380f, 0.197700f, 0.198060f, 0.198470f, 0.198930f, 0.199430f, 0.199980f, 0.200560f, 0.201170f, 0.201800f,
        0.202450f, 0.203110f, 0.203780f, 0.204440f, 0.205100f, 0.205750f, 0.206380f, 0.206980f, 0.207560f, 0.208100f,
        0.208600f, 0.209050f, 0.209460f, 0.209800f, 0.210090f, 0.210320f, 0.210500f, 0.210630f, 0.210710f, 0.210760f,
        0.210770f, 0.210740f, 0.210700f, 0.210630f, 0.210540f, 0.210440f, 0.210320f, 0.210180f, 0.210030f, 0.209870f,
        0.209690f, 0.209500f, 0.209300f, 0.209080f, 0.208860f, 0.208620f, 0.208370f, 0.208110f, 0.207840f, 0.207560f,
        0.207260f, 0.206960f, 0.206630f, 0.206300f, 0.205950f, 0.205590f, 0.205220f, 0.204830f, 0.204420f, 0.204000f,
        0.203570f, 0.203120f, 0.202660f, 0.202180f, 0.202020f, 0.201640f, 0.201270f, 0.200890f, 0.200510f, 0.200130f,
        0.199760f, 0.199380f, 0.199000f, 0.198630f, 0.198250f, 0.197880f, 0.197510f, 0.197130f, 0.196760f, 0.196390f,
        0.196020f, 0.195650f, 0.195280f, 0.194910f, 0.194540f, 0.194170f, 0.193800f, 0.193430f, 0.193070f, 0.192700f,
        0.192330f, 0.191970f, 0.191600f, 0.191240f, 0.190880f, 0.190510f, 0.190150f, 0.189790f, 0.189420f, 0.189060f,
        0.188700f, 0.188340f, 0.187980f, 0.187620f, 0.187260f, 0.186910f, 0.186550f, 0.186190f, 0.185840f, 0.185480f,
        0.185120f, 0.184770f, 0.184410f, 0.184060f, 0.183710f, 0.183350f, 0.183000f, 0.182650f, 0.182300f, 0.181950f,
        0.181600f, 0.181250f, 0.180900f, 0.180550f, 0.180200f, 0.179850f, 0.179510f, 0.179160f, 0.178820f, 0.178470f,
        0.178120f, 0.177780f, 0.177440f, 0.177090f, 0.176750f, 0.176410f, 0.176070f, 0.175720f, 0.175380f, 0.175040f,
        0.174700f, 0.174360f, 0.174030f});
    static constexpr Spectrum CES70({0.317940f, 0.323480f, 0.329060f, 0.334680f, 0.340360f, 0.346080f, 0.351850f,
        0.357660f, 0.363520f, 0.369410f, 0.375340f, 0.381310f, 0.387320f, 0.393360f, 0.399440f, 0.405540f, 0.411680f,
        0.417840f, 0.424030f, 0.430240f, 0.435370f, 0.442130f, 0.448810f, 0.455430f, 0.461970f, 0.468450f, 0.474860f,
        0.481210f, 0.487490f, 0.493710f, 0.499880f, 0.505990f, 0.512040f, 0.518040f, 0.523990f, 0.529900f, 0.535750f,
        0.541560f, 0.547320f, 0.553040f, 0.558720f, 0.564360f, 0.569960f, 0.575510f, 0.581000f, 0.586440f, 0.591820f,
        0.597130f, 0.602370f, 0.607530f, 0.612610f, 0.617610f, 0.622520f, 0.627340f, 0.632070f, 0.636720f, 0.641270f,
        0.645740f, 0.650110f, 0.654390f, 0.658570f, 0.662660f, 0.666660f, 0.670560f, 0.674350f, 0.678050f, 0.681640f,
        0.685130f, 0.688510f, 0.691790f, 0.694950f, 0.698000f, 0.700930f, 0.703720f, 0.706380f, 0.708900f, 0.711260f,
        0.713450f, 0.715470f, 0.717310f, 0.718960f, 0.720410f, 0.721670f, 0.722750f, 0.723640f, 0.724380f, 0.724940f,
        0.725360f, 0.725630f, 0.725760f, 0.725760f, 0.725630f, 0.725390f, 0.725020f, 0.724530f, 0.723920f, 0.723200f,
        0.722360f, 0.721410f, 0.720340f, 0.719160f, 0.717870f, 0.716470f, 0.714960f, 0.713340f, 0.711610f, 0.709760f,
        0.707800f, 0.705720f, 0.703540f, 0.701230f, 0.698820f, 0.696290f, 0.693640f, 0.690880f, 0.688000f, 0.685010f,
        0.681900f, 0.678670f, 0.675330f, 0.671870f, 0.668290f, 0.664600f, 0.660800f, 0.656890f, 0.652880f, 0.648780f,
        0.644580f, 0.640300f, 0.635920f, 0.631470f, 0.626940f, 0.622330f, 0.617630f, 0.612840f, 0.607960f, 0.602980f,
        0.597900f, 0.592710f, 0.587420f, 0.582010f, 0.576480f, 0.570850f, 0.565130f, 0.559320f, 0.553430f, 0.547480f,
        0.541480f, 0.535440f, 0.529370f, 0.523270f, 0.517170f, 0.511060f, 0.504970f, 0.498880f, 0.492830f, 0.486800f,
        0.480820f, 0.474880f, 0.469010f, 0.463200f, 0.457460f, 0.451790f, 0.446200f, 0.440680f, 0.435230f, 0.429850f,
        0.424540f, 0.419290f, 0.414110f, 0.408990f, 0.403940f, 0.398960f, 0.394070f, 0.389270f, 0.384580f, 0.380000f,
        0.375540f, 0.371230f, 0.367050f, 0.363030f, 0.359180f, 0.355480f, 0.351940f, 0.348540f, 0.345280f, 0.342150f,
        0.339140f, 0.336250f, 0.333460f, 0.330780f, 0.328190f, 0.325700f, 0.323310f, 0.321010f, 0.318800f, 0.316700f,
        0.314690f, 0.312770f, 0.310960f, 0.309240f, 0.307630f, 0.306100f, 0.304670f, 0.303330f, 0.302080f, 0.300910f,
        0.299810f, 0.298800f, 0.297850f, 0.296980f, 0.296170f, 0.295430f, 0.294750f, 0.294120f, 0.293550f, 0.293030f,
        0.292560f, 0.292140f, 0.291760f, 0.291420f, 0.291110f, 0.290850f, 0.290620f, 0.290420f, 0.290270f, 0.290140f,
        0.290050f, 0.289990f, 0.289970f, 0.289970f, 0.290010f, 0.290080f, 0.290170f, 0.290280f, 0.290410f, 0.290550f,
        0.290710f, 0.290870f, 0.291040f, 0.291210f, 0.291380f, 0.291550f, 0.291710f, 0.291880f, 0.292060f, 0.292240f,
        0.292430f, 0.292630f, 0.292840f, 0.293070f, 0.293300f, 0.293560f, 0.293840f, 0.294150f, 0.294500f, 0.294880f,
        0.295300f, 0.295770f, 0.296300f, 0.296880f, 0.297520f, 0.298210f, 0.298960f, 0.299740f, 0.300570f, 0.301420f,
        0.302300f, 0.303190f, 0.304100f, 0.305020f, 0.305930f, 0.306840f, 0.307730f, 0.308590f, 0.309410f, 0.310180f,
        0.310910f, 0.311560f, 0.312140f, 0.312640f, 0.313050f, 0.313380f, 0.313630f, 0.313800f, 0.313920f, 0.313970f,
        0.313980f, 0.313940f, 0.313880f, 0.313780f, 0.313660f, 0.313520f, 0.313360f, 0.313190f, 0.312990f, 0.312790f,
        0.312570f, 0.312330f, 0.312080f, 0.311820f, 0.311550f, 0.311260f, 0.310970f, 0.310660f, 0.310340f, 0.310010f,
        0.309660f, 0.309300f, 0.308920f, 0.308520f, 0.308110f, 0.307680f, 0.307240f, 0.306770f, 0.306290f, 0.305790f,
        0.305270f, 0.304720f, 0.304160f, 0.303580f, 0.303410f, 0.302960f, 0.302500f, 0.302050f, 0.301600f, 0.301140f,
        0.300690f, 0.300240f, 0.299780f, 0.299330f, 0.298880f, 0.298430f, 0.297980f, 0.297530f, 0.297080f, 0.296630f,
        0.296180f, 0.295730f, 0.295290f, 0.294840f, 0.294390f, 0.293940f, 0.293500f, 0.293050f, 0.292600f, 0.292160f,
        0.291710f, 0.291270f, 0.290830f, 0.290380f, 0.289940f, 0.289500f, 0.289050f, 0.288610f, 0.288170f, 0.287730f,
        0.287290f, 0.286850f, 0.286410f, 0.285970f, 0.285530f, 0.285090f, 0.284650f, 0.284210f, 0.283770f, 0.283340f,
        0.282900f, 0.282460f, 0.282030f, 0.281590f, 0.281160f, 0.280720f, 0.280290f, 0.279850f, 0.279420f, 0.278990f,
        0.278550f, 0.278120f, 0.277690f, 0.277260f, 0.276830f, 0.276400f, 0.275970f, 0.275540f, 0.275110f, 0.274680f,
        0.274250f, 0.273820f, 0.273390f, 0.272960f, 0.272540f, 0.272110f, 0.271680f, 0.271260f, 0.270830f, 0.270410f,
        0.269980f, 0.269560f, 0.269140f});
    static constexpr Spectrum CES71({0.036458f, 0.037052f, 0.037656f, 0.038269f, 0.038891f, 0.039523f, 0.040165f,
        0.040817f, 0.041479f, 0.042151f, 0.042834f, 0.043527f, 0.044231f, 0.044946f, 0.045672f, 0.046409f, 0.047157f,
        0.047917f, 0.048689f, 0.049472f, 0.050400f, 0.051000f, 0.051761f, 0.052667f, 0.053702f, 0.054850f, 0.056093f,
        0.057416f, 0.058803f, 0.060236f, 0.061700f, 0.063182f, 0.064688f, 0.066225f, 0.067802f, 0.069429f, 0.071113f,
        0.072863f, 0.074689f, 0.076598f, 0.078600f, 0.080701f, 0.082898f, 0.085187f, 0.087563f, 0.090022f, 0.092560f,
        0.095171f, 0.097851f, 0.100600f, 0.103400f, 0.106260f, 0.109170f, 0.112120f, 0.115120f, 0.118140f, 0.121200f,
        0.124280f, 0.127380f, 0.130480f, 0.133600f, 0.136720f, 0.139840f, 0.142970f, 0.146110f, 0.149250f, 0.152410f,
        0.155580f, 0.158770f, 0.161970f, 0.165200f, 0.168450f, 0.171700f, 0.174960f, 0.178200f, 0.181410f, 0.184590f,
        0.187720f, 0.190790f, 0.193790f, 0.196700f, 0.199510f, 0.202220f, 0.204800f, 0.207240f, 0.209530f, 0.211650f,
        0.213590f, 0.215340f, 0.216880f, 0.218200f, 0.219290f, 0.220140f, 0.220780f, 0.221190f, 0.221400f, 0.221400f,
        0.221210f, 0.220820f, 0.220250f, 0.219500f, 0.218580f, 0.217490f, 0.216250f, 0.214860f, 0.213330f, 0.211680f,
        0.209890f, 0.208000f, 0.206000f, 0.203900f, 0.201710f, 0.199430f, 0.197060f, 0.194610f, 0.192060f, 0.189430f,
        0.186700f, 0.183890f, 0.180990f, 0.178000f, 0.174920f, 0.171770f, 0.168550f, 0.165280f, 0.161970f, 0.158630f,
        0.155270f, 0.151900f, 0.148540f, 0.145200f, 0.141890f, 0.138620f, 0.135400f, 0.132240f, 0.129160f, 0.126170f,
        0.123280f, 0.120490f, 0.117830f, 0.115300f, 0.112910f, 0.110660f, 0.108540f, 0.106550f, 0.104670f, 0.102910f,
        0.101260f, 0.099714f, 0.098262f, 0.096900f, 0.095621f, 0.094409f, 0.093246f, 0.092115f, 0.090998f, 0.089876f,
        0.088734f, 0.087552f, 0.086313f, 0.085000f, 0.083599f, 0.082110f, 0.080540f, 0.078893f, 0.077176f, 0.075392f,
        0.073547f, 0.071647f, 0.069696f, 0.067700f, 0.065665f, 0.063601f, 0.061518f, 0.059427f, 0.057338f, 0.055262f,
        0.053210f, 0.051191f, 0.049218f, 0.047300f, 0.045449f, 0.043679f, 0.042005f, 0.040444f, 0.039011f, 0.037721f,
        0.036589f, 0.035632f, 0.034863f, 0.034300f, 0.033954f, 0.033829f, 0.033925f, 0.034242f, 0.034780f, 0.035540f,
        0.036522f, 0.037726f, 0.039151f, 0.040800f, 0.042663f, 0.044699f, 0.046857f, 0.049088f, 0.051342f, 0.053568f,
        0.055717f, 0.057739f, 0.059583f, 0.061200f, 0.062549f, 0.063633f, 0.064462f, 0.065047f, 0.065402f, 0.065536f,
        0.065461f, 0.065189f, 0.064732f, 0.064100f, 0.063308f, 0.062378f, 0.061336f, 0.060206f, 0.059014f, 0.057785f,
        0.056544f, 0.055316f, 0.054126f, 0.053000f, 0.051957f, 0.050999f, 0.050122f, 0.049321f, 0.048593f, 0.047933f,
        0.047338f, 0.046803f, 0.046325f, 0.045900f, 0.045523f, 0.045192f, 0.044901f, 0.044648f, 0.044428f, 0.044238f,
        0.044074f, 0.043932f, 0.043809f, 0.043700f, 0.043603f, 0.043516f, 0.043438f, 0.043370f, 0.043309f, 0.043256f,
        0.043209f, 0.043168f, 0.043132f, 0.043100f, 0.043073f, 0.043057f, 0.043057f, 0.043082f, 0.043136f, 0.043228f,
        0.043364f, 0.043550f, 0.043793f, 0.044100f, 0.044477f, 0.044931f, 0.045467f, 0.046091f, 0.046808f, 0.047624f,
        0.048546f, 0.049579f, 0.050728f, 0.052000f, 0.053398f, 0.054919f, 0.056556f, 0.058304f, 0.060157f, 0.062110f,
        0.064156f, 0.066291f, 0.068507f, 0.070800f, 0.073163f, 0.075589f, 0.078069f, 0.080597f, 0.083163f, 0.085760f,
        0.088381f, 0.091016f, 0.093658f, 0.096300f, 0.098936f, 0.101580f, 0.104230f, 0.106910f, 0.109630f, 0.112400f,
        0.115230f, 0.118130f, 0.121120f, 0.124200f, 0.127380f, 0.130610f, 0.133840f, 0.137000f, 0.140040f, 0.142900f,
        0.145520f, 0.147850f, 0.149830f, 0.151400f, 0.153990f, 0.156220f, 0.158470f, 0.160740f, 0.163050f, 0.165380f,
        0.167730f, 0.170120f, 0.172520f, 0.174960f, 0.177430f, 0.179920f, 0.182430f, 0.184980f, 0.187550f, 0.190150f,
        0.192780f, 0.195440f, 0.198120f, 0.200830f, 0.203570f, 0.206330f, 0.209130f, 0.211950f, 0.214800f, 0.217670f,
        0.220580f, 0.223510f, 0.226470f, 0.229460f, 0.232470f, 0.235510f, 0.238580f, 0.241680f, 0.244800f, 0.247960f,
        0.251130f, 0.254340f, 0.257570f, 0.260830f, 0.264120f, 0.267430f, 0.270770f, 0.274130f, 0.277520f, 0.280940f,
        0.284380f, 0.287850f, 0.291340f, 0.294860f, 0.298400f, 0.301960f, 0.305550f, 0.309170f, 0.312800f, 0.316460f,
        0.320150f, 0.323850f, 0.327580f, 0.331330f, 0.335100f, 0.338890f, 0.342710f, 0.346540f, 0.350400f, 0.354270f,
        0.358160f, 0.362070f, 0.366000f, 0.369950f, 0.373910f, 0.377900f, 0.381900f, 0.385910f, 0.389940f, 0.393990f,
        0.398050f, 0.402120f, 0.406210f});
    static constexpr Spectrum CES72({0.377560f, 0.380670f, 0.383790f, 0.386910f, 0.390050f, 0.393190f, 0.396350f,
        0.399510f, 0.402680f, 0.405860f, 0.409050f, 0.412240f, 0.415440f, 0.418650f, 0.421870f, 0.425090f, 0.428320f,
        0.431560f, 0.434800f, 0.438040f, 0.440220f, 0.443920f, 0.447570f, 0.451160f, 0.454690f, 0.458150f, 0.461550f,
        0.464880f, 0.468130f, 0.471320f, 0.474430f, 0.477460f, 0.480410f, 0.483280f, 0.486070f, 0.488770f, 0.491380f,
        0.493910f, 0.496340f, 0.498670f, 0.500910f, 0.503050f, 0.505100f, 0.507060f, 0.508930f, 0.510730f, 0.512450f,
        0.514110f, 0.515710f, 0.517250f, 0.518740f, 0.520180f, 0.521580f, 0.522950f, 0.524270f, 0.525560f, 0.526820f,
        0.528060f, 0.529260f, 0.530450f, 0.531620f, 0.532770f, 0.533900f, 0.535010f, 0.536100f, 0.537160f, 0.538200f,
        0.539210f, 0.540180f, 0.541120f, 0.542030f, 0.542890f, 0.543720f, 0.544510f, 0.545270f, 0.545980f, 0.546670f,
        0.547310f, 0.547920f, 0.548490f, 0.549030f, 0.549540f, 0.550010f, 0.550450f, 0.550850f, 0.551220f, 0.551550f,
        0.551850f, 0.552110f, 0.552340f, 0.552540f, 0.552700f, 0.552820f, 0.552920f, 0.552970f, 0.553000f, 0.552990f,
        0.552950f, 0.552880f, 0.552780f, 0.552640f, 0.552470f, 0.552270f, 0.552040f, 0.551780f, 0.551490f, 0.551170f,
        0.550810f, 0.550420f, 0.550000f, 0.549550f, 0.549060f, 0.548550f, 0.548000f, 0.547420f, 0.546810f, 0.546160f,
        0.545480f, 0.544780f, 0.544030f, 0.543260f, 0.542460f, 0.541620f, 0.540760f, 0.539860f, 0.538940f, 0.537980f,
        0.537000f, 0.535990f, 0.534950f, 0.533890f, 0.532790f, 0.531670f, 0.530520f, 0.529330f, 0.528090f, 0.526820f,
        0.525500f, 0.524120f, 0.522700f, 0.521210f, 0.519670f, 0.518060f, 0.516410f, 0.514720f, 0.512980f, 0.511210f,
        0.509410f, 0.507590f, 0.505750f, 0.503900f, 0.502040f, 0.500180f, 0.498320f, 0.496460f, 0.494610f, 0.492780f,
        0.490960f, 0.489170f, 0.487400f, 0.485660f, 0.483950f, 0.482280f, 0.480630f, 0.479010f, 0.477420f, 0.475850f,
        0.474310f, 0.472780f, 0.471280f, 0.469790f, 0.468320f, 0.466870f, 0.465430f, 0.464020f, 0.462630f, 0.461260f,
        0.459920f, 0.458610f, 0.457330f, 0.456090f, 0.454870f, 0.453690f, 0.452540f, 0.451430f, 0.450340f, 0.449290f,
        0.448260f, 0.447270f, 0.446310f, 0.445370f, 0.444460f, 0.443580f, 0.442730f, 0.441920f, 0.441150f, 0.440420f,
        0.439730f, 0.439090f, 0.438490f, 0.437950f, 0.437460f, 0.437020f, 0.436630f, 0.436280f, 0.435980f, 0.435710f,
        0.435480f, 0.435280f, 0.435110f, 0.434960f, 0.434840f, 0.434740f, 0.434670f, 0.434620f, 0.434600f, 0.434610f,
        0.434650f, 0.434720f, 0.434820f, 0.434960f, 0.435130f, 0.435330f, 0.435560f, 0.435820f, 0.436100f, 0.436400f,
        0.436710f, 0.437050f, 0.437390f, 0.437740f, 0.438100f, 0.438470f, 0.438830f, 0.439190f, 0.439550f, 0.439900f,
        0.440230f, 0.440550f, 0.440860f, 0.441140f, 0.441410f, 0.441650f, 0.441870f, 0.442080f, 0.442270f, 0.442450f,
        0.442620f, 0.442790f, 0.442940f, 0.443100f, 0.443260f, 0.443420f, 0.443590f, 0.443770f, 0.443970f, 0.444190f,
        0.444430f, 0.444710f, 0.445020f, 0.445370f, 0.445760f, 0.446190f, 0.446660f, 0.447150f, 0.447670f, 0.448220f,
        0.448780f, 0.449350f, 0.449930f, 0.450520f, 0.451110f, 0.451690f, 0.452250f, 0.452800f, 0.453320f, 0.453800f,
        0.454240f, 0.454640f, 0.454980f, 0.455260f, 0.455470f, 0.455630f, 0.455720f, 0.455770f, 0.455770f, 0.455740f,
        0.455670f, 0.455590f, 0.455480f, 0.455360f, 0.455240f, 0.455120f, 0.454990f, 0.454870f, 0.454750f, 0.454640f,
        0.454540f, 0.454460f, 0.454390f, 0.454330f, 0.454300f, 0.454290f, 0.454300f, 0.454320f, 0.454360f, 0.454410f,
        0.454480f, 0.454560f, 0.454650f, 0.454750f, 0.454850f, 0.454970f, 0.455090f, 0.455210f, 0.455340f, 0.455470f,
        0.455600f, 0.455730f, 0.455860f, 0.455980f, 0.456050f, 0.456170f, 0.456280f, 0.456400f, 0.456520f, 0.456630f,
        0.456750f, 0.456860f, 0.456980f, 0.457090f, 0.457210f, 0.457320f, 0.457440f, 0.457560f, 0.457670f, 0.457790f,
        0.457900f, 0.458020f, 0.458130f, 0.458250f, 0.458360f, 0.458480f, 0.458600f, 0.458710f, 0.458830f, 0.458940f,
        0.459060f, 0.459170f, 0.459290f, 0.459410f, 0.459520f, 0.459640f, 0.459750f, 0.459870f, 0.459980f, 0.460100f,
        0.460220f, 0.460330f, 0.460450f, 0.460560f, 0.460680f, 0.460790f, 0.460910f, 0.461030f, 0.461140f, 0.461260f,
        0.461370f, 0.461490f, 0.461600f, 0.461720f, 0.461840f, 0.461950f, 0.462070f, 0.462180f, 0.462300f, 0.462410f,
        0.462530f, 0.462650f, 0.462760f, 0.462880f, 0.462990f, 0.463110f, 0.463230f, 0.463340f, 0.463460f, 0.463570f,
        0.463690f, 0.463800f, 0.463920f, 0.464040f, 0.464150f, 0.464270f, 0.464380f, 0.464500f, 0.464620f, 0.464730f,
        0.464850f, 0.464960f, 0.465080f});
    static constexpr Spectrum CES73({0.161600f, 0.165180f, 0.168830f, 0.172540f, 0.176310f, 0.180160f, 0.184060f,
        0.188030f, 0.192070f, 0.196170f, 0.200340f, 0.204570f, 0.208870f, 0.213240f, 0.217670f, 0.222170f, 0.226730f,
        0.231360f, 0.236060f, 0.240820f, 0.239980f, 0.247840f, 0.255080f, 0.261770f, 0.267950f, 0.273700f, 0.279080f,
        0.284140f, 0.288960f, 0.293580f, 0.298070f, 0.302490f, 0.306910f, 0.311380f, 0.315960f, 0.320720f, 0.325720f,
        0.331010f, 0.336670f, 0.342740f, 0.349300f, 0.356380f, 0.363970f, 0.372010f, 0.380460f, 0.389270f, 0.398420f,
        0.407840f, 0.417510f, 0.427360f, 0.437370f, 0.447480f, 0.457620f, 0.467710f, 0.477680f, 0.487440f, 0.496930f,
        0.506050f, 0.514740f, 0.522920f, 0.530510f, 0.537440f, 0.543740f, 0.549420f, 0.554510f, 0.559040f, 0.563040f,
        0.566520f, 0.569530f, 0.572080f, 0.574190f, 0.575910f, 0.577240f, 0.578210f, 0.578840f, 0.579160f, 0.579170f,
        0.578910f, 0.578390f, 0.577640f, 0.576680f, 0.575520f, 0.574170f, 0.572640f, 0.570940f, 0.569070f, 0.567040f,
        0.564850f, 0.562520f, 0.560030f, 0.557410f, 0.554660f, 0.551800f, 0.548830f, 0.545780f, 0.542670f, 0.539500f,
        0.536290f, 0.533070f, 0.529850f, 0.526630f, 0.523450f, 0.520290f, 0.517150f, 0.514040f, 0.510950f, 0.507880f,
        0.504830f, 0.501790f, 0.498770f, 0.495750f, 0.492750f, 0.489740f, 0.486720f, 0.483680f, 0.480590f, 0.477460f,
        0.474260f, 0.470990f, 0.467630f, 0.464180f, 0.460620f, 0.456960f, 0.453200f, 0.449360f, 0.445440f, 0.441440f,
        0.437370f, 0.433250f, 0.429080f, 0.424860f, 0.420610f, 0.416310f, 0.411990f, 0.407630f, 0.403230f, 0.398800f,
        0.394340f, 0.389850f, 0.385330f, 0.380780f, 0.376200f, 0.371590f, 0.366970f, 0.362330f, 0.357690f, 0.353040f,
        0.348390f, 0.343750f, 0.339120f, 0.334510f, 0.329920f, 0.325330f, 0.320760f, 0.316170f, 0.311570f, 0.306940f,
        0.302270f, 0.297550f, 0.292780f, 0.287940f, 0.283030f, 0.278050f, 0.273010f, 0.267930f, 0.262810f, 0.257660f,
        0.252500f, 0.247330f, 0.242160f, 0.237010f, 0.231880f, 0.226780f, 0.221740f, 0.216770f, 0.211870f, 0.207060f,
        0.202360f, 0.197790f, 0.193340f, 0.189050f, 0.184910f, 0.180940f, 0.177120f, 0.173450f, 0.169940f, 0.166580f,
        0.163370f, 0.160300f, 0.157380f, 0.154590f, 0.151950f, 0.149450f, 0.147090f, 0.144870f, 0.142790f, 0.140870f,
        0.139080f, 0.137450f, 0.135970f, 0.134640f, 0.133460f, 0.132420f, 0.131530f, 0.130760f, 0.130120f, 0.129600f,
        0.129200f, 0.128900f, 0.128690f, 0.128580f, 0.128550f, 0.128600f, 0.128720f, 0.128910f, 0.129140f, 0.129420f,
        0.129750f, 0.130100f, 0.130470f, 0.130860f, 0.131260f, 0.131660f, 0.132050f, 0.132410f, 0.132750f, 0.133040f,
        0.133290f, 0.133480f, 0.133600f, 0.133640f, 0.133600f, 0.133480f, 0.133290f, 0.133030f, 0.132710f, 0.132340f,
        0.131920f, 0.131470f, 0.130980f, 0.130470f, 0.129940f, 0.129390f, 0.128850f, 0.128310f, 0.127770f, 0.127260f,
        0.126760f, 0.126300f, 0.125880f, 0.125500f, 0.125180f, 0.124910f, 0.124690f, 0.124550f, 0.124470f, 0.124460f,
        0.124520f, 0.124670f, 0.124890f, 0.125200f, 0.125610f, 0.126100f, 0.126670f, 0.127330f, 0.128070f, 0.128900f,
        0.129800f, 0.130770f, 0.131830f, 0.132950f, 0.134140f, 0.135400f, 0.136720f, 0.138080f, 0.139490f, 0.140940f,
        0.142410f, 0.143910f, 0.145420f, 0.146950f, 0.148480f, 0.150010f, 0.151540f, 0.153080f, 0.154620f, 0.156160f,
        0.157700f, 0.159250f, 0.160790f, 0.162340f, 0.163890f, 0.165490f, 0.167160f, 0.168950f, 0.170890f, 0.173030f,
        0.175400f, 0.178040f, 0.180990f, 0.184280f, 0.187950f, 0.192000f, 0.196390f, 0.201120f, 0.206170f, 0.211520f,
        0.217150f, 0.223060f, 0.229220f, 0.235610f, 0.242230f, 0.249050f, 0.256050f, 0.263220f, 0.270550f, 0.278010f,
        0.285590f, 0.293280f, 0.301050f, 0.308890f, 0.316320f, 0.324310f, 0.332410f, 0.340600f, 0.348900f, 0.357280f,
        0.365750f, 0.374310f, 0.382950f, 0.391660f, 0.400440f, 0.409280f, 0.418180f, 0.427140f, 0.436150f, 0.445190f,
        0.454280f, 0.463390f, 0.472530f, 0.481690f, 0.490860f, 0.500030f, 0.509210f, 0.518380f, 0.527540f, 0.536680f,
        0.545790f, 0.554880f, 0.563920f, 0.572930f, 0.581880f, 0.590790f, 0.599630f, 0.608410f, 0.617120f, 0.625760f,
        0.634310f, 0.642780f, 0.651170f, 0.659460f, 0.667650f, 0.675750f, 0.683740f, 0.691620f, 0.699400f, 0.707060f,
        0.714600f, 0.722030f, 0.729340f, 0.736520f, 0.743580f, 0.750520f, 0.757330f, 0.764010f, 0.770570f, 0.776990f,
        0.783290f, 0.789450f, 0.795490f, 0.801400f, 0.807170f, 0.812820f, 0.818340f, 0.823740f, 0.829000f, 0.834150f,
        0.839160f, 0.844050f, 0.848830f, 0.853480f, 0.858010f, 0.862420f, 0.866720f, 0.870900f, 0.874970f, 0.878930f,
        0.882790f, 0.886530f, 0.890170f});
    static constexpr Spectrum CES74({0.033875f, 0.034064f, 0.034254f, 0.034445f, 0.034637f, 0.034830f, 0.035025f,
        0.035220f, 0.035416f, 0.035613f, 0.035812f, 0.036011f, 0.036212f, 0.036413f, 0.036616f, 0.036820f, 0.037024f,
        0.037230f, 0.037437f, 0.037645f, 0.037827f, 0.038052f, 0.038275f, 0.038496f, 0.038715f, 0.038932f, 0.039148f,
        0.039364f, 0.039579f, 0.039795f, 0.040011f, 0.040228f, 0.040447f, 0.040667f, 0.040890f, 0.041115f, 0.041343f,
        0.041575f, 0.041811f, 0.042050f, 0.042295f, 0.042544f, 0.042799f, 0.043059f, 0.043326f, 0.043599f, 0.043879f,
        0.044166f, 0.044461f, 0.044764f, 0.045075f, 0.045394f, 0.045724f, 0.046064f, 0.046416f, 0.046782f, 0.047161f,
        0.047557f, 0.047969f, 0.048399f, 0.048847f, 0.049316f, 0.049801f, 0.050299f, 0.050809f, 0.051326f, 0.051847f,
        0.052369f, 0.052889f, 0.053404f, 0.053911f, 0.054406f, 0.054887f, 0.055353f, 0.055800f, 0.056227f, 0.056631f,
        0.057011f, 0.057364f, 0.057689f, 0.057981f, 0.058242f, 0.058470f, 0.058669f, 0.058840f, 0.058984f, 0.059104f,
        0.059201f, 0.059276f, 0.059333f, 0.059371f, 0.059394f, 0.059404f, 0.059402f, 0.059392f, 0.059374f, 0.059353f,
        0.059330f, 0.059307f, 0.059287f, 0.059272f, 0.059264f, 0.059259f, 0.059253f, 0.059244f, 0.059226f, 0.059195f,
        0.059149f, 0.059083f, 0.058993f, 0.058875f, 0.058726f, 0.058548f, 0.058342f, 0.058110f, 0.057854f, 0.057577f,
        0.057281f, 0.056966f, 0.056637f, 0.056294f, 0.055939f, 0.055577f, 0.055211f, 0.054843f, 0.054478f, 0.054120f,
        0.053771f, 0.053436f, 0.053117f, 0.052819f, 0.052543f, 0.052290f, 0.052057f, 0.051842f, 0.051643f, 0.051460f,
        0.051289f, 0.051128f, 0.050977f, 0.050833f, 0.050694f, 0.050557f, 0.050419f, 0.050276f, 0.050126f, 0.049965f,
        0.049789f, 0.049597f, 0.049383f, 0.049145f, 0.048882f, 0.048595f, 0.048289f, 0.047969f, 0.047638f, 0.047300f,
        0.046960f, 0.046622f, 0.046290f, 0.045968f, 0.045660f, 0.045370f, 0.045099f, 0.044850f, 0.044627f, 0.044431f,
        0.044266f, 0.044135f, 0.044039f, 0.043982f, 0.043967f, 0.043991f, 0.044054f, 0.044154f, 0.044291f, 0.044463f,
        0.044668f, 0.044905f, 0.045174f, 0.045472f, 0.045797f, 0.046144f, 0.046504f, 0.046873f, 0.047241f, 0.047603f,
        0.047952f, 0.048280f, 0.048581f, 0.048847f, 0.049074f, 0.049259f, 0.049404f, 0.049510f, 0.049577f, 0.049605f,
        0.049595f, 0.049548f, 0.049464f, 0.049344f, 0.049189f, 0.049002f, 0.048788f, 0.048551f, 0.048294f, 0.048022f,
        0.047738f, 0.047448f, 0.047154f, 0.046862f, 0.046574f, 0.046293f, 0.046020f, 0.045757f, 0.045505f, 0.045267f,
        0.045044f, 0.044837f, 0.044648f, 0.044479f, 0.044331f, 0.044203f, 0.044094f, 0.044004f, 0.043931f, 0.043874f,
        0.043831f, 0.043803f, 0.043788f, 0.043784f, 0.043791f, 0.043807f, 0.043831f, 0.043860f, 0.043895f, 0.043932f,
        0.043971f, 0.044010f, 0.044047f, 0.044082f, 0.044112f, 0.044138f, 0.044161f, 0.044181f, 0.044199f, 0.044215f,
        0.044231f, 0.044247f, 0.044263f, 0.044280f, 0.044300f, 0.044322f, 0.044349f, 0.044382f, 0.044422f, 0.044471f,
        0.044529f, 0.044599f, 0.044681f, 0.044777f, 0.044888f, 0.045018f, 0.045169f, 0.045344f, 0.045545f, 0.045777f,
        0.046041f, 0.046341f, 0.046680f, 0.047060f, 0.047485f, 0.047955f, 0.048472f, 0.049037f, 0.049652f, 0.050318f,
        0.051036f, 0.051807f, 0.052632f, 0.053514f, 0.054453f, 0.055455f, 0.056525f, 0.057670f, 0.058895f, 0.060205f,
        0.061607f, 0.063107f, 0.064709f, 0.066421f, 0.068246f, 0.070189f, 0.072253f, 0.074440f, 0.076755f, 0.079200f,
        0.081779f, 0.084495f, 0.087350f, 0.090348f, 0.093492f, 0.096787f, 0.100240f, 0.103840f, 0.107620f, 0.111550f,
        0.115660f, 0.119950f, 0.124420f, 0.129070f, 0.133910f, 0.138940f, 0.144170f, 0.149600f, 0.155240f, 0.161080f,
        0.167140f, 0.173420f, 0.179920f, 0.186650f, 0.192590f, 0.199380f, 0.206350f, 0.213500f, 0.220830f, 0.228340f,
        0.236020f, 0.243880f, 0.251920f, 0.260130f, 0.268520f, 0.277070f, 0.285790f, 0.294670f, 0.303710f, 0.312900f,
        0.322240f, 0.331730f, 0.341350f, 0.351110f, 0.361000f, 0.371000f, 0.381120f, 0.391340f, 0.401650f, 0.412060f,
        0.422540f, 0.433090f, 0.443710f, 0.454370f, 0.465080f, 0.475820f, 0.486580f, 0.497360f, 0.508130f, 0.518900f,
        0.529660f, 0.540380f, 0.551070f, 0.561710f, 0.572290f, 0.582810f, 0.593260f, 0.603620f, 0.613880f, 0.624050f,
        0.634110f, 0.644050f, 0.653880f, 0.663570f, 0.673120f, 0.682540f, 0.691800f, 0.700920f, 0.709880f, 0.718680f,
        0.727310f, 0.735780f, 0.744070f, 0.752200f, 0.760140f, 0.767920f, 0.775510f, 0.782930f, 0.790160f, 0.797220f,
        0.804100f, 0.810810f, 0.817330f, 0.823680f, 0.829850f, 0.835860f, 0.841690f, 0.847350f, 0.852840f, 0.858170f,
        0.863340f, 0.868340f, 0.873190f});
    static constexpr Spectrum CES75({0.029193f, 0.030018f, 0.030865f, 0.031736f, 0.032631f, 0.033549f, 0.034493f,
        0.035463f, 0.036458f, 0.037481f, 0.038531f, 0.039609f, 0.040716f, 0.041853f, 0.043020f, 0.044218f, 0.045447f,
        0.046710f, 0.048005f, 0.049335f, 0.050900f, 0.051980f, 0.053339f, 0.054933f, 0.056720f, 0.058657f, 0.060701f,
        0.062810f, 0.064942f, 0.067052f, 0.069100f, 0.071052f, 0.072915f, 0.074706f, 0.076444f, 0.078145f, 0.079826f,
        0.081505f, 0.083198f, 0.084924f, 0.086700f, 0.088539f, 0.090442f, 0.092405f, 0.094426f, 0.096502f, 0.098630f,
        0.100810f, 0.103030f, 0.105290f, 0.107600f, 0.109940f, 0.112320f, 0.114720f, 0.117160f, 0.119620f, 0.122100f,
        0.124610f, 0.127120f, 0.129660f, 0.132200f, 0.134750f, 0.137320f, 0.139900f, 0.142490f, 0.145120f, 0.147760f,
        0.150440f, 0.153160f, 0.155910f, 0.158700f, 0.161540f, 0.164400f, 0.167290f, 0.170180f, 0.173070f, 0.175930f,
        0.178760f, 0.181540f, 0.184260f, 0.186900f, 0.189450f, 0.191900f, 0.194230f, 0.196410f, 0.198430f, 0.200270f,
        0.201920f, 0.203350f, 0.204550f, 0.205500f, 0.206190f, 0.206610f, 0.206770f, 0.206690f, 0.206360f, 0.205780f,
        0.204980f, 0.203940f, 0.202680f, 0.201200f, 0.199510f, 0.197620f, 0.195560f, 0.193320f, 0.190940f, 0.188430f,
        0.185810f, 0.183080f, 0.180270f, 0.177400f, 0.174470f, 0.171500f, 0.168490f, 0.165430f, 0.162340f, 0.159220f,
        0.156070f, 0.152900f, 0.149710f, 0.146500f, 0.143280f, 0.140050f, 0.136830f, 0.133610f, 0.130400f, 0.127220f,
        0.124060f, 0.120930f, 0.117840f, 0.114800f, 0.111810f, 0.108880f, 0.106020f, 0.103240f, 0.100540f, 0.097933f,
        0.095428f, 0.093032f, 0.090753f, 0.088600f, 0.086578f, 0.084685f, 0.082915f, 0.081265f, 0.079729f, 0.078303f,
        0.076982f, 0.075761f, 0.074635f, 0.073600f, 0.072649f, 0.071773f, 0.070958f, 0.070194f, 0.069468f, 0.068769f,
        0.068084f, 0.067402f, 0.066712f, 0.066000f, 0.065258f, 0.064484f, 0.063678f, 0.062841f, 0.061974f, 0.061076f,
        0.060150f, 0.059195f, 0.058211f, 0.057200f, 0.056163f, 0.055105f, 0.054034f, 0.052955f, 0.051875f, 0.050802f,
        0.049741f, 0.048699f, 0.047683f, 0.046700f, 0.045756f, 0.044860f, 0.044021f, 0.043248f, 0.042550f, 0.041937f,
        0.041416f, 0.040996f, 0.040688f, 0.040500f, 0.040439f, 0.040505f, 0.040696f, 0.041011f, 0.041448f, 0.042005f,
        0.042680f, 0.043473f, 0.044380f, 0.045400f, 0.046527f, 0.047736f, 0.048997f, 0.050281f, 0.051558f, 0.052798f,
        0.053972f, 0.055050f, 0.056003f, 0.056800f, 0.057419f, 0.057862f, 0.058140f, 0.058260f, 0.058232f, 0.058066f,
        0.057770f, 0.057354f, 0.056828f, 0.056200f, 0.055481f, 0.054684f, 0.053824f, 0.052917f, 0.051977f, 0.051018f,
        0.050057f, 0.049106f, 0.048183f, 0.047300f, 0.046471f, 0.045696f, 0.044975f, 0.044305f, 0.043686f, 0.043116f,
        0.042594f, 0.042118f, 0.041687f, 0.041300f, 0.040955f, 0.040650f, 0.040385f, 0.040157f, 0.039965f, 0.039809f,
        0.039685f, 0.039594f, 0.039532f, 0.039500f, 0.039495f, 0.039515f, 0.039557f, 0.039620f, 0.039702f, 0.039799f,
        0.039910f, 0.040031f, 0.040162f, 0.040300f, 0.040443f, 0.040594f, 0.040758f, 0.040939f, 0.041140f, 0.041366f,
        0.041621f, 0.041909f, 0.042234f, 0.042600f, 0.043011f, 0.043471f, 0.043982f, 0.044549f, 0.045175f, 0.045862f,
        0.046615f, 0.047437f, 0.048331f, 0.049300f, 0.050347f, 0.051469f, 0.052663f, 0.053924f, 0.055249f, 0.056634f,
        0.058075f, 0.059569f, 0.061112f, 0.062700f, 0.064329f, 0.065994f, 0.067691f, 0.069413f, 0.071156f, 0.072914f,
        0.074683f, 0.076457f, 0.078231f, 0.080000f, 0.081761f, 0.083519f, 0.085282f, 0.087058f, 0.088853f, 0.090677f,
        0.092536f, 0.094437f, 0.096390f, 0.098400f, 0.100470f, 0.102570f, 0.104660f, 0.106700f, 0.108670f, 0.110520f,
        0.112210f, 0.113710f, 0.114990f, 0.116000f, 0.117660f, 0.119090f, 0.120530f, 0.121990f, 0.123460f, 0.124950f,
        0.126450f, 0.127970f, 0.129510f, 0.131060f, 0.132620f, 0.134200f, 0.135800f, 0.137410f, 0.139040f, 0.140690f,
        0.142350f, 0.144030f, 0.145720f, 0.147430f, 0.149160f, 0.150900f, 0.152660f, 0.154440f, 0.156230f, 0.158050f,
        0.159870f, 0.161720f, 0.163580f, 0.165460f, 0.167360f, 0.169270f, 0.171200f, 0.173150f, 0.175120f, 0.177100f,
        0.179110f, 0.181130f, 0.183160f, 0.185220f, 0.187290f, 0.189380f, 0.191490f, 0.193620f, 0.195760f, 0.197920f,
        0.200100f, 0.202300f, 0.204520f, 0.206750f, 0.209000f, 0.211270f, 0.213560f, 0.215870f, 0.218190f, 0.220530f,
        0.222890f, 0.225270f, 0.227660f, 0.230080f, 0.232510f, 0.234960f, 0.237430f, 0.239910f, 0.242410f, 0.244930f,
        0.247470f, 0.250030f, 0.252600f, 0.255190f, 0.257800f, 0.260420f, 0.263070f, 0.265730f, 0.268400f, 0.271100f,
        0.273810f, 0.276540f, 0.279280f});
    static constexpr Spectrum CES76({0.211990f, 0.215410f, 0.218880f, 0.222380f, 0.225930f, 0.229510f, 0.233130f,
        0.236800f, 0.240500f, 0.244240f, 0.248030f, 0.251850f, 0.255710f, 0.259610f, 0.263540f, 0.267520f, 0.271530f,
        0.275580f, 0.279670f, 0.283790f, 0.288600f, 0.291780f, 0.295760f, 0.300400f, 0.305560f, 0.311110f, 0.316890f,
        0.322770f, 0.328610f, 0.334270f, 0.339600f, 0.344500f, 0.348990f, 0.353140f, 0.357010f, 0.360660f, 0.364150f,
        0.367530f, 0.370880f, 0.374250f, 0.377700f, 0.381280f, 0.385000f, 0.388820f, 0.392750f, 0.396750f, 0.400830f,
        0.404960f, 0.409120f, 0.413310f, 0.417500f, 0.421680f, 0.425850f, 0.429980f, 0.434060f, 0.438080f, 0.442020f,
        0.445880f, 0.449640f, 0.453280f, 0.456800f, 0.460180f, 0.463400f, 0.466440f, 0.469290f, 0.471920f, 0.474320f,
        0.476460f, 0.478340f, 0.479920f, 0.481200f, 0.482150f, 0.482790f, 0.483130f, 0.483170f, 0.482930f, 0.482420f,
        0.481650f, 0.480630f, 0.479380f, 0.477900f, 0.476210f, 0.474320f, 0.472260f, 0.470030f, 0.467660f, 0.465170f,
        0.462570f, 0.459880f, 0.457120f, 0.454300f, 0.451440f, 0.448550f, 0.445610f, 0.442630f, 0.439590f, 0.436510f,
        0.433370f, 0.430180f, 0.426920f, 0.423600f, 0.420220f, 0.416760f, 0.413250f, 0.409670f, 0.406030f, 0.402320f,
        0.398560f, 0.394730f, 0.390850f, 0.386900f, 0.382900f, 0.378830f, 0.374690f, 0.370490f, 0.366210f, 0.361840f,
        0.357400f, 0.352860f, 0.348230f, 0.343500f, 0.338670f, 0.333740f, 0.328730f, 0.323640f, 0.318480f, 0.313260f,
        0.307980f, 0.302650f, 0.297290f, 0.291900f, 0.286490f, 0.281070f, 0.275680f, 0.270330f, 0.265030f, 0.259820f,
        0.254710f, 0.249730f, 0.244880f, 0.240200f, 0.235700f, 0.231370f, 0.227200f, 0.223180f, 0.219290f, 0.215540f,
        0.211890f, 0.208340f, 0.204880f, 0.201500f, 0.198180f, 0.194910f, 0.191670f, 0.188440f, 0.185210f, 0.181960f,
        0.178670f, 0.175320f, 0.171910f, 0.168400f, 0.164790f, 0.161080f, 0.157270f, 0.153350f, 0.149350f, 0.145250f,
        0.141070f, 0.136790f, 0.132440f, 0.128000f, 0.123490f, 0.118930f, 0.114330f, 0.109730f, 0.105160f, 0.100620f,
        0.096153f, 0.091778f, 0.087519f, 0.083400f, 0.079445f, 0.075672f, 0.072100f, 0.068746f, 0.065627f, 0.062763f,
        0.060170f, 0.057867f, 0.055871f, 0.054200f, 0.052870f, 0.051885f, 0.051249f, 0.050965f, 0.051035f, 0.051462f,
        0.052249f, 0.053399f, 0.054915f, 0.056800f, 0.059047f, 0.061612f, 0.064441f, 0.067483f, 0.070683f, 0.073989f,
        0.077347f, 0.080703f, 0.084005f, 0.087200f, 0.090241f, 0.093109f, 0.095795f, 0.098285f, 0.100570f, 0.102640f,
        0.104470f, 0.106070f, 0.107420f, 0.108500f, 0.109320f, 0.109890f, 0.110250f, 0.110430f, 0.110450f, 0.110360f,
        0.110180f, 0.109940f, 0.109670f, 0.109400f, 0.109160f, 0.108960f, 0.108790f, 0.108660f, 0.108570f, 0.108510f,
        0.108500f, 0.108520f, 0.108590f, 0.108700f, 0.108850f, 0.109050f, 0.109280f, 0.109550f, 0.109860f, 0.110210f,
        0.110580f, 0.110990f, 0.111430f, 0.111900f, 0.112400f, 0.112930f, 0.113500f, 0.114120f, 0.114810f, 0.115560f,
        0.116390f, 0.117300f, 0.118300f, 0.119400f, 0.120610f, 0.121920f, 0.123340f, 0.124850f, 0.126460f, 0.128160f,
        0.129940f, 0.131820f, 0.133770f, 0.135800f, 0.137900f, 0.140060f, 0.142270f, 0.144490f, 0.146720f, 0.148950f,
        0.151140f, 0.153300f, 0.155390f, 0.157400f, 0.159320f, 0.161130f, 0.162820f, 0.164380f, 0.165790f, 0.167050f,
        0.168130f, 0.169020f, 0.169720f, 0.170200f, 0.170470f, 0.170520f, 0.170390f, 0.170090f, 0.169640f, 0.169050f,
        0.168340f, 0.167540f, 0.166650f, 0.165700f, 0.164700f, 0.163640f, 0.162500f, 0.161270f, 0.159940f, 0.158490f,
        0.156910f, 0.155170f, 0.153280f, 0.151200f, 0.148940f, 0.146550f, 0.144090f, 0.141610f, 0.139170f, 0.136840f,
        0.134660f, 0.132710f, 0.131040f, 0.129700f, 0.127700f, 0.125980f, 0.124280f, 0.122600f, 0.120950f, 0.119310f,
        0.117690f, 0.116090f, 0.114500f, 0.112940f, 0.111400f, 0.109870f, 0.108360f, 0.106870f, 0.105400f, 0.103950f,
        0.102510f, 0.101100f, 0.099695f, 0.098312f, 0.096946f, 0.095596f, 0.094264f, 0.092948f, 0.091649f, 0.090366f,
        0.089100f, 0.087849f, 0.086614f, 0.085395f, 0.084192f, 0.083004f, 0.081831f, 0.080673f, 0.079531f, 0.078403f,
        0.077290f, 0.076191f, 0.075106f, 0.074036f, 0.072980f, 0.071938f, 0.070909f, 0.069895f, 0.068893f, 0.067905f,
        0.066930f, 0.065968f, 0.065019f, 0.064082f, 0.063159f, 0.062247f, 0.061348f, 0.060461f, 0.059586f, 0.058723f,
        0.057872f, 0.057032f, 0.056204f, 0.055387f, 0.054581f, 0.053786f, 0.053002f, 0.052229f, 0.051467f, 0.050715f,
        0.049974f, 0.049243f, 0.048522f, 0.047811f, 0.047110f, 0.046418f, 0.045737f, 0.045065f, 0.044402f, 0.043749f,
        0.043105f, 0.042469f, 0.041843f});
    static constexpr Spectrum CES77({0.028048f, 0.029636f, 0.030644f, 0.031716f, 0.033095f, 0.035374f, 0.038409f,
        0.041975f, 0.045803f, 0.049918f, 0.054221f, 0.058601f, 0.063164f, 0.067928f, 0.072997f, 0.078300f, 0.083754f,
        0.089247f, 0.094851f, 0.100429f, 0.105929f, 0.110653f, 0.115484f, 0.120416f, 0.125438f, 0.130460f, 0.135439f,
        0.140408f, 0.145471f, 0.150663f, 0.155909f, 0.161138f, 0.166309f, 0.171404f, 0.176431f, 0.181404f, 0.186392f,
        0.191397f, 0.196365f, 0.201203f, 0.205860f, 0.210364f, 0.214749f, 0.219038f, 0.223237f, 0.227350f, 0.231406f,
        0.235411f, 0.239248f, 0.242941f, 0.246491f, 0.249887f, 0.253170f, 0.256342f, 0.259408f, 0.262339f, 0.265047f,
        0.267509f, 0.269804f, 0.272095f, 0.274254f, 0.276187f, 0.277968f, 0.279676f, 0.281415f, 0.283126f, 0.284764f,
        0.286389f, 0.287955f, 0.289360f, 0.290615f, 0.291725f, 0.292714f, 0.293622f, 0.294353f, 0.294800f, 0.295042f,
        0.295204f, 0.295292f, 0.295278f, 0.295171f, 0.294927f, 0.294596f, 0.294248f, 0.293806f, 0.293265f, 0.292636f,
        0.291842f, 0.290888f, 0.289879f, 0.288751f, 0.287494f, 0.286152f, 0.284753f, 0.283292f, 0.281756f, 0.280152f,
        0.278521f, 0.276920f, 0.275235f, 0.273344f, 0.271420f, 0.269572f, 0.267697f, 0.265678f, 0.263461f, 0.261109f,
        0.258736f, 0.256338f, 0.253889f, 0.251406f, 0.248796f, 0.246147f, 0.243468f, 0.240740f, 0.238093f, 0.235597f,
        0.233172f, 0.230644f, 0.228009f, 0.225347f, 0.222715f, 0.220238f, 0.217823f, 0.215384f, 0.212966f, 0.210551f,
        0.208133f, 0.205734f, 0.203358f, 0.200979f, 0.198636f, 0.196339f, 0.194047f, 0.191716f, 0.189335f, 0.186942f,
        0.184556f, 0.182180f, 0.179816f, 0.177495f, 0.175235f, 0.172945f, 0.170624f, 0.168293f, 0.165978f, 0.163678f,
        0.161394f, 0.159147f, 0.156940f, 0.154781f, 0.152678f, 0.150625f, 0.148603f, 0.146599f, 0.144639f, 0.142720f,
        0.140838f, 0.138990f, 0.137160f, 0.135341f, 0.133527f, 0.131706f, 0.129894f, 0.128098f, 0.126321f, 0.124534f,
        0.122739f, 0.120943f, 0.119139f, 0.117338f, 0.115547f, 0.113774f, 0.112044f, 0.110347f, 0.108700f, 0.107081f,
        0.105515f, 0.103996f, 0.102540f, 0.101133f, 0.099788f, 0.098497f, 0.097258f, 0.096071f, 0.094913f, 0.093783f,
        0.092708f, 0.091654f, 0.090658f, 0.089671f, 0.088729f, 0.087813f, 0.086936f, 0.086091f, 0.085258f, 0.084467f,
        0.083713f, 0.082988f, 0.082331f, 0.081693f, 0.081119f, 0.080567f, 0.080059f, 0.079579f, 0.079136f, 0.078750f,
        0.078403f, 0.078109f, 0.077860f, 0.077643f, 0.077463f, 0.077294f, 0.077190f, 0.077107f, 0.077073f, 0.077070f,
        0.077102f, 0.077154f, 0.077216f, 0.077304f, 0.077422f, 0.077565f, 0.077766f, 0.077987f, 0.078251f, 0.078553f,
        0.078919f, 0.079302f, 0.079750f, 0.080234f, 0.080781f, 0.081388f, 0.082056f, 0.082779f, 0.083572f, 0.084460f,
        0.085425f, 0.086461f, 0.087624f, 0.088856f, 0.090208f, 0.091657f, 0.093221f, 0.094951f, 0.096814f, 0.098804f,
        0.100895f, 0.103096f, 0.105410f, 0.107842f, 0.110433f, 0.113158f, 0.116034f, 0.119055f, 0.122160f, 0.125408f,
        0.128794f, 0.132327f, 0.136042f, 0.139918f, 0.143925f, 0.148079f, 0.152382f, 0.156798f, 0.161328f, 0.166002f,
        0.170765f, 0.175649f, 0.180657f, 0.185741f, 0.190912f, 0.196131f, 0.201379f, 0.206639f, 0.211931f, 0.217258f,
        0.222593f, 0.227972f, 0.233377f, 0.238768f, 0.244108f, 0.249364f, 0.254606f, 0.259816f, 0.264961f, 0.269996f,
        0.274886f, 0.279626f, 0.284172f, 0.288517f, 0.292634f, 0.296597f, 0.300432f, 0.304092f, 0.307370f, 0.310359f,
        0.313562f, 0.316648f, 0.319661f, 0.322488f, 0.324876f, 0.327081f, 0.329079f, 0.330870f, 0.332473f, 0.333919f,
        0.335155f, 0.336239f, 0.337125f, 0.337817f, 0.338346f, 0.338703f, 0.338919f, 0.339003f, 0.338975f, 0.338805f,
        0.338483f, 0.338036f, 0.337425f, 0.336668f, 0.335771f, 0.334793f, 0.333760f, 0.332598f, 0.331333f, 0.329946f,
        0.328513f, 0.327074f, 0.325634f, 0.324248f, 0.322923f, 0.321682f, 0.320432f, 0.319149f, 0.317879f, 0.316795f,
        0.315878f, 0.314556f, 0.313161f, 0.311652f, 0.310112f, 0.308814f, 0.307519f, 0.306227f, 0.304938f, 0.303651f,
        0.302369f, 0.301089f, 0.299812f, 0.298538f, 0.297267f, 0.296000f, 0.294735f, 0.293474f, 0.292216f, 0.290961f,
        0.289710f, 0.288461f, 0.287216f, 0.285973f, 0.284735f, 0.283499f, 0.282267f, 0.281038f, 0.279812f, 0.278589f,
        0.277370f, 0.276154f, 0.274941f, 0.273732f, 0.272526f, 0.271323f, 0.270124f, 0.268928f, 0.267735f, 0.266546f,
        0.265360f, 0.264177f, 0.262998f, 0.261822f, 0.260627f, 0.259436f, 0.258250f, 0.257068f, 0.255891f, 0.254718f,
        0.253550f, 0.252386f, 0.251227f, 0.250073f, 0.248923f, 0.247778f, 0.246637f, 0.245501f, 0.244369f, 0.243242f,
        0.242120f, 0.241002f, 0.239890f});
    static constexpr Spectrum CES78({0.129470f, 0.131840f, 0.134240f, 0.136680f, 0.139150f, 0.141670f, 0.144220f,
        0.146810f, 0.149440f, 0.152100f, 0.154810f, 0.157560f, 0.160340f, 0.163160f, 0.166030f, 0.168930f, 0.171880f,
        0.174860f, 0.177890f, 0.180950f, 0.181200f, 0.185940f, 0.190320f, 0.194380f, 0.198170f, 0.201720f, 0.205090f,
        0.208310f, 0.211430f, 0.214500f, 0.217540f, 0.220620f, 0.223770f, 0.227040f, 0.230460f, 0.234080f, 0.237960f,
        0.242120f, 0.246610f, 0.251480f, 0.256760f, 0.262500f, 0.268670f, 0.275230f, 0.282150f, 0.289400f, 0.296950f,
        0.304760f, 0.312810f, 0.321040f, 0.329440f, 0.337970f, 0.346560f, 0.355150f, 0.363670f, 0.372060f, 0.380260f,
        0.388190f, 0.395800f, 0.403010f, 0.409770f, 0.416010f, 0.421730f, 0.426910f, 0.431560f, 0.435670f, 0.439220f,
        0.442230f, 0.444680f, 0.446570f, 0.447900f, 0.448660f, 0.448870f, 0.448570f, 0.447790f, 0.446540f, 0.444860f,
        0.442780f, 0.440330f, 0.437520f, 0.434390f, 0.430970f, 0.427280f, 0.423320f, 0.419130f, 0.414720f, 0.410110f,
        0.405320f, 0.400350f, 0.395250f, 0.390010f, 0.384660f, 0.379230f, 0.373740f, 0.368200f, 0.362650f, 0.357110f,
        0.351600f, 0.346150f, 0.340770f, 0.335500f, 0.330350f, 0.325330f, 0.320430f, 0.315650f, 0.310990f, 0.306450f,
        0.302030f, 0.297720f, 0.293520f, 0.289430f, 0.285450f, 0.281560f, 0.277740f, 0.273990f, 0.270270f, 0.266580f,
        0.262890f, 0.259190f, 0.255470f, 0.251700f, 0.247870f, 0.244000f, 0.240090f, 0.236140f, 0.232190f, 0.228230f,
        0.224280f, 0.220340f, 0.216440f, 0.212580f, 0.208770f, 0.205030f, 0.201350f, 0.197740f, 0.194210f, 0.190770f,
        0.187430f, 0.184190f, 0.181050f, 0.178030f, 0.175120f, 0.172330f, 0.169650f, 0.167070f, 0.164580f, 0.162180f,
        0.159850f, 0.157610f, 0.155420f, 0.153300f, 0.151230f, 0.149210f, 0.147200f, 0.145220f, 0.143230f, 0.141230f,
        0.139210f, 0.137140f, 0.135030f, 0.132850f, 0.130600f, 0.128280f, 0.125890f, 0.123460f, 0.120980f, 0.118470f,
        0.115930f, 0.113370f, 0.110800f, 0.108230f, 0.105660f, 0.103120f, 0.100600f, 0.098122f, 0.095695f, 0.093329f,
        0.091035f, 0.088824f, 0.086707f, 0.084694f, 0.082796f, 0.081011f, 0.079340f, 0.077780f, 0.076331f, 0.074990f,
        0.073758f, 0.072632f, 0.071611f, 0.070694f, 0.069880f, 0.069169f, 0.068561f, 0.068057f, 0.067656f, 0.067359f,
        0.067167f, 0.067080f, 0.067097f, 0.067219f, 0.067446f, 0.067775f, 0.068199f, 0.068715f, 0.069318f, 0.070003f,
        0.070766f, 0.071602f, 0.072506f, 0.073475f, 0.074501f, 0.075577f, 0.076692f, 0.077838f, 0.079003f, 0.080180f,
        0.081356f, 0.082524f, 0.083673f, 0.084794f, 0.085876f, 0.086912f, 0.087893f, 0.088812f, 0.089659f, 0.090427f,
        0.091107f, 0.091691f, 0.092171f, 0.092538f, 0.092788f, 0.092925f, 0.092958f, 0.092894f, 0.092743f, 0.092511f,
        0.092208f, 0.091841f, 0.091419f, 0.090950f, 0.090441f, 0.089903f, 0.089343f, 0.088771f, 0.088196f, 0.087627f,
        0.087072f, 0.086542f, 0.086044f, 0.085588f, 0.085182f, 0.084831f, 0.084539f, 0.084310f, 0.084148f, 0.084058f,
        0.084044f, 0.084109f, 0.084259f, 0.084496f, 0.084824f, 0.085242f, 0.085748f, 0.086339f, 0.087013f, 0.087767f,
        0.088600f, 0.089509f, 0.090491f, 0.091545f, 0.092668f, 0.093854f, 0.095098f, 0.096393f, 0.097736f, 0.099119f,
        0.100540f, 0.101990f, 0.103460f, 0.104950f, 0.106450f, 0.107970f, 0.109490f, 0.111020f, 0.112560f, 0.114100f,
        0.115640f, 0.117170f, 0.118710f, 0.120240f, 0.121770f, 0.123330f, 0.124940f, 0.126660f, 0.128500f, 0.130500f,
        0.132690f, 0.135120f, 0.137810f, 0.140790f, 0.144100f, 0.147730f, 0.151660f, 0.155880f, 0.160390f, 0.165170f,
        0.170210f, 0.175500f, 0.181020f, 0.186760f, 0.192720f, 0.198880f, 0.205220f, 0.211740f, 0.218420f, 0.225260f,
        0.232230f, 0.239340f, 0.246560f, 0.253880f, 0.261040f, 0.268630f, 0.276350f, 0.284210f, 0.292200f, 0.300320f,
        0.308570f, 0.316940f, 0.325440f, 0.334050f, 0.342770f, 0.351600f, 0.360530f, 0.369560f, 0.378680f, 0.387890f,
        0.397180f, 0.406550f, 0.415980f, 0.425470f, 0.435020f, 0.444620f, 0.454260f, 0.463940f, 0.473640f, 0.483360f,
        0.493090f, 0.502830f, 0.512570f, 0.522300f, 0.532010f, 0.541700f, 0.551350f, 0.560970f, 0.570540f, 0.580060f,
        0.589520f, 0.598910f, 0.608230f, 0.617470f, 0.626630f, 0.635700f, 0.644680f, 0.653550f, 0.662320f, 0.670980f,
        0.679520f, 0.687940f, 0.696240f, 0.704420f, 0.712470f, 0.720380f, 0.728160f, 0.735800f, 0.743310f, 0.750670f,
        0.757890f, 0.764970f, 0.771900f, 0.778690f, 0.785330f, 0.791820f, 0.798170f, 0.804380f, 0.810430f, 0.816350f,
        0.822120f, 0.827740f, 0.833230f, 0.838570f, 0.843770f, 0.848840f, 0.853770f, 0.858570f, 0.863230f, 0.867770f,
        0.872180f, 0.876460f, 0.880610f});
    static constexpr Spectrum CES79({0.567490f, 0.569710f, 0.571930f, 0.574150f, 0.576360f, 0.578570f, 0.580770f,
        0.582970f, 0.585170f, 0.587370f, 0.589560f, 0.591750f, 0.593930f, 0.596110f, 0.598290f, 0.600460f, 0.602630f,
        0.604800f, 0.606960f, 0.609120f, 0.606620f, 0.610970f, 0.614880f, 0.618370f, 0.621470f, 0.624220f, 0.626630f,
        0.628740f, 0.630570f, 0.632160f, 0.633530f, 0.634700f, 0.635710f, 0.636590f, 0.637360f, 0.638050f, 0.638680f,
        0.639300f, 0.639910f, 0.640560f, 0.641270f, 0.642060f, 0.642920f, 0.643830f, 0.644790f, 0.645760f, 0.646730f,
        0.647690f, 0.648610f, 0.649490f, 0.650310f, 0.651040f, 0.651680f, 0.652210f, 0.652640f, 0.652950f, 0.653130f,
        0.653180f, 0.653070f, 0.652810f, 0.652390f, 0.651800f, 0.651070f, 0.650230f, 0.649320f, 0.648370f, 0.647410f,
        0.646470f, 0.645600f, 0.644810f, 0.644150f, 0.643640f, 0.643280f, 0.643050f, 0.642960f, 0.642970f, 0.643100f,
        0.643320f, 0.643620f, 0.644000f, 0.644450f, 0.644940f, 0.645460f, 0.645970f, 0.646440f, 0.646840f, 0.647150f,
        0.647320f, 0.647320f, 0.647140f, 0.646730f, 0.646080f, 0.645180f, 0.644060f, 0.642720f, 0.641190f, 0.639470f,
        0.637580f, 0.635540f, 0.633360f, 0.631040f, 0.628620f, 0.626070f, 0.623420f, 0.620660f, 0.617780f, 0.614790f,
        0.611700f, 0.608490f, 0.605180f, 0.601760f, 0.598230f, 0.594610f, 0.590900f, 0.587130f, 0.583300f, 0.579420f,
        0.575520f, 0.571590f, 0.567660f, 0.563730f, 0.559820f, 0.555930f, 0.552070f, 0.548230f, 0.544430f, 0.540660f,
        0.536930f, 0.533240f, 0.529600f, 0.526000f, 0.522460f, 0.518960f, 0.515490f, 0.512050f, 0.508630f, 0.505220f,
        0.501800f, 0.498370f, 0.494930f, 0.491450f, 0.487940f, 0.484410f, 0.480850f, 0.477280f, 0.473700f, 0.470120f,
        0.466560f, 0.463010f, 0.459490f, 0.456010f, 0.452560f, 0.449160f, 0.445810f, 0.442520f, 0.439280f, 0.436090f,
        0.432980f, 0.429930f, 0.426950f, 0.424040f, 0.421210f, 0.418460f, 0.415780f, 0.413180f, 0.410660f, 0.408220f,
        0.405860f, 0.403570f, 0.401350f, 0.399220f, 0.397160f, 0.395170f, 0.393250f, 0.391400f, 0.389610f, 0.387880f,
        0.386200f, 0.384570f, 0.382990f, 0.381450f, 0.379940f, 0.378470f, 0.377030f, 0.375610f, 0.374200f, 0.372800f,
        0.371400f, 0.370000f, 0.368580f, 0.367150f, 0.365700f, 0.364240f, 0.362790f, 0.361350f, 0.359950f, 0.358600f,
        0.357310f, 0.356090f, 0.354970f, 0.353940f, 0.353040f, 0.352250f, 0.351580f, 0.351030f, 0.350600f, 0.350280f,
        0.350080f, 0.350000f, 0.350030f, 0.350170f, 0.350430f, 0.350800f, 0.351270f, 0.351830f, 0.352470f, 0.353200f,
        0.353990f, 0.354850f, 0.355760f, 0.356720f, 0.357730f, 0.358770f, 0.359840f, 0.360940f, 0.362070f, 0.363220f,
        0.364390f, 0.365560f, 0.366750f, 0.367940f, 0.369140f, 0.370340f, 0.371560f, 0.372800f, 0.374080f, 0.375390f,
        0.376740f, 0.378150f, 0.379610f, 0.381150f, 0.382760f, 0.384460f, 0.386280f, 0.388240f, 0.390350f, 0.392640f,
        0.395130f, 0.397840f, 0.400780f, 0.403980f, 0.407460f, 0.411220f, 0.415260f, 0.419570f, 0.424150f, 0.429010f,
        0.434140f, 0.439540f, 0.445210f, 0.451140f, 0.457340f, 0.463800f, 0.470500f, 0.477440f, 0.484600f, 0.491980f,
        0.499560f, 0.507330f, 0.515290f, 0.523420f, 0.531710f, 0.540130f, 0.548650f, 0.557230f, 0.565830f, 0.574430f,
        0.583000f, 0.591480f, 0.599870f, 0.608110f, 0.616180f, 0.624070f, 0.631760f, 0.639240f, 0.646500f, 0.653530f,
        0.660320f, 0.666850f, 0.673110f, 0.679100f, 0.684800f, 0.690220f, 0.695380f, 0.700290f, 0.704960f, 0.709400f,
        0.713630f, 0.717660f, 0.721500f, 0.725160f, 0.728670f, 0.732020f, 0.735240f, 0.738340f, 0.741340f, 0.744250f,
        0.747080f, 0.749840f, 0.752560f, 0.755250f, 0.757920f, 0.760580f, 0.763250f, 0.765950f, 0.768680f, 0.771470f,
        0.774330f, 0.777270f, 0.780300f, 0.783440f, 0.785090f, 0.787670f, 0.790220f, 0.792750f, 0.795250f, 0.797740f,
        0.800200f, 0.802640f, 0.805050f, 0.807450f, 0.809820f, 0.812170f, 0.814500f, 0.816800f, 0.819080f, 0.821340f,
        0.823580f, 0.825790f, 0.827990f, 0.830160f, 0.832310f, 0.834440f, 0.836540f, 0.838630f, 0.840690f, 0.842730f,
        0.844750f, 0.846750f, 0.848730f, 0.850690f, 0.852620f, 0.854540f, 0.856430f, 0.858310f, 0.860160f, 0.861990f,
        0.863810f, 0.865600f, 0.867370f, 0.869130f, 0.870860f, 0.872570f, 0.874270f, 0.875940f, 0.877600f, 0.879230f,
        0.880850f, 0.882450f, 0.884030f, 0.885590f, 0.887140f, 0.888660f, 0.890170f, 0.891660f, 0.893130f, 0.894580f,
        0.896020f, 0.897440f, 0.898840f, 0.900230f, 0.901590f, 0.902950f, 0.904280f, 0.905600f, 0.906900f, 0.908190f,
        0.909460f, 0.910710f, 0.911950f, 0.913170f, 0.914380f, 0.915570f, 0.916750f, 0.917910f, 0.919060f, 0.920190f,
        0.921310f, 0.922410f, 0.923500f});
    static constexpr Spectrum CES80({0.316840f, 0.322680f, 0.328320f, 0.333790f, 0.339120f, 0.344330f, 0.349450f,
        0.354530f, 0.359580f, 0.364630f, 0.369720f, 0.374870f, 0.380080f, 0.385340f, 0.390640f, 0.395990f, 0.401380f,
        0.406790f, 0.412220f, 0.417680f, 0.423140f, 0.428610f, 0.434070f, 0.439500f, 0.444900f, 0.450250f, 0.455530f,
        0.460740f, 0.465850f, 0.470860f, 0.475750f, 0.480510f, 0.485130f, 0.489630f, 0.493990f, 0.498230f, 0.502340f,
        0.506320f, 0.510170f, 0.513900f, 0.517510f, 0.521000f, 0.524360f, 0.527610f, 0.530750f, 0.533770f, 0.536690f,
        0.539500f, 0.542210f, 0.544820f, 0.547330f, 0.549750f, 0.552070f, 0.554300f, 0.556440f, 0.558470f, 0.560410f,
        0.562260f, 0.564000f, 0.565650f, 0.567200f, 0.568650f, 0.569990f, 0.571240f, 0.572370f, 0.573400f, 0.574320f,
        0.575120f, 0.575810f, 0.576380f, 0.576830f, 0.577160f, 0.577360f, 0.577450f, 0.577400f, 0.577230f, 0.576930f,
        0.576500f, 0.575930f, 0.575240f, 0.574410f, 0.573450f, 0.572340f, 0.571090f, 0.569700f, 0.568150f, 0.566450f,
        0.564590f, 0.562560f, 0.560370f, 0.558010f, 0.555470f, 0.552760f, 0.549880f, 0.546820f, 0.543600f, 0.540200f,
        0.536640f, 0.532920f, 0.529020f, 0.524970f, 0.520760f, 0.516390f, 0.511900f, 0.507280f, 0.502560f, 0.497740f,
        0.492840f, 0.487870f, 0.482860f, 0.477800f, 0.472720f, 0.467620f, 0.462510f, 0.457400f, 0.452290f, 0.447210f,
        0.442150f, 0.437120f, 0.432140f, 0.427200f, 0.422320f, 0.417500f, 0.412750f, 0.408060f, 0.403440f, 0.398880f,
        0.394400f, 0.389990f, 0.385660f, 0.381400f, 0.377220f, 0.373120f, 0.369100f, 0.365150f, 0.361270f, 0.357470f,
        0.353730f, 0.350050f, 0.346440f, 0.342890f, 0.339400f, 0.335970f, 0.332620f, 0.329330f, 0.326120f, 0.322990f,
        0.319950f, 0.317000f, 0.314140f, 0.311380f, 0.308720f, 0.306160f, 0.303690f, 0.301320f, 0.299030f, 0.296830f,
        0.294700f, 0.292660f, 0.290680f, 0.288770f, 0.286930f, 0.285150f, 0.283430f, 0.281790f, 0.280210f, 0.278690f,
        0.277250f, 0.275870f, 0.274570f, 0.273330f, 0.272160f, 0.271070f, 0.270040f, 0.269090f, 0.268200f, 0.267380f,
        0.266630f, 0.265950f, 0.265340f, 0.264790f, 0.264310f, 0.263900f, 0.263550f, 0.263280f, 0.263070f, 0.262930f,
        0.262850f, 0.262850f, 0.262910f, 0.263040f, 0.263240f, 0.263500f, 0.263820f, 0.264200f, 0.264620f, 0.265080f,
        0.265590f, 0.266130f, 0.266690f, 0.267280f, 0.267890f, 0.268520f, 0.269160f, 0.269830f, 0.270520f, 0.271240f,
        0.271980f, 0.272750f, 0.273550f, 0.274370f, 0.275220f, 0.276110f, 0.277010f, 0.277940f, 0.278900f, 0.279870f,
        0.280860f, 0.281860f, 0.282870f, 0.283900f, 0.284940f, 0.285980f, 0.287040f, 0.288120f, 0.289210f, 0.290330f,
        0.291480f, 0.292650f, 0.293850f, 0.295080f, 0.296350f, 0.297650f, 0.298990f, 0.300370f, 0.301780f, 0.303230f,
        0.304720f, 0.306230f, 0.307790f, 0.309380f, 0.311010f, 0.312670f, 0.314390f, 0.316150f, 0.317960f, 0.319840f,
        0.321780f, 0.323780f, 0.325860f, 0.328020f, 0.330260f, 0.332580f, 0.334980f, 0.337470f, 0.340040f, 0.342700f,
        0.345440f, 0.348260f, 0.351170f, 0.354160f, 0.357240f, 0.360410f, 0.363670f, 0.367020f, 0.370480f, 0.374050f,
        0.377720f, 0.381510f, 0.385410f, 0.389430f, 0.393580f, 0.397840f, 0.402230f, 0.406720f, 0.411330f, 0.416040f,
        0.420850f, 0.425760f, 0.430770f, 0.435860f, 0.441040f, 0.446300f, 0.451650f, 0.457080f, 0.462600f, 0.468190f,
        0.473860f, 0.479610f, 0.485440f, 0.491340f, 0.497320f, 0.503370f, 0.509490f, 0.515690f, 0.521960f, 0.528300f,
        0.534710f, 0.541190f, 0.547730f, 0.554350f, 0.561030f, 0.567760f, 0.574540f, 0.581340f, 0.588150f, 0.594960f,
        0.601750f, 0.608520f, 0.615240f, 0.621910f, 0.628510f, 0.635040f, 0.641490f, 0.647880f, 0.654190f, 0.660430f,
        0.666590f, 0.672680f, 0.678690f, 0.684630f, 0.690490f, 0.696260f, 0.701940f, 0.707520f, 0.712990f, 0.718350f,
        0.723570f, 0.728670f, 0.733620f, 0.738430f, 0.743080f, 0.747590f, 0.751950f, 0.756180f, 0.760290f, 0.764280f,
        0.768170f, 0.771960f, 0.775650f, 0.779270f, 0.782810f, 0.786270f, 0.789660f, 0.792950f, 0.796150f, 0.799250f,
        0.802260f, 0.805150f, 0.807940f, 0.810610f, 0.813160f, 0.815610f, 0.817950f, 0.820220f, 0.822410f, 0.824530f,
        0.826610f, 0.828650f, 0.830670f, 0.832670f, 0.834660f, 0.836630f, 0.838540f, 0.840380f, 0.842110f, 0.843720f,
        0.845180f, 0.846460f, 0.847540f, 0.848390f, 0.849000f, 0.849390f, 0.849600f, 0.849640f, 0.849570f, 0.849390f,
        0.849160f, 0.848890f, 0.848630f, 0.848390f, 0.848210f, 0.848090f, 0.848020f, 0.847990f, 0.848000f, 0.848040f,
        0.848100f, 0.848190f, 0.848290f, 0.848390f, 0.848490f, 0.848590f, 0.848680f, 0.848740f, 0.848780f, 0.848790f,
        0.848760f, 0.848690f, 0.848570f});
    static constexpr Spectrum CES81({0.120670f, 0.131180f, 0.141180f, 0.152590f, 0.164170f, 0.176020f, 0.187220f,
        0.199450f, 0.211100f, 0.222000f, 0.232980f, 0.242680f, 0.251120f, 0.259720f, 0.267250f, 0.274470f, 0.281400f,
        0.288430f, 0.294010f, 0.299480f, 0.304420f, 0.308950f, 0.313560f, 0.317470f, 0.320710f, 0.323570f, 0.326750f,
        0.329860f, 0.332420f, 0.335020f, 0.336820f, 0.338770f, 0.340590f, 0.341930f, 0.343520f, 0.345010f, 0.346320f,
        0.346650f, 0.346590f, 0.347020f, 0.347630f, 0.348340f, 0.348400f, 0.348590f, 0.349460f, 0.350690f, 0.352120f,
        0.352870f, 0.353690f, 0.354780f, 0.355000f, 0.354680f, 0.354280f, 0.354340f, 0.354210f, 0.352710f, 0.351000f,
        0.349390f, 0.348480f, 0.347200f, 0.345190f, 0.344070f, 0.342980f, 0.341420f, 0.339610f, 0.338140f, 0.337690f,
        0.337120f, 0.336220f, 0.334680f, 0.333350f, 0.332490f, 0.330730f, 0.328930f, 0.327460f, 0.325850f, 0.323600f,
        0.321110f, 0.318460f, 0.315730f, 0.313430f, 0.311260f, 0.309130f, 0.307180f, 0.305510f, 0.302960f, 0.301300f,
        0.299910f, 0.297990f, 0.295800f, 0.293430f, 0.290730f, 0.287910f, 0.285030f, 0.282340f, 0.279070f, 0.276630f,
        0.272930f, 0.268700f, 0.265050f, 0.261790f, 0.258710f, 0.255830f, 0.252790f, 0.249640f, 0.246440f, 0.243440f,
        0.239990f, 0.236830f, 0.233390f, 0.230070f, 0.227110f, 0.224360f, 0.221760f, 0.219180f, 0.216520f, 0.213820f,
        0.210880f, 0.208180f, 0.205540f, 0.203560f, 0.201230f, 0.198760f, 0.196060f, 0.193030f, 0.189940f, 0.186950f,
        0.184030f, 0.181060f, 0.177910f, 0.174810f, 0.171830f, 0.169010f, 0.165340f, 0.161710f, 0.158380f, 0.155090f,
        0.151550f, 0.147710f, 0.144160f, 0.140460f, 0.136570f, 0.132730f, 0.129020f, 0.125990f, 0.123150f, 0.120250f,
        0.117300f, 0.114570f, 0.112480f, 0.110690f, 0.109470f, 0.108550f, 0.107620f, 0.106900f, 0.106180f, 0.105350f,
        0.104460f, 0.104080f, 0.103800f, 0.103100f, 0.102510f, 0.102070f, 0.101670f, 0.100970f, 0.099973f, 0.099195f,
        0.098641f, 0.098719f, 0.098215f, 0.097486f, 0.096614f, 0.095262f, 0.093774f, 0.092362f, 0.091213f, 0.089899f,
        0.088341f, 0.086642f, 0.084249f, 0.082113f, 0.080282f, 0.078914f, 0.077949f, 0.076426f, 0.074395f, 0.072563f,
        0.071100f, 0.069971f, 0.068740f, 0.067994f, 0.067429f, 0.067089f, 0.066726f, 0.066693f, 0.067378f, 0.068152f,
        0.068785f, 0.069414f, 0.070249f, 0.071705f, 0.073490f, 0.075620f, 0.077706f, 0.079242f, 0.081241f, 0.083152f,
        0.084934f, 0.086570f, 0.088421f, 0.090241f, 0.091636f, 0.092684f, 0.093529f, 0.094517f, 0.096300f, 0.097743f,
        0.098657f, 0.099700f, 0.101020f, 0.102330f, 0.103420f, 0.104230f, 0.104550f, 0.105120f, 0.105170f, 0.104520f,
        0.103220f, 0.101900f, 0.101140f, 0.100100f, 0.098970f, 0.097916f, 0.097323f, 0.096479f, 0.095355f, 0.094920f,
        0.095497f, 0.096479f, 0.098077f, 0.099781f, 0.101520f, 0.104490f, 0.107120f, 0.109660f, 0.113680f, 0.116940f,
        0.120710f, 0.124770f, 0.128830f, 0.132970f, 0.137500f, 0.144020f, 0.149500f, 0.155350f, 0.161670f, 0.167590f,
        0.174540f, 0.181700f, 0.188800f, 0.196160f, 0.204360f, 0.212650f, 0.219980f, 0.228610f, 0.238980f, 0.249540f,
        0.260170f, 0.271000f, 0.281850f, 0.291730f, 0.302970f, 0.314030f, 0.324840f, 0.335000f, 0.343450f, 0.350870f,
        0.358250f, 0.365260f, 0.371920f, 0.378350f, 0.386170f, 0.392150f, 0.396790f, 0.400290f, 0.404030f, 0.408800f,
        0.413240f, 0.417090f, 0.420990f, 0.424720f, 0.427290f, 0.428390f, 0.429690f, 0.432540f, 0.436000f, 0.438770f,
        0.440960f, 0.442970f, 0.444770f, 0.446740f, 0.448850f, 0.450910f, 0.453060f, 0.454990f, 0.455910f, 0.456710f,
        0.458110f, 0.459470f, 0.460260f, 0.460760f, 0.461060f, 0.461630f, 0.462620f, 0.463320f, 0.461650f, 0.463520f,
        0.469320f, 0.470430f, 0.465310f, 0.465310f, 0.468510f, 0.469240f, 0.469980f, 0.470710f, 0.471450f, 0.472180f,
        0.472920f, 0.473660f, 0.474390f, 0.475130f, 0.475870f, 0.476600f, 0.477340f, 0.478070f, 0.478810f, 0.479550f,
        0.480280f, 0.481020f, 0.481760f, 0.482500f, 0.483230f, 0.483970f, 0.484710f, 0.485450f, 0.486180f, 0.486920f,
        0.487660f, 0.488400f, 0.489130f, 0.489870f, 0.490610f, 0.491350f, 0.492080f, 0.492820f, 0.493560f, 0.494300f,
        0.495040f, 0.495770f, 0.496510f, 0.497250f, 0.497990f, 0.498730f, 0.499460f, 0.500200f, 0.500940f, 0.501680f,
        0.502420f, 0.503150f, 0.503890f, 0.504630f, 0.505370f, 0.506110f, 0.506840f, 0.507580f, 0.508320f, 0.509060f,
        0.509800f, 0.510530f, 0.511270f, 0.512010f, 0.512750f, 0.513480f, 0.514220f, 0.514960f, 0.515700f, 0.516430f,
        0.517170f, 0.517910f, 0.518650f, 0.519380f, 0.520120f, 0.520860f, 0.521590f, 0.522330f, 0.523070f, 0.523800f,
        0.524540f, 0.525280f, 0.526010f});
    static constexpr Spectrum CES82({0.707720f, 0.709190f, 0.710650f, 0.712120f, 0.713570f, 0.715030f, 0.716470f,
        0.717920f, 0.719360f, 0.720790f, 0.722220f, 0.723640f, 0.725070f, 0.726480f, 0.727890f, 0.729300f, 0.730700f,
        0.732100f, 0.733500f, 0.734880f, 0.736500f, 0.737520f, 0.738800f, 0.740300f, 0.741970f, 0.743770f, 0.745650f,
        0.747580f, 0.749510f, 0.751400f, 0.753200f, 0.754880f, 0.756440f, 0.757910f, 0.759280f, 0.760590f, 0.761850f,
        0.763060f, 0.764240f, 0.765420f, 0.766600f, 0.767800f, 0.769000f, 0.770210f, 0.771410f, 0.772590f, 0.773740f,
        0.774860f, 0.775930f, 0.776950f, 0.777900f, 0.778780f, 0.779600f, 0.780350f, 0.781060f, 0.781710f, 0.782330f,
        0.782900f, 0.783460f, 0.783990f, 0.784500f, 0.785000f, 0.785490f, 0.785960f, 0.786410f, 0.786820f, 0.787190f,
        0.787520f, 0.787810f, 0.788030f, 0.788200f, 0.788300f, 0.788340f, 0.788310f, 0.788210f, 0.788040f, 0.787810f,
        0.787510f, 0.787140f, 0.786700f, 0.786200f, 0.785630f, 0.784980f, 0.784270f, 0.783480f, 0.782610f, 0.781670f,
        0.780650f, 0.779550f, 0.778370f, 0.777100f, 0.775750f, 0.774320f, 0.772820f, 0.771260f, 0.769640f, 0.767980f,
        0.766270f, 0.764540f, 0.762780f, 0.761000f, 0.759210f, 0.757410f, 0.755590f, 0.753740f, 0.751860f, 0.749950f,
        0.747990f, 0.745980f, 0.743920f, 0.741800f, 0.739620f, 0.737380f, 0.735090f, 0.732780f, 0.730440f, 0.728090f,
        0.725740f, 0.723400f, 0.721090f, 0.718800f, 0.716550f, 0.714350f, 0.712180f, 0.710050f, 0.707960f, 0.705890f,
        0.703860f, 0.701850f, 0.699860f, 0.697900f, 0.695960f, 0.694030f, 0.692120f, 0.690230f, 0.688340f, 0.686470f,
        0.684590f, 0.682730f, 0.680860f, 0.679000f, 0.677130f, 0.675270f, 0.673430f, 0.671600f, 0.669800f, 0.668040f,
        0.666320f, 0.664650f, 0.663040f, 0.661500f, 0.660030f, 0.658630f, 0.657310f, 0.656040f, 0.654840f, 0.653710f,
        0.652630f, 0.651600f, 0.650620f, 0.649700f, 0.648820f, 0.648000f, 0.647230f, 0.646520f, 0.645870f, 0.645300f,
        0.644800f, 0.644380f, 0.644050f, 0.643800f, 0.643650f, 0.643580f, 0.643590f, 0.643670f, 0.643810f, 0.644000f,
        0.644230f, 0.644500f, 0.644790f, 0.645100f, 0.645420f, 0.645750f, 0.646080f, 0.646410f, 0.646740f, 0.647070f,
        0.647390f, 0.647710f, 0.648010f, 0.648300f, 0.648570f, 0.648840f, 0.649090f, 0.649330f, 0.649560f, 0.649790f,
        0.650020f, 0.650240f, 0.650470f, 0.650700f, 0.650940f, 0.651180f, 0.651440f, 0.651710f, 0.652000f, 0.652310f,
        0.652640f, 0.653000f, 0.653390f, 0.653800f, 0.654250f, 0.654720f, 0.655220f, 0.655750f, 0.656290f, 0.656840f,
        0.657400f, 0.657970f, 0.658530f, 0.659100f, 0.659660f, 0.660210f, 0.660740f, 0.661270f, 0.661770f, 0.662250f,
        0.662710f, 0.663140f, 0.663530f, 0.663900f, 0.664230f, 0.664530f, 0.664800f, 0.665050f, 0.665280f, 0.665500f,
        0.665700f, 0.665900f, 0.666100f, 0.666300f, 0.666510f, 0.666720f, 0.666940f, 0.667160f, 0.667370f, 0.667590f,
        0.667800f, 0.668010f, 0.668210f, 0.668400f, 0.668580f, 0.668780f, 0.669010f, 0.669290f, 0.669660f, 0.670130f,
        0.670730f, 0.671470f, 0.672390f, 0.673500f, 0.674820f, 0.676350f, 0.678070f, 0.679960f, 0.682010f, 0.684210f,
        0.686540f, 0.689000f, 0.691550f, 0.694200f, 0.696930f, 0.699740f, 0.702630f, 0.705620f, 0.708710f, 0.711900f,
        0.715200f, 0.718610f, 0.722140f, 0.725800f, 0.729580f, 0.733490f, 0.737510f, 0.741650f, 0.745880f, 0.750220f,
        0.754640f, 0.759150f, 0.763740f, 0.768400f, 0.773120f, 0.777880f, 0.782650f, 0.787410f, 0.792130f, 0.796780f,
        0.801340f, 0.805780f, 0.810080f, 0.814200f, 0.818130f, 0.821900f, 0.825540f, 0.829070f, 0.832530f, 0.835950f,
        0.839370f, 0.842810f, 0.846310f, 0.849900f, 0.853590f, 0.857340f, 0.861080f, 0.864740f, 0.868260f, 0.871570f,
        0.874600f, 0.877300f, 0.879590f, 0.881400f, 0.884130f, 0.886440f, 0.888720f, 0.890950f, 0.893150f, 0.895300f,
        0.897420f, 0.899500f, 0.901540f, 0.903550f, 0.905520f, 0.907450f, 0.909350f, 0.911210f, 0.913040f, 0.914830f,
        0.916590f, 0.918320f, 0.920010f, 0.921670f, 0.923300f, 0.924900f, 0.926470f, 0.928010f, 0.929520f, 0.930990f,
        0.932440f, 0.933870f, 0.935260f, 0.936630f, 0.937970f, 0.939280f, 0.940570f, 0.941830f, 0.943070f, 0.944280f,
        0.945470f, 0.946630f, 0.947770f, 0.948890f, 0.949980f, 0.951050f, 0.952100f, 0.953130f, 0.954140f, 0.955130f,
        0.956090f, 0.957040f, 0.957970f, 0.958880f, 0.959770f, 0.960640f, 0.961490f, 0.962330f, 0.963140f, 0.963950f,
        0.964730f, 0.965500f, 0.966250f, 0.966980f, 0.967700f, 0.968410f, 0.969100f, 0.969770f, 0.970440f, 0.971080f,
        0.971720f, 0.972340f, 0.972940f, 0.973540f, 0.974120f, 0.974690f, 0.975240f, 0.975790f, 0.976320f, 0.976840f,
        0.977350f, 0.977850f, 0.978340f});
    static constexpr Spectrum CES83({0.563490f, 0.565110f, 0.566720f, 0.568330f, 0.569940f, 0.571540f, 0.573150f,
        0.574750f, 0.576360f, 0.577960f, 0.579560f, 0.581150f, 0.582750f, 0.584350f, 0.585940f, 0.587530f, 0.589120f,
        0.590700f, 0.592290f, 0.593870f, 0.595700f, 0.596890f, 0.598380f, 0.600090f, 0.601980f, 0.603960f, 0.605980f,
        0.607980f, 0.609890f, 0.611650f, 0.613200f, 0.614490f, 0.615520f, 0.616330f, 0.616950f, 0.617390f, 0.617690f,
        0.617880f, 0.617970f, 0.618010f, 0.618000f, 0.617980f, 0.617940f, 0.617870f, 0.617770f, 0.617620f, 0.617420f,
        0.617150f, 0.616820f, 0.616400f, 0.615900f, 0.615300f, 0.614620f, 0.613860f, 0.613030f, 0.612130f, 0.611190f,
        0.610210f, 0.609190f, 0.608150f, 0.607100f, 0.606040f, 0.604960f, 0.603860f, 0.602710f, 0.601520f, 0.600270f,
        0.598940f, 0.597520f, 0.596020f, 0.594400f, 0.592670f, 0.590830f, 0.588900f, 0.586870f, 0.584780f, 0.582620f,
        0.580420f, 0.578170f, 0.575890f, 0.573600f, 0.571300f, 0.568970f, 0.566610f, 0.564200f, 0.561730f, 0.559180f,
        0.556530f, 0.553780f, 0.550910f, 0.547900f, 0.544750f, 0.541470f, 0.538090f, 0.534630f, 0.531110f, 0.527570f,
        0.524030f, 0.520500f, 0.517020f, 0.513600f, 0.510270f, 0.507020f, 0.503850f, 0.500730f, 0.497660f, 0.494630f,
        0.491630f, 0.488650f, 0.485680f, 0.482700f, 0.479710f, 0.476710f, 0.473710f, 0.470710f, 0.467710f, 0.464730f,
        0.461760f, 0.458810f, 0.455890f, 0.453000f, 0.450140f, 0.447330f, 0.444560f, 0.441830f, 0.439160f, 0.436540f,
        0.433980f, 0.431490f, 0.429060f, 0.426700f, 0.424410f, 0.422200f, 0.420050f, 0.417960f, 0.415920f, 0.413940f,
        0.412000f, 0.410100f, 0.408230f, 0.406400f, 0.404590f, 0.402820f, 0.401070f, 0.399360f, 0.397690f, 0.396050f,
        0.394470f, 0.392930f, 0.391440f, 0.390000f, 0.388620f, 0.387300f, 0.386050f, 0.384880f, 0.383780f, 0.382770f,
        0.381850f, 0.381030f, 0.380310f, 0.379700f, 0.379200f, 0.378810f, 0.378510f, 0.378300f, 0.378170f, 0.378100f,
        0.378090f, 0.378120f, 0.378200f, 0.378300f, 0.378420f, 0.378580f, 0.378760f, 0.379000f, 0.379280f, 0.379630f,
        0.380050f, 0.380540f, 0.381120f, 0.381800f, 0.382580f, 0.383440f, 0.384390f, 0.385390f, 0.386430f, 0.387510f,
        0.388600f, 0.389680f, 0.390760f, 0.391800f, 0.392800f, 0.393760f, 0.394670f, 0.395550f, 0.396390f, 0.397200f,
        0.397970f, 0.398710f, 0.399420f, 0.400100f, 0.400760f, 0.401390f, 0.401990f, 0.402550f, 0.403090f, 0.403590f,
        0.404050f, 0.404480f, 0.404860f, 0.405200f, 0.405500f, 0.405750f, 0.405970f, 0.406150f, 0.406300f, 0.406420f,
        0.406520f, 0.406600f, 0.406660f, 0.406700f, 0.406730f, 0.406750f, 0.406760f, 0.406750f, 0.406720f, 0.406660f,
        0.406590f, 0.406490f, 0.406360f, 0.406200f, 0.406010f, 0.405790f, 0.405550f, 0.405280f, 0.404990f, 0.404680f,
        0.404350f, 0.404010f, 0.403660f, 0.403300f, 0.402930f, 0.402580f, 0.402260f, 0.401980f, 0.401770f, 0.401650f,
        0.401630f, 0.401740f, 0.401990f, 0.402400f, 0.402990f, 0.403750f, 0.404680f, 0.405770f, 0.407030f, 0.408450f,
        0.410020f, 0.411730f, 0.413600f, 0.415600f, 0.417740f, 0.420020f, 0.422450f, 0.425040f, 0.427790f, 0.430710f,
        0.433800f, 0.437080f, 0.440540f, 0.444200f, 0.448060f, 0.452110f, 0.456330f, 0.460720f, 0.465260f, 0.469940f,
        0.474750f, 0.479670f, 0.484690f, 0.489800f, 0.494990f, 0.500270f, 0.505630f, 0.511080f, 0.516640f, 0.522300f,
        0.528070f, 0.533950f, 0.539960f, 0.546100f, 0.552360f, 0.558710f, 0.565100f, 0.571480f, 0.577810f, 0.584030f,
        0.590110f, 0.596000f, 0.601640f, 0.607000f, 0.612040f, 0.616820f, 0.621400f, 0.625840f, 0.630230f, 0.634620f,
        0.639090f, 0.643690f, 0.648510f, 0.653600f, 0.659010f, 0.664630f, 0.670360f, 0.676070f, 0.681620f, 0.686910f,
        0.691810f, 0.696180f, 0.699920f, 0.702900f, 0.707460f, 0.711360f, 0.715230f, 0.719060f, 0.722870f, 0.726640f,
        0.730390f, 0.734090f, 0.737770f, 0.741410f, 0.745030f, 0.748600f, 0.752150f, 0.755660f, 0.759130f, 0.762570f,
        0.765980f, 0.769360f, 0.772700f, 0.776000f, 0.779270f, 0.782510f, 0.785710f, 0.788880f, 0.792020f, 0.795110f,
        0.798180f, 0.801210f, 0.804210f, 0.807170f, 0.810100f, 0.812990f, 0.815850f, 0.818670f, 0.821470f, 0.824220f,
        0.826950f, 0.829640f, 0.832300f, 0.834920f, 0.837510f, 0.840070f, 0.842590f, 0.845090f, 0.847550f, 0.849980f,
        0.852370f, 0.854740f, 0.857070f, 0.859370f, 0.861640f, 0.863880f, 0.866090f, 0.868270f, 0.870420f, 0.872540f,
        0.874630f, 0.876690f, 0.878720f, 0.880720f, 0.882690f, 0.884640f, 0.886560f, 0.888440f, 0.890310f, 0.892140f,
        0.893950f, 0.895730f, 0.897480f, 0.899210f, 0.900910f, 0.902590f, 0.904240f, 0.905860f, 0.907460f, 0.909040f,
        0.910590f, 0.912120f, 0.913630f});
    static constexpr Spectrum CES84({0.033754f, 0.034106f, 0.034463f, 0.034823f, 0.035186f, 0.035553f, 0.035924f,
        0.036299f, 0.036677f, 0.037059f, 0.037445f, 0.037835f, 0.038229f, 0.038626f, 0.039028f, 0.039434f, 0.039843f,
        0.040257f, 0.040675f, 0.041097f, 0.041383f, 0.041884f, 0.042374f, 0.042854f, 0.043326f, 0.043791f, 0.044249f,
        0.044703f, 0.045153f, 0.045601f, 0.046047f, 0.046494f, 0.046942f, 0.047392f, 0.047846f, 0.048304f, 0.048769f,
        0.049241f, 0.049721f, 0.050211f, 0.050712f, 0.051224f, 0.051746f, 0.052275f, 0.052808f, 0.053344f, 0.053880f,
        0.054413f, 0.054941f, 0.055461f, 0.055971f, 0.056469f, 0.056952f, 0.057419f, 0.057867f, 0.058293f, 0.058697f,
        0.059076f, 0.059427f, 0.059750f, 0.060040f, 0.060298f, 0.060525f, 0.060724f, 0.060896f, 0.061045f, 0.061173f,
        0.061283f, 0.061378f, 0.061459f, 0.061529f, 0.061590f, 0.061641f, 0.061680f, 0.061704f, 0.061712f, 0.061701f,
        0.061669f, 0.061615f, 0.061536f, 0.061430f, 0.061295f, 0.061135f, 0.060953f, 0.060752f, 0.060536f, 0.060308f,
        0.060072f, 0.059830f, 0.059587f, 0.059346f, 0.059110f, 0.058885f, 0.058677f, 0.058491f, 0.058333f, 0.058208f,
        0.058122f, 0.058081f, 0.058090f, 0.058155f, 0.058279f, 0.058458f, 0.058685f, 0.058953f, 0.059255f, 0.059583f,
        0.059933f, 0.060295f, 0.060664f, 0.061033f, 0.061394f, 0.061744f, 0.062076f, 0.062387f, 0.062670f, 0.062922f,
        0.063137f, 0.063311f, 0.063438f, 0.063514f, 0.063535f, 0.063504f, 0.063426f, 0.063305f, 0.063146f, 0.062952f,
        0.062729f, 0.062480f, 0.062211f, 0.061926f, 0.061627f, 0.061315f, 0.060986f, 0.060639f, 0.060272f, 0.059882f,
        0.059468f, 0.059027f, 0.058556f, 0.058055f, 0.057522f, 0.056958f, 0.056367f, 0.055752f, 0.055115f, 0.054459f,
        0.053787f, 0.053102f, 0.052407f, 0.051704f, 0.050997f, 0.050286f, 0.049572f, 0.048856f, 0.048138f, 0.047420f,
        0.046702f, 0.045986f, 0.045271f, 0.044559f, 0.043850f, 0.043146f, 0.042447f, 0.041756f, 0.041072f, 0.040397f,
        0.039732f, 0.039079f, 0.038438f, 0.037810f, 0.037198f, 0.036603f, 0.036029f, 0.035479f, 0.034956f, 0.034462f,
        0.034002f, 0.033577f, 0.033192f, 0.032848f, 0.032550f, 0.032294f, 0.032079f, 0.031904f, 0.031765f, 0.031662f,
        0.031591f, 0.031552f, 0.031542f, 0.031558f, 0.031601f, 0.031674f, 0.031782f, 0.031929f, 0.032122f, 0.032365f,
        0.032664f, 0.033022f, 0.033446f, 0.033940f, 0.034508f, 0.035145f, 0.035847f, 0.036607f, 0.037419f, 0.038279f,
        0.039181f, 0.040118f, 0.041086f, 0.042078f, 0.043089f, 0.044110f, 0.045135f, 0.046154f, 0.047160f, 0.048144f,
        0.049099f, 0.050016f, 0.050887f, 0.051704f, 0.052461f, 0.053161f, 0.053809f, 0.054408f, 0.054964f, 0.055482f,
        0.055966f, 0.056420f, 0.056851f, 0.057261f, 0.057656f, 0.058035f, 0.058398f, 0.058745f, 0.059074f, 0.059385f,
        0.059677f, 0.059950f, 0.060204f, 0.060437f, 0.060650f, 0.060843f, 0.061017f, 0.061175f, 0.061316f, 0.061443f,
        0.061556f, 0.061657f, 0.061747f, 0.061827f, 0.061898f, 0.061963f, 0.062024f, 0.062083f, 0.062142f, 0.062204f,
        0.062271f, 0.062344f, 0.062427f, 0.062521f, 0.062628f, 0.062747f, 0.062876f, 0.063014f, 0.063158f, 0.063306f,
        0.063457f, 0.063610f, 0.063761f, 0.063911f, 0.064056f, 0.064195f, 0.064329f, 0.064455f, 0.064573f, 0.064681f,
        0.064779f, 0.064866f, 0.064941f, 0.065002f, 0.065050f, 0.065084f, 0.065105f, 0.065116f, 0.065117f, 0.065108f,
        0.065091f, 0.065067f, 0.065037f, 0.065002f, 0.064963f, 0.064921f, 0.064875f, 0.064826f, 0.064776f, 0.064724f,
        0.064670f, 0.064616f, 0.064561f, 0.064506f, 0.064452f, 0.064397f, 0.064342f, 0.064287f, 0.064230f, 0.064171f,
        0.064111f, 0.064047f, 0.063981f, 0.063911f, 0.063837f, 0.063758f, 0.063675f, 0.063586f, 0.063492f, 0.063391f,
        0.063284f, 0.063170f, 0.063048f, 0.062918f, 0.062921f, 0.062833f, 0.062746f, 0.062659f, 0.062572f, 0.062485f,
        0.062398f, 0.062311f, 0.062224f, 0.062138f, 0.062051f, 0.061965f, 0.061879f, 0.061793f, 0.061707f, 0.061621f,
        0.061535f, 0.061450f, 0.061364f, 0.061279f, 0.061193f, 0.061108f, 0.061023f, 0.060938f, 0.060853f, 0.060769f,
        0.060684f, 0.060599f, 0.060515f, 0.060431f, 0.060347f, 0.060262f, 0.060178f, 0.060095f, 0.060011f, 0.059927f,
        0.059844f, 0.059760f, 0.059677f, 0.059594f, 0.059511f, 0.059428f, 0.059345f, 0.059262f, 0.059179f, 0.059097f,
        0.059014f, 0.058932f, 0.058850f, 0.058768f, 0.058686f, 0.058604f, 0.058522f, 0.058440f, 0.058358f, 0.058277f,
        0.058196f, 0.058114f, 0.058033f, 0.057952f, 0.057871f, 0.057790f, 0.057710f, 0.057629f, 0.057548f, 0.057468f,
        0.057388f, 0.057307f, 0.057227f, 0.057147f, 0.057067f, 0.056988f, 0.056908f, 0.056828f, 0.056749f, 0.056670f,
        0.056590f, 0.056511f, 0.056432f});
    static constexpr Spectrum CES85({0.014249f, 0.013601f, 0.011656f, 0.011680f, 0.012398f, 0.012915f, 0.012881f,
        0.013581f, 0.014289f, 0.015239f, 0.017314f, 0.018985f, 0.019990f, 0.022170f, 0.023181f, 0.024437f, 0.026441f,
        0.029210f, 0.032328f, 0.035656f, 0.038310f, 0.039964f, 0.041976f, 0.044240f, 0.046645f, 0.049395f, 0.052798f,
        0.056377f, 0.059012f, 0.061242f, 0.063243f, 0.065841f, 0.069235f, 0.072588f, 0.076166f, 0.079552f, 0.082672f,
        0.085384f, 0.088418f, 0.091860f, 0.095230f, 0.098943f, 0.101950f, 0.104700f, 0.108480f, 0.112330f, 0.116270f,
        0.120570f, 0.124480f, 0.128220f, 0.131660f, 0.135460f, 0.139270f, 0.143670f, 0.147560f, 0.149500f, 0.151390f,
        0.153380f, 0.155190f, 0.156720f, 0.158200f, 0.160290f, 0.161460f, 0.161590f, 0.161550f, 0.162750f, 0.165140f,
        0.166730f, 0.168370f, 0.169610f, 0.171070f, 0.171860f, 0.172030f, 0.173020f, 0.174320f, 0.174960f, 0.174800f,
        0.174970f, 0.175520f, 0.175000f, 0.173950f, 0.173150f, 0.173480f, 0.174220f, 0.174400f, 0.173700f, 0.174110f,
        0.174470f, 0.173730f, 0.173010f, 0.172160f, 0.171710f, 0.171270f, 0.169930f, 0.168450f, 0.167270f, 0.167120f,
        0.165620f, 0.163550f, 0.161780f, 0.160070f, 0.159260f, 0.158500f, 0.157780f, 0.156920f, 0.155760f, 0.154480f,
        0.152970f, 0.152000f, 0.150700f, 0.149270f, 0.147950f, 0.146550f, 0.145620f, 0.144570f, 0.143210f, 0.141730f,
        0.139730f, 0.137650f, 0.135450f, 0.134360f, 0.133270f, 0.132070f, 0.130440f, 0.128370f, 0.126440f, 0.125070f,
        0.123830f, 0.122800f, 0.121920f, 0.121230f, 0.120180f, 0.118940f, 0.117330f, 0.116430f, 0.115990f, 0.115200f,
        0.113980f, 0.112750f, 0.111520f, 0.110040f, 0.108310f, 0.106840f, 0.105930f, 0.105250f, 0.104000f, 0.102410f,
        0.100770f, 0.099384f, 0.098325f, 0.097159f, 0.096076f, 0.095405f, 0.094822f, 0.094267f, 0.093593f, 0.093181f,
        0.092920f, 0.092830f, 0.092639f, 0.092134f, 0.092305f, 0.092678f, 0.093017f, 0.093029f, 0.092614f, 0.092157f,
        0.091279f, 0.090637f, 0.090438f, 0.090657f, 0.090810f, 0.090329f, 0.089886f, 0.089573f, 0.089016f, 0.088773f,
        0.088818f, 0.089194f, 0.089247f, 0.089068f, 0.088876f, 0.088609f, 0.088624f, 0.088360f, 0.087918f, 0.088044f,
        0.088265f, 0.088420f, 0.088144f, 0.087941f, 0.087402f, 0.086601f, 0.086218f, 0.085813f, 0.085875f, 0.085788f,
        0.085625f, 0.085083f, 0.084757f, 0.084629f, 0.084470f, 0.084794f, 0.085620f, 0.086012f, 0.087252f, 0.088404f,
        0.089045f, 0.089310f, 0.089975f, 0.090875f, 0.092511f, 0.094109f, 0.095129f, 0.095883f, 0.097581f, 0.099109f,
        0.099961f, 0.100600f, 0.101520f, 0.102740f, 0.104080f, 0.105060f, 0.106250f, 0.107710f, 0.108740f, 0.108860f,
        0.107990f, 0.107040f, 0.108160f, 0.109590f, 0.110490f, 0.110570f, 0.110690f, 0.111220f, 0.111710f, 0.111890f,
        0.112480f, 0.113780f, 0.115750f, 0.116030f, 0.116060f, 0.117560f, 0.119660f, 0.121610f, 0.122760f, 0.124320f,
        0.125600f, 0.126890f, 0.127910f, 0.129160f, 0.131510f, 0.135090f, 0.136550f, 0.137460f, 0.138260f, 0.139280f,
        0.140280f, 0.141920f, 0.143620f, 0.144640f, 0.145510f, 0.146010f, 0.146020f, 0.147030f, 0.150490f, 0.153660f,
        0.156800f, 0.160010f, 0.163940f, 0.167850f, 0.172910f, 0.178570f, 0.184360f, 0.189100f, 0.192650f, 0.194600f,
        0.196900f, 0.200190f, 0.203910f, 0.207840f, 0.212290f, 0.216580f, 0.220300f, 0.223870f, 0.227940f, 0.233120f,
        0.238120f, 0.243820f, 0.248490f, 0.252550f, 0.256260f, 0.259640f, 0.263050f, 0.266580f, 0.269740f, 0.273160f,
        0.277130f, 0.280980f, 0.283950f, 0.287370f, 0.290550f, 0.294450f, 0.298320f, 0.301660f, 0.304810f, 0.307440f,
        0.309660f, 0.311930f, 0.314130f, 0.316300f, 0.318190f, 0.320520f, 0.322430f, 0.324190f, 0.324170f, 0.324540f,
        0.328790f, 0.335870f, 0.338780f, 0.338780f, 0.340040f, 0.342290f, 0.344550f, 0.346810f, 0.349080f, 0.351360f,
        0.353640f, 0.355930f, 0.358230f, 0.360540f, 0.362850f, 0.365160f, 0.367490f, 0.369820f, 0.372150f, 0.374490f,
        0.376840f, 0.379190f, 0.381550f, 0.383920f, 0.386290f, 0.388660f, 0.391040f, 0.393430f, 0.395820f, 0.398220f,
        0.400620f, 0.403020f, 0.405440f, 0.407850f, 0.410270f, 0.412690f, 0.415120f, 0.417550f, 0.419990f, 0.422430f,
        0.424870f, 0.427320f, 0.429770f, 0.432230f, 0.434690f, 0.437150f, 0.439610f, 0.442080f, 0.444550f, 0.447020f,
        0.449500f, 0.451970f, 0.454450f, 0.456940f, 0.459420f, 0.461910f, 0.464400f, 0.466890f, 0.469380f, 0.471870f,
        0.474370f, 0.476860f, 0.479360f, 0.481860f, 0.484360f, 0.486860f, 0.489360f, 0.491860f, 0.494360f, 0.496870f,
        0.499370f, 0.501870f, 0.504370f, 0.506870f, 0.509380f, 0.511880f, 0.514380f, 0.516880f, 0.519380f, 0.521870f,
        0.524370f, 0.526870f, 0.529360f});
    static constexpr Spectrum CES86({0.502520f, 0.504820f, 0.507110f, 0.509410f, 0.511700f, 0.514000f, 0.516290f,
        0.518580f, 0.520870f, 0.523170f, 0.525460f, 0.527740f, 0.530030f, 0.532320f, 0.534600f, 0.536890f, 0.539170f,
        0.541450f, 0.543730f, 0.546010f, 0.546360f, 0.549520f, 0.552510f, 0.555340f, 0.558020f, 0.560560f, 0.562960f,
        0.565240f, 0.567390f, 0.569430f, 0.571380f, 0.573220f, 0.574980f, 0.576660f, 0.578260f, 0.579810f, 0.581300f,
        0.582740f, 0.584140f, 0.585510f, 0.586860f, 0.588190f, 0.589480f, 0.590700f, 0.591830f, 0.592840f, 0.593700f,
        0.594400f, 0.594900f, 0.595170f, 0.595200f, 0.594970f, 0.594480f, 0.593770f, 0.592850f, 0.591760f, 0.590510f,
        0.589140f, 0.587660f, 0.586100f, 0.584480f, 0.582830f, 0.581130f, 0.579390f, 0.577600f, 0.575750f, 0.573840f,
        0.571860f, 0.569800f, 0.567650f, 0.565420f, 0.563090f, 0.560670f, 0.558160f, 0.555560f, 0.552870f, 0.550100f,
        0.547240f, 0.544310f, 0.541300f, 0.538210f, 0.535060f, 0.531830f, 0.528530f, 0.525170f, 0.521740f, 0.518240f,
        0.514690f, 0.511070f, 0.507400f, 0.503660f, 0.499880f, 0.496050f, 0.492180f, 0.488280f, 0.484370f, 0.480450f,
        0.476520f, 0.472610f, 0.468710f, 0.464840f, 0.461010f, 0.457210f, 0.453440f, 0.449700f, 0.445970f, 0.442270f,
        0.438590f, 0.434920f, 0.431260f, 0.427610f, 0.423970f, 0.420320f, 0.416670f, 0.413010f, 0.409340f, 0.405650f,
        0.401950f, 0.398220f, 0.394470f, 0.390680f, 0.386860f, 0.383020f, 0.379160f, 0.375300f, 0.371450f, 0.367610f,
        0.363810f, 0.360040f, 0.356320f, 0.352650f, 0.349060f, 0.345540f, 0.342100f, 0.338750f, 0.335490f, 0.332330f,
        0.329280f, 0.326330f, 0.323500f, 0.320780f, 0.318200f, 0.315730f, 0.313360f, 0.311090f, 0.308900f, 0.306780f,
        0.304730f, 0.302730f, 0.300770f, 0.298840f, 0.296930f, 0.295030f, 0.293130f, 0.291210f, 0.289270f, 0.287290f,
        0.285260f, 0.283180f, 0.281020f, 0.278790f, 0.276470f, 0.274070f, 0.271630f, 0.269140f, 0.266640f, 0.264150f,
        0.261670f, 0.259240f, 0.256860f, 0.254560f, 0.252360f, 0.250270f, 0.248310f, 0.246490f, 0.244830f, 0.243340f,
        0.242050f, 0.240960f, 0.240100f, 0.239470f, 0.239090f, 0.238940f, 0.239010f, 0.239270f, 0.239690f, 0.240270f,
        0.240970f, 0.241780f, 0.242680f, 0.243640f, 0.244650f, 0.245680f, 0.246720f, 0.247750f, 0.248750f, 0.249690f,
        0.250560f, 0.251350f, 0.252030f, 0.252580f, 0.252980f, 0.253250f, 0.253390f, 0.253410f, 0.253310f, 0.253100f,
        0.252800f, 0.252410f, 0.251930f, 0.251390f, 0.250780f, 0.250130f, 0.249490f, 0.248880f, 0.248330f, 0.247880f,
        0.247560f, 0.247410f, 0.247440f, 0.247710f, 0.248240f, 0.249050f, 0.250160f, 0.251600f, 0.253390f, 0.255550f,
        0.258110f, 0.261080f, 0.264490f, 0.268360f, 0.272710f, 0.277540f, 0.282850f, 0.288640f, 0.294910f, 0.301660f,
        0.308880f, 0.316580f, 0.324750f, 0.333390f, 0.342510f, 0.352070f, 0.362060f, 0.372460f, 0.383250f, 0.394400f,
        0.405890f, 0.417710f, 0.429820f, 0.442210f, 0.454850f, 0.467680f, 0.480640f, 0.493670f, 0.506700f, 0.519670f,
        0.532510f, 0.545180f, 0.557590f, 0.569690f, 0.581420f, 0.592770f, 0.603720f, 0.614260f, 0.624380f, 0.634080f,
        0.643340f, 0.652150f, 0.660500f, 0.668370f, 0.675780f, 0.682710f, 0.689190f, 0.695240f, 0.700870f, 0.706080f,
        0.710910f, 0.715360f, 0.719440f, 0.723180f, 0.726580f, 0.729680f, 0.732490f, 0.735030f, 0.737330f, 0.739410f,
        0.741300f, 0.743010f, 0.744580f, 0.746010f, 0.747350f, 0.748580f, 0.749740f, 0.750820f, 0.751850f, 0.752820f,
        0.753760f, 0.754670f, 0.755560f, 0.756440f, 0.757330f, 0.758220f, 0.759130f, 0.760070f, 0.761020f, 0.762010f,
        0.763030f, 0.764100f, 0.765210f, 0.766370f, 0.767590f, 0.768870f, 0.770210f, 0.771630f, 0.773130f, 0.774710f,
        0.776370f, 0.778130f, 0.779990f, 0.781950f, 0.782060f, 0.783420f, 0.784770f, 0.786120f, 0.787470f, 0.788800f,
        0.790130f, 0.791460f, 0.792780f, 0.794090f, 0.795390f, 0.796690f, 0.797990f, 0.799270f, 0.800560f, 0.801830f,
        0.803100f, 0.804360f, 0.805620f, 0.806870f, 0.808110f, 0.809350f, 0.810580f, 0.811810f, 0.813030f, 0.814240f,
        0.815450f, 0.816650f, 0.817850f, 0.819040f, 0.820220f, 0.821400f, 0.822570f, 0.823740f, 0.824900f, 0.826050f,
        0.827200f, 0.828340f, 0.829470f, 0.830600f, 0.831730f, 0.832840f, 0.833960f, 0.835060f, 0.836160f, 0.837250f,
        0.838340f, 0.839420f, 0.840500f, 0.841570f, 0.842640f, 0.843690f, 0.844750f, 0.845790f, 0.846830f, 0.847870f,
        0.848900f, 0.849920f, 0.850940f, 0.851950f, 0.852960f, 0.853960f, 0.854960f, 0.855950f, 0.856930f, 0.857910f,
        0.858880f, 0.859850f, 0.860810f, 0.861770f, 0.862720f, 0.863670f, 0.864610f, 0.865540f, 0.866470f, 0.867390f,
        0.868310f, 0.869220f, 0.870130f});
    static constexpr Spectrum CES87({0.106130f, 0.108020f, 0.109930f, 0.111880f, 0.113850f, 0.115860f, 0.117890f,
        0.119960f, 0.122060f, 0.124190f, 0.126350f, 0.128540f, 0.130760f, 0.133020f, 0.135310f, 0.137640f, 0.139990f,
        0.142390f, 0.144810f, 0.147270f, 0.145960f, 0.150310f, 0.154350f, 0.158090f, 0.161550f, 0.164760f, 0.167720f,
        0.170470f, 0.173000f, 0.175350f, 0.177530f, 0.179560f, 0.181460f, 0.183240f, 0.184920f, 0.186530f, 0.188070f,
        0.189580f, 0.191060f, 0.192530f, 0.194010f, 0.195520f, 0.197070f, 0.198650f, 0.200270f, 0.201930f, 0.203640f,
        0.205400f, 0.207200f, 0.209070f, 0.210990f, 0.212970f, 0.214990f, 0.217040f, 0.219100f, 0.221140f, 0.223160f,
        0.225120f, 0.227020f, 0.228840f, 0.230550f, 0.232140f, 0.233600f, 0.234900f, 0.236040f, 0.237000f, 0.237750f,
        0.238300f, 0.238610f, 0.238680f, 0.238490f, 0.238040f, 0.237320f, 0.236370f, 0.235210f, 0.233840f, 0.232290f,
        0.230590f, 0.228740f, 0.226770f, 0.224690f, 0.222530f, 0.220300f, 0.217990f, 0.215630f, 0.213220f, 0.210770f,
        0.208290f, 0.205790f, 0.203280f, 0.200760f, 0.198250f, 0.195750f, 0.193250f, 0.190770f, 0.188290f, 0.185830f,
        0.183380f, 0.180950f, 0.178540f, 0.176140f, 0.173760f, 0.171410f, 0.169080f, 0.166770f, 0.164480f, 0.162210f,
        0.159970f, 0.157760f, 0.155570f, 0.153400f, 0.151260f, 0.149150f, 0.147050f, 0.144970f, 0.142900f, 0.140830f,
        0.138770f, 0.136700f, 0.134630f, 0.132550f, 0.130460f, 0.128350f, 0.126240f, 0.124110f, 0.121970f, 0.119830f,
        0.117680f, 0.115520f, 0.113370f, 0.111200f, 0.109040f, 0.106890f, 0.104750f, 0.102630f, 0.100540f, 0.098496f,
        0.096494f, 0.094547f, 0.092663f, 0.090850f, 0.089115f, 0.087458f, 0.085880f, 0.084381f, 0.082961f, 0.081620f,
        0.080358f, 0.079176f, 0.078073f, 0.077049f, 0.076103f, 0.075224f, 0.074401f, 0.073623f, 0.072876f, 0.072149f,
        0.071430f, 0.070708f, 0.069970f, 0.069205f, 0.068403f, 0.067562f, 0.066684f, 0.065770f, 0.064821f, 0.063836f,
        0.062819f, 0.061768f, 0.060687f, 0.059574f, 0.058433f, 0.057271f, 0.056096f, 0.054918f, 0.053745f, 0.052585f,
        0.051448f, 0.050341f, 0.049274f, 0.048255f, 0.047292f, 0.046387f, 0.045541f, 0.044756f, 0.044034f, 0.043375f,
        0.042782f, 0.042254f, 0.041795f, 0.041404f, 0.041088f, 0.040872f, 0.040783f, 0.040852f, 0.041106f, 0.041575f,
        0.042286f, 0.043270f, 0.044555f, 0.046170f, 0.048135f, 0.050444f, 0.053079f, 0.056024f, 0.059265f, 0.062784f,
        0.066567f, 0.070596f, 0.074857f, 0.079333f, 0.084001f, 0.088811f, 0.093705f, 0.098627f, 0.103520f, 0.108320f,
        0.112980f, 0.117430f, 0.121630f, 0.125500f, 0.129020f, 0.132180f, 0.135030f, 0.137570f, 0.139840f, 0.141870f,
        0.143680f, 0.145290f, 0.146740f, 0.148040f, 0.149230f, 0.150310f, 0.151290f, 0.152190f, 0.153010f, 0.153770f,
        0.154480f, 0.155140f, 0.155770f, 0.156380f, 0.156980f, 0.157560f, 0.158140f, 0.158700f, 0.159270f, 0.159820f,
        0.160380f, 0.160930f, 0.161490f, 0.162040f, 0.162600f, 0.163160f, 0.163740f, 0.164320f, 0.164920f, 0.165530f,
        0.166160f, 0.166820f, 0.167490f, 0.168200f, 0.168930f, 0.169690f, 0.170480f, 0.171300f, 0.172150f, 0.173030f,
        0.173930f, 0.174870f, 0.175840f, 0.176840f, 0.177860f, 0.178900f, 0.179950f, 0.180990f, 0.182000f, 0.182980f,
        0.183910f, 0.184770f, 0.185560f, 0.186270f, 0.186870f, 0.187380f, 0.187780f, 0.188090f, 0.188300f, 0.188410f,
        0.188440f, 0.188370f, 0.188210f, 0.187960f, 0.187620f, 0.187210f, 0.186750f, 0.186240f, 0.185710f, 0.185160f,
        0.184620f, 0.184100f, 0.183620f, 0.183190f, 0.182820f, 0.182510f, 0.182250f, 0.182040f, 0.181850f, 0.181700f,
        0.181570f, 0.181450f, 0.181330f, 0.181200f, 0.181070f, 0.180920f, 0.180730f, 0.180520f, 0.180260f, 0.179940f,
        0.179570f, 0.179130f, 0.178620f, 0.178030f, 0.178470f, 0.178230f, 0.177990f, 0.177750f, 0.177510f, 0.177270f,
        0.177030f, 0.176790f, 0.176560f, 0.176320f, 0.176080f, 0.175840f, 0.175600f, 0.175370f, 0.175130f, 0.174890f,
        0.174660f, 0.174420f, 0.174190f, 0.173950f, 0.173710f, 0.173480f, 0.173240f, 0.173010f, 0.172780f, 0.172540f,
        0.172310f, 0.172070f, 0.171840f, 0.171610f, 0.171370f, 0.171140f, 0.170910f, 0.170680f, 0.170450f, 0.170210f,
        0.169980f, 0.169750f, 0.169520f, 0.169290f, 0.169060f, 0.168830f, 0.168600f, 0.168370f, 0.168140f, 0.167910f,
        0.167680f, 0.167450f, 0.167230f, 0.167000f, 0.166770f, 0.166540f, 0.166310f, 0.166090f, 0.165860f, 0.165630f,
        0.165410f, 0.165180f, 0.164960f, 0.164730f, 0.164510f, 0.164280f, 0.164060f, 0.163830f, 0.163610f, 0.163380f,
        0.163160f, 0.162930f, 0.162710f, 0.162490f, 0.162260f, 0.162040f, 0.161820f, 0.161600f, 0.161380f, 0.161150f,
        0.160930f, 0.160710f, 0.160490f});
    static constexpr Spectrum CES88({0.426400f, 0.429170f, 0.431960f, 0.434740f, 0.437530f, 0.440330f, 0.443130f,
        0.445930f, 0.448730f, 0.451540f, 0.454350f, 0.457170f, 0.459990f, 0.462810f, 0.465630f, 0.468450f, 0.471280f,
        0.474110f, 0.476940f, 0.479770f, 0.474900f, 0.481290f, 0.487020f, 0.492110f, 0.496610f, 0.500540f, 0.503940f,
        0.506830f, 0.509260f, 0.511250f, 0.512830f, 0.514050f, 0.514920f, 0.515490f, 0.515790f, 0.515840f, 0.515690f,
        0.515350f, 0.514880f, 0.514290f, 0.513630f, 0.512910f, 0.512130f, 0.511290f, 0.510380f, 0.509380f, 0.508290f,
        0.507090f, 0.505790f, 0.504360f, 0.502800f, 0.501110f, 0.499290f, 0.497340f, 0.495280f, 0.493130f, 0.490880f,
        0.488550f, 0.486140f, 0.483680f, 0.481160f, 0.478600f, 0.476000f, 0.473370f, 0.470710f, 0.468040f, 0.465360f,
        0.462670f, 0.459990f, 0.457310f, 0.454650f, 0.452010f, 0.449390f, 0.446790f, 0.444210f, 0.441640f, 0.439090f,
        0.436560f, 0.434040f, 0.431530f, 0.429030f, 0.426540f, 0.424060f, 0.421590f, 0.419110f, 0.416620f, 0.414130f,
        0.411630f, 0.409110f, 0.406570f, 0.404010f, 0.401420f, 0.398800f, 0.396140f, 0.393440f, 0.390680f, 0.387870f,
        0.384990f, 0.382040f, 0.379020f, 0.375910f, 0.372720f, 0.369440f, 0.366070f, 0.362630f, 0.359120f, 0.355530f,
        0.351880f, 0.348170f, 0.344390f, 0.340560f, 0.336690f, 0.332780f, 0.328850f, 0.324930f, 0.321030f, 0.317180f,
        0.313380f, 0.309660f, 0.306040f, 0.302540f, 0.299160f, 0.295900f, 0.292730f, 0.289650f, 0.286620f, 0.283640f,
        0.280680f, 0.277720f, 0.274760f, 0.271760f, 0.268710f, 0.265610f, 0.262460f, 0.259250f, 0.255980f, 0.252660f,
        0.249270f, 0.245810f, 0.242290f, 0.238690f, 0.235030f, 0.231340f, 0.227650f, 0.224000f, 0.220420f, 0.216950f,
        0.213630f, 0.210500f, 0.207590f, 0.204930f, 0.202570f, 0.200480f, 0.198670f, 0.197110f, 0.195800f, 0.194730f,
        0.193870f, 0.193230f, 0.192780f, 0.192520f, 0.192430f, 0.192480f, 0.192640f, 0.192870f, 0.193150f, 0.193430f,
        0.193690f, 0.193900f, 0.194020f, 0.194010f, 0.193860f, 0.193570f, 0.193130f, 0.192570f, 0.191900f, 0.191110f,
        0.190230f, 0.189250f, 0.188190f, 0.187060f, 0.185880f, 0.184680f, 0.183510f, 0.182440f, 0.181510f, 0.180770f,
        0.180260f, 0.180050f, 0.180180f, 0.180710f, 0.181670f, 0.183060f, 0.184890f, 0.187150f, 0.189830f, 0.192930f,
        0.196440f, 0.200360f, 0.204680f, 0.209400f, 0.214510f, 0.219960f, 0.225710f, 0.231720f, 0.237950f, 0.244350f,
        0.250870f, 0.257480f, 0.264120f, 0.270760f, 0.277360f, 0.283890f, 0.290350f, 0.296740f, 0.303020f, 0.309210f,
        0.315280f, 0.321230f, 0.327050f, 0.332720f, 0.338240f, 0.343610f, 0.348820f, 0.353900f, 0.358830f, 0.363620f,
        0.368280f, 0.372810f, 0.377200f, 0.381470f, 0.385620f, 0.389630f, 0.393500f, 0.397210f, 0.400760f, 0.404130f,
        0.407310f, 0.410300f, 0.413070f, 0.415630f, 0.417960f, 0.420060f, 0.421940f, 0.423610f, 0.425080f, 0.426350f,
        0.427420f, 0.428300f, 0.429000f, 0.429530f, 0.429890f, 0.430090f, 0.430150f, 0.430080f, 0.429900f, 0.429610f,
        0.429240f, 0.428800f, 0.428290f, 0.427740f, 0.427160f, 0.426570f, 0.425990f, 0.425450f, 0.424970f, 0.424580f,
        0.424280f, 0.424120f, 0.424110f, 0.424270f, 0.424620f, 0.425190f, 0.425990f, 0.427030f, 0.428330f, 0.429900f,
        0.431770f, 0.433940f, 0.436430f, 0.439260f, 0.442440f, 0.445980f, 0.449880f, 0.454140f, 0.458770f, 0.463780f,
        0.469160f, 0.474920f, 0.481070f, 0.487610f, 0.494540f, 0.501830f, 0.509460f, 0.517400f, 0.525620f, 0.534090f,
        0.542780f, 0.551670f, 0.560730f, 0.569920f, 0.579230f, 0.588620f, 0.598060f, 0.607510f, 0.616950f, 0.626350f,
        0.635670f, 0.644880f, 0.653950f, 0.662860f, 0.671560f, 0.680040f, 0.688250f, 0.696160f, 0.703750f, 0.710990f,
        0.717840f, 0.724260f, 0.730240f, 0.735740f, 0.745970f, 0.752850f, 0.759610f, 0.766250f, 0.772750f, 0.779130f,
        0.785380f, 0.791490f, 0.797480f, 0.803340f, 0.809070f, 0.814670f, 0.820140f, 0.825490f, 0.830710f, 0.835810f,
        0.840780f, 0.845630f, 0.850350f, 0.854960f, 0.859450f, 0.863820f, 0.868080f, 0.872220f, 0.876250f, 0.880180f,
        0.883990f, 0.887700f, 0.891300f, 0.894810f, 0.898210f, 0.901510f, 0.904720f, 0.907840f, 0.910860f, 0.913790f,
        0.916640f, 0.919400f, 0.922080f, 0.924670f, 0.927190f, 0.929620f, 0.931980f, 0.934270f, 0.936490f, 0.938640f,
        0.940720f, 0.942730f, 0.944680f, 0.946560f, 0.948390f, 0.950150f, 0.951860f, 0.953510f, 0.955110f, 0.956660f,
        0.958160f, 0.959600f, 0.961000f, 0.962350f, 0.963660f, 0.964920f, 0.966140f, 0.967320f, 0.968460f, 0.969560f,
        0.970620f, 0.971650f, 0.972650f, 0.973600f, 0.974530f, 0.975430f, 0.976290f, 0.977120f, 0.977930f, 0.978710f,
        0.979460f, 0.980180f, 0.980880f});
    static constexpr Spectrum CES89({0.072745f, 0.074530f, 0.076356f, 0.078222f, 0.080130f, 0.082080f, 0.084073f,
        0.086111f, 0.088192f, 0.090320f, 0.092493f, 0.094713f, 0.096981f, 0.099297f, 0.101660f, 0.104080f, 0.106540f,
        0.109060f, 0.111630f, 0.114250f, 0.115460f, 0.118710f, 0.121990f, 0.125300f, 0.128610f, 0.131910f, 0.135180f,
        0.138400f, 0.141550f, 0.144610f, 0.147570f, 0.150400f, 0.153090f, 0.155620f, 0.157960f, 0.160110f, 0.162050f,
        0.163740f, 0.165190f, 0.166360f, 0.167240f, 0.167820f, 0.168150f, 0.168300f, 0.168320f, 0.168280f, 0.168240f,
        0.168250f, 0.168390f, 0.168710f, 0.169270f, 0.170120f, 0.171220f, 0.172500f, 0.173900f, 0.175370f, 0.176850f,
        0.178270f, 0.179580f, 0.180710f, 0.181600f, 0.182220f, 0.182570f, 0.182680f, 0.182580f, 0.182300f, 0.181850f,
        0.181280f, 0.180600f, 0.179850f, 0.179060f, 0.178230f, 0.177380f, 0.176500f, 0.175560f, 0.174580f, 0.173530f,
        0.172410f, 0.171220f, 0.169940f, 0.168560f, 0.167080f, 0.165500f, 0.163820f, 0.162050f, 0.160170f, 0.158200f,
        0.156140f, 0.153980f, 0.151740f, 0.149400f, 0.146980f, 0.144490f, 0.141950f, 0.139380f, 0.136790f, 0.134200f,
        0.131630f, 0.129100f, 0.126630f, 0.124230f, 0.121910f, 0.119680f, 0.117520f, 0.115410f, 0.113360f, 0.111340f,
        0.109350f, 0.107380f, 0.105410f, 0.103440f, 0.101460f, 0.099471f, 0.097490f, 0.095523f, 0.093581f, 0.091674f,
        0.089811f, 0.088002f, 0.086257f, 0.084586f, 0.082996f, 0.081487f, 0.080056f, 0.078699f, 0.077414f, 0.076196f,
        0.075044f, 0.073954f, 0.072924f, 0.071949f, 0.071027f, 0.070152f, 0.069319f, 0.068523f, 0.067757f, 0.067017f,
        0.066298f, 0.065592f, 0.064896f, 0.064204f, 0.063511f, 0.062819f, 0.062131f, 0.061451f, 0.060780f, 0.060122f,
        0.059481f, 0.058858f, 0.058257f, 0.057681f, 0.057133f, 0.056615f, 0.056129f, 0.055678f, 0.055263f, 0.054886f,
        0.054550f, 0.054257f, 0.054009f, 0.053809f, 0.053657f, 0.053557f, 0.053510f, 0.053518f, 0.053583f, 0.053707f,
        0.053891f, 0.054138f, 0.054450f, 0.054828f, 0.055273f, 0.055781f, 0.056346f, 0.056964f, 0.057627f, 0.058332f,
        0.059073f, 0.059843f, 0.060638f, 0.061452f, 0.062280f, 0.063117f, 0.063959f, 0.064801f, 0.065639f, 0.066468f,
        0.067283f, 0.068081f, 0.068856f, 0.069605f, 0.070324f, 0.071020f, 0.071699f, 0.072370f, 0.073040f, 0.073716f,
        0.074406f, 0.075118f, 0.075859f, 0.076637f, 0.077458f, 0.078323f, 0.079234f, 0.080191f, 0.081195f, 0.082245f,
        0.083344f, 0.084490f, 0.085685f, 0.086930f, 0.088224f, 0.089566f, 0.090954f, 0.092386f, 0.093861f, 0.095377f,
        0.096931f, 0.098523f, 0.100150f, 0.101810f, 0.103500f, 0.105210f, 0.106950f, 0.108700f, 0.110450f, 0.112210f,
        0.113960f, 0.115710f, 0.117430f, 0.119130f, 0.120810f, 0.122450f, 0.124050f, 0.125620f, 0.127140f, 0.128610f,
        0.130020f, 0.131380f, 0.132680f, 0.133910f, 0.135080f, 0.136180f, 0.137230f, 0.138240f, 0.139220f, 0.140170f,
        0.141100f, 0.142030f, 0.142960f, 0.143900f, 0.144860f, 0.145850f, 0.146870f, 0.147950f, 0.149090f, 0.150300f,
        0.151580f, 0.152960f, 0.154440f, 0.156030f, 0.157730f, 0.159560f, 0.161520f, 0.163610f, 0.165830f, 0.168190f,
        0.170680f, 0.173330f, 0.176120f, 0.179060f, 0.182150f, 0.185400f, 0.188800f, 0.192360f, 0.196070f, 0.199940f,
        0.203960f, 0.208140f, 0.212480f, 0.216970f, 0.221620f, 0.226420f, 0.231370f, 0.236470f, 0.241700f, 0.247080f,
        0.252580f, 0.258220f, 0.263980f, 0.269860f, 0.275860f, 0.281980f, 0.288220f, 0.294580f, 0.301070f, 0.307690f,
        0.314440f, 0.321320f, 0.328330f, 0.335490f, 0.342780f, 0.350200f, 0.357730f, 0.365370f, 0.373090f, 0.380880f,
        0.388740f, 0.396640f, 0.404580f, 0.412530f, 0.420500f, 0.428460f, 0.436410f, 0.444320f, 0.452190f, 0.460000f,
        0.467730f, 0.475390f, 0.482950f, 0.490390f, 0.499550f, 0.507560f, 0.515570f, 0.523570f, 0.531550f, 0.539530f,
        0.547480f, 0.555400f, 0.563300f, 0.571170f, 0.579000f, 0.586790f, 0.594540f, 0.602240f, 0.609890f, 0.617490f,
        0.625030f, 0.632510f, 0.639930f, 0.647280f, 0.654560f, 0.661770f, 0.668910f, 0.675960f, 0.682940f, 0.689840f,
        0.696660f, 0.703390f, 0.710030f, 0.716580f, 0.723050f, 0.729420f, 0.735700f, 0.741880f, 0.747970f, 0.753960f,
        0.759860f, 0.765660f, 0.771360f, 0.776960f, 0.782470f, 0.787870f, 0.793180f, 0.798390f, 0.803500f, 0.808510f,
        0.813420f, 0.818230f, 0.822950f, 0.827570f, 0.832100f, 0.836530f, 0.840860f, 0.845110f, 0.849250f, 0.853310f,
        0.857280f, 0.861150f, 0.864940f, 0.868640f, 0.872260f, 0.875780f, 0.879230f, 0.882590f, 0.885870f, 0.889070f,
        0.892190f, 0.895240f, 0.898200f, 0.901100f, 0.903920f, 0.906660f, 0.909340f, 0.911950f, 0.914490f, 0.916960f,
        0.919370f, 0.921710f, 0.923990f});
    static constexpr Spectrum CES90({0.221780f, 0.223650f, 0.225540f, 0.227440f, 0.229340f, 0.231270f, 0.233200f,
        0.235140f, 0.237090f, 0.239060f, 0.241030f, 0.243020f, 0.245020f, 0.247030f, 0.249050f, 0.251080f, 0.253120f,
        0.255180f, 0.257240f, 0.259320f, 0.261700f, 0.263320f, 0.265330f, 0.267620f, 0.270120f, 0.272710f, 0.275300f,
        0.277800f, 0.280120f, 0.282150f, 0.283800f, 0.285000f, 0.285770f, 0.286160f, 0.286210f, 0.285980f, 0.285520f,
        0.284860f, 0.284050f, 0.283150f, 0.282200f, 0.281240f, 0.280260f, 0.279260f, 0.278210f, 0.277110f, 0.275940f,
        0.274680f, 0.273330f, 0.271880f, 0.270300f, 0.268590f, 0.266760f, 0.264810f, 0.262770f, 0.260630f, 0.258410f,
        0.256130f, 0.253800f, 0.251410f, 0.249000f, 0.246560f, 0.244110f, 0.241630f, 0.239140f, 0.236630f, 0.234110f,
        0.231570f, 0.229020f, 0.226470f, 0.223900f, 0.221330f, 0.218750f, 0.216180f, 0.213620f, 0.211070f, 0.208540f,
        0.206030f, 0.203550f, 0.201110f, 0.198700f, 0.196340f, 0.194010f, 0.191710f, 0.189440f, 0.187180f, 0.184930f,
        0.182690f, 0.180440f, 0.178180f, 0.175900f, 0.173600f, 0.171270f, 0.168920f, 0.166560f, 0.164170f, 0.161760f,
        0.159340f, 0.156910f, 0.154460f, 0.152000f, 0.149530f, 0.147070f, 0.144610f, 0.142180f, 0.139790f, 0.137430f,
        0.135140f, 0.132910f, 0.130760f, 0.128700f, 0.126740f, 0.124870f, 0.123090f, 0.121390f, 0.119770f, 0.118210f,
        0.116730f, 0.115300f, 0.113930f, 0.112600f, 0.111310f, 0.110060f, 0.108830f, 0.107620f, 0.106420f, 0.105210f,
        0.103990f, 0.102760f, 0.101500f, 0.100200f, 0.098861f, 0.097485f, 0.096079f, 0.094651f, 0.093207f, 0.091754f,
        0.090300f, 0.088851f, 0.087416f, 0.086000f, 0.084611f, 0.083256f, 0.081939f, 0.080667f, 0.079445f, 0.078280f,
        0.077177f, 0.076142f, 0.075181f, 0.074300f, 0.073503f, 0.072791f, 0.072160f, 0.071610f, 0.071139f, 0.070745f,
        0.070425f, 0.070179f, 0.070005f, 0.069900f, 0.069862f, 0.069882f, 0.069949f, 0.070054f, 0.070187f, 0.070336f,
        0.070492f, 0.070645f, 0.070784f, 0.070900f, 0.070984f, 0.071040f, 0.071073f, 0.071088f, 0.071090f, 0.071085f,
        0.071078f, 0.071075f, 0.071080f, 0.071100f, 0.071140f, 0.071213f, 0.071331f, 0.071507f, 0.071754f, 0.072085f,
        0.072512f, 0.073048f, 0.073706f, 0.074500f, 0.075441f, 0.076537f, 0.077798f, 0.079231f, 0.080844f, 0.082646f,
        0.084644f, 0.086847f, 0.089263f, 0.091900f, 0.094761f, 0.097832f, 0.101090f, 0.104520f, 0.108110f, 0.111820f,
        0.115650f, 0.119560f, 0.123560f, 0.127600f, 0.131680f, 0.135770f, 0.139850f, 0.143910f, 0.147910f, 0.151850f,
        0.155700f, 0.159440f, 0.163040f, 0.166500f, 0.169790f, 0.172910f, 0.175860f, 0.178660f, 0.181300f, 0.183780f,
        0.186120f, 0.188320f, 0.190380f, 0.192300f, 0.194090f, 0.195760f, 0.197310f, 0.198760f, 0.200110f, 0.201380f,
        0.202560f, 0.203670f, 0.204710f, 0.205700f, 0.206640f, 0.207540f, 0.208390f, 0.209210f, 0.209990f, 0.210730f,
        0.211440f, 0.212120f, 0.212770f, 0.213400f, 0.214010f, 0.214590f, 0.215160f, 0.215710f, 0.216260f, 0.216790f,
        0.217320f, 0.217850f, 0.218370f, 0.218900f, 0.219430f, 0.219970f, 0.220520f, 0.221080f, 0.221640f, 0.222210f,
        0.222790f, 0.223380f, 0.223990f, 0.224600f, 0.225220f, 0.225850f, 0.226480f, 0.227090f, 0.227680f, 0.228250f,
        0.228780f, 0.229270f, 0.229720f, 0.230100f, 0.230420f, 0.230700f, 0.230930f, 0.231150f, 0.231370f, 0.231590f,
        0.231840f, 0.232130f, 0.232480f, 0.232900f, 0.233400f, 0.233980f, 0.234620f, 0.235320f, 0.236060f, 0.236840f,
        0.237640f, 0.238460f, 0.239280f, 0.240100f, 0.240900f, 0.241690f, 0.242450f, 0.243200f, 0.243930f, 0.244630f,
        0.245310f, 0.245970f, 0.246600f, 0.247200f, 0.247770f, 0.248320f, 0.248830f, 0.249300f, 0.249730f, 0.250120f,
        0.250470f, 0.250760f, 0.251010f, 0.251200f, 0.251520f, 0.251790f, 0.252060f, 0.252340f, 0.252610f, 0.252880f,
        0.253150f, 0.253420f, 0.253690f, 0.253970f, 0.254240f, 0.254510f, 0.254780f, 0.255060f, 0.255330f, 0.255600f,
        0.255880f, 0.256150f, 0.256420f, 0.256700f, 0.256970f, 0.257250f, 0.257520f, 0.257800f, 0.258070f, 0.258350f,
        0.258620f, 0.258900f, 0.259180f, 0.259450f, 0.259730f, 0.260000f, 0.260280f, 0.260560f, 0.260830f, 0.261110f,
        0.261390f, 0.261670f, 0.261940f, 0.262220f, 0.262500f, 0.262780f, 0.263060f, 0.263340f, 0.263620f, 0.263890f,
        0.264170f, 0.264450f, 0.264730f, 0.265010f, 0.265290f, 0.265570f, 0.265850f, 0.266130f, 0.266420f, 0.266700f,
        0.266980f, 0.267260f, 0.267540f, 0.267820f, 0.268100f, 0.268390f, 0.268670f, 0.268950f, 0.269230f, 0.269520f,
        0.269800f, 0.270080f, 0.270370f, 0.270650f, 0.270930f, 0.271220f, 0.271500f, 0.271790f, 0.272070f, 0.272360f,
        0.272640f, 0.272930f, 0.273210f});
    static constexpr Spectrum CES91({0.029659f, 0.029735f, 0.030120f, 0.030910f, 0.032140f, 0.033792f, 0.035833f,
        0.038217f, 0.040943f, 0.044085f, 0.047769f, 0.052163f, 0.057429f, 0.063686f, 0.070954f, 0.079168f, 0.088215f,
        0.097973f, 0.108346f, 0.119276f, 0.130745f, 0.139272f, 0.148185f, 0.157455f, 0.167044f, 0.176905f, 0.186995f,
        0.197276f, 0.207725f, 0.218323f, 0.229051f, 0.239898f, 0.250861f, 0.261941f, 0.273141f, 0.284463f, 0.295915f,
        0.307493f, 0.319186f, 0.330962f, 0.342785f, 0.354613f, 0.366402f, 0.378107f, 0.389684f, 0.401096f, 0.412311f,
        0.423307f, 0.434062f, 0.444560f, 0.454787f, 0.464732f, 0.474390f, 0.483749f, 0.492788f, 0.501498f, 0.509870f,
        0.517893f, 0.525540f, 0.532791f, 0.539635f, 0.546061f, 0.552061f, 0.557638f, 0.562815f, 0.567622f, 0.572087f,
        0.576241f, 0.580111f, 0.583716f, 0.587077f, 0.590206f, 0.593118f, 0.595822f, 0.598321f, 0.600611f, 0.602682f,
        0.604519f, 0.606103f, 0.607427f, 0.608495f, 0.609324f, 0.609932f, 0.610343f, 0.610573f, 0.610639f, 0.610551f,
        0.610318f, 0.609942f, 0.609428f, 0.608792f, 0.608047f, 0.607200f, 0.606247f, 0.605195f, 0.604056f, 0.602836f,
        0.601536f, 0.600170f, 0.598766f, 0.597343f, 0.595904f, 0.594453f, 0.592999f, 0.591547f, 0.590092f, 0.588628f,
        0.587161f, 0.585694f, 0.584227f, 0.582752f, 0.581265f, 0.579765f, 0.578244f, 0.576693f, 0.575109f, 0.573493f,
        0.571853f, 0.570199f, 0.568549f, 0.566923f, 0.565338f, 0.563800f, 0.562313f, 0.560876f, 0.559486f, 0.558136f,
        0.556820f, 0.555528f, 0.554248f, 0.552966f, 0.551660f, 0.550314f, 0.548916f, 0.547467f, 0.545970f, 0.544429f,
        0.542851f, 0.541242f, 0.539602f, 0.537934f, 0.536239f, 0.534526f, 0.532809f, 0.531106f, 0.529435f, 0.527813f,
        0.526257f, 0.524776f, 0.523378f, 0.522066f, 0.520841f, 0.519699f, 0.518632f, 0.517633f, 0.516693f, 0.515799f,
        0.514938f, 0.514099f, 0.513271f, 0.512442f, 0.511602f, 0.510745f, 0.509870f, 0.508976f, 0.508068f, 0.507156f,
        0.506250f, 0.505365f, 0.504512f, 0.503709f, 0.502973f, 0.502319f, 0.501760f, 0.501308f, 0.500973f, 0.500760f,
        0.500670f, 0.500699f, 0.500843f, 0.501095f, 0.501439f, 0.501861f, 0.502344f, 0.502869f, 0.503417f, 0.503970f,
        0.504515f, 0.505041f, 0.505540f, 0.506010f, 0.506452f, 0.506871f, 0.507278f, 0.507680f, 0.508091f, 0.508523f,
        0.508990f, 0.509504f, 0.510076f, 0.510715f, 0.511431f, 0.512228f, 0.513109f, 0.514076f, 0.515127f, 0.516260f,
        0.517473f, 0.518762f, 0.520122f, 0.521545f, 0.523022f, 0.524544f, 0.526099f, 0.527673f, 0.529253f, 0.530831f,
        0.532400f, 0.533959f, 0.535507f, 0.537053f, 0.538605f, 0.540175f, 0.541775f, 0.543416f, 0.545110f, 0.546866f,
        0.548691f, 0.550587f, 0.552553f, 0.554584f, 0.556670f, 0.558802f, 0.560968f, 0.563159f, 0.565367f, 0.567590f,
        0.569826f, 0.572076f, 0.574339f, 0.576610f, 0.578875f, 0.581119f, 0.583318f, 0.585450f, 0.587491f, 0.589422f,
        0.591231f, 0.592910f, 0.594454f, 0.595864f, 0.597150f, 0.598327f, 0.599418f, 0.600450f, 0.601456f, 0.602469f,
        0.603514f, 0.604605f, 0.605749f, 0.606944f, 0.608180f, 0.609443f, 0.610717f, 0.611989f, 0.613245f, 0.614470f,
        0.615654f, 0.616783f, 0.617852f, 0.618857f, 0.619801f, 0.620699f, 0.621571f, 0.622438f, 0.623315f, 0.624212f,
        0.625132f, 0.626066f, 0.626997f, 0.627900f, 0.628765f, 0.629575f, 0.630313f, 0.630959f, 0.631500f, 0.631929f,
        0.632241f, 0.632442f, 0.632546f, 0.632582f, 0.632576f, 0.632549f, 0.632514f, 0.632478f, 0.632336f, 0.632312f,
        0.632364f, 0.632431f, 0.632434f, 0.632349f, 0.631869f, 0.631270f, 0.630558f, 0.629746f, 0.628840f, 0.627849f,
        0.626780f, 0.625644f, 0.624451f, 0.623210f, 0.621928f, 0.620607f, 0.619242f, 0.617821f, 0.616327f, 0.614749f,
        0.613078f, 0.611308f, 0.609442f, 0.607494f, 0.605480f, 0.603413f, 0.601315f, 0.599215f, 0.597151f, 0.595158f,
        0.593262f, 0.591485f, 0.589834f, 0.588300f, 0.586859f, 0.585481f, 0.584140f, 0.582809f, 0.581463f, 0.580184f,
        0.578724f, 0.577099f, 0.575349f, 0.573542f, 0.571696f, 0.570118f, 0.568539f, 0.566958f, 0.565376f, 0.563793f,
        0.562208f, 0.560622f, 0.559035f, 0.557447f, 0.555858f, 0.554267f, 0.552676f, 0.551083f, 0.549489f, 0.547894f,
        0.546299f, 0.544702f, 0.543104f, 0.541505f, 0.539906f, 0.538306f, 0.536704f, 0.535102f, 0.533500f, 0.531897f,
        0.530293f, 0.528688f, 0.527083f, 0.525478f, 0.523871f, 0.522265f, 0.520657f, 0.519050f, 0.517442f, 0.515833f,
        0.514225f, 0.512616f, 0.511007f, 0.509397f, 0.507789f, 0.506180f, 0.504570f, 0.502960f, 0.501349f, 0.499738f,
        0.498126f, 0.496514f, 0.494902f, 0.493290f, 0.491678f, 0.490066f, 0.488455f, 0.486843f, 0.485232f, 0.483620f,
        0.482010f, 0.480400f, 0.478790f});
    static constexpr Spectrum CES92({0.005260f, 0.008165f, 0.011513f, 0.014950f, 0.018986f, 0.024483f, 0.030356f,
        0.036139f, 0.041619f, 0.046919f, 0.051650f, 0.055890f, 0.060028f, 0.064418f, 0.068980f, 0.073734f, 0.078563f,
        0.083261f, 0.087780f, 0.092207f, 0.096701f, 0.101293f, 0.105973f, 0.110737f, 0.115562f, 0.120482f, 0.125563f,
        0.130789f, 0.136153f, 0.141616f, 0.147052f, 0.152431f, 0.157753f, 0.163019f, 0.168272f, 0.173554f, 0.178816f,
        0.183998f, 0.189103f, 0.194102f, 0.199025f, 0.203942f, 0.208860f, 0.213757f, 0.218606f, 0.223370f, 0.227988f,
        0.232458f, 0.236784f, 0.240956f, 0.244973f, 0.248849f, 0.252592f, 0.256212f, 0.259747f, 0.263213f, 0.266606f,
        0.269933f, 0.273181f, 0.276317f, 0.279320f, 0.282178f, 0.284872f, 0.287406f, 0.289791f, 0.292027f, 0.294074f,
        0.295895f, 0.297472f, 0.298812f, 0.299902f, 0.300780f, 0.301492f, 0.302067f, 0.302473f, 0.302696f, 0.302736f,
        0.302617f, 0.302362f, 0.301998f, 0.301543f, 0.301041f, 0.300470f, 0.299808f, 0.299062f, 0.298263f, 0.297410f,
        0.296528f, 0.295636f, 0.294738f, 0.293857f, 0.292966f, 0.292029f, 0.291046f, 0.290022f, 0.288921f, 0.287750f,
        0.286534f, 0.285250f, 0.283890f, 0.282468f, 0.280970f, 0.279409f, 0.277801f, 0.276152f, 0.274454f, 0.272715f,
        0.270918f, 0.269055f, 0.267130f, 0.265144f, 0.263112f, 0.261050f, 0.258977f, 0.256915f, 0.254876f, 0.252869f,
        0.250905f, 0.248995f, 0.247135f, 0.245322f, 0.243570f, 0.241857f, 0.240168f, 0.238485f, 0.236787f, 0.235054f,
        0.233275f, 0.231433f, 0.229528f, 0.227558f, 0.225512f, 0.223409f, 0.221278f, 0.219123f, 0.216957f, 0.214807f,
        0.212666f, 0.210529f, 0.208415f, 0.206344f, 0.204319f, 0.202366f, 0.200462f, 0.198596f, 0.196761f, 0.194964f,
        0.193191f, 0.191466f, 0.189792f, 0.188162f, 0.186562f, 0.184995f, 0.183449f, 0.181934f, 0.180448f, 0.178996f,
        0.177569f, 0.176163f, 0.174766f, 0.173397f, 0.172069f, 0.170787f, 0.169551f, 0.168376f, 0.167254f, 0.166183f,
        0.165177f, 0.164244f, 0.163390f, 0.162626f, 0.161949f, 0.161348f, 0.160827f, 0.160383f, 0.160004f, 0.159691f,
        0.159439f, 0.159231f, 0.159063f, 0.158930f, 0.158824f, 0.158743f, 0.158690f, 0.158671f, 0.158684f, 0.158729f,
        0.158797f, 0.158896f, 0.159018f, 0.159176f, 0.159361f, 0.159580f, 0.159841f, 0.160151f, 0.160484f, 0.160858f,
        0.161274f, 0.161740f, 0.162250f, 0.162830f, 0.163457f, 0.164150f, 0.164911f, 0.165732f, 0.166612f, 0.167577f,
        0.168614f, 0.169720f, 0.170917f, 0.172219f, 0.173622f, 0.175137f, 0.176774f, 0.178530f, 0.180408f, 0.182434f,
        0.184621f, 0.186987f, 0.189553f, 0.192328f, 0.195308f, 0.198500f, 0.201927f, 0.205595f, 0.209527f, 0.213749f,
        0.218264f, 0.223068f, 0.228148f, 0.233499f, 0.239095f, 0.244944f, 0.251031f, 0.257342f, 0.263842f, 0.270529f,
        0.277365f, 0.284320f, 0.291371f, 0.298504f, 0.305699f, 0.312958f, 0.320279f, 0.327681f, 0.335178f, 0.342766f,
        0.350425f, 0.358153f, 0.365920f, 0.373727f, 0.381617f, 0.389588f, 0.397642f, 0.405786f, 0.413987f, 0.422183f,
        0.430360f, 0.438497f, 0.446565f, 0.454565f, 0.462458f, 0.470204f, 0.477791f, 0.485187f, 0.492361f, 0.499291f,
        0.505980f, 0.512397f, 0.518531f, 0.524388f, 0.529981f, 0.535327f, 0.540436f, 0.545325f, 0.550008f, 0.554483f,
        0.558757f, 0.562840f, 0.566758f, 0.570493f, 0.574074f, 0.577567f, 0.581031f, 0.584494f, 0.587992f, 0.591541f,
        0.595091f, 0.598621f, 0.602144f, 0.605699f, 0.609305f, 0.612904f, 0.616428f, 0.619817f, 0.622869f, 0.625653f,
        0.628230f, 0.630680f, 0.633010f, 0.635228f, 0.637651f, 0.640030f, 0.642406f, 0.644786f, 0.647172f, 0.649577f,
        0.652004f, 0.654420f, 0.656850f, 0.659312f, 0.661817f, 0.664354f, 0.666920f, 0.669509f, 0.672117f, 0.674741f,
        0.677410f, 0.680132f, 0.682914f, 0.685752f, 0.688678f, 0.691680f, 0.694701f, 0.697691f, 0.700595f, 0.703355f,
        0.705942f, 0.708402f, 0.710759f, 0.713026f, 0.715189f, 0.717219f, 0.719117f, 0.720946f, 0.722755f, 0.724735f,
        0.726844f, 0.729082f, 0.731366f, 0.733696f, 0.736072f, 0.738168f, 0.740253f, 0.742328f, 0.744391f, 0.746444f,
        0.748485f, 0.750516f, 0.752536f, 0.754545f, 0.756543f, 0.758530f, 0.760506f, 0.762471f, 0.764424f, 0.766367f,
        0.768299f, 0.770219f, 0.772128f, 0.774026f, 0.775912f, 0.777788f, 0.779652f, 0.781505f, 0.783347f, 0.785178f,
        0.786997f, 0.788806f, 0.790602f, 0.792388f, 0.794163f, 0.795926f, 0.797678f, 0.799419f, 0.801149f, 0.802867f,
        0.804574f, 0.806270f, 0.807955f, 0.809629f, 0.811366f, 0.813088f, 0.814795f, 0.816487f, 0.818164f, 0.819826f,
        0.821472f, 0.823104f, 0.824721f, 0.826323f, 0.827911f, 0.829483f, 0.831041f, 0.832584f, 0.834113f, 0.835626f,
        0.837126f, 0.838610f, 0.840080f});
    static constexpr Spectrum CES93({0.230770f, 0.233310f, 0.235850f, 0.238400f, 0.240960f, 0.243520f, 0.246080f,
        0.248650f, 0.251220f, 0.253790f, 0.256360f, 0.258930f, 0.261500f, 0.264070f, 0.266640f, 0.269200f, 0.271760f,
        0.274320f, 0.276870f, 0.279410f, 0.281950f, 0.284480f, 0.287010f, 0.289530f, 0.292060f, 0.294600f, 0.297150f,
        0.299720f, 0.302300f, 0.304910f, 0.307540f, 0.310200f, 0.312910f, 0.315650f, 0.318450f, 0.321300f, 0.324210f,
        0.327190f, 0.330250f, 0.333390f, 0.336620f, 0.339940f, 0.343320f, 0.346740f, 0.350170f, 0.353570f, 0.356930f,
        0.360200f, 0.363360f, 0.366390f, 0.369250f, 0.371920f, 0.374410f, 0.376750f, 0.378960f, 0.381050f, 0.383050f,
        0.384980f, 0.386860f, 0.388700f, 0.390540f, 0.392380f, 0.394230f, 0.396080f, 0.397920f, 0.399750f, 0.401550f,
        0.403320f, 0.405050f, 0.406750f, 0.408390f, 0.409970f, 0.411480f, 0.412900f, 0.414210f, 0.415390f, 0.416430f,
        0.417300f, 0.418000f, 0.418490f, 0.418770f, 0.418820f, 0.418670f, 0.418330f, 0.417830f, 0.417190f, 0.416440f,
        0.415610f, 0.414710f, 0.413780f, 0.412830f, 0.411880f, 0.410930f, 0.409930f, 0.408860f, 0.407710f, 0.406440f,
        0.405030f, 0.403460f, 0.401690f, 0.399710f, 0.397500f, 0.395070f, 0.392460f, 0.389700f, 0.386820f, 0.383840f,
        0.380790f, 0.377720f, 0.374630f, 0.371570f, 0.368560f, 0.365590f, 0.362680f, 0.359820f, 0.357000f, 0.354230f,
        0.351500f, 0.348810f, 0.346170f, 0.343560f, 0.340990f, 0.338460f, 0.335960f, 0.333500f, 0.331080f, 0.328700f,
        0.326350f, 0.324040f, 0.321760f, 0.319520f, 0.317320f, 0.315150f, 0.313040f, 0.310980f, 0.308970f, 0.307030f,
        0.305160f, 0.303370f, 0.301660f, 0.300030f, 0.298490f, 0.297040f, 0.295680f, 0.294380f, 0.293160f, 0.292000f,
        0.290900f, 0.289850f, 0.288840f, 0.287880f, 0.286950f, 0.286060f, 0.285220f, 0.284410f, 0.283650f, 0.282940f,
        0.282290f, 0.281680f, 0.281140f, 0.280660f, 0.280240f, 0.279880f, 0.279580f, 0.279340f, 0.279140f, 0.278990f,
        0.278890f, 0.278820f, 0.278790f, 0.278800f, 0.278840f, 0.278910f, 0.279010f, 0.279160f, 0.279350f, 0.279590f,
        0.279880f, 0.280220f, 0.280610f, 0.281070f, 0.281590f, 0.282160f, 0.282790f, 0.283470f, 0.284190f, 0.284950f,
        0.285740f, 0.286560f, 0.287410f, 0.288270f, 0.289150f, 0.290050f, 0.290980f, 0.291940f, 0.292950f, 0.294000f,
        0.295120f, 0.296300f, 0.297560f, 0.298890f, 0.300310f, 0.301820f, 0.303410f, 0.305080f, 0.306830f, 0.308660f,
        0.310560f, 0.312540f, 0.314580f, 0.316700f, 0.318880f, 0.321130f, 0.323460f, 0.325850f, 0.328320f, 0.330870f,
        0.333500f, 0.336210f, 0.339010f, 0.341890f, 0.344860f, 0.347920f, 0.351070f, 0.354310f, 0.357640f, 0.361050f,
        0.364550f, 0.368140f, 0.371820f, 0.375580f, 0.379430f, 0.383360f, 0.387380f, 0.391480f, 0.395670f, 0.399930f,
        0.404290f, 0.408720f, 0.413230f, 0.417830f, 0.422510f, 0.427260f, 0.432110f, 0.437030f, 0.442040f, 0.447140f,
        0.452320f, 0.457590f, 0.462950f, 0.468400f, 0.473940f, 0.479560f, 0.485250f, 0.491000f, 0.496790f, 0.502620f,
        0.508470f, 0.514340f, 0.520200f, 0.526060f, 0.531900f, 0.537710f, 0.543500f, 0.549270f, 0.555020f, 0.560740f,
        0.566430f, 0.572100f, 0.577750f, 0.583370f, 0.588960f, 0.594520f, 0.600030f, 0.605490f, 0.610880f, 0.616210f,
        0.621460f, 0.626610f, 0.631680f, 0.636630f, 0.641470f, 0.646200f, 0.650810f, 0.655300f, 0.659680f, 0.663950f,
        0.668090f, 0.672130f, 0.676040f, 0.679840f, 0.683520f, 0.687090f, 0.690540f, 0.693870f, 0.697090f, 0.700200f,
        0.703190f, 0.706070f, 0.708830f, 0.711480f, 0.714020f, 0.716450f, 0.718790f, 0.721040f, 0.723210f, 0.725300f,
        0.727330f, 0.729300f, 0.731210f, 0.733090f, 0.734930f, 0.736730f, 0.738490f, 0.740220f, 0.741890f, 0.743530f,
        0.745110f, 0.746650f, 0.748140f, 0.749570f, 0.750950f, 0.752260f, 0.753520f, 0.754700f, 0.755810f, 0.756850f,
        0.757800f, 0.758670f, 0.759450f, 0.760130f, 0.760720f, 0.761210f, 0.761630f, 0.761970f, 0.762260f, 0.762490f,
        0.762670f, 0.762820f, 0.762950f, 0.763060f, 0.763160f, 0.763270f, 0.763390f, 0.763520f, 0.763690f, 0.763900f,
        0.764150f, 0.764470f, 0.764840f, 0.765300f, 0.764990f, 0.765170f, 0.765350f, 0.765520f, 0.765700f, 0.765880f,
        0.766050f, 0.766230f, 0.766400f, 0.766580f, 0.766750f, 0.766930f, 0.767110f, 0.767280f, 0.767460f, 0.767630f,
        0.767810f, 0.767980f, 0.768160f, 0.768330f, 0.768500f, 0.768680f, 0.768850f, 0.769030f, 0.769200f, 0.769380f,
        0.769550f, 0.769720f, 0.769900f, 0.770070f, 0.770250f, 0.770420f, 0.770590f, 0.770770f, 0.770940f, 0.771110f,
        0.771290f, 0.771460f, 0.771630f, 0.771810f, 0.771980f, 0.772150f, 0.772320f, 0.772500f, 0.772670f, 0.772840f,
        0.773010f, 0.773190f, 0.773360f});
    static constexpr Spectrum CES94({0.061017f, 0.062174f, 0.063351f, 0.064549f, 0.065767f, 0.067008f, 0.068269f,
        0.069553f, 0.070859f, 0.072188f, 0.073540f, 0.074915f, 0.076314f, 0.077736f, 0.079183f, 0.080654f, 0.082150f,
        0.083672f, 0.085219f, 0.086792f, 0.088100f, 0.089822f, 0.091571f, 0.093344f, 0.095136f, 0.096944f, 0.098762f,
        0.100590f, 0.102420f, 0.104240f, 0.106070f, 0.107880f, 0.109680f, 0.111460f, 0.113230f, 0.114960f, 0.116670f,
        0.118340f, 0.119980f, 0.121570f, 0.123120f, 0.124610f, 0.126070f, 0.127470f, 0.128840f, 0.130160f, 0.131450f,
        0.132700f, 0.133920f, 0.135120f, 0.136290f, 0.137430f, 0.138540f, 0.139610f, 0.140640f, 0.141620f, 0.142540f,
        0.143390f, 0.144170f, 0.144860f, 0.145470f, 0.145990f, 0.146410f, 0.146750f, 0.147010f, 0.147180f, 0.147280f,
        0.147310f, 0.147270f, 0.147170f, 0.147000f, 0.146780f, 0.146500f, 0.146160f, 0.145760f, 0.145280f, 0.144740f,
        0.144120f, 0.143430f, 0.142650f, 0.141800f, 0.140860f, 0.139840f, 0.138750f, 0.137590f, 0.136370f, 0.135090f,
        0.133770f, 0.132400f, 0.130990f, 0.129550f, 0.128080f, 0.126590f, 0.125060f, 0.123520f, 0.121940f, 0.120330f,
        0.118700f, 0.117040f, 0.115340f, 0.113620f, 0.111870f, 0.110100f, 0.108320f, 0.106530f, 0.104760f, 0.103000f,
        0.101280f, 0.099590f, 0.097951f, 0.096369f, 0.094854f, 0.093401f, 0.092006f, 0.090664f, 0.089368f, 0.088113f,
        0.086895f, 0.085708f, 0.084546f, 0.083405f, 0.082279f, 0.081170f, 0.080077f, 0.079003f, 0.077947f, 0.076912f,
        0.075898f, 0.074905f, 0.073936f, 0.072992f, 0.072072f, 0.071177f, 0.070308f, 0.069464f, 0.068645f, 0.067851f,
        0.067082f, 0.066339f, 0.065620f, 0.064927f, 0.064259f, 0.063616f, 0.063001f, 0.062413f, 0.061854f, 0.061324f,
        0.060824f, 0.060356f, 0.059920f, 0.059516f, 0.059146f, 0.058810f, 0.058507f, 0.058236f, 0.057997f, 0.057790f,
        0.057614f, 0.057469f, 0.057355f, 0.057270f, 0.057216f, 0.057190f, 0.057191f, 0.057219f, 0.057272f, 0.057350f,
        0.057451f, 0.057574f, 0.057718f, 0.057883f, 0.058066f, 0.058268f, 0.058485f, 0.058717f, 0.058962f, 0.059218f,
        0.059485f, 0.059761f, 0.060044f, 0.060333f, 0.060627f, 0.060926f, 0.061233f, 0.061548f, 0.061874f, 0.062211f,
        0.062561f, 0.062925f, 0.063305f, 0.063702f, 0.064118f, 0.064553f, 0.065010f, 0.065489f, 0.065992f, 0.066520f,
        0.067074f, 0.067656f, 0.068267f, 0.068908f, 0.069581f, 0.070288f, 0.071032f, 0.071814f, 0.072639f, 0.073509f,
        0.074425f, 0.075391f, 0.076410f, 0.077484f, 0.078615f, 0.079809f, 0.081070f, 0.082403f, 0.083813f, 0.085305f,
        0.086883f, 0.088552f, 0.090317f, 0.092184f, 0.094155f, 0.096234f, 0.098419f, 0.100710f, 0.103110f, 0.105630f,
        0.108250f, 0.110980f, 0.113830f, 0.116790f, 0.119860f, 0.123050f, 0.126370f, 0.129820f, 0.133410f, 0.137140f,
        0.141030f, 0.145070f, 0.149270f, 0.153640f, 0.158180f, 0.162890f, 0.167760f, 0.172770f, 0.177910f, 0.183180f,
        0.188560f, 0.194050f, 0.199630f, 0.205300f, 0.211040f, 0.216850f, 0.222750f, 0.228720f, 0.234770f, 0.240900f,
        0.247110f, 0.253400f, 0.259780f, 0.266240f, 0.272790f, 0.279420f, 0.286130f, 0.292910f, 0.299750f, 0.306670f,
        0.313630f, 0.320660f, 0.327730f, 0.334840f, 0.342000f, 0.349190f, 0.356420f, 0.363670f, 0.370960f, 0.378270f,
        0.385610f, 0.392970f, 0.400340f, 0.407730f, 0.415140f, 0.422550f, 0.429960f, 0.437360f, 0.444750f, 0.452130f,
        0.459480f, 0.466810f, 0.474090f, 0.481340f, 0.488530f, 0.495680f, 0.502790f, 0.509850f, 0.516860f, 0.523830f,
        0.530750f, 0.537630f, 0.544470f, 0.551270f, 0.558020f, 0.564730f, 0.571390f, 0.578000f, 0.584550f, 0.591040f,
        0.597460f, 0.603810f, 0.610090f, 0.616300f, 0.622420f, 0.628450f, 0.634400f, 0.640250f, 0.646000f, 0.651650f,
        0.657200f, 0.662630f, 0.667950f, 0.673160f, 0.679550f, 0.685040f, 0.690490f, 0.695880f, 0.701210f, 0.706500f,
        0.711720f, 0.716900f, 0.722010f, 0.727070f, 0.732070f, 0.737010f, 0.741900f, 0.746720f, 0.751480f, 0.756190f,
        0.760830f, 0.765410f, 0.769930f, 0.774390f, 0.778790f, 0.783120f, 0.787390f, 0.791610f, 0.795760f, 0.799850f,
        0.803870f, 0.807840f, 0.811740f, 0.815580f, 0.819370f, 0.823090f, 0.826750f, 0.830350f, 0.833890f, 0.837370f,
        0.840790f, 0.844150f, 0.847460f, 0.850710f, 0.853900f, 0.857030f, 0.860110f, 0.863130f, 0.866100f, 0.869010f,
        0.871870f, 0.874670f, 0.877420f, 0.880120f, 0.882770f, 0.885370f, 0.887920f, 0.890410f, 0.892860f, 0.895260f,
        0.897620f, 0.899920f, 0.902180f, 0.904400f, 0.906570f, 0.908690f, 0.910770f, 0.912810f, 0.914810f, 0.916760f,
        0.918670f, 0.920550f, 0.922380f, 0.924180f, 0.925930f, 0.927650f, 0.929340f, 0.930980f, 0.932590f, 0.934170f,
        0.935710f, 0.937220f, 0.938690f});
    static constexpr Spectrum CES95({0.073718f, 0.076853f, 0.080110f, 0.083493f, 0.087005f, 0.090650f, 0.094431f,
        0.098354f, 0.102420f, 0.106640f, 0.111000f, 0.115530f, 0.120210f, 0.125050f, 0.130070f, 0.135250f, 0.140600f,
        0.146130f, 0.151840f, 0.157740f, 0.164280f, 0.170140f, 0.176340f, 0.182870f, 0.189690f, 0.196780f, 0.204120f,
        0.211680f, 0.219430f, 0.227350f, 0.235410f, 0.243590f, 0.251860f, 0.260190f, 0.268560f, 0.276950f, 0.285320f,
        0.293650f, 0.301920f, 0.310100f, 0.318160f, 0.326090f, 0.333880f, 0.341530f, 0.349040f, 0.356430f, 0.363680f,
        0.370800f, 0.377800f, 0.384680f, 0.391440f, 0.398070f, 0.404550f, 0.410830f, 0.416890f, 0.422670f, 0.428150f,
        0.433290f, 0.438040f, 0.442380f, 0.446270f, 0.449670f, 0.452610f, 0.455130f, 0.457270f, 0.459050f, 0.460510f,
        0.461700f, 0.462630f, 0.463350f, 0.463900f, 0.464300f, 0.464560f, 0.464700f, 0.464720f, 0.464640f, 0.464450f,
        0.464180f, 0.463810f, 0.463380f, 0.462880f, 0.462310f, 0.461680f, 0.460970f, 0.460170f, 0.459260f, 0.458240f,
        0.457090f, 0.455810f, 0.454380f, 0.452790f, 0.451030f, 0.449130f, 0.447090f, 0.444940f, 0.442690f, 0.440370f,
        0.437990f, 0.435560f, 0.433120f, 0.430670f, 0.428240f, 0.425820f, 0.423410f, 0.421010f, 0.418620f, 0.416220f,
        0.413830f, 0.411440f, 0.409030f, 0.406620f, 0.404200f, 0.401770f, 0.399320f, 0.396860f, 0.394370f, 0.391870f,
        0.389350f, 0.386810f, 0.384250f, 0.381660f, 0.379040f, 0.376400f, 0.373730f, 0.371040f, 0.368320f, 0.365590f,
        0.362830f, 0.360060f, 0.357260f, 0.354450f, 0.351610f, 0.348780f, 0.345940f, 0.343120f, 0.340320f, 0.337560f,
        0.334840f, 0.332170f, 0.329560f, 0.327030f, 0.324580f, 0.322220f, 0.319930f, 0.317730f, 0.315620f, 0.313580f,
        0.311620f, 0.309750f, 0.307960f, 0.306240f, 0.304610f, 0.303050f, 0.301560f, 0.300150f, 0.298800f, 0.297520f,
        0.296300f, 0.295140f, 0.294040f, 0.292990f, 0.292000f, 0.291050f, 0.290150f, 0.289290f, 0.288480f, 0.287700f,
        0.286970f, 0.286260f, 0.285590f, 0.284940f, 0.284330f, 0.283740f, 0.283200f, 0.282710f, 0.282270f, 0.281900f,
        0.281590f, 0.281360f, 0.281220f, 0.281170f, 0.281220f, 0.281360f, 0.281610f, 0.281950f, 0.282380f, 0.282910f,
        0.283530f, 0.284250f, 0.285060f, 0.285960f, 0.286950f, 0.288050f, 0.289260f, 0.290590f, 0.292060f, 0.293680f,
        0.295460f, 0.297400f, 0.299530f, 0.301860f, 0.304390f, 0.307120f, 0.310070f, 0.313240f, 0.316620f, 0.320220f,
        0.324050f, 0.328110f, 0.332390f, 0.336920f, 0.341670f, 0.346660f, 0.351870f, 0.357300f, 0.362930f, 0.368770f,
        0.374800f, 0.381020f, 0.387420f, 0.393990f, 0.400720f, 0.407610f, 0.414650f, 0.421840f, 0.429150f, 0.436600f,
        0.444160f, 0.451830f, 0.459600f, 0.467460f, 0.475410f, 0.483420f, 0.491460f, 0.499490f, 0.507500f, 0.515450f,
        0.523320f, 0.531080f, 0.538690f, 0.546140f, 0.553390f, 0.560460f, 0.567360f, 0.574090f, 0.580680f, 0.587140f,
        0.593470f, 0.599690f, 0.605820f, 0.611870f, 0.617850f, 0.623760f, 0.629620f, 0.635430f, 0.641190f, 0.646930f,
        0.652630f, 0.658320f, 0.663990f, 0.669650f, 0.675320f, 0.680980f, 0.686640f, 0.692300f, 0.697960f, 0.703610f,
        0.709250f, 0.714880f, 0.720500f, 0.726110f, 0.731710f, 0.737280f, 0.742800f, 0.748270f, 0.753670f, 0.758970f,
        0.764170f, 0.769260f, 0.774200f, 0.779000f, 0.783640f, 0.788110f, 0.792410f, 0.796530f, 0.800480f, 0.804240f,
        0.807820f, 0.811220f, 0.814420f, 0.817420f, 0.820240f, 0.822870f, 0.825340f, 0.827670f, 0.829880f, 0.831990f,
        0.834010f, 0.835960f, 0.837870f, 0.839740f, 0.841610f, 0.843450f, 0.845270f, 0.847060f, 0.848810f, 0.850500f,
        0.852150f, 0.853720f, 0.855230f, 0.856660f, 0.858000f, 0.859250f, 0.860400f, 0.861430f, 0.862350f, 0.863140f,
        0.863810f, 0.864330f, 0.864700f, 0.864920f, 0.867590f, 0.868580f, 0.869560f, 0.870540f, 0.871510f, 0.872470f,
        0.873430f, 0.874390f, 0.875330f, 0.876270f, 0.877210f, 0.878140f, 0.879060f, 0.879980f, 0.880890f, 0.881790f,
        0.882690f, 0.883580f, 0.884470f, 0.885350f, 0.886220f, 0.887090f, 0.887960f, 0.888810f, 0.889660f, 0.890510f,
        0.891350f, 0.892190f, 0.893010f, 0.893840f, 0.894660f, 0.895470f, 0.896270f, 0.897080f, 0.897870f, 0.898660f,
        0.899450f, 0.900230f, 0.901000f, 0.901770f, 0.902530f, 0.903290f, 0.904040f, 0.904790f, 0.905530f, 0.906270f,
        0.907000f, 0.907730f, 0.908450f, 0.909170f, 0.909880f, 0.910580f, 0.911290f, 0.911980f, 0.912670f, 0.913360f,
        0.914040f, 0.914720f, 0.915390f, 0.916060f, 0.916720f, 0.917380f, 0.918030f, 0.918680f, 0.919320f, 0.919960f,
        0.920600f, 0.921230f, 0.921850f, 0.922470f, 0.923090f, 0.923700f, 0.924310f, 0.924910f, 0.925510f, 0.926100f,
        0.926690f, 0.927280f, 0.927860f});
    static constexpr Spectrum CES96({0.338510f, 0.343040f, 0.347600f, 0.352190f, 0.356810f, 0.361450f, 0.366120f,
        0.370810f, 0.375530f, 0.380270f, 0.385040f, 0.389830f, 0.394630f, 0.399460f, 0.404310f, 0.409180f, 0.414070f,
        0.418970f, 0.423890f, 0.428820f, 0.433940f, 0.438690f, 0.443550f, 0.448500f, 0.453520f, 0.458590f, 0.463690f,
        0.468810f, 0.473930f, 0.479020f, 0.484080f, 0.489070f, 0.493990f, 0.498810f, 0.503510f, 0.508080f, 0.512500f,
        0.516760f, 0.520820f, 0.524670f, 0.528300f, 0.531700f, 0.534880f, 0.537870f, 0.540700f, 0.543400f, 0.546000f,
        0.548520f, 0.551000f, 0.553450f, 0.555920f, 0.558420f, 0.560920f, 0.563390f, 0.565790f, 0.568100f, 0.570280f,
        0.572300f, 0.574110f, 0.575700f, 0.577020f, 0.578050f, 0.578800f, 0.579290f, 0.579550f, 0.579600f, 0.579450f,
        0.579120f, 0.578650f, 0.578040f, 0.577320f, 0.576520f, 0.575620f, 0.574640f, 0.573590f, 0.572450f, 0.571240f,
        0.569960f, 0.568610f, 0.567190f, 0.565710f, 0.564160f, 0.562550f, 0.560880f, 0.559130f, 0.557310f, 0.555410f,
        0.553430f, 0.551360f, 0.549200f, 0.546950f, 0.544610f, 0.542180f, 0.539660f, 0.537080f, 0.534430f, 0.531720f,
        0.528960f, 0.526160f, 0.523320f, 0.520460f, 0.517580f, 0.514680f, 0.511770f, 0.508850f, 0.505930f, 0.503010f,
        0.500100f, 0.497190f, 0.494290f, 0.491410f, 0.488550f, 0.485710f, 0.482880f, 0.480080f, 0.477290f, 0.474520f,
        0.471760f, 0.469020f, 0.466300f, 0.463590f, 0.460900f, 0.458220f, 0.455550f, 0.452890f, 0.450250f, 0.447620f,
        0.445000f, 0.442390f, 0.439790f, 0.437200f, 0.434620f, 0.432050f, 0.429520f, 0.427030f, 0.424590f, 0.422220f,
        0.419910f, 0.417690f, 0.415570f, 0.413550f, 0.411650f, 0.409860f, 0.408190f, 0.406630f, 0.405170f, 0.403830f,
        0.402590f, 0.401450f, 0.400420f, 0.399490f, 0.398660f, 0.397920f, 0.397270f, 0.396720f, 0.396250f, 0.395870f,
        0.395560f, 0.395340f, 0.395190f, 0.395110f, 0.395100f, 0.395160f, 0.395290f, 0.395490f, 0.395750f, 0.396090f,
        0.396500f, 0.396990f, 0.397540f, 0.398160f, 0.398860f, 0.399640f, 0.400490f, 0.401420f, 0.402440f, 0.403560f,
        0.404760f, 0.406060f, 0.407460f, 0.408970f, 0.410580f, 0.412290f, 0.414100f, 0.416010f, 0.418000f, 0.420080f,
        0.422240f, 0.424470f, 0.426780f, 0.429150f, 0.431580f, 0.434080f, 0.436650f, 0.439310f, 0.442060f, 0.444900f,
        0.447850f, 0.450910f, 0.454080f, 0.457370f, 0.460800f, 0.464350f, 0.468030f, 0.471830f, 0.475740f, 0.479760f,
        0.483900f, 0.488140f, 0.492480f, 0.496920f, 0.501450f, 0.506080f, 0.510800f, 0.515620f, 0.520530f, 0.525540f,
        0.530650f, 0.535850f, 0.541150f, 0.546550f, 0.552040f, 0.557630f, 0.563320f, 0.569110f, 0.575000f, 0.580990f,
        0.587080f, 0.593270f, 0.599560f, 0.605960f, 0.612460f, 0.619040f, 0.625680f, 0.632370f, 0.639090f, 0.645820f,
        0.652530f, 0.659210f, 0.665840f, 0.672410f, 0.678880f, 0.685270f, 0.691560f, 0.697770f, 0.703870f, 0.709880f,
        0.715800f, 0.721610f, 0.727330f, 0.732940f, 0.738450f, 0.743860f, 0.749160f, 0.754350f, 0.759430f, 0.764410f,
        0.769270f, 0.774020f, 0.778660f, 0.783180f, 0.787590f, 0.791880f, 0.796050f, 0.800110f, 0.804050f, 0.807870f,
        0.811580f, 0.815170f, 0.818650f, 0.822010f, 0.825260f, 0.828390f, 0.831390f, 0.834280f, 0.837040f, 0.839680f,
        0.842190f, 0.844560f, 0.846810f, 0.848920f, 0.850890f, 0.852730f, 0.854440f, 0.856030f, 0.857510f, 0.858880f,
        0.860150f, 0.861320f, 0.862400f, 0.863390f, 0.864300f, 0.865150f, 0.865950f, 0.866730f, 0.867490f, 0.868260f,
        0.869050f, 0.869890f, 0.870780f, 0.871740f, 0.872790f, 0.873920f, 0.875110f, 0.876350f, 0.877610f, 0.878900f,
        0.880180f, 0.881440f, 0.882680f, 0.883870f, 0.885000f, 0.886050f, 0.887010f, 0.887870f, 0.888600f, 0.889190f,
        0.889630f, 0.889900f, 0.889990f, 0.889880f, 0.892460f, 0.893240f, 0.894010f, 0.894780f, 0.895540f, 0.896300f,
        0.897050f, 0.897800f, 0.898540f, 0.899280f, 0.900010f, 0.900740f, 0.901460f, 0.902180f, 0.902900f, 0.903610f,
        0.904310f, 0.905010f, 0.905710f, 0.906400f, 0.907090f, 0.907770f, 0.908450f, 0.909120f, 0.909790f, 0.910450f,
        0.911110f, 0.911770f, 0.912420f, 0.913070f, 0.913710f, 0.914350f, 0.914980f, 0.915610f, 0.916240f, 0.916860f,
        0.917480f, 0.918090f, 0.918700f, 0.919300f, 0.919900f, 0.920500f, 0.921090f, 0.921680f, 0.922260f, 0.922840f,
        0.923420f, 0.923990f, 0.924560f, 0.925130f, 0.925690f, 0.926240f, 0.926800f, 0.927350f, 0.927890f, 0.928430f,
        0.928970f, 0.929500f, 0.930040f, 0.930560f, 0.931080f, 0.931600f, 0.932120f, 0.932630f, 0.933140f, 0.933650f,
        0.934150f, 0.934640f, 0.935140f, 0.935630f, 0.936120f, 0.936600f, 0.937080f, 0.937560f, 0.938030f, 0.938500f,
        0.938970f, 0.939430f, 0.939890f});
    static constexpr Spectrum CES97({0.119310f, 0.121730f, 0.124190f, 0.126690f, 0.129240f, 0.131830f, 0.134460f,
        0.137130f, 0.139860f, 0.142620f, 0.145440f, 0.148290f, 0.151200f, 0.154150f, 0.157150f, 0.160190f, 0.163280f,
        0.166420f, 0.169610f, 0.172850f, 0.176600f, 0.179190f, 0.182410f, 0.186120f, 0.190200f, 0.194510f, 0.198910f,
        0.203260f, 0.207440f, 0.211290f, 0.214700f, 0.217550f, 0.219880f, 0.221750f, 0.223230f, 0.224380f, 0.225270f,
        0.225960f, 0.226520f, 0.227010f, 0.227500f, 0.228040f, 0.228630f, 0.229270f, 0.229940f, 0.230630f, 0.231350f,
        0.232070f, 0.232790f, 0.233510f, 0.234200f, 0.234870f, 0.235500f, 0.236080f, 0.236610f, 0.237070f, 0.237460f,
        0.237770f, 0.237980f, 0.238100f, 0.238100f, 0.237980f, 0.237750f, 0.237380f, 0.236880f, 0.236250f, 0.235470f,
        0.234550f, 0.233490f, 0.232270f, 0.230900f, 0.229370f, 0.227690f, 0.225870f, 0.223910f, 0.221820f, 0.219620f,
        0.217310f, 0.214900f, 0.212390f, 0.209800f, 0.207130f, 0.204410f, 0.201620f, 0.198800f, 0.195950f, 0.193090f,
        0.190220f, 0.187350f, 0.184510f, 0.181700f, 0.178930f, 0.176200f, 0.173490f, 0.170810f, 0.168140f, 0.165490f,
        0.162830f, 0.160170f, 0.157490f, 0.154800f, 0.152080f, 0.149330f, 0.146560f, 0.143760f, 0.140930f, 0.138080f,
        0.135200f, 0.132290f, 0.129360f, 0.126400f, 0.123420f, 0.120430f, 0.117450f, 0.114520f, 0.111630f, 0.108830f,
        0.106120f, 0.103540f, 0.101090f, 0.098800f, 0.096689f, 0.094752f, 0.092984f, 0.091378f, 0.089928f, 0.088628f,
        0.087471f, 0.086452f, 0.085564f, 0.084800f, 0.084153f, 0.083603f, 0.083132f, 0.082718f, 0.082341f, 0.081981f,
        0.081617f, 0.081229f, 0.080797f, 0.080300f, 0.079721f, 0.079060f, 0.078317f, 0.077494f, 0.076594f, 0.075619f,
        0.074570f, 0.073449f, 0.072258f, 0.071000f, 0.069679f, 0.068312f, 0.066920f, 0.065524f, 0.064144f, 0.062800f,
        0.061513f, 0.060304f, 0.059193f, 0.058200f, 0.057350f, 0.056684f, 0.056246f, 0.056080f, 0.056231f, 0.056744f,
        0.057663f, 0.059032f, 0.060896f, 0.063300f, 0.066280f, 0.069839f, 0.073974f, 0.078681f, 0.083956f, 0.089794f,
        0.096192f, 0.103140f, 0.110650f, 0.118700f, 0.127280f, 0.136300f, 0.145660f, 0.155250f, 0.164980f, 0.174750f,
        0.184440f, 0.193970f, 0.203220f, 0.212100f, 0.220520f, 0.228450f, 0.235890f, 0.242830f, 0.249270f, 0.255180f,
        0.260570f, 0.265420f, 0.269740f, 0.273500f, 0.276710f, 0.279400f, 0.281600f, 0.283340f, 0.284660f, 0.285590f,
        0.286160f, 0.286420f, 0.286380f, 0.286100f, 0.285590f, 0.284890f, 0.284010f, 0.282980f, 0.281810f, 0.280530f,
        0.279160f, 0.277710f, 0.276220f, 0.274700f, 0.273170f, 0.271650f, 0.270140f, 0.268650f, 0.267210f, 0.265810f,
        0.264470f, 0.263190f, 0.262000f, 0.260900f, 0.259900f, 0.258990f, 0.258190f, 0.257470f, 0.256850f, 0.256320f,
        0.255870f, 0.255500f, 0.255210f, 0.255000f, 0.254860f, 0.254800f, 0.254790f, 0.254850f, 0.254960f, 0.255130f,
        0.255340f, 0.255590f, 0.255880f, 0.256200f, 0.256550f, 0.256930f, 0.257360f, 0.257830f, 0.258350f, 0.258940f,
        0.259590f, 0.260310f, 0.261110f, 0.262000f, 0.262980f, 0.264050f, 0.265200f, 0.266450f, 0.267780f, 0.269200f,
        0.270700f, 0.272290f, 0.273950f, 0.275700f, 0.277520f, 0.279410f, 0.281340f, 0.283290f, 0.285250f, 0.287190f,
        0.289110f, 0.290980f, 0.292780f, 0.294500f, 0.296120f, 0.297620f, 0.299010f, 0.300280f, 0.301420f, 0.302420f,
        0.303270f, 0.303970f, 0.304520f, 0.304900f, 0.305110f, 0.305170f, 0.305070f, 0.304850f, 0.304500f, 0.304040f,
        0.303490f, 0.302860f, 0.302160f, 0.301400f, 0.300590f, 0.299730f, 0.298810f, 0.297810f, 0.296730f, 0.295550f,
        0.294270f, 0.292870f, 0.291350f, 0.289700f, 0.287910f, 0.286020f, 0.284080f, 0.282130f, 0.280210f, 0.278380f,
        0.276680f, 0.275150f, 0.273840f, 0.272800f, 0.271190f, 0.269800f, 0.268420f, 0.267050f, 0.265680f, 0.264310f,
        0.262950f, 0.261590f, 0.260240f, 0.258890f, 0.257550f, 0.256210f, 0.254880f, 0.253550f, 0.252220f, 0.250900f,
        0.249580f, 0.248270f, 0.246970f, 0.245660f, 0.244370f, 0.243070f, 0.241790f, 0.240500f, 0.239220f, 0.237950f,
        0.236680f, 0.235420f, 0.234150f, 0.232900f, 0.231650f, 0.230400f, 0.229160f, 0.227920f, 0.226690f, 0.225470f,
        0.224240f, 0.223030f, 0.221810f, 0.220600f, 0.219400f, 0.218200f, 0.217010f, 0.215820f, 0.214630f, 0.213450f,
        0.212280f, 0.211110f, 0.209940f, 0.208780f, 0.207630f, 0.206470f, 0.205330f, 0.204180f, 0.203050f, 0.201910f,
        0.200790f, 0.199660f, 0.198550f, 0.197430f, 0.196320f, 0.195220f, 0.194120f, 0.193020f, 0.191930f, 0.190850f,
        0.189770f, 0.188690f, 0.187620f, 0.186550f, 0.185490f, 0.184430f, 0.183380f, 0.182330f, 0.181290f, 0.180250f,
        0.179220f, 0.178190f, 0.177160f});
    static constexpr Spectrum CES98({0.274250f, 0.275980f, 0.277720f, 0.279460f, 0.281210f, 0.282970f, 0.284730f,
        0.286500f, 0.288280f, 0.290060f, 0.291850f, 0.293650f, 0.295450f, 0.297260f, 0.299080f, 0.300900f, 0.302730f,
        0.304560f, 0.306400f, 0.308250f, 0.310400f, 0.311790f, 0.313540f, 0.315580f, 0.317850f, 0.320280f, 0.322820f,
        0.325400f, 0.327960f, 0.330440f, 0.332770f, 0.334910f, 0.336870f, 0.338660f, 0.340310f, 0.341840f, 0.343260f,
        0.344610f, 0.345890f, 0.347130f, 0.348360f, 0.349580f, 0.350790f, 0.351980f, 0.353130f, 0.354240f, 0.355270f,
        0.356240f, 0.357110f, 0.357890f, 0.358550f, 0.359080f, 0.359480f, 0.359730f, 0.359820f, 0.359750f, 0.359490f,
        0.359050f, 0.358410f, 0.357570f, 0.356500f, 0.355210f, 0.353700f, 0.351970f, 0.350030f, 0.347890f, 0.345560f,
        0.343040f, 0.340330f, 0.337450f, 0.334400f, 0.331190f, 0.327830f, 0.324330f, 0.320720f, 0.316990f, 0.313180f,
        0.309290f, 0.305330f, 0.301320f, 0.297280f, 0.293210f, 0.289130f, 0.285040f, 0.280950f, 0.276870f, 0.272790f,
        0.268730f, 0.264700f, 0.260700f, 0.256730f, 0.252810f, 0.248940f, 0.245120f, 0.241370f, 0.237680f, 0.234070f,
        0.230540f, 0.227100f, 0.223740f, 0.220490f, 0.217340f, 0.214280f, 0.211310f, 0.208430f, 0.205610f, 0.202860f,
        0.200160f, 0.197500f, 0.194890f, 0.192300f, 0.189740f, 0.187190f, 0.184670f, 0.182180f, 0.179700f, 0.177240f,
        0.174810f, 0.172400f, 0.170010f, 0.167650f, 0.165310f, 0.163010f, 0.160750f, 0.158550f, 0.156420f, 0.154380f,
        0.152430f, 0.150580f, 0.148860f, 0.147270f, 0.145820f, 0.144510f, 0.143330f, 0.142280f, 0.141360f, 0.140560f,
        0.139880f, 0.139310f, 0.138850f, 0.138490f, 0.138230f, 0.138060f, 0.137950f, 0.137900f, 0.137880f, 0.137890f,
        0.137900f, 0.137900f, 0.137870f, 0.137810f, 0.137690f, 0.137520f, 0.137310f, 0.137070f, 0.136790f, 0.136490f,
        0.136170f, 0.135840f, 0.135510f, 0.135180f, 0.134860f, 0.134570f, 0.134330f, 0.134150f, 0.134050f, 0.134050f,
        0.134180f, 0.134430f, 0.134840f, 0.135430f, 0.136210f, 0.137190f, 0.138410f, 0.139890f, 0.141630f, 0.143670f,
        0.146030f, 0.148720f, 0.151770f, 0.155190f, 0.159010f, 0.163230f, 0.167860f, 0.172900f, 0.178350f, 0.184240f,
        0.190550f, 0.197290f, 0.204470f, 0.212100f, 0.220170f, 0.228660f, 0.237530f, 0.246750f, 0.256290f, 0.266120f,
        0.276200f, 0.286510f, 0.297000f, 0.307650f, 0.318420f, 0.329280f, 0.340200f, 0.351130f, 0.362050f, 0.372920f,
        0.383710f, 0.394370f, 0.404880f, 0.415200f, 0.425300f, 0.435140f, 0.444680f, 0.453900f, 0.462760f, 0.471230f,
        0.479270f, 0.486840f, 0.493920f, 0.500470f, 0.506470f, 0.511930f, 0.516880f, 0.521370f, 0.525410f, 0.529040f,
        0.532290f, 0.535190f, 0.537770f, 0.540060f, 0.542090f, 0.543870f, 0.545440f, 0.546790f, 0.547960f, 0.548950f,
        0.549790f, 0.550480f, 0.551060f, 0.551540f, 0.551930f, 0.552260f, 0.552530f, 0.552790f, 0.553030f, 0.553290f,
        0.553580f, 0.553920f, 0.554330f, 0.554840f, 0.555450f, 0.556180f, 0.557040f, 0.558010f, 0.559120f, 0.560360f,
        0.561740f, 0.563270f, 0.564940f, 0.566770f, 0.568760f, 0.570910f, 0.573220f, 0.575710f, 0.578380f, 0.581230f,
        0.584270f, 0.587490f, 0.590920f, 0.594540f, 0.598370f, 0.602390f, 0.606590f, 0.610950f, 0.615480f, 0.620140f,
        0.624930f, 0.629840f, 0.634860f, 0.639960f, 0.645140f, 0.650400f, 0.655720f, 0.661110f, 0.666550f, 0.672040f,
        0.677580f, 0.683150f, 0.688760f, 0.694390f, 0.700040f, 0.705700f, 0.711350f, 0.716970f, 0.722550f, 0.728070f,
        0.733520f, 0.738880f, 0.744120f, 0.749250f, 0.754240f, 0.759110f, 0.763850f, 0.768500f, 0.773060f, 0.777530f,
        0.781940f, 0.786300f, 0.790610f, 0.794890f, 0.799140f, 0.803320f, 0.807390f, 0.811280f, 0.814940f, 0.818330f,
        0.821400f, 0.824090f, 0.826360f, 0.828140f, 0.830940f, 0.833300f, 0.835630f, 0.837940f, 0.840220f, 0.842480f,
        0.844710f, 0.846910f, 0.849090f, 0.851250f, 0.853370f, 0.855480f, 0.857550f, 0.859610f, 0.861630f, 0.863640f,
        0.865620f, 0.867570f, 0.869500f, 0.871410f, 0.873290f, 0.875150f, 0.876980f, 0.878790f, 0.880580f, 0.882350f,
        0.884090f, 0.885810f, 0.887510f, 0.889190f, 0.890840f, 0.892480f, 0.894090f, 0.895680f, 0.897250f, 0.898800f,
        0.900320f, 0.901830f, 0.903320f, 0.904780f, 0.906230f, 0.907660f, 0.909060f, 0.910450f, 0.911820f, 0.913170f,
        0.914500f, 0.915810f, 0.917110f, 0.918380f, 0.919640f, 0.920880f, 0.922100f, 0.923310f, 0.924500f, 0.925670f,
        0.926820f, 0.927960f, 0.929080f, 0.930190f, 0.931280f, 0.932350f, 0.933410f, 0.934460f, 0.935480f, 0.936500f,
        0.937490f, 0.938480f, 0.939450f, 0.940400f, 0.941340f, 0.942270f, 0.943180f, 0.944080f, 0.944960f, 0.945840f,
        0.946700f, 0.947540f, 0.948380f});
    static constexpr Spectrum CES99({0.169580f, 0.171400f, 0.173230f, 0.175080f, 0.176950f, 0.178820f, 0.180720f,
        0.182630f, 0.184560f, 0.186500f, 0.188460f, 0.190430f, 0.192420f, 0.194420f, 0.196440f, 0.198480f, 0.200530f,
        0.202600f, 0.204680f, 0.206780f, 0.209200f, 0.210850f, 0.212900f, 0.215240f, 0.217790f, 0.220440f, 0.223100f,
        0.225670f, 0.228060f, 0.230170f, 0.231900f, 0.233180f, 0.234040f, 0.234520f, 0.234680f, 0.234570f, 0.234220f,
        0.233710f, 0.233060f, 0.232350f, 0.231600f, 0.230870f, 0.230150f, 0.229420f, 0.228690f, 0.227940f, 0.227140f,
        0.226310f, 0.225410f, 0.224450f, 0.223400f, 0.222260f, 0.221040f, 0.219750f, 0.218380f, 0.216950f, 0.215480f,
        0.213950f, 0.212390f, 0.210810f, 0.209200f, 0.207580f, 0.205950f, 0.204310f, 0.202660f, 0.201000f, 0.199330f,
        0.197650f, 0.195970f, 0.194290f, 0.192600f, 0.190910f, 0.189210f, 0.187520f, 0.185830f, 0.184140f, 0.182450f,
        0.180770f, 0.179110f, 0.177450f, 0.175800f, 0.174170f, 0.172540f, 0.170930f, 0.169330f, 0.167730f, 0.166140f,
        0.164550f, 0.162970f, 0.161390f, 0.159800f, 0.158210f, 0.156620f, 0.155030f, 0.153430f, 0.151830f, 0.150230f,
        0.148630f, 0.147020f, 0.145410f, 0.143800f, 0.142190f, 0.140580f, 0.138970f, 0.137380f, 0.135810f, 0.134260f,
        0.132740f, 0.131250f, 0.129810f, 0.128400f, 0.127040f, 0.125730f, 0.124460f, 0.123230f, 0.122040f, 0.120880f,
        0.119740f, 0.118640f, 0.117560f, 0.116500f, 0.115460f, 0.114430f, 0.113400f, 0.112390f, 0.111370f, 0.110340f,
        0.109310f, 0.108260f, 0.107190f, 0.106100f, 0.104980f, 0.103840f, 0.102680f, 0.101500f, 0.100320f, 0.099129f,
        0.097937f, 0.096748f, 0.095568f, 0.094400f, 0.093249f, 0.092122f, 0.091024f, 0.089962f, 0.088943f, 0.087972f,
        0.087056f, 0.086201f, 0.085413f, 0.084700f, 0.084065f, 0.083508f, 0.083025f, 0.082614f, 0.082272f, 0.081996f,
        0.081784f, 0.081632f, 0.081539f, 0.081500f, 0.081513f, 0.081570f, 0.081664f, 0.081788f, 0.081932f, 0.082090f,
        0.082253f, 0.082415f, 0.082566f, 0.082700f, 0.082810f, 0.082900f, 0.082973f, 0.083034f, 0.083088f, 0.083138f,
        0.083190f, 0.083248f, 0.083317f, 0.083400f, 0.083504f, 0.083643f, 0.083829f, 0.084078f, 0.084405f, 0.084823f,
        0.085347f, 0.085992f, 0.086771f, 0.087700f, 0.088793f, 0.090063f, 0.091525f, 0.093193f, 0.095080f, 0.097201f,
        0.099570f, 0.102200f, 0.105110f, 0.108300f, 0.111800f, 0.115620f, 0.119770f, 0.124270f, 0.129140f, 0.134380f,
        0.140030f, 0.146090f, 0.152570f, 0.159500f, 0.166880f, 0.174680f, 0.182860f, 0.191380f, 0.200200f, 0.209290f,
        0.218600f, 0.228100f, 0.237750f, 0.247500f, 0.257330f, 0.267220f, 0.277160f, 0.287160f, 0.297190f, 0.307250f,
        0.317350f, 0.327450f, 0.337580f, 0.347700f, 0.357820f, 0.367910f, 0.377960f, 0.387930f, 0.397820f, 0.407590f,
        0.417230f, 0.426710f, 0.436000f, 0.445100f, 0.453980f, 0.462620f, 0.471030f, 0.479200f, 0.487120f, 0.494790f,
        0.502190f, 0.509340f, 0.516210f, 0.522800f, 0.529110f, 0.535140f, 0.540900f, 0.546380f, 0.551590f, 0.556540f,
        0.561220f, 0.565630f, 0.569790f, 0.573700f, 0.577360f, 0.580770f, 0.583970f, 0.586950f, 0.589750f, 0.592360f,
        0.594810f, 0.597100f, 0.599260f, 0.601300f, 0.603230f, 0.605050f, 0.606770f, 0.608390f, 0.609920f, 0.611350f,
        0.612690f, 0.613950f, 0.615110f, 0.616200f, 0.617210f, 0.618140f, 0.619010f, 0.619820f, 0.620570f, 0.621280f,
        0.621950f, 0.622590f, 0.623200f, 0.623800f, 0.624380f, 0.624950f, 0.625500f, 0.626040f, 0.626560f, 0.627050f,
        0.627530f, 0.627980f, 0.628400f, 0.628800f, 0.629170f, 0.629510f, 0.629830f, 0.630130f, 0.630410f, 0.630680f,
        0.630940f, 0.631190f, 0.631450f, 0.631700f, 0.631960f, 0.632210f, 0.632460f, 0.632710f, 0.632940f, 0.633160f,
        0.633360f, 0.633530f, 0.633680f, 0.633800f, 0.633990f, 0.634150f, 0.634310f, 0.634470f, 0.634630f, 0.634790f,
        0.634940f, 0.635100f, 0.635260f, 0.635420f, 0.635580f, 0.635740f, 0.635900f, 0.636060f, 0.636220f, 0.636380f,
        0.636540f, 0.636700f, 0.636860f, 0.637020f, 0.637180f, 0.637340f, 0.637500f, 0.637660f, 0.637820f, 0.637970f,
        0.638130f, 0.638290f, 0.638450f, 0.638610f, 0.638770f, 0.638930f, 0.639090f, 0.639250f, 0.639410f, 0.639560f,
        0.639720f, 0.639880f, 0.640040f, 0.640200f, 0.640360f, 0.640520f, 0.640680f, 0.640830f, 0.640990f, 0.641150f,
        0.641310f, 0.641470f, 0.641630f, 0.641790f, 0.641940f, 0.642100f, 0.642260f, 0.642420f, 0.642580f, 0.642740f,
        0.642890f, 0.643050f, 0.643210f, 0.643370f, 0.643530f, 0.643680f, 0.643840f, 0.644000f, 0.644160f, 0.644320f,
        0.644470f, 0.644630f, 0.644790f, 0.644950f, 0.645110f, 0.645260f, 0.645420f, 0.645580f, 0.645740f, 0.645890f,
        0.646050f, 0.646210f, 0.646370f});

    static constexpr std::array<Spectrum, 99> sample{CES01, CES02, CES03, CES04, CES05, CES06, CES07, CES08, CES09,
        CES10, CES11, CES12, CES13, CES14, CES15, CES16, CES17, CES18, CES19, CES20, CES21, CES22, CES23, CES24, CES25,
        CES26, CES27, CES28, CES29, CES30, CES31, CES32, CES33, CES34, CES35, CES36, CES37, CES38, CES39, CES40, CES41,
        CES42, CES43, CES44, CES45, CES46, CES47, CES48, CES49, CES50, CES51, CES52, CES53, CES54, CES55, CES56, CES57,
        CES58, CES59, CES60, CES61, CES62, CES63, CES64, CES65, CES66, CES67, CES68, CES69, CES70, CES71, CES72, CES73,
        CES74, CES75, CES76, CES77, CES78, CES79, CES80, CES81, CES82, CES83, CES84, CES85, CES86, CES87, CES88, CES89,
        CES90, CES91, CES92, CES93, CES94, CES95, CES96, CES97, CES98, CES99};

} // namespace TM_30_15

// ------------------- implements macbeth chart.

namespace Macbeth
{
    static constexpr Spectrum Macbeth01({0.054745777f, 0.055103879f, 0.055461981f, 0.055820083f, 0.056178185f,
        0.056536288f, 0.05689439f, 0.057252492f, 0.057610594f, 0.057968696f, 0.058326798f, 0.058610495f, 0.058894191f,
        0.059177887f, 0.059461584f, 0.05974528f, 0.060028976f, 0.060312673f, 0.060596369f, 0.060880065f, 0.061163762f,
        0.061285862f, 0.061407963f, 0.061530064f, 0.061652164f, 0.061774265f, 0.061896366f, 0.062018466f, 0.062140567f,
        0.062262668f, 0.062384768f, 0.062377367f, 0.062369965f, 0.062362563f, 0.062355162f, 0.06234776f, 0.062340358f,
        0.062332957f, 0.062325555f, 0.062318153f, 0.062310752f, 0.062286679f, 0.062262607f, 0.062238535f, 0.062214462f,
        0.06219039f, 0.062166318f, 0.062142245f, 0.062118173f, 0.062094101f, 0.062070028f, 0.062045879f, 0.06202173f,
        0.06199758f, 0.061973431f, 0.061949282f, 0.061925132f, 0.061900983f, 0.061876834f, 0.061852684f, 0.061828535f,
        0.061804253f, 0.061779972f, 0.06175569f, 0.061731408f, 0.061707127f, 0.061682845f, 0.061658563f, 0.061634282f,
        0.06161f, 0.061585718f, 0.061581288f, 0.061576857f, 0.061572426f, 0.061567995f, 0.061563564f, 0.061559133f,
        0.061554703f, 0.061550272f, 0.061545841f, 0.06154141f, 0.061549581f, 0.061557751f, 0.061565922f, 0.061574092f,
        0.061582263f, 0.061590433f, 0.061598604f, 0.061606774f, 0.061614945f, 0.061623115f, 0.061663561f, 0.061704008f,
        0.061744454f, 0.0617849f, 0.061825347f, 0.061865793f, 0.061906239f, 0.061946686f, 0.061987132f, 0.062027578f,
        0.062121243f, 0.062214907f, 0.062308571f, 0.062402235f, 0.062495899f, 0.062589563f, 0.062683228f, 0.062776892f,
        0.062870556f, 0.06296422f, 0.063185879f, 0.063407537f, 0.063629196f, 0.063850855f, 0.064072513f, 0.064294172f,
        0.064515831f, 0.064737489f, 0.064959148f, 0.065180807f, 0.065689872f, 0.066198937f, 0.066708003f, 0.067217068f,
        0.067726133f, 0.068235199f, 0.068744264f, 0.069253329f, 0.069762395f, 0.07027146f, 0.070883902f, 0.071496344f,
        0.072108786f, 0.072721227f, 0.073333669f, 0.073946111f, 0.074558553f, 0.075170995f, 0.075783437f, 0.076395878f,
        0.076705358f, 0.077014837f, 0.077324316f, 0.077633795f, 0.077943274f, 0.078252753f, 0.078562233f, 0.078871712f,
        0.079181191f, 0.07949067f, 0.079669186f, 0.079847702f, 0.080026218f, 0.080204734f, 0.08038325f, 0.080561766f,
        0.080740282f, 0.080918798f, 0.081097314f, 0.08127583f, 0.081577204f, 0.081878577f, 0.082179951f, 0.082481324f,
        0.082782698f, 0.083084071f, 0.083385445f, 0.083686818f, 0.083988192f, 0.084289565f, 0.084918979f, 0.085548393f,
        0.086177807f, 0.086807221f, 0.087436635f, 0.088066049f, 0.088695463f, 0.089324877f, 0.089954291f, 0.090583705f,
        0.091815509f, 0.093047312f, 0.094279116f, 0.09551092f, 0.096742723f, 0.097974527f, 0.099206331f, 0.100438134f,
        0.101669938f, 0.102901742f, 0.104516701f, 0.10613166f, 0.107746619f, 0.109361578f, 0.110976538f, 0.112591497f,
        0.114206456f, 0.115821415f, 0.117436374f, 0.119051333f, 0.120572658f, 0.122093982f, 0.123615306f, 0.12513663f,
        0.126657954f, 0.128179278f, 0.129700603f, 0.131221927f, 0.132743251f, 0.134264575f, 0.135157913f, 0.13605125f,
        0.136944588f, 0.137837925f, 0.138731263f, 0.1396246f, 0.140517938f, 0.141411275f, 0.142304613f, 0.14319795f,
        0.143566098f, 0.143934245f, 0.144302393f, 0.14467054f, 0.145038688f, 0.145406835f, 0.145774983f, 0.14614313f,
        0.146511278f, 0.146879425f, 0.147268983f, 0.14765854f, 0.148048098f, 0.148437655f, 0.148827213f, 0.14921677f,
        0.149606328f, 0.149995885f, 0.150385443f, 0.150775f, 0.151507093f, 0.152239185f, 0.152971278f, 0.15370337f,
        0.154435463f, 0.155167555f, 0.155899648f, 0.15663174f, 0.157363833f, 0.158095925f, 0.159105434f, 0.160114943f,
        0.161124453f, 0.162133962f, 0.163143471f, 0.16415298f, 0.165162489f, 0.166171998f, 0.167181508f, 0.168191017f,
        0.169262342f, 0.170333667f, 0.171404992f, 0.172476317f, 0.173547642f, 0.174618967f, 0.175690292f, 0.176761617f,
        0.177832942f, 0.178904267f, 0.179768972f, 0.180633677f, 0.181498382f, 0.182363087f, 0.183227792f, 0.184092497f,
        0.184957202f, 0.185821907f, 0.186686612f, 0.187551317f, 0.187759782f, 0.187968247f, 0.188176712f, 0.188385177f,
        0.188593642f, 0.188802107f, 0.189010572f, 0.189219037f, 0.189427502f, 0.189635967f, 0.189249621f, 0.188863275f,
        0.188476929f, 0.188090583f, 0.187704238f, 0.187317892f, 0.186931546f, 0.1865452f, 0.186158854f, 0.185772508f,
        0.185344593f, 0.184916678f, 0.184488763f, 0.184060848f, 0.183632933f, 0.183205018f, 0.182777103f, 0.182349188f,
        0.181921273f, 0.181493358f, 0.181504526f, 0.181515693f, 0.181526861f, 0.181538028f, 0.181549196f, 0.181560363f,
        0.181571531f, 0.181582698f, 0.181593866f, 0.181605033f, 0.18216564f, 0.182726247f, 0.183286853f, 0.18384746f,
        0.184408067f, 0.184968673f, 0.18552928f, 0.186089887f, 0.186650493f, 0.1872111f, 0.188095342f, 0.188979583f,
        0.189863825f, 0.190748067f, 0.191632308f, 0.19251655f, 0.193400792f, 0.194285033f, 0.195169275f, 0.196053517f,
        0.197397435f, 0.198741353f, 0.200085272f, 0.20142919f, 0.202773108f, 0.204117027f, 0.205460945f, 0.206804863f,
        0.208148782f, 0.2094927f, 0.18854343f, 0.16759416f, 0.14664489f, 0.12569562f, 0.10474635f, 0.08379708f,
        0.06284781f, 0.04189854f, 0.02094927f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth02({0.117132943f, 0.119765069f, 0.122397195f, 0.12502932f, 0.127661446f,
        0.130293572f, 0.132925697f, 0.135557823f, 0.138189949f, 0.140822074f, 0.1434542f, 0.146561686f, 0.149669172f,
        0.152776658f, 0.155884143f, 0.158991629f, 0.162099115f, 0.165206601f, 0.168314087f, 0.171421573f, 0.174529058f,
        0.176169633f, 0.177810208f, 0.179450783f, 0.181091358f, 0.182731933f, 0.184372508f, 0.186013083f, 0.187653658f,
        0.189294233f, 0.190934808f, 0.19140154f, 0.191868272f, 0.192335003f, 0.192801735f, 0.193268467f, 0.193735198f,
        0.19420193f, 0.194668662f, 0.195135393f, 0.195602125f, 0.195941709f, 0.196281293f, 0.196620878f, 0.196960462f,
        0.197300046f, 0.19763963f, 0.197979214f, 0.198318798f, 0.198658383f, 0.198997967f, 0.199521663f, 0.20004536f,
        0.200569057f, 0.201092753f, 0.20161645f, 0.202140147f, 0.202663843f, 0.20318754f, 0.203711237f, 0.204234933f,
        0.20512993f, 0.206024927f, 0.206919923f, 0.20781492f, 0.208709917f, 0.209604913f, 0.21049991f, 0.211394907f,
        0.212289903f, 0.2131849f, 0.214708742f, 0.216232583f, 0.217756425f, 0.219280267f, 0.220804108f, 0.22232795f,
        0.223851792f, 0.225375633f, 0.226899475f, 0.228423317f, 0.230708312f, 0.232993307f, 0.235278302f, 0.237563297f,
        0.239848292f, 0.242133287f, 0.244418282f, 0.246703277f, 0.248988272f, 0.251273267f, 0.254150625f, 0.257027983f,
        0.259905342f, 0.2627827f, 0.265660058f, 0.268537417f, 0.271414775f, 0.274292133f, 0.277169492f, 0.28004685f,
        0.282920274f, 0.285793698f, 0.288667123f, 0.291540547f, 0.294413971f, 0.297287395f, 0.300160819f, 0.303034243f,
        0.305907668f, 0.308781092f, 0.310848142f, 0.312915192f, 0.314982242f, 0.317049292f, 0.319116342f, 0.321183392f,
        0.323250442f, 0.325317492f, 0.327384542f, 0.329451592f, 0.329842878f, 0.330234163f, 0.330625449f, 0.331016735f,
        0.331408021f, 0.331799307f, 0.332190593f, 0.332581878f, 0.332973164f, 0.33336445f, 0.331487665f, 0.32961088f,
        0.327734095f, 0.32585731f, 0.323980525f, 0.32210374f, 0.320226955f, 0.31835017f, 0.316473385f, 0.3145966f,
        0.31176484f, 0.30893308f, 0.30610132f, 0.30326956f, 0.3004378f, 0.29760604f, 0.29477428f, 0.29194252f,
        0.28911076f, 0.286279f, 0.284999998f, 0.283720997f, 0.282441995f, 0.281162993f, 0.279883992f, 0.27860499f,
        0.277325988f, 0.276046987f, 0.274767985f, 0.273488983f, 0.27378576f, 0.274082537f, 0.274379313f, 0.27467609f,
        0.274972867f, 0.275269643f, 0.27556642f, 0.275863197f, 0.276159973f, 0.27645675f, 0.276531368f, 0.276605987f,
        0.276680605f, 0.276755223f, 0.276829842f, 0.27690446f, 0.276979078f, 0.277053697f, 0.277128315f, 0.277202933f,
        0.278412386f, 0.279621838f, 0.280831291f, 0.282040743f, 0.283250196f, 0.284459648f, 0.285669101f, 0.286878553f,
        0.288088006f, 0.289297458f, 0.294305343f, 0.299313227f, 0.304321111f, 0.309328995f, 0.314336879f, 0.319344763f,
        0.324352648f, 0.329360532f, 0.334368416f, 0.3393763f, 0.347460277f, 0.355544253f, 0.36362823f, 0.371712207f,
        0.379796183f, 0.38788016f, 0.395964137f, 0.404048113f, 0.41213209f, 0.420216067f, 0.426973407f, 0.433730747f,
        0.440488087f, 0.447245427f, 0.454002767f, 0.460760107f, 0.467517447f, 0.474274787f, 0.481032127f, 0.487789467f,
        0.491521339f, 0.495253212f, 0.498985084f, 0.502716957f, 0.506448829f, 0.510180702f, 0.513912574f, 0.517644447f,
        0.521376319f, 0.525108192f, 0.527171099f, 0.529234007f, 0.531296914f, 0.533359822f, 0.535422729f, 0.537485637f,
        0.539548544f, 0.541611452f, 0.543674359f, 0.545737267f, 0.547319515f, 0.548901763f, 0.550484012f, 0.55206626f,
        0.553648508f, 0.555230757f, 0.556813005f, 0.558395253f, 0.559977502f, 0.56155975f, 0.563192083f, 0.564824417f,
        0.56645675f, 0.568089083f, 0.569721417f, 0.57135375f, 0.572986083f, 0.574618417f, 0.57625075f, 0.577883083f,
        0.579591896f, 0.581300708f, 0.583009521f, 0.584718333f, 0.586427146f, 0.588135958f, 0.589844771f, 0.591553583f,
        0.593262396f, 0.594971208f, 0.596653813f, 0.598336417f, 0.600019021f, 0.601701625f, 0.603384229f, 0.605066833f,
        0.606749438f, 0.608432042f, 0.610114646f, 0.61179725f, 0.613092223f, 0.614387195f, 0.615682168f, 0.61697714f,
        0.618272113f, 0.619567085f, 0.620862058f, 0.62215703f, 0.623452003f, 0.624746975f, 0.626082417f, 0.627417858f,
        0.6287533f, 0.630088742f, 0.631424183f, 0.632759625f, 0.634095067f, 0.635430508f, 0.63676595f, 0.638101392f,
        0.639887074f, 0.641672757f, 0.643458439f, 0.645244122f, 0.647029804f, 0.648815487f, 0.650601169f, 0.652386852f,
        0.654172534f, 0.655958217f, 0.65818396f, 0.660409703f, 0.662635447f, 0.66486119f, 0.667086933f, 0.669312677f,
        0.67153842f, 0.673764163f, 0.675989907f, 0.67821565f, 0.680352183f, 0.682488717f, 0.68462525f, 0.686761783f,
        0.688898317f, 0.69103485f, 0.693171383f, 0.695307917f, 0.69744445f, 0.699580983f, 0.701331758f, 0.703082532f,
        0.704833306f, 0.70658408f, 0.708334854f, 0.710085628f, 0.711836403f, 0.713587177f, 0.715337951f, 0.717088725f,
        0.718762279f, 0.720435833f, 0.722109388f, 0.723782942f, 0.725456496f, 0.72713005f, 0.728803604f, 0.730477158f,
        0.732150713f, 0.733824267f, 0.66044184f, 0.587059413f, 0.513676987f, 0.44029456f, 0.366912133f, 0.293529707f,
        0.22014728f, 0.146764853f, 0.073382427f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth03({0.130362525f, 0.135033678f, 0.13970483f, 0.144375983f, 0.149047135f,
        0.153718288f, 0.15838944f, 0.163060593f, 0.167731745f, 0.172402898f, 0.17707405f, 0.184467636f, 0.191861222f,
        0.199254808f, 0.206648393f, 0.214041979f, 0.221435565f, 0.228829151f, 0.236222737f, 0.243616323f, 0.251009908f,
        0.256533426f, 0.262056943f, 0.267580461f, 0.273103978f, 0.278627496f, 0.284151013f, 0.289674531f, 0.295198048f,
        0.300721566f, 0.306245083f, 0.308012125f, 0.309779167f, 0.311546208f, 0.31331325f, 0.315080292f, 0.316847333f,
        0.318614375f, 0.320381417f, 0.322148458f, 0.3239155f, 0.324517105f, 0.32511871f, 0.325720315f, 0.32632192f,
        0.326923525f, 0.32752513f, 0.328126735f, 0.32872834f, 0.329329945f, 0.32993155f, 0.330221306f, 0.330511062f,
        0.330800818f, 0.331090573f, 0.331380329f, 0.331670085f, 0.331959841f, 0.332249597f, 0.332539353f, 0.332829108f,
        0.332643691f, 0.332458273f, 0.332272856f, 0.332087438f, 0.331902021f, 0.331716603f, 0.331531186f, 0.331345768f,
        0.331160351f, 0.330974933f, 0.330219834f, 0.329464735f, 0.328709636f, 0.327954537f, 0.327199438f, 0.326444338f,
        0.325689239f, 0.32493414f, 0.324179041f, 0.323423942f, 0.322215453f, 0.321006963f, 0.319798474f, 0.318589985f,
        0.317381496f, 0.316173007f, 0.314964518f, 0.313756028f, 0.312547539f, 0.31133905f, 0.310027743f, 0.308716437f,
        0.30740513f, 0.306093823f, 0.304782517f, 0.30347121f, 0.302159903f, 0.300848597f, 0.29953729f, 0.298225983f,
        0.296936445f, 0.295646907f, 0.294357368f, 0.29306783f, 0.291778292f, 0.290488753f, 0.289199215f, 0.287909677f,
        0.286620138f, 0.2853306f, 0.283740298f, 0.282149995f, 0.280559693f, 0.27896939f, 0.277379088f, 0.275788785f,
        0.274198483f, 0.27260818f, 0.271017878f, 0.269427575f, 0.26752225f, 0.265616925f, 0.2637116f, 0.261806275f,
        0.25990095f, 0.257995625f, 0.2560903f, 0.254184975f, 0.25227965f, 0.250374325f, 0.248481174f, 0.246588023f,
        0.244694873f, 0.242801722f, 0.240908571f, 0.23901542f, 0.237122269f, 0.235229118f, 0.233335968f, 0.231442817f,
        0.229724124f, 0.228005432f, 0.226286739f, 0.224568047f, 0.222849354f, 0.221130662f, 0.219411969f, 0.217693277f,
        0.215974584f, 0.214255892f, 0.212772542f, 0.211289192f, 0.209805842f, 0.208322492f, 0.206839142f, 0.205355792f,
        0.203872442f, 0.202389092f, 0.200905742f, 0.199422392f, 0.197931317f, 0.196440242f, 0.194949167f, 0.193458092f,
        0.191967017f, 0.190475942f, 0.188984867f, 0.187493792f, 0.186002717f, 0.184511642f, 0.182998908f, 0.181486173f,
        0.179973439f, 0.178460705f, 0.176947971f, 0.175435237f, 0.173922503f, 0.172409768f, 0.170897034f, 0.1693843f,
        0.168174663f, 0.166965025f, 0.165755388f, 0.16454575f, 0.163336113f, 0.162126475f, 0.160916838f, 0.1597072f,
        0.158497563f, 0.157287925f, 0.156469936f, 0.155651947f, 0.154833958f, 0.154015968f, 0.153197979f, 0.15237999f,
        0.151562001f, 0.150744012f, 0.149926023f, 0.149108033f, 0.14867956f, 0.148251087f, 0.147822613f, 0.14739414f,
        0.146965667f, 0.146537193f, 0.14610872f, 0.145680247f, 0.145251773f, 0.1448233f, 0.144526793f, 0.144230287f,
        0.14393378f, 0.143637273f, 0.143340767f, 0.14304426f, 0.142747753f, 0.142451247f, 0.14215474f, 0.141858233f,
        0.141729483f, 0.141600733f, 0.141471983f, 0.141343233f, 0.141214483f, 0.141085733f, 0.140956983f, 0.140828233f,
        0.140699483f, 0.140570733f, 0.140581125f, 0.140591517f, 0.140601908f, 0.1406123f, 0.140622692f, 0.140633083f,
        0.140643475f, 0.140653867f, 0.140664258f, 0.14067465f, 0.140716511f, 0.140758372f, 0.140800233f, 0.140842093f,
        0.140883954f, 0.140925815f, 0.140967676f, 0.141009537f, 0.141051398f, 0.141093258f, 0.141240748f, 0.141388238f,
        0.141535728f, 0.141683218f, 0.141830708f, 0.141978198f, 0.142125688f, 0.142273178f, 0.142420668f, 0.142568158f,
        0.142965302f, 0.143362445f, 0.143759588f, 0.144156732f, 0.144553875f, 0.144951018f, 0.145348162f, 0.145745305f,
        0.146142448f, 0.146539592f, 0.147069393f, 0.147599195f, 0.148128997f, 0.148658798f, 0.1491886f, 0.149718402f,
        0.150248203f, 0.150778005f, 0.151307807f, 0.151837608f, 0.152005322f, 0.152173035f, 0.152340748f, 0.152508462f,
        0.152676175f, 0.152843888f, 0.153011602f, 0.153179315f, 0.153347028f, 0.153514742f, 0.153172203f, 0.152829663f,
        0.152487124f, 0.152144585f, 0.151802046f, 0.151459507f, 0.151116968f, 0.150774428f, 0.150431889f, 0.15008935f,
        0.149475664f, 0.148861978f, 0.148248293f, 0.147634607f, 0.147020921f, 0.146407235f, 0.145793549f, 0.145179863f,
        0.144566178f, 0.143952492f, 0.143196713f, 0.142440933f, 0.141685154f, 0.140929375f, 0.140173596f, 0.139417817f,
        0.138662038f, 0.137906258f, 0.137150479f, 0.1363947f, 0.135990375f, 0.135586049f, 0.135181724f, 0.134777398f,
        0.134373073f, 0.133968747f, 0.133564422f, 0.133160096f, 0.132755771f, 0.132351445f, 0.132612045f, 0.132872644f,
        0.133133244f, 0.133393844f, 0.133654443f, 0.133915043f, 0.134175643f, 0.134436242f, 0.134696842f, 0.134957442f,
        0.136135052f, 0.137312662f, 0.138490272f, 0.139667882f, 0.140845493f, 0.142023103f, 0.143200713f, 0.144378323f,
        0.145555933f, 0.146733543f, 0.132060189f, 0.117386835f, 0.10271348f, 0.088040126f, 0.073366772f, 0.058693417f,
        0.044020063f, 0.029346709f, 0.014673354f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth04({0.05123929f, 0.051538762f, 0.051838234f, 0.052137706f, 0.052437178f,
        0.05273665f, 0.053036122f, 0.053335594f, 0.053635066f, 0.053934538f, 0.05423401f, 0.054409636f, 0.054585261f,
        0.054760887f, 0.054936512f, 0.055112138f, 0.055287763f, 0.055463389f, 0.055639014f, 0.05581464f, 0.055990265f,
        0.056095125f, 0.056199985f, 0.056304845f, 0.056409705f, 0.056514565f, 0.056619425f, 0.056724285f, 0.056829145f,
        0.056934005f, 0.057038865f, 0.0571212f, 0.057203536f, 0.057285871f, 0.057368206f, 0.057450542f, 0.057532877f,
        0.057615212f, 0.057697548f, 0.057779883f, 0.057862218f, 0.05797104f, 0.058079861f, 0.058188682f, 0.058297503f,
        0.058406324f, 0.058515145f, 0.058623967f, 0.058732788f, 0.058841609f, 0.05895043f, 0.059085396f, 0.059220361f,
        0.059355327f, 0.059490292f, 0.059625258f, 0.059760223f, 0.059895189f, 0.060030154f, 0.06016512f, 0.060300085f,
        0.060401482f, 0.060502878f, 0.060604275f, 0.060705671f, 0.060807068f, 0.060908464f, 0.061009861f, 0.061111257f,
        0.061212654f, 0.06131405f, 0.0614103f, 0.06150655f, 0.0616028f, 0.06169905f, 0.0617953f, 0.06189155f,
        0.0619878f, 0.06208405f, 0.0621803f, 0.06227655f, 0.062373402f, 0.062470254f, 0.062567106f, 0.062663958f,
        0.06276081f, 0.062857662f, 0.062954514f, 0.063051366f, 0.063148218f, 0.06324507f, 0.063398907f, 0.063552743f,
        0.06370658f, 0.063860416f, 0.064014253f, 0.064168089f, 0.064321926f, 0.064475762f, 0.064629599f, 0.064783435f,
        0.065042848f, 0.065302261f, 0.065561674f, 0.065821086f, 0.066080499f, 0.066339912f, 0.066599325f, 0.066858738f,
        0.067118151f, 0.067377563f, 0.068171087f, 0.068964611f, 0.069758135f, 0.070551659f, 0.071345183f, 0.072138707f,
        0.072932231f, 0.073725755f, 0.074519279f, 0.075312803f, 0.077901374f, 0.080489944f, 0.083078515f, 0.085667085f,
        0.088255656f, 0.090844226f, 0.093432797f, 0.096021367f, 0.098609938f, 0.101198508f, 0.105614732f, 0.110030955f,
        0.114447178f, 0.118863402f, 0.123279625f, 0.127695848f, 0.132112072f, 0.136528295f, 0.140944518f, 0.145360742f,
        0.148651089f, 0.151941437f, 0.155231784f, 0.158522132f, 0.161812479f, 0.165102827f, 0.168393174f, 0.171683522f,
        0.174973869f, 0.178264217f, 0.178832204f, 0.179400192f, 0.179968179f, 0.180536167f, 0.181104154f, 0.181672142f,
        0.182240129f, 0.182808117f, 0.183376104f, 0.183944092f, 0.182560788f, 0.181177485f, 0.179794182f, 0.178410878f,
        0.177027575f, 0.175644272f, 0.174260968f, 0.172877665f, 0.171494362f, 0.170111058f, 0.168038043f, 0.165965027f,
        0.163892011f, 0.161818995f, 0.159745979f, 0.157672963f, 0.155599948f, 0.153526932f, 0.151453916f, 0.1493809f,
        0.147717223f, 0.146053547f, 0.14438987f, 0.142726193f, 0.141062517f, 0.13939884f, 0.137735163f, 0.136071487f,
        0.13440781f, 0.132744133f, 0.131655826f, 0.130567518f, 0.129479211f, 0.128390903f, 0.127302596f, 0.126214288f,
        0.125125981f, 0.124037673f, 0.122949366f, 0.121861058f, 0.121192263f, 0.120523468f, 0.119854673f, 0.119185878f,
        0.118517083f, 0.117848288f, 0.117179493f, 0.116510698f, 0.115841903f, 0.115173108f, 0.114603934f, 0.11403476f,
        0.113465586f, 0.112896412f, 0.112327238f, 0.111758063f, 0.111188889f, 0.110619715f, 0.110050541f, 0.109481367f,
        0.109069129f, 0.108656892f, 0.108244654f, 0.107832417f, 0.107420179f, 0.107007942f, 0.106595704f, 0.106183467f,
        0.105771229f, 0.105358992f, 0.105257307f, 0.105155623f, 0.105053938f, 0.104952254f, 0.104850569f, 0.104748885f,
        0.1046472f, 0.104545516f, 0.104443831f, 0.104342147f, 0.104506853f, 0.104671559f, 0.104836265f, 0.105000971f,
        0.105165678f, 0.105330384f, 0.10549509f, 0.105659796f, 0.105824502f, 0.105989208f, 0.106280923f, 0.106572637f,
        0.106864351f, 0.107156065f, 0.107447779f, 0.107739493f, 0.108031208f, 0.108322922f, 0.108614636f, 0.10890635f,
        0.109205153f, 0.109503956f, 0.10980276f, 0.110101563f, 0.110400366f, 0.110699169f, 0.110997972f, 0.111296775f,
        0.111595579f, 0.111894382f, 0.112111166f, 0.11232795f, 0.112544735f, 0.112761519f, 0.112978303f, 0.113195088f,
        0.113411872f, 0.113628656f, 0.113845441f, 0.114062225f, 0.114051304f, 0.114040383f, 0.114029463f, 0.114018542f,
        0.114007621f, 0.1139967f, 0.113985779f, 0.113974858f, 0.113963938f, 0.113953017f, 0.113797518f, 0.11364202f,
        0.113486522f, 0.113331023f, 0.113175525f, 0.113020027f, 0.112864528f, 0.11270903f, 0.112553532f, 0.112398033f,
        0.112373548f, 0.112349063f, 0.112324578f, 0.112300093f, 0.112275608f, 0.112251123f, 0.112226638f, 0.112202153f,
        0.112177668f, 0.112153183f, 0.112419765f, 0.112686347f, 0.112952928f, 0.11321951f, 0.113486092f, 0.113752673f,
        0.114019255f, 0.114285837f, 0.114552418f, 0.114819f, 0.115313643f, 0.115808285f, 0.116302928f, 0.11679757f,
        0.117292213f, 0.117786855f, 0.118281498f, 0.11877614f, 0.119270783f, 0.119765425f, 0.120248358f, 0.120731292f,
        0.121214225f, 0.121697158f, 0.122180092f, 0.122663025f, 0.123145958f, 0.123628892f, 0.124111825f, 0.124594758f,
        0.125165738f, 0.125736717f, 0.126307696f, 0.126878675f, 0.127449654f, 0.128020633f, 0.128591613f, 0.129162592f,
        0.129733571f, 0.13030455f, 0.117274095f, 0.10424364f, 0.091213185f, 0.07818273f, 0.065152275f, 0.05212182f,
        0.039091365f, 0.02606091f, 0.013030455f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth05({0.144234675f, 0.149638121f, 0.155041567f, 0.160445013f, 0.165848458f,
        0.171251904f, 0.17665535f, 0.182058796f, 0.187462242f, 0.192865688f, 0.198269133f, 0.207885608f, 0.217502082f,
        0.227118556f, 0.23673503f, 0.246351504f, 0.255967978f, 0.265584453f, 0.275200927f, 0.284817401f, 0.294433875f,
        0.30253424f, 0.310634605f, 0.31873497f, 0.326835335f, 0.3349357f, 0.343036065f, 0.35113643f, 0.359236795f,
        0.36733716f, 0.375437525f, 0.378731121f, 0.382024717f, 0.385318313f, 0.388611908f, 0.391905504f, 0.3951991f,
        0.398492696f, 0.401786292f, 0.405079888f, 0.408373483f, 0.409631083f, 0.410888683f, 0.412146283f, 0.413403883f,
        0.414661483f, 0.415919083f, 0.417176683f, 0.418434283f, 0.419691883f, 0.420949483f, 0.421472757f, 0.42199603f,
        0.422519303f, 0.423042577f, 0.42356585f, 0.424089123f, 0.424612397f, 0.42513567f, 0.425658943f, 0.426182217f,
        0.426172657f, 0.426163097f, 0.426153537f, 0.426143977f, 0.426134417f, 0.426124857f, 0.426115297f, 0.426105737f,
        0.426096177f, 0.426086617f, 0.425410258f, 0.424733898f, 0.424057539f, 0.42338118f, 0.422704821f, 0.422028462f,
        0.421352103f, 0.420675743f, 0.419999384f, 0.419323025f, 0.417734222f, 0.416145418f, 0.414556615f, 0.412967812f,
        0.411379008f, 0.409790205f, 0.408201402f, 0.406612598f, 0.405023795f, 0.403434992f, 0.401018862f, 0.398602732f,
        0.396186602f, 0.393770472f, 0.391354342f, 0.388938212f, 0.386522082f, 0.384105952f, 0.381689822f, 0.379273692f,
        0.375982226f, 0.37269076f, 0.369399294f, 0.366107828f, 0.362816363f, 0.359524897f, 0.356233431f, 0.352941965f,
        0.349650499f, 0.346359033f, 0.342835295f, 0.339311557f, 0.335787818f, 0.33226408f, 0.328740342f, 0.325216603f,
        0.321692865f, 0.318169127f, 0.314645388f, 0.31112165f, 0.308132998f, 0.305144345f, 0.302155693f, 0.29916704f,
        0.296178388f, 0.293189735f, 0.290201083f, 0.28721243f, 0.284223778f, 0.281235125f, 0.278499214f, 0.275763303f,
        0.273027393f, 0.270291482f, 0.267555571f, 0.26481966f, 0.262083749f, 0.259347838f, 0.256611928f, 0.253876017f,
        0.251377398f, 0.248878778f, 0.246380159f, 0.24388154f, 0.241382921f, 0.238884302f, 0.236385683f, 0.233887063f,
        0.231388444f, 0.228889825f, 0.227421307f, 0.225952788f, 0.22448427f, 0.223015752f, 0.221547233f, 0.220078715f,
        0.218610197f, 0.217141678f, 0.21567316f, 0.214204642f, 0.213619655f, 0.213034668f, 0.212449682f, 0.211864695f,
        0.211279708f, 0.210694722f, 0.210109735f, 0.209524748f, 0.208939762f, 0.208354775f, 0.207681531f, 0.207008287f,
        0.206335043f, 0.205661798f, 0.204988554f, 0.20431531f, 0.203642066f, 0.202968822f, 0.202295578f, 0.201622333f,
        0.20089973f, 0.200177127f, 0.199454523f, 0.19873192f, 0.198009317f, 0.197286713f, 0.19656411f, 0.195841507f,
        0.195118903f, 0.1943963f, 0.194213338f, 0.194030375f, 0.193847413f, 0.19366445f, 0.193481488f, 0.193298525f,
        0.193115563f, 0.1929326f, 0.192749638f, 0.192566675f, 0.193328157f, 0.194089638f, 0.19485112f, 0.195612602f,
        0.196374083f, 0.197135565f, 0.197897047f, 0.198658528f, 0.19942001f, 0.200181492f, 0.201604331f, 0.20302717f,
        0.204450009f, 0.205872848f, 0.207295688f, 0.208718527f, 0.210141366f, 0.211564205f, 0.212987044f, 0.214409883f,
        0.215920689f, 0.217431495f, 0.218942301f, 0.220453107f, 0.221963913f, 0.223474718f, 0.224985524f, 0.22649633f,
        0.228007136f, 0.229517942f, 0.230623808f, 0.231729675f, 0.232835542f, 0.233941408f, 0.235047275f, 0.236153142f,
        0.237259008f, 0.238364875f, 0.239470742f, 0.240576608f, 0.241915044f, 0.24325348f, 0.244591916f, 0.245930352f,
        0.247268788f, 0.248607223f, 0.249945659f, 0.251284095f, 0.252622531f, 0.253960967f, 0.256415994f, 0.258871022f,
        0.261326049f, 0.263781077f, 0.266236104f, 0.268691132f, 0.271146159f, 0.273601187f, 0.276056214f, 0.278511242f,
        0.281982214f, 0.285453187f, 0.288924159f, 0.292395132f, 0.295866104f, 0.299337077f, 0.302808049f, 0.306279022f,
        0.309749994f, 0.313220967f, 0.316678026f, 0.320135085f, 0.323592144f, 0.327049203f, 0.330506263f, 0.333963322f,
        0.337420381f, 0.34087744f, 0.344334499f, 0.347791558f, 0.349599475f, 0.351407392f, 0.353215308f, 0.355023225f,
        0.356831142f, 0.358639058f, 0.360446975f, 0.362254892f, 0.364062808f, 0.365870725f, 0.36586278f, 0.365854835f,
        0.36584689f, 0.365838945f, 0.365831f, 0.365823055f, 0.36581511f, 0.365807165f, 0.36579922f, 0.365791275f,
        0.36515456f, 0.364517845f, 0.36388113f, 0.363244415f, 0.3626077f, 0.361970985f, 0.36133427f, 0.360697555f,
        0.36006084f, 0.359424125f, 0.359281038f, 0.35913795f, 0.358994863f, 0.358851775f, 0.358708688f, 0.3585656f,
        0.358422513f, 0.358279425f, 0.358136338f, 0.35799325f, 0.35868661f, 0.35937997f, 0.36007333f, 0.36076669f,
        0.36146005f, 0.36215341f, 0.36284677f, 0.36354013f, 0.36423349f, 0.36492685f, 0.366156758f, 0.367386667f,
        0.368616575f, 0.369846483f, 0.371076392f, 0.3723063f, 0.373536208f, 0.374766117f, 0.375996025f, 0.377225933f,
        0.379285967f, 0.381346f, 0.383406033f, 0.385466067f, 0.3875261f, 0.389586133f, 0.391646167f, 0.3937062f,
        0.395766233f, 0.397826267f, 0.35804364f, 0.318261013f, 0.278478387f, 0.23869576f, 0.198913133f, 0.159130507f,
        0.11934788f, 0.079565253f, 0.039782627f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth06({0.136268951f, 0.140587868f, 0.144906785f, 0.149225702f, 0.153544619f,
        0.157863536f, 0.162182452f, 0.166501369f, 0.170820286f, 0.175139203f, 0.17945812f, 0.186201515f, 0.19294491f,
        0.199688305f, 0.2064317f, 0.213175095f, 0.21991849f, 0.226661885f, 0.23340528f, 0.240148675f, 0.24689207f,
        0.251884762f, 0.256877454f, 0.261870146f, 0.266862839f, 0.271855531f, 0.276848223f, 0.281840915f, 0.286833607f,
        0.291826299f, 0.296818991f, 0.299165242f, 0.301511493f, 0.303857744f, 0.306203995f, 0.308550246f, 0.310896496f,
        0.313242747f, 0.315588998f, 0.317935249f, 0.3202815f, 0.321961795f, 0.323642089f, 0.325322384f, 0.327002679f,
        0.328682974f, 0.330363268f, 0.332043563f, 0.333723858f, 0.335404153f, 0.337084447f, 0.338925824f, 0.3407672f,
        0.342608576f, 0.344449953f, 0.346291329f, 0.348132705f, 0.349974082f, 0.351815458f, 0.353656834f, 0.355498211f,
        0.358067014f, 0.360635818f, 0.363204621f, 0.365773425f, 0.368342228f, 0.370911032f, 0.373479835f, 0.376048639f,
        0.378617442f, 0.381186246f, 0.38498025f, 0.388774254f, 0.392568259f, 0.396362263f, 0.400156268f, 0.403950272f,
        0.407744276f, 0.411538281f, 0.415332285f, 0.419126289f, 0.423809552f, 0.428492814f, 0.433176076f, 0.437859339f,
        0.442542601f, 0.447225863f, 0.451909125f, 0.456592388f, 0.46127565f, 0.465958912f, 0.470411019f, 0.474863126f,
        0.479315233f, 0.48376734f, 0.488219447f, 0.492671554f, 0.497123661f, 0.501575768f, 0.506027875f, 0.510479982f,
        0.514012641f, 0.5175453f, 0.521077959f, 0.524610618f, 0.528143276f, 0.531675935f, 0.535208594f, 0.538741253f,
        0.542273911f, 0.54580657f, 0.547945156f, 0.550083742f, 0.552222328f, 0.554360914f, 0.5564995f, 0.558638086f,
        0.560776672f, 0.562915258f, 0.565053844f, 0.56719243f, 0.567899658f, 0.568606886f, 0.569314114f, 0.570021342f,
        0.57072857f, 0.571435798f, 0.572143026f, 0.572850254f, 0.573557482f, 0.574264711f, 0.573746582f, 0.573228454f,
        0.572710326f, 0.572192198f, 0.57167407f, 0.571155942f, 0.570637814f, 0.570119686f, 0.569601558f, 0.56908343f,
        0.567243478f, 0.565403526f, 0.563563575f, 0.561723623f, 0.559883671f, 0.558043719f, 0.556203768f, 0.554363816f,
        0.552523864f, 0.550683912f, 0.547966075f, 0.545248237f, 0.542530399f, 0.539812561f, 0.537094724f, 0.534376886f,
        0.531659048f, 0.528941211f, 0.526223373f, 0.523505535f, 0.519998335f, 0.516491135f, 0.512983935f, 0.509476735f,
        0.505969535f, 0.502462335f, 0.498955135f, 0.495447935f, 0.491940735f, 0.488433535f, 0.48411129f, 0.479789046f,
        0.475466801f, 0.471144556f, 0.466822311f, 0.462500067f, 0.458177822f, 0.453855577f, 0.449533332f, 0.445211088f,
        0.440677259f, 0.43614343f, 0.431609601f, 0.427075772f, 0.422541943f, 0.418008114f, 0.413474285f, 0.408940456f,
        0.404406627f, 0.399872798f, 0.394928727f, 0.389984656f, 0.385040585f, 0.380096514f, 0.375152443f, 0.370208372f,
        0.365264301f, 0.36032023f, 0.355376159f, 0.350432088f, 0.345327729f, 0.34022337f, 0.335119011f, 0.330014653f,
        0.324910294f, 0.319805935f, 0.314701576f, 0.309597218f, 0.304492859f, 0.2993885f, 0.294693059f, 0.289997618f,
        0.285302176f, 0.280606735f, 0.275911294f, 0.271215853f, 0.266520411f, 0.26182497f, 0.257129529f, 0.252434088f,
        0.249286877f, 0.246139667f, 0.242992456f, 0.239845246f, 0.236698035f, 0.233550825f, 0.230403614f, 0.227256404f,
        0.224109193f, 0.220961982f, 0.219297194f, 0.217632405f, 0.215967617f, 0.214302828f, 0.212638039f, 0.210973251f,
        0.209308462f, 0.207643674f, 0.205978885f, 0.204314096f, 0.20346194f, 0.202609784f, 0.201757628f, 0.200905472f,
        0.200053316f, 0.19920116f, 0.198349004f, 0.197496847f, 0.196644691f, 0.195792535f, 0.195300989f, 0.194809442f,
        0.194317896f, 0.193826349f, 0.193334803f, 0.192843256f, 0.19235171f, 0.191860163f, 0.191368617f, 0.19087707f,
        0.190612275f, 0.190347479f, 0.190082683f, 0.189817888f, 0.189553092f, 0.189288296f, 0.189023501f, 0.188758705f,
        0.18849391f, 0.188229114f, 0.188478167f, 0.188727219f, 0.188976272f, 0.189225325f, 0.189474377f, 0.18972343f,
        0.189972482f, 0.190221535f, 0.190470588f, 0.19071964f, 0.191589889f, 0.192460137f, 0.193330385f, 0.194200633f,
        0.195070882f, 0.19594113f, 0.196811378f, 0.197681626f, 0.198551875f, 0.199422123f, 0.200638499f, 0.201854875f,
        0.203071252f, 0.204287628f, 0.205504004f, 0.206720381f, 0.207936757f, 0.209153133f, 0.21036951f, 0.211585886f,
        0.212737445f, 0.213889004f, 0.215040562f, 0.216192121f, 0.21734368f, 0.218495239f, 0.219646797f, 0.220798356f,
        0.221949915f, 0.223101474f, 0.223955268f, 0.224809061f, 0.225662855f, 0.226516649f, 0.227370443f, 0.228224237f,
        0.229078031f, 0.229931825f, 0.230785618f, 0.231639412f, 0.231807176f, 0.23197494f, 0.232142705f, 0.232310469f,
        0.232478233f, 0.232645997f, 0.232813761f, 0.232981525f, 0.233149289f, 0.233317053f, 0.232926626f, 0.232536198f,
        0.232145771f, 0.231755343f, 0.231364916f, 0.230974488f, 0.23058406f, 0.230193633f, 0.229803205f, 0.229412778f,
        0.229406829f, 0.22940088f, 0.229394931f, 0.229388982f, 0.229383033f, 0.229377084f, 0.229371136f, 0.229365187f,
        0.229359238f, 0.229353289f, 0.20641796f, 0.183482631f, 0.160547302f, 0.137611973f, 0.114676644f, 0.091741316f,
        0.068805987f, 0.045870658f, 0.022935329f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth11({0.053807993f, 0.053796048f, 0.053784102f, 0.053772156f, 0.053760211f,
        0.053748265f, 0.053736319f, 0.053724374f, 0.053712428f, 0.053700482f, 0.053688537f, 0.053645508f, 0.053602479f,
        0.053559451f, 0.053516422f, 0.053473393f, 0.053430365f, 0.053387336f, 0.053344307f, 0.053301279f, 0.05325825f,
        0.053302142f, 0.053346034f, 0.053389926f, 0.053433817f, 0.053477709f, 0.053521601f, 0.053565493f, 0.053609385f,
        0.053653277f, 0.053697168f, 0.053729406f, 0.053761644f, 0.053793882f, 0.05382612f, 0.053858358f, 0.053890596f,
        0.053922834f, 0.053955072f, 0.05398731f, 0.054019548f, 0.054069433f, 0.054119318f, 0.054169203f, 0.054219088f,
        0.054268973f, 0.054318858f, 0.054368743f, 0.054418628f, 0.054468513f, 0.054518398f, 0.054561244f, 0.05460409f,
        0.054646936f, 0.054689782f, 0.054732628f, 0.054775473f, 0.054818319f, 0.054861165f, 0.054904011f, 0.054946857f,
        0.054967891f, 0.054988926f, 0.055009961f, 0.055030995f, 0.05505203f, 0.055073065f, 0.055094099f, 0.055115134f,
        0.055136169f, 0.055157203f, 0.055209092f, 0.05526098f, 0.055312868f, 0.055364756f, 0.055416644f, 0.055468532f,
        0.055520421f, 0.055572309f, 0.055624197f, 0.055676085f, 0.0557723f, 0.055868515f, 0.055964731f, 0.056060946f,
        0.056157161f, 0.056253376f, 0.056349591f, 0.056445806f, 0.056542022f, 0.056638237f, 0.056814154f, 0.056990071f,
        0.057165988f, 0.057341905f, 0.057517822f, 0.057693739f, 0.057869656f, 0.058045573f, 0.05822149f, 0.058397407f,
        0.058679266f, 0.058961125f, 0.059242984f, 0.059524843f, 0.059806703f, 0.060088562f, 0.060370421f, 0.06065228f,
        0.060934139f, 0.061215998f, 0.061916158f, 0.062616317f, 0.063316477f, 0.064016636f, 0.064716796f, 0.065416955f,
        0.066117115f, 0.066817274f, 0.067517434f, 0.068217593f, 0.070337773f, 0.072457953f, 0.074578133f, 0.076698313f,
        0.078818493f, 0.080938673f, 0.083058853f, 0.085179033f, 0.087299213f, 0.089419393f, 0.092938903f, 0.096458412f,
        0.099977921f, 0.10349743f, 0.107016939f, 0.110536448f, 0.114055958f, 0.117575467f, 0.121094976f, 0.124614485f,
        0.12750334f, 0.130392195f, 0.13328105f, 0.136169904f, 0.139058759f, 0.141947614f, 0.144836469f, 0.147725324f,
        0.150614179f, 0.153503033f, 0.155532077f, 0.15756112f, 0.159590163f, 0.161619207f, 0.16364825f, 0.165677293f,
        0.167706337f, 0.16973538f, 0.171764423f, 0.173793467f, 0.176358407f, 0.178923347f, 0.181488287f, 0.184053227f,
        0.186618167f, 0.189183107f, 0.191748047f, 0.194312987f, 0.196877927f, 0.199442867f, 0.204325613f, 0.20920836f,
        0.214091107f, 0.218973853f, 0.2238566f, 0.228739347f, 0.233622093f, 0.23850484f, 0.243387587f, 0.248270333f,
        0.256985634f, 0.265700935f, 0.274416236f, 0.283131537f, 0.291846838f, 0.300562138f, 0.309277439f, 0.31799274f,
        0.326708041f, 0.335423342f, 0.34628038f, 0.357137418f, 0.367994457f, 0.378851495f, 0.389708533f, 0.400565572f,
        0.41142261f, 0.422279648f, 0.433136687f, 0.443993725f, 0.453441677f, 0.462889628f, 0.47233758f, 0.481785532f,
        0.491233483f, 0.500681435f, 0.510129387f, 0.519577338f, 0.52902529f, 0.538473242f, 0.543293256f, 0.54811327f,
        0.552933284f, 0.557753298f, 0.562573313f, 0.567393327f, 0.572213341f, 0.577033355f, 0.581853369f, 0.586673383f,
        0.587489945f, 0.588306507f, 0.589123068f, 0.58993963f, 0.590756192f, 0.591572753f, 0.592389315f, 0.593205877f,
        0.594022438f, 0.594839f, 0.594414013f, 0.593989027f, 0.59356404f, 0.593139053f, 0.592714067f, 0.59228908f,
        0.591864093f, 0.591439107f, 0.59101412f, 0.590589133f, 0.590191936f, 0.589794738f, 0.589397541f, 0.589000343f,
        0.588603146f, 0.588205948f, 0.587808751f, 0.587411553f, 0.587014356f, 0.586617158f, 0.586372875f, 0.586128592f,
        0.585884308f, 0.585640025f, 0.585395742f, 0.585151458f, 0.584907175f, 0.584662892f, 0.584418608f, 0.584174325f,
        0.584142693f, 0.58411106f, 0.584079428f, 0.584047795f, 0.584016163f, 0.58398453f, 0.583952898f, 0.583921265f,
        0.583889633f, 0.583858f, 0.584446924f, 0.585035848f, 0.585624773f, 0.586213697f, 0.586802621f, 0.587391545f,
        0.587980469f, 0.588569393f, 0.589158318f, 0.589747242f, 0.591023664f, 0.592300087f, 0.593576509f, 0.594852932f,
        0.596129354f, 0.597405777f, 0.598682199f, 0.599958622f, 0.601235044f, 0.602511467f, 0.604299346f, 0.606087225f,
        0.607875104f, 0.609662983f, 0.611450863f, 0.613238742f, 0.615026621f, 0.6168145f, 0.618602379f, 0.620390258f,
        0.622230919f, 0.62407158f, 0.625912241f, 0.627752902f, 0.629593563f, 0.631434223f, 0.633274884f, 0.635115545f,
        0.636956206f, 0.638796867f, 0.640398172f, 0.641999477f, 0.643600782f, 0.645202087f, 0.646803392f, 0.648404697f,
        0.650006002f, 0.651607307f, 0.653208612f, 0.654809917f, 0.655583534f, 0.656357152f, 0.657130769f, 0.657904387f,
        0.658678004f, 0.659451622f, 0.660225239f, 0.660998857f, 0.661772474f, 0.662546092f, 0.662546326f, 0.66254656f,
        0.662546794f, 0.662547028f, 0.662547263f, 0.662547497f, 0.662547731f, 0.662547965f, 0.662548199f, 0.662548433f,
        0.662975086f, 0.663401738f, 0.663828391f, 0.664255043f, 0.664681696f, 0.665108348f, 0.665535001f, 0.665961653f,
        0.666388306f, 0.666814958f, 0.600133463f, 0.533451967f, 0.466770471f, 0.400088975f, 0.333407479f, 0.266725983f,
        0.200044488f, 0.133362992f, 0.066681496f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth12({0.122364925f, 0.126576671f, 0.130788417f, 0.135000163f, 0.139211908f,
        0.143423654f, 0.1476354f, 0.151847146f, 0.156058892f, 0.160270638f, 0.164482383f, 0.170884575f, 0.177286767f,
        0.183688958f, 0.19009115f, 0.196493342f, 0.202895533f, 0.209297725f, 0.215699917f, 0.222102108f, 0.2285043f,
        0.234261801f, 0.240019302f, 0.245776803f, 0.251534303f, 0.257291804f, 0.263049305f, 0.268806806f, 0.274564307f,
        0.280321808f, 0.286079308f, 0.290201841f, 0.294324373f, 0.298446906f, 0.302569438f, 0.306691971f, 0.310814503f,
        0.314937036f, 0.319059568f, 0.323182101f, 0.327304633f, 0.330682438f, 0.334060242f, 0.337438046f, 0.34081585f,
        0.344193654f, 0.347571458f, 0.350949263f, 0.354327067f, 0.357704871f, 0.361082675f, 0.363731076f, 0.366379477f,
        0.369027878f, 0.371676278f, 0.374324679f, 0.37697308f, 0.379621481f, 0.382269882f, 0.384918283f, 0.387566683f,
        0.388773437f, 0.38998019f, 0.391186943f, 0.392393697f, 0.39360045f, 0.394807203f, 0.396013957f, 0.39722071f,
        0.398427463f, 0.399634217f, 0.398827703f, 0.398021188f, 0.397214674f, 0.39640816f, 0.395601646f, 0.394795132f,
        0.393988618f, 0.393182103f, 0.392375589f, 0.391569075f, 0.388654937f, 0.385740798f, 0.38282666f, 0.379912522f,
        0.376998383f, 0.374084245f, 0.371170107f, 0.368255968f, 0.36534183f, 0.362427692f, 0.357796759f, 0.353165827f,
        0.348534894f, 0.343903962f, 0.339273029f, 0.334642097f, 0.330011164f, 0.325380232f, 0.320749299f, 0.316118367f,
        0.310530894f, 0.304943422f, 0.299355949f, 0.293768477f, 0.288181004f, 0.282593532f, 0.277006059f, 0.271418587f,
        0.265831114f, 0.260243642f, 0.255077422f, 0.249911202f, 0.244744982f, 0.239578762f, 0.234412542f, 0.229246322f,
        0.224080102f, 0.218913882f, 0.213747662f, 0.208581442f, 0.204554178f, 0.200526913f, 0.196499649f, 0.192472385f,
        0.188445121f, 0.184417857f, 0.180390593f, 0.176363328f, 0.172336064f, 0.1683088f, 0.165246354f, 0.162183908f,
        0.159121463f, 0.156059017f, 0.152996571f, 0.149934125f, 0.146871679f, 0.143809233f, 0.140746788f, 0.137684342f,
        0.135571868f, 0.133459393f, 0.131346919f, 0.129234445f, 0.127121971f, 0.125009497f, 0.122897023f, 0.120784548f,
        0.118672074f, 0.1165596f, 0.115328577f, 0.114097553f, 0.11286653f, 0.111635507f, 0.110404483f, 0.10917346f,
        0.107942437f, 0.106711413f, 0.10548039f, 0.104249367f, 0.103460951f, 0.102672535f, 0.10188412f, 0.101095704f,
        0.100307288f, 0.099518873f, 0.098730457f, 0.097942041f, 0.097153626f, 0.09636521f, 0.095708316f, 0.095051423f,
        0.094394529f, 0.093737635f, 0.093080742f, 0.092423848f, 0.091766954f, 0.091110061f, 0.090453167f, 0.089796273f,
        0.08936739f, 0.088938507f, 0.088509623f, 0.08808074f, 0.087651857f, 0.087222973f, 0.08679409f, 0.086365207f,
        0.085936323f, 0.08550744f, 0.085329032f, 0.085150623f, 0.084972215f, 0.084793807f, 0.084615398f, 0.08443699f,
        0.084258582f, 0.084080173f, 0.083901765f, 0.083723357f, 0.083747411f, 0.083771464f, 0.083795518f, 0.083819572f,
        0.083843626f, 0.08386768f, 0.083891734f, 0.083915787f, 0.083939841f, 0.083963895f, 0.083999012f, 0.084034128f,
        0.084069245f, 0.084104362f, 0.084139478f, 0.084174595f, 0.084209712f, 0.084244828f, 0.084279945f, 0.084315062f,
        0.084294631f, 0.084274199f, 0.084253768f, 0.084233337f, 0.084212906f, 0.084192475f, 0.084172044f, 0.084151612f,
        0.084131181f, 0.08411075f, 0.084085745f, 0.08406074f, 0.084035736f, 0.084010731f, 0.083985726f, 0.083960721f,
        0.083935716f, 0.083910711f, 0.083885707f, 0.083860702f, 0.083991644f, 0.084122585f, 0.084253527f, 0.084384469f,
        0.084515411f, 0.084646353f, 0.084777295f, 0.084908236f, 0.085039178f, 0.08517012f, 0.085629799f, 0.086089477f,
        0.086549156f, 0.087008834f, 0.087468513f, 0.087928191f, 0.08838787f, 0.088847548f, 0.089307227f, 0.089766905f,
        0.090575043f, 0.091383181f, 0.092191319f, 0.092999456f, 0.093807594f, 0.094615732f, 0.09542387f, 0.096232008f,
        0.097040146f, 0.097848283f, 0.098975833f, 0.100103383f, 0.101230933f, 0.102358483f, 0.103486033f, 0.104613583f,
        0.105741133f, 0.106868683f, 0.107996233f, 0.109123783f, 0.110557826f, 0.111991868f, 0.113425911f, 0.114859953f,
        0.116293996f, 0.117728038f, 0.119162081f, 0.120596123f, 0.122030166f, 0.123464208f, 0.125387027f, 0.127309845f,
        0.129232663f, 0.131155482f, 0.1330783f, 0.135001118f, 0.136923937f, 0.138846755f, 0.140769573f, 0.142692392f,
        0.145353644f, 0.148014897f, 0.150676149f, 0.153337402f, 0.155998654f, 0.158659907f, 0.161321159f, 0.163982412f,
        0.166643664f, 0.169304917f, 0.172839598f, 0.176374278f, 0.179908959f, 0.18344364f, 0.186978321f, 0.190513002f,
        0.194047683f, 0.197582363f, 0.201117044f, 0.204651725f, 0.208581847f, 0.212511968f, 0.21644209f, 0.220372212f,
        0.224302333f, 0.228232455f, 0.232162577f, 0.236092698f, 0.24002282f, 0.243952942f, 0.248276378f, 0.252599813f,
        0.256923249f, 0.261246685f, 0.265570121f, 0.269893557f, 0.274216993f, 0.278540428f, 0.282863864f, 0.2871873f,
        0.291717082f, 0.296246863f, 0.300776645f, 0.305306427f, 0.309836208f, 0.31436599f, 0.318895772f, 0.323425553f,
        0.327955335f, 0.332485117f, 0.299236605f, 0.265988093f, 0.232739582f, 0.19949107f, 0.166242558f, 0.132994047f,
        0.099745535f, 0.066497023f, 0.033248512f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth13({0.096004948f, 0.097870687f, 0.099736425f, 0.101602164f, 0.103467902f,
        0.105333641f, 0.107199379f, 0.109065118f, 0.110930856f, 0.112796595f, 0.114662333f, 0.116254343f, 0.117846353f,
        0.119438363f, 0.121030373f, 0.122622383f, 0.124214393f, 0.125806403f, 0.127398413f, 0.128990423f, 0.130582433f,
        0.131032368f, 0.131482303f, 0.131932238f, 0.132382173f, 0.132832108f, 0.133282043f, 0.133731978f, 0.134181913f,
        0.134631848f, 0.135081783f, 0.134918414f, 0.134755045f, 0.134591676f, 0.134428307f, 0.134264938f, 0.134101568f,
        0.133938199f, 0.13377483f, 0.133611461f, 0.133448092f, 0.13326246f, 0.133076828f, 0.132891197f, 0.132705565f,
        0.132519933f, 0.132334302f, 0.13214867f, 0.131963038f, 0.131777407f, 0.131591775f, 0.131453653f, 0.13131553f,
        0.131177408f, 0.131039285f, 0.130901163f, 0.13076304f, 0.130624918f, 0.130486795f, 0.130348673f, 0.13021055f,
        0.130000184f, 0.129789818f, 0.129579453f, 0.129369087f, 0.129158721f, 0.128948355f, 0.128737989f, 0.128527623f,
        0.128317258f, 0.128106892f, 0.127801283f, 0.127495673f, 0.127190064f, 0.126884455f, 0.126578846f, 0.126273237f,
        0.125967628f, 0.125662018f, 0.125356409f, 0.1250508f, 0.124593303f, 0.124135807f, 0.12367831f, 0.123220813f,
        0.122763317f, 0.12230582f, 0.121848323f, 0.121390827f, 0.12093333f, 0.120475833f, 0.119940628f, 0.119405423f,
        0.118870218f, 0.118335013f, 0.117799808f, 0.117264603f, 0.116729398f, 0.116194193f, 0.115658988f, 0.115123783f,
        0.11459622f, 0.114068657f, 0.113541093f, 0.11301353f, 0.112485967f, 0.111958403f, 0.11143084f, 0.110903277f,
        0.110375713f, 0.10984815f, 0.109357321f, 0.108866492f, 0.108375663f, 0.107884833f, 0.107394004f, 0.106903175f,
        0.106412346f, 0.105921517f, 0.105430688f, 0.104939858f, 0.104427923f, 0.103915987f, 0.103404051f, 0.102892116f,
        0.10238018f, 0.101868244f, 0.101356309f, 0.100844373f, 0.100332437f, 0.099820502f, 0.099354759f, 0.098889016f,
        0.098423273f, 0.09795753f, 0.097491787f, 0.097026044f, 0.096560301f, 0.096094558f, 0.095628815f, 0.095163072f,
        0.094912235f, 0.094661398f, 0.094410561f, 0.094159724f, 0.093908888f, 0.093658051f, 0.093407214f, 0.093156377f,
        0.09290554f, 0.092654703f, 0.09263596f, 0.092617217f, 0.092598474f, 0.092579731f, 0.092560988f, 0.092542244f,
        0.092523501f, 0.092504758f, 0.092486015f, 0.092467272f, 0.092539468f, 0.092611665f, 0.092683862f, 0.092756058f,
        0.092828255f, 0.092900452f, 0.092972648f, 0.093044845f, 0.093117042f, 0.093189238f, 0.093490886f, 0.093792533f,
        0.094094181f, 0.094395828f, 0.094697476f, 0.094999123f, 0.095300771f, 0.095602418f, 0.095904066f, 0.096205713f,
        0.097396776f, 0.098587839f, 0.099778902f, 0.100969965f, 0.102161028f, 0.10335209f, 0.104543153f, 0.105734216f,
        0.106925279f, 0.108116342f, 0.112861428f, 0.117606515f, 0.122351602f, 0.127096688f, 0.131841775f, 0.136586862f,
        0.141331948f, 0.146077035f, 0.150822122f, 0.155567208f, 0.166549495f, 0.177531782f, 0.188514068f, 0.199496355f,
        0.210478642f, 0.221460928f, 0.232443215f, 0.243425502f, 0.254407788f, 0.265390075f, 0.278722221f, 0.292054367f,
        0.305386513f, 0.318718658f, 0.332050804f, 0.34538295f, 0.358715096f, 0.372047242f, 0.385379388f, 0.398711533f,
        0.408848528f, 0.418985522f, 0.429122516f, 0.43925951f, 0.449396504f, 0.459533498f, 0.469670493f, 0.479807487f,
        0.489944481f, 0.500081475f, 0.505705238f, 0.511329f, 0.516952763f, 0.522576525f, 0.528200288f, 0.53382405f,
        0.539447813f, 0.545071575f, 0.550695338f, 0.5563191f, 0.558632255f, 0.56094541f, 0.563258565f, 0.56557172f,
        0.567884875f, 0.57019803f, 0.572511185f, 0.57482434f, 0.577137495f, 0.57945065f, 0.580278085f, 0.58110552f,
        0.581932955f, 0.58276039f, 0.583587825f, 0.58441526f, 0.585242695f, 0.58607013f, 0.586897565f, 0.587725f,
        0.588015169f, 0.588305338f, 0.588595508f, 0.588885677f, 0.589175846f, 0.589466015f, 0.589756184f, 0.590046353f,
        0.590336523f, 0.590626692f, 0.590815022f, 0.591003352f, 0.591191682f, 0.591380012f, 0.591568342f, 0.591756672f,
        0.591945002f, 0.592133332f, 0.592321662f, 0.592509992f, 0.592704069f, 0.592898147f, 0.593092224f, 0.593286302f,
        0.593480379f, 0.593674457f, 0.593868534f, 0.594062612f, 0.594256689f, 0.594450767f, 0.594790606f, 0.595130445f,
        0.595470284f, 0.595810123f, 0.596149963f, 0.596489802f, 0.596829641f, 0.59716948f, 0.597509319f, 0.597849158f,
        0.598283257f, 0.598717355f, 0.599151453f, 0.599585552f, 0.60001965f, 0.600453748f, 0.600887847f, 0.601321945f,
        0.601756043f, 0.602190142f, 0.602661484f, 0.603132827f, 0.603604169f, 0.604075512f, 0.604546854f, 0.605018197f,
        0.605489539f, 0.605960882f, 0.606432224f, 0.606903567f, 0.607138388f, 0.607373208f, 0.607608029f, 0.60784285f,
        0.608077671f, 0.608312492f, 0.608547313f, 0.608782133f, 0.609016954f, 0.609251775f, 0.609222884f, 0.609193993f,
        0.609165103f, 0.609136212f, 0.609107321f, 0.60907843f, 0.609049539f, 0.609020648f, 0.608991758f, 0.608962867f,
        0.609090122f, 0.609217377f, 0.609344632f, 0.609471887f, 0.609599142f, 0.609726397f, 0.609853652f, 0.609980907f,
        0.610108162f, 0.610235417f, 0.549211875f, 0.488188333f, 0.427164792f, 0.36614125f, 0.305117708f, 0.244094167f,
        0.183070625f, 0.122047083f, 0.061023542f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    static constexpr Spectrum Macbeth14({0.091993695f, 0.094395723f, 0.096797751f, 0.099199779f, 0.101601807f,
        0.104003835f, 0.106405863f, 0.108807891f, 0.111209919f, 0.113611947f, 0.116013975f, 0.118973427f, 0.121932878f,
        0.12489233f, 0.127851782f, 0.130811233f, 0.133770685f, 0.136730137f, 0.139689588f, 0.14264904f, 0.145608492f,
        0.147901087f, 0.150193682f, 0.152486277f, 0.154778872f, 0.157071467f, 0.159364062f, 0.161656657f, 0.163949252f,
        0.166241847f, 0.168534442f, 0.169527616f, 0.17052079f, 0.171513964f, 0.172507138f, 0.173500313f, 0.174493487f,
        0.175486661f, 0.176479835f, 0.177473009f, 0.178466183f, 0.177920935f, 0.177375687f, 0.176830438f, 0.17628519f,
        0.175739942f, 0.175194693f, 0.174649445f, 0.174104197f, 0.173558948f, 0.1730137f, 0.171508923f, 0.170004147f,
        0.16849937f, 0.166994593f, 0.165489817f, 0.16398504f, 0.162480263f, 0.160975487f, 0.15947071f, 0.157965933f,
        0.156047798f, 0.154129663f, 0.152211528f, 0.150293393f, 0.148375258f, 0.146457123f, 0.144538988f, 0.142620853f,
        0.140702718f, 0.138784583f, 0.136819271f, 0.134853958f, 0.132888646f, 0.130923333f, 0.128958021f, 0.126992708f,
        0.125027396f, 0.123062083f, 0.121096771f, 0.119131458f, 0.117357856f, 0.115584253f, 0.113810651f, 0.112037048f,
        0.110263446f, 0.108489843f, 0.106716241f, 0.104942638f, 0.103169036f, 0.101395433f, 0.099950526f, 0.098505619f,
        0.097060712f, 0.095615805f, 0.094170898f, 0.09272599f, 0.091281083f, 0.089836176f, 0.088391269f, 0.086946362f,
        0.085769857f, 0.084593351f, 0.083416846f, 0.082240341f, 0.081063836f, 0.079887331f, 0.078710826f, 0.07753432f,
        0.076357815f, 0.07518131f, 0.074272676f, 0.073364041f, 0.072455407f, 0.071546773f, 0.070638138f, 0.069729504f,
        0.06882087f, 0.067912235f, 0.067003601f, 0.066094967f, 0.065517321f, 0.064939675f, 0.064362029f, 0.063784383f,
        0.063206737f, 0.062629091f, 0.062051445f, 0.061473799f, 0.060896153f, 0.060318507f, 0.059932656f, 0.059546804f,
        0.059160953f, 0.058775102f, 0.058389251f, 0.0580034f, 0.057617549f, 0.057231697f, 0.056845846f, 0.056459995f,
        0.056125823f, 0.055791651f, 0.05545748f, 0.055123308f, 0.054789136f, 0.054454964f, 0.054120792f, 0.05378662f,
        0.053452449f, 0.053118277f, 0.052927354f, 0.052736432f, 0.05254551f, 0.052354587f, 0.052163665f, 0.051972743f,
        0.05178182f, 0.051590898f, 0.051399976f, 0.051209053f, 0.05121224f, 0.051215426f, 0.051218612f, 0.051221798f,
        0.051224984f, 0.05122817f, 0.051231357f, 0.051234543f, 0.051237729f, 0.051240915f, 0.051311842f, 0.051382769f,
        0.051453697f, 0.051524624f, 0.051595551f, 0.051666478f, 0.051737405f, 0.051808332f, 0.05187926f, 0.051950187f,
        0.05194242f, 0.051934653f, 0.051926886f, 0.051919119f, 0.051911352f, 0.051903585f, 0.051895818f, 0.051888051f,
        0.051880284f, 0.051872517f, 0.051805078f, 0.051737639f, 0.051670201f, 0.051602762f, 0.051535323f, 0.051467885f,
        0.051400446f, 0.051333007f, 0.051265569f, 0.05119813f, 0.051320318f, 0.051442506f, 0.051564695f, 0.051686883f,
        0.051809071f, 0.051931259f, 0.052053447f, 0.052175635f, 0.052297824f, 0.052420012f, 0.053019264f, 0.053618516f,
        0.054217768f, 0.05481702f, 0.055416273f, 0.056015525f, 0.056614777f, 0.057214029f, 0.057813281f, 0.058412533f,
        0.0598892f, 0.061365867f, 0.062842534f, 0.064319201f, 0.065795868f, 0.067272535f, 0.068749202f, 0.070225869f,
        0.071702536f, 0.073179203f, 0.075412842f, 0.077646481f, 0.07988012f, 0.082113759f, 0.084347398f, 0.086581037f,
        0.088814676f, 0.091048315f, 0.093281954f, 0.095515593f, 0.097856567f, 0.100197541f, 0.102538515f, 0.104879489f,
        0.107220463f, 0.109561437f, 0.111902411f, 0.114243385f, 0.116584359f, 0.118925333f, 0.121172248f, 0.123419163f,
        0.125666078f, 0.127912993f, 0.130159908f, 0.132406823f, 0.134653738f, 0.136900653f, 0.139147568f, 0.141394483f,
        0.143809122f, 0.14622376f, 0.148638398f, 0.151053037f, 0.153467675f, 0.155882313f, 0.158296952f, 0.16071159f,
        0.163126228f, 0.165540867f, 0.168391362f, 0.171241857f, 0.174092352f, 0.176942847f, 0.179793342f, 0.182643837f,
        0.185494332f, 0.188344827f, 0.191195322f, 0.194045817f, 0.197347446f, 0.200649075f, 0.203950704f, 0.207252333f,
        0.210553963f, 0.213855592f, 0.217157221f, 0.22045885f, 0.223760479f, 0.227062108f, 0.230895335f, 0.234728562f,
        0.238561788f, 0.242395015f, 0.246228242f, 0.250061468f, 0.253894695f, 0.257727922f, 0.261561148f, 0.265394375f,
        0.269747375f, 0.274100375f, 0.278453375f, 0.282806375f, 0.287159375f, 0.291512375f, 0.295865375f, 0.300218375f,
        0.304571375f, 0.308924375f, 0.313487192f, 0.318050008f, 0.322612825f, 0.327175642f, 0.331738458f, 0.336301275f,
        0.340864092f, 0.345426908f, 0.349989725f, 0.354552542f, 0.358673831f, 0.36279512f, 0.366916409f, 0.371037698f,
        0.375158988f, 0.379280277f, 0.383401566f, 0.387522855f, 0.391644144f, 0.395765433f, 0.399772599f, 0.403779765f,
        0.407786931f, 0.411794097f, 0.415801263f, 0.419808428f, 0.423815594f, 0.42782276f, 0.431829926f, 0.435837092f,
        0.440100286f, 0.44436348f, 0.448626674f, 0.452889868f, 0.457153063f, 0.461416257f, 0.465679451f, 0.469942645f,
        0.474205839f, 0.478469033f, 0.43062213f, 0.382775227f, 0.334928323f, 0.28708142f, 0.239234517f, 0.191387613f,
        0.14354071f, 0.095693807f, 0.047846903f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth15({0.061031217f, 0.061052864f, 0.061074512f, 0.061096159f, 0.061117807f,
        0.061139454f, 0.061161102f, 0.061182749f, 0.061204397f, 0.061226044f, 0.061247692f, 0.061315304f, 0.061382916f,
        0.061450528f, 0.06151814f, 0.061585752f, 0.061653364f, 0.061720976f, 0.061788588f, 0.0618562f, 0.061923812f,
        0.062022906f, 0.062122f, 0.062221094f, 0.062320188f, 0.062419282f, 0.062518376f, 0.06261747f, 0.062716564f,
        0.062815658f, 0.062914752f, 0.063020304f, 0.063125856f, 0.063231408f, 0.06333696f, 0.063442512f, 0.063548064f,
        0.063653616f, 0.063759168f, 0.06386472f, 0.063970272f, 0.064166072f, 0.064361873f, 0.064557674f, 0.064753474f,
        0.064949275f, 0.065145076f, 0.065340876f, 0.065536677f, 0.065732478f, 0.065928278f, 0.066256456f, 0.066584634f,
        0.066912811f, 0.067240989f, 0.067569167f, 0.067897344f, 0.068225522f, 0.0685537f, 0.068881877f, 0.069210055f,
        0.069761582f, 0.07031311f, 0.070864637f, 0.071416164f, 0.071967692f, 0.072519219f, 0.073070746f, 0.073622274f,
        0.074173801f, 0.074725328f, 0.075801767f, 0.076878206f, 0.077954645f, 0.079031084f, 0.080107523f, 0.081183961f,
        0.0822604f, 0.083336839f, 0.084413278f, 0.085489717f, 0.087446972f, 0.089404227f, 0.091361482f, 0.093318737f,
        0.095275992f, 0.097233247f, 0.099190502f, 0.101147757f, 0.103105012f, 0.105062267f, 0.108422853f, 0.111783438f,
        0.115144024f, 0.11850461f, 0.121865196f, 0.125225782f, 0.128586368f, 0.131946953f, 0.135307539f, 0.138668125f,
        0.144010023f, 0.14935192f, 0.154693818f, 0.160035715f, 0.165377613f, 0.17071951f, 0.176061408f, 0.181403305f,
        0.186745203f, 0.1920871f, 0.199951233f, 0.207815367f, 0.2156795f, 0.223543633f, 0.231407767f, 0.2392719f,
        0.247136033f, 0.255000167f, 0.2628643f, 0.270728433f, 0.281266263f, 0.291804092f, 0.302341921f, 0.31287975f,
        0.323417579f, 0.333955408f, 0.344493238f, 0.355031067f, 0.365568896f, 0.376106725f, 0.386073688f, 0.39604065f,
        0.406007613f, 0.415974575f, 0.425941538f, 0.4359085f, 0.445875463f, 0.455842425f, 0.465809388f, 0.47577635f,
        0.481321102f, 0.486865853f, 0.492410605f, 0.497955357f, 0.503500108f, 0.50904486f, 0.514589612f, 0.520134363f,
        0.525679115f, 0.531223867f, 0.533017229f, 0.534810592f, 0.536603954f, 0.538397317f, 0.540190679f, 0.541984042f,
        0.543777404f, 0.545570767f, 0.547364129f, 0.549157492f, 0.548813237f, 0.548468982f, 0.548124727f, 0.547780472f,
        0.547436217f, 0.547091962f, 0.546747707f, 0.546403452f, 0.546059197f, 0.545714942f, 0.543950689f, 0.542186437f,
        0.540422184f, 0.538657932f, 0.536893679f, 0.535129427f, 0.533365174f, 0.531600922f, 0.529836669f, 0.528072417f,
        0.525710778f, 0.523349138f, 0.520987499f, 0.51862586f, 0.516264221f, 0.513902582f, 0.511540943f, 0.509179303f,
        0.506817664f, 0.504456025f, 0.501062431f, 0.497668837f, 0.494275243f, 0.490881648f, 0.487488054f, 0.48409446f,
        0.480700866f, 0.477307272f, 0.473913678f, 0.470520083f, 0.466231866f, 0.461943648f, 0.457655431f, 0.453367213f,
        0.449078996f, 0.444790778f, 0.440502561f, 0.436214343f, 0.431926126f, 0.427637908f, 0.422998625f, 0.418359342f,
        0.413720058f, 0.409080775f, 0.404441492f, 0.399802208f, 0.395162925f, 0.390523642f, 0.385884358f, 0.381245075f,
        0.377800438f, 0.374355802f, 0.370911165f, 0.367466528f, 0.364021892f, 0.360577255f, 0.357132618f, 0.353687982f,
        0.350243345f, 0.346798708f, 0.344863173f, 0.342927637f, 0.340992101f, 0.339056565f, 0.337121029f, 0.335185493f,
        0.333249958f, 0.331314422f, 0.329378886f, 0.32744335f, 0.326469735f, 0.32549612f, 0.324522505f, 0.32354889f,
        0.322575275f, 0.32160166f, 0.320628045f, 0.31965443f, 0.318680815f, 0.3177072f, 0.317183179f, 0.316659158f,
        0.316135138f, 0.315611117f, 0.315087096f, 0.314563075f, 0.314039054f, 0.313515033f, 0.312991013f, 0.312466992f,
        0.312213849f, 0.311960707f, 0.311707564f, 0.311454422f, 0.311201279f, 0.310948137f, 0.310694994f, 0.310441852f,
        0.310188709f, 0.309935567f, 0.31038308f, 0.310830593f, 0.311278107f, 0.31172562f, 0.312173133f, 0.312620647f,
        0.31306816f, 0.313515673f, 0.313963187f, 0.3144107f, 0.315710298f, 0.317009897f, 0.318309495f, 0.319609093f,
        0.320908692f, 0.32220829f, 0.323507888f, 0.324807487f, 0.326107085f, 0.327406683f, 0.329188613f, 0.330970543f,
        0.332752473f, 0.334534403f, 0.336316333f, 0.338098263f, 0.339880193f, 0.341662123f, 0.343444053f, 0.345225983f,
        0.346958421f, 0.348690858f, 0.350423296f, 0.352155733f, 0.353888171f, 0.355620608f, 0.357353046f, 0.359085483f,
        0.360817921f, 0.362550358f, 0.363916925f, 0.365283492f, 0.366650058f, 0.368016625f, 0.369383192f, 0.370749758f,
        0.372116325f, 0.373482892f, 0.374849458f, 0.376216025f, 0.376648712f, 0.377081398f, 0.377514085f, 0.377946772f,
        0.378379458f, 0.378812145f, 0.379244832f, 0.379677518f, 0.380110205f, 0.380542892f, 0.380256029f, 0.379969167f,
        0.379682304f, 0.379395442f, 0.379108579f, 0.378821717f, 0.378534854f, 0.378247992f, 0.377961129f, 0.377674267f,
        0.377847467f, 0.378020667f, 0.378193867f, 0.378367067f, 0.378540267f, 0.378713467f, 0.378886667f, 0.379059867f,
        0.379233067f, 0.379406267f, 0.34146564f, 0.303525013f, 0.265584387f, 0.22764376f, 0.189703133f, 0.151762507f,
        0.11382188f, 0.075881253f, 0.037940627f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth16({0.06281745f, 0.062819804f, 0.062822157f, 0.062824511f, 0.062826865f,
        0.062829218f, 0.062831572f, 0.062833926f, 0.062836279f, 0.062838633f, 0.062840987f, 0.062890569f, 0.062940151f,
        0.062989733f, 0.063039315f, 0.063088897f, 0.063138479f, 0.063188061f, 0.063237643f, 0.063287225f, 0.063336807f,
        0.063357351f, 0.063377895f, 0.063398439f, 0.063418983f, 0.063439528f, 0.063460072f, 0.063480616f, 0.06350116f,
        0.063521704f, 0.063542248f, 0.063558928f, 0.063575608f, 0.063592288f, 0.063608968f, 0.063625648f, 0.063642327f,
        0.063659007f, 0.063675687f, 0.063692367f, 0.063709047f, 0.063779971f, 0.063850896f, 0.063921821f, 0.063992745f,
        0.06406367f, 0.064134595f, 0.064205519f, 0.064276444f, 0.064347369f, 0.064418293f, 0.064512504f, 0.064606714f,
        0.064700924f, 0.064795135f, 0.064889345f, 0.064983555f, 0.065077766f, 0.065171976f, 0.065266186f, 0.065360397f,
        0.065423373f, 0.065486349f, 0.065549325f, 0.065612301f, 0.065675277f, 0.065738253f, 0.065801229f, 0.065864205f,
        0.065927181f, 0.065990157f, 0.066084808f, 0.066179459f, 0.066274111f, 0.066368762f, 0.066463413f, 0.066558065f,
        0.066652716f, 0.066747367f, 0.066842019f, 0.06693667f, 0.06708427f, 0.06723187f, 0.06737947f, 0.06752707f,
        0.06767467f, 0.06782227f, 0.06796987f, 0.06811747f, 0.06826507f, 0.06841267f, 0.068699645f, 0.068986619f,
        0.069273594f, 0.069560569f, 0.069847543f, 0.070134518f, 0.070421493f, 0.070708467f, 0.070995442f, 0.071282417f,
        0.071724769f, 0.072167121f, 0.072609473f, 0.073051825f, 0.073494177f, 0.073936529f, 0.074378881f, 0.074821233f,
        0.075263585f, 0.075705937f, 0.076857519f, 0.078009102f, 0.079160684f, 0.080312267f, 0.081463849f, 0.082615432f,
        0.083767014f, 0.084918597f, 0.086070179f, 0.087221762f, 0.091030431f, 0.094839101f, 0.098647771f, 0.10245644f,
        0.10626511f, 0.11007378f, 0.113882449f, 0.117691119f, 0.121499789f, 0.125308458f, 0.133360902f, 0.141413345f,
        0.149465788f, 0.157518232f, 0.165570675f, 0.173623118f, 0.181675562f, 0.189728005f, 0.197780448f, 0.205832892f,
        0.215775898f, 0.225718903f, 0.235661909f, 0.245604915f, 0.255547921f, 0.265490927f, 0.275433933f, 0.285376938f,
        0.295319944f, 0.30526295f, 0.31305211f, 0.32084127f, 0.32863043f, 0.33641959f, 0.34420875f, 0.35199791f,
        0.35978707f, 0.36757623f, 0.37536539f, 0.38315455f, 0.387933383f, 0.392712215f, 0.397491048f, 0.40226988f,
        0.407048713f, 0.411827545f, 0.416606378f, 0.42138521f, 0.426164043f, 0.430942875f, 0.434763645f, 0.438584415f,
        0.442405185f, 0.446225955f, 0.450046725f, 0.453867495f, 0.457688265f, 0.461509035f, 0.465329805f, 0.469150575f,
        0.474024545f, 0.478898515f, 0.483772485f, 0.488646455f, 0.493520425f, 0.498394395f, 0.503268365f, 0.508142335f,
        0.513016305f, 0.517890275f, 0.522894058f, 0.52789784f, 0.532901623f, 0.537905405f, 0.542909188f, 0.54791297f,
        0.552916753f, 0.557920535f, 0.562924318f, 0.5679281f, 0.571823276f, 0.575718452f, 0.579613628f, 0.583508803f,
        0.587403979f, 0.591299155f, 0.595194331f, 0.599089507f, 0.602984683f, 0.606879858f, 0.608996518f, 0.611113177f,
        0.613229836f, 0.615346495f, 0.617463154f, 0.619579813f, 0.621696473f, 0.623813132f, 0.625929791f, 0.62804645f,
        0.628944887f, 0.629843323f, 0.63074176f, 0.631640197f, 0.632538633f, 0.63343707f, 0.634335507f, 0.635233943f,
        0.63613238f, 0.637030817f, 0.637326863f, 0.637622908f, 0.637918954f, 0.638215f, 0.638511046f, 0.638807092f,
        0.639103138f, 0.639399183f, 0.639695229f, 0.639991275f, 0.640190068f, 0.640388862f, 0.640587655f, 0.640786448f,
        0.640985242f, 0.641184035f, 0.641382828f, 0.641581622f, 0.641780415f, 0.641979208f, 0.642326169f, 0.64267313f,
        0.643020091f, 0.643367052f, 0.643714013f, 0.644060973f, 0.644407934f, 0.644754895f, 0.645101856f, 0.645448817f,
        0.645727897f, 0.646006977f, 0.646286057f, 0.646565137f, 0.646844217f, 0.647123297f, 0.647402377f, 0.647681457f,
        0.647960537f, 0.648239617f, 0.648517497f, 0.648795377f, 0.649073257f, 0.649351137f, 0.649629017f, 0.649906897f,
        0.650184777f, 0.650462657f, 0.650740537f, 0.651018417f, 0.651223242f, 0.651428067f, 0.651632892f, 0.651837717f,
        0.652042542f, 0.652247367f, 0.652452192f, 0.652657017f, 0.652861842f, 0.653066667f, 0.65349558f, 0.653924493f,
        0.654353407f, 0.65478232f, 0.655211233f, 0.655640147f, 0.65606906f, 0.656497973f, 0.656926887f, 0.6573558f,
        0.658023494f, 0.658691188f, 0.659358883f, 0.660026577f, 0.660694271f, 0.661361965f, 0.662029659f, 0.662697353f,
        0.663365048f, 0.664032742f, 0.66489469f, 0.665756638f, 0.666618587f, 0.667480535f, 0.668342483f, 0.669204432f,
        0.67006638f, 0.670928328f, 0.671790277f, 0.672652225f, 0.67335715f, 0.674062075f, 0.674767f, 0.675471925f,
        0.67617685f, 0.676881775f, 0.6775867f, 0.678291625f, 0.67899655f, 0.679701475f, 0.680107636f, 0.680513797f,
        0.680919958f, 0.681326118f, 0.681732279f, 0.68213844f, 0.682544601f, 0.682950762f, 0.683356923f, 0.683763083f,
        0.684215498f, 0.684667913f, 0.685120328f, 0.685572743f, 0.686025158f, 0.686477573f, 0.686929988f, 0.687382403f,
        0.687834818f, 0.688287233f, 0.61945851f, 0.550629787f, 0.481801063f, 0.41297234f, 0.344143617f, 0.275314893f,
        0.20648617f, 0.137657447f, 0.068828723f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth21({0.066244793f, 0.067483929f, 0.068723065f, 0.069962201f, 0.071201337f,
        0.072440473f, 0.073679608f, 0.074918744f, 0.07615788f, 0.077397016f, 0.078636152f, 0.080931439f, 0.083226725f,
        0.085522012f, 0.087817299f, 0.090112586f, 0.092407873f, 0.09470316f, 0.096998446f, 0.099293733f, 0.10158902f,
        0.105984264f, 0.110379508f, 0.114774752f, 0.119169995f, 0.123565239f, 0.127960483f, 0.132355727f, 0.136750971f,
        0.141146215f, 0.145541458f, 0.150938738f, 0.156336017f, 0.161733296f, 0.167130575f, 0.172527854f, 0.177925133f,
        0.183322413f, 0.188719692f, 0.194116971f, 0.19951425f, 0.204002765f, 0.20849128f, 0.212979795f, 0.21746831f,
        0.221956825f, 0.22644534f, 0.230933855f, 0.23542237f, 0.239910885f, 0.2443994f, 0.248209173f, 0.252018945f,
        0.255828718f, 0.25963849f, 0.263448263f, 0.267258035f, 0.271067808f, 0.27487758f, 0.278687353f, 0.282497125f,
        0.28518319f, 0.287869255f, 0.29055532f, 0.293241385f, 0.29592745f, 0.298613515f, 0.30129958f, 0.303985645f,
        0.30667171f, 0.309357775f, 0.309180931f, 0.309004087f, 0.308827243f, 0.308650398f, 0.308473554f, 0.30829671f,
        0.308119866f, 0.307943022f, 0.307766178f, 0.307589333f, 0.304611171f, 0.301633008f, 0.298654846f, 0.295676683f,
        0.292698521f, 0.289720358f, 0.286742196f, 0.283764033f, 0.280785871f, 0.277807708f, 0.273113663f, 0.268419618f,
        0.263725573f, 0.259031528f, 0.254337483f, 0.249643438f, 0.244949393f, 0.240255348f, 0.235561303f, 0.230867258f,
        0.225534513f, 0.220201767f, 0.214869021f, 0.209536275f, 0.204203529f, 0.198870783f, 0.193538038f, 0.188205292f,
        0.182872546f, 0.1775398f, 0.172757168f, 0.167974537f, 0.163191905f, 0.158409273f, 0.153626642f, 0.14884401f,
        0.144061378f, 0.139278747f, 0.134496115f, 0.129713483f, 0.126170058f, 0.122626632f, 0.119083206f, 0.11553978f,
        0.111996354f, 0.108452928f, 0.104909503f, 0.101366077f, 0.097822651f, 0.094279225f, 0.091805435f, 0.089331646f,
        0.086857856f, 0.084384066f, 0.081910277f, 0.079436487f, 0.076962697f, 0.074488908f, 0.072015118f, 0.069541328f,
        0.067986135f, 0.066430942f, 0.064875749f, 0.063320556f, 0.061765363f, 0.06021017f, 0.058654977f, 0.057099784f,
        0.055544591f, 0.053989398f, 0.053172617f, 0.052355836f, 0.051539055f, 0.050722274f, 0.049905493f, 0.049088711f,
        0.04827193f, 0.047455149f, 0.046638368f, 0.045821587f, 0.04540613f, 0.044990673f, 0.044575217f, 0.04415976f,
        0.043744303f, 0.043328847f, 0.04291339f, 0.042497933f, 0.042082477f, 0.04166702f, 0.041444174f, 0.041221328f,
        0.040998483f, 0.040775637f, 0.040552791f, 0.040329945f, 0.040107099f, 0.039884253f, 0.039661408f, 0.039438562f,
        0.03932571f, 0.039212859f, 0.039100007f, 0.038987156f, 0.038874304f, 0.038761453f, 0.038648601f, 0.03853575f,
        0.038422898f, 0.038310047f, 0.038254032f, 0.038198017f, 0.038142002f, 0.038085987f, 0.038029973f, 0.037973958f,
        0.037917943f, 0.037861928f, 0.037805913f, 0.037749898f, 0.037750194f, 0.03775049f, 0.037750786f, 0.037751082f,
        0.037751378f, 0.037751674f, 0.03775197f, 0.037752266f, 0.037752562f, 0.037752858f, 0.037774517f, 0.037796176f,
        0.037817835f, 0.037839494f, 0.037861153f, 0.037882811f, 0.03790447f, 0.037926129f, 0.037947788f, 0.037969447f,
        0.038024552f, 0.038079657f, 0.038134763f, 0.038189868f, 0.038244973f, 0.038300079f, 0.038355184f, 0.038410289f,
        0.038465395f, 0.0385205f, 0.038593144f, 0.038665789f, 0.038738433f, 0.038811077f, 0.038883722f, 0.038956366f,
        0.03902901f, 0.039101655f, 0.039174299f, 0.039246943f, 0.039308816f, 0.039370688f, 0.039432561f, 0.039494433f,
        0.039556306f, 0.039618178f, 0.039680051f, 0.039741923f, 0.039803796f, 0.039865668f, 0.039958673f, 0.040051677f,
        0.040144682f, 0.040237686f, 0.040330691f, 0.040423695f, 0.0405167f, 0.040609704f, 0.040702709f, 0.040795713f,
        0.040948515f, 0.041101316f, 0.041254117f, 0.041406919f, 0.04155972f, 0.041712521f, 0.041865323f, 0.042018124f,
        0.042170925f, 0.042323727f, 0.042509255f, 0.042694784f, 0.042880312f, 0.043065841f, 0.043251369f, 0.043436898f,
        0.043622426f, 0.043807955f, 0.043993483f, 0.044179012f, 0.04430846f, 0.044437908f, 0.044567356f, 0.044696804f,
        0.044826252f, 0.0449557f, 0.045085148f, 0.045214596f, 0.045344044f, 0.045473492f, 0.045510572f, 0.045547653f,
        0.045584734f, 0.045621814f, 0.045658895f, 0.045695976f, 0.045733056f, 0.045770137f, 0.045807218f, 0.045844298f,
        0.045901646f, 0.045958994f, 0.046016342f, 0.04607369f, 0.046131038f, 0.046188385f, 0.046245733f, 0.046303081f,
        0.046360429f, 0.046417777f, 0.046612773f, 0.04680777f, 0.047002766f, 0.047197763f, 0.047392759f, 0.047587756f,
        0.047782752f, 0.047977749f, 0.048172745f, 0.048367742f, 0.048748155f, 0.049128568f, 0.049508982f, 0.049889395f,
        0.050269808f, 0.050650222f, 0.051030635f, 0.051411048f, 0.051791462f, 0.052171875f, 0.05268579f, 0.053199704f,
        0.053713619f, 0.054227534f, 0.054741448f, 0.055255363f, 0.055769278f, 0.056283192f, 0.056797107f, 0.057311022f,
        0.058078237f, 0.058845453f, 0.059612669f, 0.060379884f, 0.0611471f, 0.061914316f, 0.062681531f, 0.063448747f,
        0.064215963f, 0.064983178f, 0.058484861f, 0.051986543f, 0.045488225f, 0.038989907f, 0.032491589f, 0.025993271f,
        0.019494954f, 0.012996636f, 0.006498318f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth22({0.051950585f, 0.052061574f, 0.052172564f, 0.052283553f, 0.052394542f,
        0.052505532f, 0.052616521f, 0.05272751f, 0.0528385f, 0.052949489f, 0.053060478f, 0.053174412f, 0.053288345f,
        0.053402279f, 0.053516212f, 0.053630146f, 0.053744079f, 0.053858013f, 0.053971946f, 0.05408588f, 0.054199813f,
        0.054324554f, 0.054449295f, 0.054574036f, 0.054698777f, 0.054823518f, 0.054948259f, 0.055073f, 0.055197741f,
        0.055322482f, 0.055447223f, 0.055579124f, 0.055711025f, 0.055842925f, 0.055974826f, 0.056106727f, 0.056238627f,
        0.056370528f, 0.056502429f, 0.056634329f, 0.05676623f, 0.056945853f, 0.057125476f, 0.057305099f, 0.057484722f,
        0.057664345f, 0.057843968f, 0.058023591f, 0.058203214f, 0.058382837f, 0.05856246f, 0.058843538f, 0.059124616f,
        0.059405694f, 0.059686771f, 0.059967849f, 0.060248927f, 0.060530005f, 0.060811083f, 0.061092161f, 0.061373238f,
        0.061811859f, 0.062250479f, 0.062689099f, 0.06312772f, 0.06356634f, 0.06400496f, 0.064443581f, 0.064882201f,
        0.065320821f, 0.065759442f, 0.066666911f, 0.06757438f, 0.068481849f, 0.069389318f, 0.070296788f, 0.071204257f,
        0.072111726f, 0.073019195f, 0.073926664f, 0.074834133f, 0.076619936f, 0.078405738f, 0.08019154f, 0.081977342f,
        0.083763144f, 0.085548946f, 0.087334749f, 0.089120551f, 0.090906353f, 0.092692155f, 0.095911162f, 0.099130169f,
        0.102349176f, 0.105568183f, 0.10878719f, 0.112006197f, 0.115225204f, 0.118444211f, 0.121663218f, 0.124882225f,
        0.130183128f, 0.135484032f, 0.140784935f, 0.146085838f, 0.151386742f, 0.156687645f, 0.161988548f, 0.167289452f,
        0.172590355f, 0.177891258f, 0.184680648f, 0.191470038f, 0.198259428f, 0.205048818f, 0.211838208f, 0.218627598f,
        0.225416988f, 0.232206378f, 0.238995768f, 0.245785158f, 0.25193172f, 0.258078282f, 0.264224843f, 0.270371405f,
        0.276517967f, 0.282664528f, 0.28881109f, 0.294957652f, 0.301104213f, 0.307250775f, 0.310241468f, 0.31323216f,
        0.316222853f, 0.319213545f, 0.322204238f, 0.32519493f, 0.328185623f, 0.331176315f, 0.334167008f, 0.3371577f,
        0.336795513f, 0.336433325f, 0.336071138f, 0.33570895f, 0.335346763f, 0.334984575f, 0.334622388f, 0.3342602f,
        0.333898013f, 0.333535825f, 0.331835113f, 0.3301344f, 0.328433688f, 0.326732975f, 0.325032263f, 0.32333155f,
        0.321630838f, 0.319930125f, 0.318229413f, 0.3165287f, 0.314175215f, 0.31182173f, 0.309468245f, 0.30711476f,
        0.304761275f, 0.30240779f, 0.300054305f, 0.29770082f, 0.295347335f, 0.29299385f, 0.289880905f, 0.28676796f,
        0.283655015f, 0.28054207f, 0.277429125f, 0.27431618f, 0.271203235f, 0.26809029f, 0.264977345f, 0.2618644f,
        0.258677209f, 0.255490018f, 0.252302828f, 0.249115637f, 0.245928446f, 0.242741255f, 0.239554064f, 0.236366873f,
        0.233179683f, 0.229992492f, 0.22675861f, 0.223524728f, 0.220290847f, 0.217056965f, 0.213823083f, 0.210589202f,
        0.20735532f, 0.204121438f, 0.200887557f, 0.197653675f, 0.194392049f, 0.191130423f, 0.187868798f, 0.184607172f,
        0.181345546f, 0.17808392f, 0.174822294f, 0.171560668f, 0.168299043f, 0.165037417f, 0.162034949f, 0.159032482f,
        0.156030014f, 0.153027547f, 0.150025079f, 0.147022612f, 0.144020144f, 0.141017677f, 0.138015209f, 0.135012742f,
        0.133001506f, 0.13099027f, 0.128979034f, 0.126967798f, 0.124956563f, 0.122945327f, 0.120934091f, 0.118922855f,
        0.116911619f, 0.114900383f, 0.113807521f, 0.112714659f, 0.111621797f, 0.110528935f, 0.109436073f, 0.108343211f,
        0.107250349f, 0.106157487f, 0.105064625f, 0.103971763f, 0.103366052f, 0.10276034f, 0.102154628f, 0.101548916f,
        0.100943204f, 0.100337492f, 0.099731781f, 0.099126069f, 0.098520357f, 0.097914645f, 0.097562141f, 0.097209637f,
        0.096857133f, 0.096504629f, 0.096152125f, 0.095799621f, 0.095447117f, 0.095094613f, 0.094742109f, 0.094389605f,
        0.094185508f, 0.093981411f, 0.093777314f, 0.093573216f, 0.093369119f, 0.093165022f, 0.092960925f, 0.092756828f,
        0.092552731f, 0.092348633f, 0.092390764f, 0.092432894f, 0.092475024f, 0.092517155f, 0.092559285f, 0.092601415f,
        0.092643546f, 0.092685676f, 0.092727806f, 0.092769937f, 0.09314561f, 0.093521283f, 0.093896957f, 0.09427263f,
        0.094648303f, 0.095023977f, 0.09539965f, 0.095775323f, 0.096150997f, 0.09652667f, 0.097114159f, 0.097701649f,
        0.098289138f, 0.098876627f, 0.099464117f, 0.100051606f, 0.100639095f, 0.101226585f, 0.101814074f, 0.102401563f,
        0.103003696f, 0.103605828f, 0.10420796f, 0.104810093f, 0.105412225f, 0.106014357f, 0.10661649f, 0.107218622f,
        0.107820754f, 0.108422887f, 0.108925907f, 0.109428927f, 0.109931948f, 0.110434968f, 0.110937988f, 0.111441009f,
        0.111944029f, 0.112447049f, 0.11295007f, 0.11345309f, 0.11364037f, 0.11382765f, 0.11401493f, 0.11420221f,
        0.11438949f, 0.11457677f, 0.11476405f, 0.11495133f, 0.11513861f, 0.11532589f, 0.115185591f, 0.115045292f,
        0.114904993f, 0.114764693f, 0.114624394f, 0.114484095f, 0.114343796f, 0.114203497f, 0.114063198f, 0.113922898f,
        0.113957618f, 0.113992338f, 0.114027058f, 0.114061778f, 0.114096498f, 0.114131218f, 0.114165938f, 0.114200658f,
        0.114235378f, 0.114270098f, 0.102843089f, 0.091416079f, 0.079989069f, 0.068562059f, 0.057135049f, 0.045708039f,
        0.03428103f, 0.02285402f, 0.01142701f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth23({0.049923225f, 0.049808127f, 0.049693028f, 0.04957793f, 0.049462832f,
        0.049347733f, 0.049232635f, 0.049117537f, 0.049002438f, 0.04888734f, 0.048772242f, 0.048654367f, 0.048536492f,
        0.048418618f, 0.048300743f, 0.048182868f, 0.048064994f, 0.047947119f, 0.047829244f, 0.04771137f, 0.047593495f,
        0.047557731f, 0.047521968f, 0.047486204f, 0.04745044f, 0.047414677f, 0.047378913f, 0.047343149f, 0.047307386f,
        0.047271622f, 0.047235858f, 0.047228584f, 0.047221309f, 0.047214034f, 0.047206759f, 0.047199484f, 0.047192209f,
        0.047184935f, 0.04717766f, 0.047170385f, 0.04716311f, 0.047182243f, 0.047201376f, 0.047220509f, 0.047239641f,
        0.047258774f, 0.047277907f, 0.04729704f, 0.047316173f, 0.047335306f, 0.047354438f, 0.047360874f, 0.047367309f,
        0.047373745f, 0.04738018f, 0.047386616f, 0.047393051f, 0.047399487f, 0.047405922f, 0.047412358f, 0.047418793f,
        0.047368927f, 0.047319061f, 0.047269195f, 0.047219329f, 0.047169463f, 0.047119597f, 0.047069731f, 0.047019865f,
        0.046969999f, 0.046920133f, 0.046835367f, 0.046750601f, 0.046665835f, 0.046581069f, 0.046496303f, 0.046411537f,
        0.046326771f, 0.046242005f, 0.046157239f, 0.046072473f, 0.045982751f, 0.045893029f, 0.045803307f, 0.045713585f,
        0.045623863f, 0.04553414f, 0.045444418f, 0.045354696f, 0.045264974f, 0.045175252f, 0.045101658f, 0.045028063f,
        0.044954469f, 0.044880875f, 0.044807281f, 0.044733687f, 0.044660093f, 0.044586498f, 0.044512904f, 0.04443931f,
        0.044424251f, 0.044409192f, 0.044394133f, 0.044379074f, 0.044364015f, 0.044348956f, 0.044333897f, 0.044318838f,
        0.044303779f, 0.04428872f, 0.04432793f, 0.04436714f, 0.04440635f, 0.04444556f, 0.04448477f, 0.04452398f,
        0.04456319f, 0.0446024f, 0.04464161f, 0.04468082f, 0.044772309f, 0.044863799f, 0.044955288f, 0.045046777f,
        0.045138267f, 0.045229756f, 0.045321245f, 0.045412735f, 0.045504224f, 0.045595713f, 0.045714163f, 0.045832613f,
        0.045951063f, 0.046069513f, 0.046187963f, 0.046306413f, 0.046424863f, 0.046543313f, 0.046661763f, 0.046780213f,
        0.046866198f, 0.046952183f, 0.047038167f, 0.047124152f, 0.047210137f, 0.047296121f, 0.047382106f, 0.047468091f,
        0.047554075f, 0.04764006f, 0.047734641f, 0.047829222f, 0.047923803f, 0.048018383f, 0.048112964f, 0.048207545f,
        0.048302126f, 0.048396707f, 0.048491288f, 0.048585868f, 0.048763807f, 0.048941745f, 0.049119684f, 0.049297622f,
        0.049475561f, 0.049653499f, 0.049831438f, 0.050009376f, 0.050187315f, 0.050365253f, 0.050713413f, 0.051061572f,
        0.051409731f, 0.051757891f, 0.05210605f, 0.052454209f, 0.052802369f, 0.053150528f, 0.053498687f, 0.053846847f,
        0.054448226f, 0.055049604f, 0.055650983f, 0.056252362f, 0.056853741f, 0.05745512f, 0.058056499f, 0.058657877f,
        0.059259256f, 0.059860635f, 0.061086313f, 0.062311991f, 0.063537669f, 0.064763347f, 0.065989025f, 0.067214703f,
        0.068440381f, 0.069666059f, 0.070891737f, 0.072117415f, 0.07526156f, 0.078405705f, 0.08154985f, 0.084693995f,
        0.08783814f, 0.090982285f, 0.09412643f, 0.097270575f, 0.10041472f, 0.103558865f, 0.110955209f, 0.118351554f,
        0.125747898f, 0.133144242f, 0.140540587f, 0.147936931f, 0.155333275f, 0.16272962f, 0.170125964f, 0.177522308f,
        0.190977247f, 0.204432185f, 0.217887123f, 0.231342062f, 0.244797f, 0.258251938f, 0.271706877f, 0.285161815f,
        0.298616753f, 0.312071692f, 0.32754781f, 0.343023928f, 0.358500047f, 0.373976165f, 0.389452283f, 0.404928402f,
        0.42040452f, 0.435880638f, 0.451356757f, 0.466832875f, 0.478232743f, 0.489632612f, 0.50103248f, 0.512432348f,
        0.523832217f, 0.535232085f, 0.546631953f, 0.558031822f, 0.56943169f, 0.580831558f, 0.587191888f, 0.593552218f,
        0.599912548f, 0.606272878f, 0.612633208f, 0.618993538f, 0.625353868f, 0.631714198f, 0.638074528f, 0.644434858f,
        0.647475183f, 0.650515507f, 0.653555831f, 0.656596155f, 0.659636479f, 0.662676803f, 0.665717128f, 0.668757452f,
        0.671797776f, 0.6748381f, 0.676371843f, 0.677905586f, 0.679439329f, 0.680973071f, 0.682506814f, 0.684040557f,
        0.6855743f, 0.687108043f, 0.688641786f, 0.690175528f, 0.690982426f, 0.691789324f, 0.692596222f, 0.69340312f,
        0.694210018f, 0.695016916f, 0.695823814f, 0.696630712f, 0.69743761f, 0.698244508f, 0.699011744f, 0.69977898f,
        0.700546216f, 0.701313452f, 0.702080688f, 0.702847923f, 0.703615159f, 0.704382395f, 0.705149631f, 0.705916867f,
        0.706820624f, 0.707724382f, 0.708628139f, 0.709531897f, 0.710435654f, 0.711339412f, 0.712243169f, 0.713146927f,
        0.714050684f, 0.714954442f, 0.715829284f, 0.716704127f, 0.717578969f, 0.718453812f, 0.719328654f, 0.720203497f,
        0.721078339f, 0.721953182f, 0.722828024f, 0.723702867f, 0.724343013f, 0.724983158f, 0.725623304f, 0.72626345f,
        0.726903596f, 0.727543742f, 0.728183888f, 0.728824033f, 0.729464179f, 0.730104325f, 0.730465293f, 0.730826262f,
        0.73118723f, 0.731548198f, 0.731909167f, 0.732270135f, 0.732631103f, 0.732992072f, 0.73335304f, 0.733714008f,
        0.734183318f, 0.734652627f, 0.735121936f, 0.735591245f, 0.736060554f, 0.736529863f, 0.736999173f, 0.737468482f,
        0.737937791f, 0.7384071f, 0.66456639f, 0.59072568f, 0.51688497f, 0.44304426f, 0.36920355f, 0.29536284f,
        0.22152213f, 0.14768142f, 0.07384071f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth24({0.057984715f, 0.057628052f, 0.057271388f, 0.056914725f, 0.056558062f,
        0.056201398f, 0.055844735f, 0.055488072f, 0.055131408f, 0.054774745f, 0.054418082f, 0.054192412f, 0.053966742f,
        0.053741073f, 0.053515403f, 0.053289733f, 0.053064064f, 0.052838394f, 0.052612724f, 0.052387055f, 0.052161385f,
        0.052143699f, 0.052126012f, 0.052108326f, 0.052090639f, 0.052072953f, 0.052055266f, 0.05203758f, 0.052019893f,
        0.052002207f, 0.05198452f, 0.052048576f, 0.052112633f, 0.052176689f, 0.052240745f, 0.052304802f, 0.052368858f,
        0.052432914f, 0.052496971f, 0.052561027f, 0.052625083f, 0.052760931f, 0.052896779f, 0.053032627f, 0.053168475f,
        0.053304323f, 0.05344017f, 0.053576018f, 0.053711866f, 0.053847714f, 0.053983562f, 0.054192932f, 0.054402303f,
        0.054611674f, 0.054821044f, 0.055030415f, 0.055239786f, 0.055449156f, 0.055658527f, 0.055867898f, 0.056077268f,
        0.056411799f, 0.05674633f, 0.057080861f, 0.057415392f, 0.057749923f, 0.058084454f, 0.058418985f, 0.058753516f,
        0.059088047f, 0.059422578f, 0.060139201f, 0.060855824f, 0.061572446f, 0.062289069f, 0.063005692f, 0.063722314f,
        0.064438937f, 0.06515556f, 0.065872182f, 0.066588805f, 0.067997994f, 0.069407183f, 0.070816372f, 0.07222556f,
        0.073634749f, 0.075043938f, 0.076453127f, 0.077862316f, 0.079271505f, 0.080680693f, 0.083301016f, 0.085921339f,
        0.088541662f, 0.091161985f, 0.093782308f, 0.09640263f, 0.099022953f, 0.101643276f, 0.104263599f, 0.106883922f,
        0.111399713f, 0.115915504f, 0.120431295f, 0.124947086f, 0.129462878f, 0.133978669f, 0.13849446f, 0.143010251f,
        0.147526042f, 0.152041833f, 0.159344925f, 0.166648017f, 0.173951108f, 0.1812542f, 0.188557292f, 0.195860383f,
        0.203163475f, 0.210466567f, 0.217769658f, 0.22507275f, 0.236118538f, 0.247164327f, 0.258210115f, 0.269255903f,
        0.280301692f, 0.29134748f, 0.302393268f, 0.313439057f, 0.324484845f, 0.335530633f, 0.348216821f, 0.360903008f,
        0.373589196f, 0.386275383f, 0.398961571f, 0.411647758f, 0.424333946f, 0.437020133f, 0.449706321f, 0.462392508f,
        0.472026513f, 0.481660517f, 0.491294521f, 0.500928525f, 0.510562529f, 0.520196533f, 0.529830538f, 0.539464542f,
        0.549098546f, 0.55873255f, 0.564431976f, 0.570131402f, 0.575830828f, 0.581530253f, 0.587229679f, 0.592929105f,
        0.598628531f, 0.604327957f, 0.610027383f, 0.615726808f, 0.619126633f, 0.622526458f, 0.625926283f, 0.629326108f,
        0.632725933f, 0.636125758f, 0.639525583f, 0.642925408f, 0.646325233f, 0.649725058f, 0.651974476f, 0.654223893f,
        0.656473311f, 0.658722728f, 0.660972146f, 0.663221563f, 0.665470981f, 0.667720398f, 0.669969816f, 0.672219233f,
        0.674384135f, 0.676549037f, 0.678713938f, 0.68087884f, 0.683043742f, 0.685208643f, 0.687373545f, 0.689538447f,
        0.691703348f, 0.69386825f, 0.695476616f, 0.697084982f, 0.698693348f, 0.700301713f, 0.701910079f, 0.703518445f,
        0.705126811f, 0.706735177f, 0.708343543f, 0.709951908f, 0.711275701f, 0.712599493f, 0.713923286f, 0.715247078f,
        0.716570871f, 0.717894663f, 0.719218456f, 0.720542248f, 0.721866041f, 0.723189833f, 0.724014358f, 0.724838883f,
        0.725663408f, 0.726487933f, 0.727312458f, 0.728136983f, 0.728961508f, 0.729786033f, 0.730610558f, 0.731435083f,
        0.732195703f, 0.732956322f, 0.733716941f, 0.73447756f, 0.735238179f, 0.735998798f, 0.736759418f, 0.737520037f,
        0.738280656f, 0.739041275f, 0.739757163f, 0.74047305f, 0.741188938f, 0.741904825f, 0.742620713f, 0.7433366f,
        0.744052488f, 0.744768375f, 0.745484263f, 0.74620015f, 0.746760517f, 0.747320883f, 0.74788125f, 0.748441617f,
        0.749001983f, 0.74956235f, 0.750122717f, 0.750683083f, 0.75124345f, 0.751803817f, 0.752439286f, 0.753074755f,
        0.753710224f, 0.754345693f, 0.754981163f, 0.755616632f, 0.756252101f, 0.75688757f, 0.757523039f, 0.758158508f,
        0.758736744f, 0.75931498f, 0.759893216f, 0.760471452f, 0.761049688f, 0.761627923f, 0.762206159f, 0.762784395f,
        0.763362631f, 0.763940867f, 0.764415688f, 0.76489051f, 0.765365332f, 0.765840153f, 0.766314975f, 0.766789797f,
        0.767264618f, 0.76773944f, 0.768214262f, 0.768689083f, 0.768918358f, 0.769147632f, 0.769376906f, 0.76960618f,
        0.769835454f, 0.770064728f, 0.770294003f, 0.770523277f, 0.770752551f, 0.770981825f, 0.771434485f, 0.771887145f,
        0.772339805f, 0.772792465f, 0.773245125f, 0.773697785f, 0.774150445f, 0.774603105f, 0.775055765f, 0.775508425f,
        0.776197791f, 0.776887157f, 0.777576523f, 0.778265888f, 0.778955254f, 0.77964462f, 0.780333986f, 0.781023352f,
        0.781712718f, 0.782402083f, 0.783179574f, 0.783957065f, 0.784734556f, 0.785512047f, 0.786289538f, 0.787067028f,
        0.787844519f, 0.78862201f, 0.789399501f, 0.790176992f, 0.79077827f, 0.791379548f, 0.791980827f, 0.792582105f,
        0.793183383f, 0.793784662f, 0.79438594f, 0.794987218f, 0.795588497f, 0.796189775f, 0.796500418f, 0.796811062f,
        0.797121705f, 0.797432348f, 0.797742992f, 0.798053635f, 0.798364278f, 0.798674922f, 0.798985565f, 0.799296208f,
        0.799732681f, 0.800169153f, 0.800605626f, 0.801042098f, 0.801478571f, 0.801915043f, 0.802351516f, 0.802787988f,
        0.803224461f, 0.803660933f, 0.72329484f, 0.642928747f, 0.562562653f, 0.48219656f, 0.401830467f, 0.321464373f,
        0.24109828f, 0.160732187f, 0.080366093f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth25({0.144547292f, 0.149603995f, 0.154660698f, 0.159717402f, 0.164774105f,
        0.169830808f, 0.174887512f, 0.179944215f, 0.185000918f, 0.190057622f, 0.195114325f, 0.20386221f, 0.212610095f,
        0.22135798f, 0.230105865f, 0.23885375f, 0.247601635f, 0.25634952f, 0.265097405f, 0.27384529f, 0.282593175f,
        0.288910735f, 0.295228295f, 0.301545855f, 0.307863415f, 0.314180975f, 0.320498535f, 0.326816095f, 0.333133655f,
        0.339451215f, 0.345768775f, 0.347374304f, 0.348979833f, 0.350585363f, 0.352190892f, 0.353796421f, 0.35540195f,
        0.357007479f, 0.358613008f, 0.360218538f, 0.361824067f, 0.361073996f, 0.360323925f, 0.359573854f, 0.358823783f,
        0.358073713f, 0.357323642f, 0.356573571f, 0.3558235f, 0.355073429f, 0.354323358f, 0.352252065f, 0.350180772f,
        0.348109478f, 0.346038185f, 0.343966892f, 0.341895598f, 0.339824305f, 0.337753012f, 0.335681718f, 0.333610425f,
        0.330820404f, 0.328030383f, 0.325240363f, 0.322450342f, 0.319660321f, 0.3168703f, 0.314080279f, 0.311290258f,
        0.308500238f, 0.305710217f, 0.302761862f, 0.299813507f, 0.296865152f, 0.293916797f, 0.290968442f, 0.288020087f,
        0.285071732f, 0.282123377f, 0.279175022f, 0.276226667f, 0.273359616f, 0.270492565f, 0.267625514f, 0.264758463f,
        0.261891413f, 0.259024362f, 0.256157311f, 0.25329026f, 0.250423209f, 0.247556158f, 0.244605828f, 0.241655498f,
        0.238705168f, 0.235754838f, 0.232804508f, 0.229854178f, 0.226903848f, 0.223953518f, 0.221003188f, 0.218052858f,
        0.215235745f, 0.212418632f, 0.209601518f, 0.206784405f, 0.203967292f, 0.201150178f, 0.198333065f, 0.195515952f,
        0.192698838f, 0.189881725f, 0.187679393f, 0.185477062f, 0.18327473f, 0.181072398f, 0.178870067f, 0.176667735f,
        0.174465403f, 0.172263072f, 0.17006074f, 0.167858408f, 0.165968076f, 0.164077743f, 0.162187411f, 0.160297078f,
        0.158406746f, 0.156516413f, 0.154626081f, 0.152735748f, 0.150845416f, 0.148955083f, 0.146756609f, 0.144558135f,
        0.142359661f, 0.140161187f, 0.137962713f, 0.135764238f, 0.133565764f, 0.13136729f, 0.129168816f, 0.126970342f,
        0.12499661f, 0.123022878f, 0.121049147f, 0.119075415f, 0.117101683f, 0.115127952f, 0.11315422f, 0.111180488f,
        0.109206757f, 0.107233025f, 0.106471288f, 0.10570955f, 0.104947813f, 0.104186075f, 0.103424338f, 0.1026626f,
        0.101900863f, 0.101139125f, 0.100377388f, 0.09961565f, 0.099842944f, 0.100070238f, 0.100297533f, 0.100524827f,
        0.100752121f, 0.100979415f, 0.101206709f, 0.101434003f, 0.101661298f, 0.101888592f, 0.102055773f, 0.102222955f,
        0.102390137f, 0.102557318f, 0.1027245f, 0.102891682f, 0.103058863f, 0.103226045f, 0.103393227f, 0.103560408f,
        0.104111032f, 0.104661655f, 0.105212278f, 0.105762902f, 0.106313525f, 0.106864148f, 0.107414772f, 0.107965395f,
        0.108516018f, 0.109066642f, 0.111840288f, 0.114613935f, 0.117387582f, 0.120161228f, 0.122934875f, 0.125708522f,
        0.128482168f, 0.131255815f, 0.134029462f, 0.136803108f, 0.14308545f, 0.149367792f, 0.155650133f, 0.161932475f,
        0.168214817f, 0.174497158f, 0.1807795f, 0.187061842f, 0.193344183f, 0.199626525f, 0.208676675f, 0.217726825f,
        0.226776975f, 0.235827125f, 0.244877275f, 0.253927425f, 0.262977575f, 0.272027725f, 0.281077875f, 0.290128025f,
        0.301121203f, 0.312114382f, 0.32310756f, 0.334100738f, 0.345093917f, 0.356087095f, 0.367080273f, 0.378073452f,
        0.38906663f, 0.400059808f, 0.411633658f, 0.423207507f, 0.434781356f, 0.446355205f, 0.457929054f, 0.469502903f,
        0.481076753f, 0.492650602f, 0.504224451f, 0.5157983f, 0.525704643f, 0.535610985f, 0.545517328f, 0.55542367f,
        0.565330013f, 0.575236355f, 0.585142698f, 0.59504904f, 0.604955383f, 0.614861725f, 0.622030803f, 0.629199882f,
        0.63636896f, 0.643538038f, 0.650707117f, 0.657876195f, 0.665045273f, 0.672214352f, 0.67938343f, 0.686552508f,
        0.69107417f, 0.695595832f, 0.700117493f, 0.704639155f, 0.709160817f, 0.713682478f, 0.71820414f, 0.722725802f,
        0.727247463f, 0.731769125f, 0.734567121f, 0.737365117f, 0.740163113f, 0.742961108f, 0.745759104f, 0.7485571f,
        0.751355096f, 0.754153092f, 0.756951088f, 0.759749083f, 0.761206923f, 0.762664762f, 0.764122601f, 0.76558044f,
        0.767038279f, 0.768496118f, 0.769953958f, 0.771411797f, 0.772869636f, 0.774327475f, 0.775209151f, 0.776090827f,
        0.776972503f, 0.777854178f, 0.778735854f, 0.77961753f, 0.780499206f, 0.781380882f, 0.782262558f, 0.783144233f,
        0.784085772f, 0.78502731f, 0.785968848f, 0.786910387f, 0.787851925f, 0.788793463f, 0.789735002f, 0.79067654f,
        0.791618078f, 0.792559617f, 0.793640925f, 0.794722233f, 0.795803542f, 0.79688485f, 0.797966158f, 0.799047467f,
        0.800128775f, 0.801210083f, 0.802291392f, 0.8033727f, 0.804190628f, 0.805008555f, 0.805826483f, 0.80664441f,
        0.807462338f, 0.808280265f, 0.809098193f, 0.80991612f, 0.810734048f, 0.811551975f, 0.812114972f, 0.812677968f,
        0.813240965f, 0.813803962f, 0.814366958f, 0.814929955f, 0.815492952f, 0.816055948f, 0.816618945f, 0.817181942f,
        0.818004698f, 0.818827453f, 0.819650209f, 0.820472965f, 0.821295721f, 0.822118477f, 0.822941233f, 0.823763988f,
        0.824586744f, 0.8254095f, 0.74286855f, 0.6603276f, 0.57778665f, 0.4952457f, 0.41270475f, 0.3301638f,
        0.24762285f, 0.1650819f, 0.08254095f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth26({0.107731417f, 0.111077022f, 0.114422627f, 0.117768232f, 0.121113837f,
        0.124459442f, 0.127805047f, 0.131150652f, 0.134496257f, 0.137841862f, 0.141187467f, 0.146315534f, 0.151443602f,
        0.156571669f, 0.161699737f, 0.166827804f, 0.171955872f, 0.177083939f, 0.182212007f, 0.187340074f, 0.192468142f,
        0.196862589f, 0.201257037f, 0.205651484f, 0.210045932f, 0.214440379f, 0.218834827f, 0.223229274f, 0.227623722f,
        0.232018169f, 0.236412617f, 0.238856745f, 0.241300873f, 0.243745002f, 0.24618913f, 0.248633258f, 0.251077387f,
        0.253521515f, 0.255965643f, 0.258409772f, 0.2608539f, 0.263318711f, 0.265783522f, 0.268248333f, 0.270713143f,
        0.273177954f, 0.275642765f, 0.278107576f, 0.280572387f, 0.283037198f, 0.285502008f, 0.2886915f, 0.291880992f,
        0.295070483f, 0.298259975f, 0.301449467f, 0.304638958f, 0.30782845f, 0.311017942f, 0.314207433f, 0.317396925f,
        0.320970598f, 0.32454427f, 0.328117943f, 0.331691615f, 0.335265288f, 0.33883896f, 0.342412633f, 0.345986305f,
        0.349559978f, 0.35313365f, 0.356843817f, 0.360553983f, 0.36426415f, 0.367974317f, 0.371684483f, 0.37539465f,
        0.379104817f, 0.382814983f, 0.38652515f, 0.390235317f, 0.393808434f, 0.397381552f, 0.400954669f, 0.404527787f,
        0.408100904f, 0.411674022f, 0.415247139f, 0.418820257f, 0.422393374f, 0.425966492f, 0.427930518f, 0.429894545f,
        0.431858572f, 0.433822598f, 0.435786625f, 0.437750652f, 0.439714678f, 0.441678705f, 0.443642732f, 0.445606758f,
        0.445469133f, 0.445331508f, 0.445193883f, 0.445056258f, 0.444918633f, 0.444781008f, 0.444643383f, 0.444505758f,
        0.444368133f, 0.444230508f, 0.442128823f, 0.440027137f, 0.437925451f, 0.435823765f, 0.433722079f, 0.431620393f,
        0.429518708f, 0.427417022f, 0.425315336f, 0.42321365f, 0.419441042f, 0.415668433f, 0.411895825f, 0.408123217f,
        0.404350608f, 0.400578f, 0.396805392f, 0.393032783f, 0.389260175f, 0.385487567f, 0.380610718f, 0.375733868f,
        0.370857019f, 0.36598017f, 0.361103321f, 0.356226472f, 0.351349623f, 0.346472773f, 0.341595924f, 0.336719075f,
        0.331320268f, 0.325921462f, 0.320522655f, 0.315123848f, 0.309725042f, 0.304326235f, 0.298927428f, 0.293528622f,
        0.288129815f, 0.282731008f, 0.277586156f, 0.272441303f, 0.267296451f, 0.262151598f, 0.257006746f, 0.251861893f,
        0.246717041f, 0.241572188f, 0.236427336f, 0.231282483f, 0.226659771f, 0.222037058f, 0.217414346f, 0.212791633f,
        0.208168921f, 0.203546208f, 0.198923496f, 0.194300783f, 0.189678071f, 0.185055358f, 0.181103803f, 0.177152248f,
        0.173200693f, 0.169249138f, 0.165297583f, 0.161346028f, 0.157394473f, 0.153442918f, 0.149491363f, 0.145539808f,
        0.142792583f, 0.140045357f, 0.137298131f, 0.134550905f, 0.131803679f, 0.129056453f, 0.126309228f, 0.123562002f,
        0.120814776f, 0.11806755f, 0.116313407f, 0.114559264f, 0.112805121f, 0.111050978f, 0.109296835f, 0.107542692f,
        0.105788549f, 0.104034406f, 0.102280263f, 0.10052612f, 0.099431749f, 0.098337378f, 0.097243008f, 0.096148637f,
        0.095054266f, 0.093959895f, 0.092865524f, 0.091771153f, 0.090676783f, 0.089582412f, 0.088780049f, 0.087977686f,
        0.087175323f, 0.08637296f, 0.085570597f, 0.084768234f, 0.083965871f, 0.083163508f, 0.082361145f, 0.081558782f,
        0.08104252f, 0.080526259f, 0.080009997f, 0.079493736f, 0.078977474f, 0.078461213f, 0.077944951f, 0.07742869f,
        0.076912428f, 0.076396167f, 0.076162454f, 0.075928742f, 0.075695029f, 0.075461317f, 0.075227604f, 0.074993892f,
        0.074760179f, 0.074526467f, 0.074292754f, 0.074059042f, 0.07395845f, 0.073857859f, 0.073757268f, 0.073656676f,
        0.073556085f, 0.073455494f, 0.073354902f, 0.073254311f, 0.07315372f, 0.073053128f, 0.073042195f, 0.073031262f,
        0.073020328f, 0.073009395f, 0.072998462f, 0.072987528f, 0.072976595f, 0.072965662f, 0.072954728f, 0.072943795f,
        0.073030288f, 0.07311678f, 0.073203273f, 0.073289766f, 0.073376258f, 0.073462751f, 0.073549244f, 0.073635736f,
        0.073722229f, 0.073808722f, 0.073986753f, 0.074164784f, 0.074342815f, 0.074520846f, 0.074698878f, 0.074876909f,
        0.07505494f, 0.075232971f, 0.075411002f, 0.075589033f, 0.075705172f, 0.075821311f, 0.07593745f, 0.076053589f,
        0.076169728f, 0.076285867f, 0.076402006f, 0.076518145f, 0.076634284f, 0.076750423f, 0.07672324f, 0.076696057f,
        0.076668874f, 0.076641691f, 0.076614508f, 0.076587324f, 0.076560141f, 0.076532958f, 0.076505775f, 0.076478592f,
        0.076329891f, 0.07618119f, 0.076032489f, 0.075883788f, 0.075735087f, 0.075586386f, 0.075437685f, 0.075288984f,
        0.075140283f, 0.074991582f, 0.074769379f, 0.074547177f, 0.074324975f, 0.074102772f, 0.07388057f, 0.073658368f,
        0.073436165f, 0.073213963f, 0.072991761f, 0.072769558f, 0.072692466f, 0.072615373f, 0.07253828f, 0.072461188f,
        0.072384095f, 0.072307002f, 0.07222991f, 0.072152817f, 0.072075724f, 0.071998632f, 0.072173265f, 0.072347897f,
        0.07252253f, 0.072697163f, 0.072871796f, 0.073046429f, 0.073221062f, 0.073395694f, 0.073570327f, 0.07374496f,
        0.074305955f, 0.074866949f, 0.075427944f, 0.075988939f, 0.076549933f, 0.077110928f, 0.077671923f, 0.078232917f,
        0.078793912f, 0.079354907f, 0.071419416f, 0.063483925f, 0.055548435f, 0.047612944f, 0.039677453f, 0.031741963f,
        0.023806472f, 0.015870981f, 0.007935491f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth31({0.189362217f, 0.195889904f, 0.202417592f, 0.208945279f, 0.215472967f,
        0.222000654f, 0.228528342f, 0.235056029f, 0.241583717f, 0.248111404f, 0.254639092f, 0.271435503f, 0.288231913f,
        0.305028324f, 0.321824735f, 0.338621146f, 0.355417557f, 0.372213968f, 0.389010378f, 0.405806789f, 0.4226032f,
        0.446363808f, 0.470124417f, 0.493885025f, 0.517645633f, 0.541406242f, 0.56516685f, 0.588927458f, 0.612688067f,
        0.636448675f, 0.660209283f, 0.675285888f, 0.690362493f, 0.705439098f, 0.720515703f, 0.735592308f, 0.750668913f,
        0.765745518f, 0.780822123f, 0.795898728f, 0.810975333f, 0.816089584f, 0.821203835f, 0.826318086f, 0.831432337f,
        0.836546588f, 0.841660838f, 0.846775089f, 0.85188934f, 0.857003591f, 0.862117842f, 0.863563629f, 0.865009417f,
        0.866455204f, 0.867900992f, 0.869346779f, 0.870792567f, 0.872238354f, 0.873684142f, 0.875129929f, 0.876575717f,
        0.877335507f, 0.878095297f, 0.878855087f, 0.879614877f, 0.880374667f, 0.881134457f, 0.881894247f, 0.882654037f,
        0.883413827f, 0.884173617f, 0.884859832f, 0.885546047f, 0.886232262f, 0.886918477f, 0.887604692f, 0.888290907f,
        0.888977122f, 0.889663337f, 0.890349552f, 0.891035767f, 0.891498596f, 0.891961425f, 0.892424254f, 0.892887083f,
        0.893349913f, 0.893812742f, 0.894275571f, 0.8947384f, 0.895201229f, 0.895664058f, 0.896029293f, 0.896394528f,
        0.896759763f, 0.897124998f, 0.897490233f, 0.897855468f, 0.898220703f, 0.898585938f, 0.898951173f, 0.899316408f,
        0.899754399f, 0.90019239f, 0.900630381f, 0.901068372f, 0.901506363f, 0.901944353f, 0.902382344f, 0.902820335f,
        0.903258326f, 0.903696317f, 0.904045051f, 0.904393785f, 0.904742519f, 0.905091253f, 0.905439988f, 0.905788722f,
        0.906137456f, 0.90648619f, 0.906834924f, 0.907183658f, 0.907373738f, 0.907563817f, 0.907753896f, 0.907943975f,
        0.908134054f, 0.908324133f, 0.908514213f, 0.908704292f, 0.908894371f, 0.90908445f, 0.90926711f, 0.90944977f,
        0.90963243f, 0.90981509f, 0.90999775f, 0.91018041f, 0.91036307f, 0.91054573f, 0.91072839f, 0.91091105f,
        0.910825224f, 0.910739398f, 0.910653573f, 0.910567747f, 0.910481921f, 0.910396095f, 0.910310269f, 0.910224443f,
        0.910138618f, 0.910052792f, 0.910169041f, 0.91028529f, 0.910401539f, 0.910517788f, 0.910634038f, 0.910750287f,
        0.910866536f, 0.910982785f, 0.911099034f, 0.911215283f, 0.911495718f, 0.911776152f, 0.912056586f, 0.91233702f,
        0.912617454f, 0.912897888f, 0.913178323f, 0.913458757f, 0.913739191f, 0.914019625f, 0.913960303f, 0.913900982f,
        0.91384166f, 0.913782338f, 0.913723017f, 0.913663695f, 0.913604373f, 0.913545052f, 0.91348573f, 0.913426408f,
        0.913685639f, 0.91394487f, 0.914204101f, 0.914463332f, 0.914722563f, 0.914981793f, 0.915241024f, 0.915500255f,
        0.915759486f, 0.916018717f, 0.915964848f, 0.915910978f, 0.915857109f, 0.91580324f, 0.915749371f, 0.915695502f,
        0.915641633f, 0.915587763f, 0.915533894f, 0.915480025f, 0.915516173f, 0.915552322f, 0.91558847f, 0.915624618f,
        0.915660767f, 0.915696915f, 0.915733063f, 0.915769212f, 0.91580536f, 0.915841508f, 0.91569081f, 0.915540112f,
        0.915389413f, 0.915238715f, 0.915088017f, 0.914937318f, 0.91478662f, 0.914635922f, 0.914485223f, 0.914334525f,
        0.914447788f, 0.91456105f, 0.914674313f, 0.914787575f, 0.914900838f, 0.9150141f, 0.915127363f, 0.915240625f,
        0.915353888f, 0.91546715f, 0.915684095f, 0.91590104f, 0.916117985f, 0.91633493f, 0.916551875f, 0.91676882f,
        0.916985765f, 0.91720271f, 0.917419655f, 0.9176366f, 0.917735608f, 0.917834617f, 0.917933625f, 0.918032633f,
        0.918131642f, 0.91823065f, 0.918329658f, 0.918428667f, 0.918527675f, 0.918626683f, 0.918865424f, 0.919104165f,
        0.919342906f, 0.919581647f, 0.919820388f, 0.920059128f, 0.920297869f, 0.92053661f, 0.920775351f, 0.921014092f,
        0.921203662f, 0.921393232f, 0.921582802f, 0.921772372f, 0.921961942f, 0.922151512f, 0.922341082f, 0.922530652f,
        0.922720222f, 0.922909792f, 0.923004354f, 0.923098917f, 0.923193479f, 0.923288042f, 0.923382604f, 0.923477167f,
        0.923571729f, 0.923666292f, 0.923760854f, 0.923855417f, 0.923668726f, 0.923482036f, 0.923295345f, 0.923108655f,
        0.922921964f, 0.922735274f, 0.922548583f, 0.922361893f, 0.922175202f, 0.921988512f, 0.922031691f, 0.922074871f,
        0.922118051f, 0.92216123f, 0.92220441f, 0.92224759f, 0.922290769f, 0.922333949f, 0.922377129f, 0.922420308f,
        0.922655382f, 0.922890455f, 0.923125528f, 0.923360602f, 0.923595675f, 0.923830748f, 0.924065822f, 0.924300895f,
        0.924535968f, 0.924771042f, 0.925042622f, 0.925314202f, 0.925585782f, 0.925857362f, 0.926128942f, 0.926400522f,
        0.926672102f, 0.926943682f, 0.927215262f, 0.927486842f, 0.927715248f, 0.927943653f, 0.928172059f, 0.928400465f,
        0.928628871f, 0.928857277f, 0.929085683f, 0.929314088f, 0.929542494f, 0.9297709f, 0.929834968f, 0.929899037f,
        0.929963105f, 0.930027173f, 0.930091242f, 0.93015531f, 0.930219378f, 0.930283447f, 0.930347515f, 0.930411583f,
        0.930699356f, 0.930987128f, 0.931274901f, 0.931562673f, 0.931850446f, 0.932138218f, 0.932425991f, 0.932713763f,
        0.933001536f, 0.933289308f, 0.839960378f, 0.746631447f, 0.653302516f, 0.559973585f, 0.466644654f, 0.373315723f,
        0.279986793f, 0.186657862f, 0.093328931f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth32({0.17084875f, 0.176969716f, 0.183090682f, 0.189211648f, 0.195332613f,
        0.201453579f, 0.207574545f, 0.213695511f, 0.219816477f, 0.225937443f, 0.232058408f, 0.245359243f, 0.258660078f,
        0.271960913f, 0.285261748f, 0.298562583f, 0.311863418f, 0.325164253f, 0.338465088f, 0.351765923f, 0.365066758f,
        0.379215755f, 0.393364752f, 0.407513748f, 0.421662745f, 0.435811742f, 0.449960738f, 0.464109735f, 0.478258732f,
        0.492407728f, 0.506556725f, 0.512649794f, 0.518742863f, 0.524835933f, 0.530929002f, 0.537022071f, 0.54311514f,
        0.549208209f, 0.555301278f, 0.561394348f, 0.567487417f, 0.569008233f, 0.57052905f, 0.572049867f, 0.573570683f,
        0.5750915f, 0.576612317f, 0.578133133f, 0.57965395f, 0.581174767f, 0.582695583f, 0.583196129f, 0.583696675f,
        0.584197221f, 0.584697767f, 0.585198313f, 0.585698858f, 0.586199404f, 0.58669995f, 0.587200496f, 0.587701042f,
        0.587940117f, 0.588179192f, 0.588418267f, 0.588657342f, 0.588896417f, 0.589135492f, 0.589374567f, 0.589613642f,
        0.589852717f, 0.590091792f, 0.590181816f, 0.59027184f, 0.590361864f, 0.590451888f, 0.590541913f, 0.590631937f,
        0.590721961f, 0.590811985f, 0.590902009f, 0.590992033f, 0.590869947f, 0.59074786f, 0.590625773f, 0.590503687f,
        0.5903816f, 0.590259513f, 0.590137427f, 0.59001534f, 0.589893253f, 0.589771167f, 0.589634701f, 0.589498235f,
        0.589361769f, 0.589225303f, 0.589088838f, 0.588952372f, 0.588815906f, 0.58867944f, 0.588542974f, 0.588406508f,
        0.588408395f, 0.588410282f, 0.588412168f, 0.588414055f, 0.588415942f, 0.588417828f, 0.588419715f, 0.588421602f,
        0.588423488f, 0.588425375f, 0.588480734f, 0.588536093f, 0.588591453f, 0.588646812f, 0.588702171f, 0.58875753f,
        0.588812889f, 0.588868248f, 0.588923608f, 0.588978967f, 0.589029458f, 0.58907995f, 0.589130442f, 0.589180933f,
        0.589231425f, 0.589281917f, 0.589332408f, 0.5893829f, 0.589433392f, 0.589483883f, 0.589594587f, 0.58970529f,
        0.589815993f, 0.589926697f, 0.5900374f, 0.590148103f, 0.590258807f, 0.59036951f, 0.590480213f, 0.590590917f,
        0.590533679f, 0.590476442f, 0.590419204f, 0.590361967f, 0.590304729f, 0.590247492f, 0.590190254f, 0.590133017f,
        0.590075779f, 0.590018542f, 0.590006389f, 0.589994237f, 0.589982084f, 0.589969932f, 0.589957779f, 0.589945627f,
        0.589933474f, 0.589921322f, 0.589909169f, 0.589897017f, 0.58993754f, 0.589978063f, 0.590018587f, 0.59005911f,
        0.590099633f, 0.590140157f, 0.59018068f, 0.590221203f, 0.590261727f, 0.59030225f, 0.590201258f, 0.590100265f,
        0.589999273f, 0.58989828f, 0.589797288f, 0.589696295f, 0.589595303f, 0.58949431f, 0.589393318f, 0.589292325f,
        0.589457419f, 0.589622513f, 0.589787608f, 0.589952702f, 0.590117796f, 0.59028289f, 0.590447984f, 0.590613078f,
        0.590778173f, 0.590943267f, 0.590879753f, 0.59081624f, 0.590752727f, 0.590689213f, 0.5906257f, 0.590562187f,
        0.590498673f, 0.59043516f, 0.590371647f, 0.590308133f, 0.590248664f, 0.590189195f, 0.590129726f, 0.590070257f,
        0.590010788f, 0.589951318f, 0.589891849f, 0.58983238f, 0.589772911f, 0.589713442f, 0.589455258f, 0.589197073f,
        0.588938889f, 0.588680705f, 0.588422521f, 0.588164337f, 0.587906153f, 0.587647968f, 0.587389784f, 0.5871316f,
        0.586933215f, 0.58673483f, 0.586536445f, 0.58633806f, 0.586139675f, 0.58594129f, 0.585742905f, 0.58554452f,
        0.585346135f, 0.58514775f, 0.584937078f, 0.584726407f, 0.584515735f, 0.584305063f, 0.584094392f, 0.58388372f,
        0.583673048f, 0.583462377f, 0.583251705f, 0.583041033f, 0.582733138f, 0.582425242f, 0.582117346f, 0.58180945f,
        0.581501554f, 0.581193658f, 0.580885763f, 0.580577867f, 0.580269971f, 0.579962075f, 0.579744815f, 0.579527555f,
        0.579310295f, 0.579093035f, 0.578875775f, 0.578658515f, 0.578441255f, 0.578223995f, 0.578006735f, 0.577789475f,
        0.577605334f, 0.577421193f, 0.577237053f, 0.577052912f, 0.576868771f, 0.57668463f, 0.576500489f, 0.576316348f,
        0.576132208f, 0.575948067f, 0.575793163f, 0.57563826f, 0.575483357f, 0.575328453f, 0.57517355f, 0.575018647f,
        0.574863743f, 0.57470884f, 0.574553937f, 0.574399033f, 0.574180169f, 0.573961305f, 0.573742441f, 0.573523577f,
        0.573304713f, 0.573085848f, 0.572866984f, 0.57264812f, 0.572429256f, 0.572210392f, 0.572050543f, 0.571890693f,
        0.571730844f, 0.571570995f, 0.571411146f, 0.571251297f, 0.571091448f, 0.570931598f, 0.570771749f, 0.5706119f,
        0.570473113f, 0.570334325f, 0.570195538f, 0.57005675f, 0.569917963f, 0.569779175f, 0.569640388f, 0.5695016f,
        0.569362813f, 0.569224025f, 0.569129288f, 0.569034552f, 0.568939815f, 0.568845078f, 0.568750342f, 0.568655605f,
        0.568560868f, 0.568466132f, 0.568371395f, 0.568276658f, 0.568245926f, 0.568215193f, 0.568184461f, 0.568153728f,
        0.568122996f, 0.568092263f, 0.568061531f, 0.568030798f, 0.568000066f, 0.567969333f, 0.567820551f, 0.567671768f,
        0.567522986f, 0.567374203f, 0.567225421f, 0.567076638f, 0.566927856f, 0.566779073f, 0.566630291f, 0.566481508f,
        0.566464265f, 0.566447022f, 0.566429778f, 0.566412535f, 0.566395292f, 0.566378048f, 0.566360805f, 0.566343562f,
        0.566326318f, 0.566309075f, 0.509678168f, 0.45304726f, 0.396416353f, 0.339785445f, 0.283154538f, 0.22652363f,
        0.169892723f, 0.113261815f, 0.056630908f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth33({0.1442082f, 0.149032892f, 0.153857583f, 0.158682275f, 0.163506967f,
        0.168331658f, 0.17315635f, 0.177981042f, 0.182805733f, 0.187630425f, 0.192455117f, 0.200393591f, 0.208332065f,
        0.216270539f, 0.224209013f, 0.232147488f, 0.240085962f, 0.248024436f, 0.25596291f, 0.263901384f, 0.271839858f,
        0.27773725f, 0.283634642f, 0.289532033f, 0.295429425f, 0.301326817f, 0.307224208f, 0.3131216f, 0.319018992f,
        0.324916383f, 0.330813775f, 0.33277435f, 0.334734925f, 0.3366955f, 0.338656075f, 0.34061665f, 0.342577225f,
        0.3445378f, 0.346498375f, 0.34845895f, 0.350419525f, 0.351069437f, 0.351719348f, 0.35236926f, 0.353019172f,
        0.353669083f, 0.354318995f, 0.354968907f, 0.355618818f, 0.35626873f, 0.356918642f, 0.357349918f, 0.357781193f,
        0.358212469f, 0.358643745f, 0.359075021f, 0.359506297f, 0.359937573f, 0.360368848f, 0.360800124f, 0.3612314f,
        0.361434438f, 0.361637477f, 0.361840515f, 0.362043553f, 0.362246592f, 0.36244963f, 0.362652668f, 0.362855707f,
        0.363058745f, 0.363261783f, 0.363232583f, 0.363203382f, 0.363174181f, 0.36314498f, 0.363115779f, 0.363086578f,
        0.363057378f, 0.363028177f, 0.362998976f, 0.362969775f, 0.3627538f, 0.362537825f, 0.36232185f, 0.362105875f,
        0.3618899f, 0.361673925f, 0.36145795f, 0.361241975f, 0.361026f, 0.360810025f, 0.360602629f, 0.360395233f,
        0.360187838f, 0.359980442f, 0.359773046f, 0.35956565f, 0.359358254f, 0.359150858f, 0.358943463f, 0.358736067f,
        0.35867327f, 0.358610473f, 0.358547677f, 0.35848488f, 0.358422083f, 0.358359287f, 0.35829649f, 0.358233693f,
        0.358170897f, 0.3581081f, 0.358142819f, 0.358177538f, 0.358212258f, 0.358246977f, 0.358281696f, 0.358316415f,
        0.358351134f, 0.358385853f, 0.358420573f, 0.358455292f, 0.358528363f, 0.358601435f, 0.358674507f, 0.358747578f,
        0.35882065f, 0.358893722f, 0.358966793f, 0.359039865f, 0.359112937f, 0.359186008f, 0.359308604f, 0.3594312f,
        0.359553796f, 0.359676392f, 0.359798988f, 0.359921583f, 0.360044179f, 0.360166775f, 0.360289371f, 0.360411967f,
        0.360416353f, 0.360420738f, 0.360425124f, 0.36042951f, 0.360433896f, 0.360438282f, 0.360442668f, 0.360447053f,
        0.360451439f, 0.360455825f, 0.360466266f, 0.360476707f, 0.360487148f, 0.360497588f, 0.360508029f, 0.36051847f,
        0.360528911f, 0.360539352f, 0.360549793f, 0.360560233f, 0.360587145f, 0.360614057f, 0.360640968f, 0.36066788f,
        0.360694792f, 0.360721703f, 0.360748615f, 0.360775527f, 0.360802438f, 0.36082935f, 0.360789041f, 0.360748732f,
        0.360708423f, 0.360668113f, 0.360627804f, 0.360587495f, 0.360547186f, 0.360506877f, 0.360466568f, 0.360426258f,
        0.360568448f, 0.360710637f, 0.360852826f, 0.360995015f, 0.361137204f, 0.361279393f, 0.361421583f, 0.361563772f,
        0.361705961f, 0.36184815f, 0.361838397f, 0.361828643f, 0.36181889f, 0.361809137f, 0.361799383f, 0.36178963f,
        0.361779877f, 0.361770123f, 0.36176037f, 0.361750617f, 0.361708134f, 0.361665652f, 0.361623169f, 0.361580687f,
        0.361538204f, 0.361495722f, 0.361453239f, 0.361410757f, 0.361368274f, 0.361325792f, 0.361125494f, 0.360925197f,
        0.360724899f, 0.360524602f, 0.360324304f, 0.360124007f, 0.359923709f, 0.359723412f, 0.359523114f, 0.359322817f,
        0.359143894f, 0.358964972f, 0.358786049f, 0.358607127f, 0.358428204f, 0.358249282f, 0.358070359f, 0.357891437f,
        0.357712514f, 0.357533592f, 0.357323623f, 0.357113653f, 0.356903684f, 0.356693715f, 0.356483746f, 0.356273777f,
        0.356063808f, 0.355853838f, 0.355643869f, 0.3554339f, 0.355132009f, 0.354830118f, 0.354528228f, 0.354226337f,
        0.353924446f, 0.353622555f, 0.353320664f, 0.353018773f, 0.352716883f, 0.352414992f, 0.352163242f, 0.351911492f,
        0.351659742f, 0.351407992f, 0.351156242f, 0.350904492f, 0.350652742f, 0.350400992f, 0.350149242f, 0.349897492f,
        0.349671781f, 0.34944607f, 0.349220359f, 0.348994648f, 0.348768938f, 0.348543227f, 0.348317516f, 0.348091805f,
        0.347866094f, 0.347640383f, 0.347423158f, 0.347205932f, 0.346988706f, 0.34677148f, 0.346554254f, 0.346337028f,
        0.346119803f, 0.345902577f, 0.345685351f, 0.345468125f, 0.345193043f, 0.344917962f, 0.34464288f, 0.344367798f,
        0.344092717f, 0.343817635f, 0.343542553f, 0.343267472f, 0.34299239f, 0.342717308f, 0.342462254f, 0.3422072f,
        0.341952146f, 0.341697092f, 0.341442038f, 0.341186983f, 0.340931929f, 0.340676875f, 0.340421821f, 0.340166767f,
        0.33991013f, 0.339653493f, 0.339396857f, 0.33914022f, 0.338883583f, 0.338626947f, 0.33837031f, 0.338113673f,
        0.337857037f, 0.3376004f, 0.337371173f, 0.337141945f, 0.336912718f, 0.33668349f, 0.336454263f, 0.336225035f,
        0.335995808f, 0.33576658f, 0.335537353f, 0.335308125f, 0.335160589f, 0.335013053f, 0.334865518f, 0.334717982f,
        0.334570446f, 0.33442291f, 0.334275374f, 0.334127838f, 0.333980303f, 0.333832767f, 0.333629798f, 0.33342683f,
        0.333223862f, 0.333020893f, 0.332817925f, 0.332614957f, 0.332411988f, 0.33220902f, 0.332006052f, 0.331803083f,
        0.331677119f, 0.331551155f, 0.331425191f, 0.331299227f, 0.331173263f, 0.331047298f, 0.330921334f, 0.33079537f,
        0.330669406f, 0.330543442f, 0.297489098f, 0.264434753f, 0.231380409f, 0.198326065f, 0.165271721f, 0.132217377f,
        0.099163033f, 0.066108688f, 0.033054344f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth34({0.105192122f, 0.107805515f, 0.110418907f, 0.1130323f, 0.115645693f,
        0.118259086f, 0.120872479f, 0.123485872f, 0.126099264f, 0.128712657f, 0.13132605f, 0.134453238f, 0.137580427f,
        0.140707615f, 0.143834803f, 0.146961992f, 0.15008918f, 0.153216368f, 0.156343557f, 0.159470745f, 0.162597933f,
        0.164354708f, 0.166111483f, 0.167868258f, 0.169625033f, 0.171381808f, 0.173138583f, 0.174895358f, 0.176652133f,
        0.178408908f, 0.180165683f, 0.180740801f, 0.181315918f, 0.181891036f, 0.182466153f, 0.183041271f, 0.183616388f,
        0.184191506f, 0.184766623f, 0.185341741f, 0.185916858f, 0.186278556f, 0.186640253f, 0.187001951f, 0.187363648f,
        0.187725346f, 0.188087043f, 0.188448741f, 0.188810438f, 0.189172136f, 0.189533833f, 0.189866533f, 0.190199233f,
        0.190531933f, 0.190864633f, 0.191197333f, 0.191530033f, 0.191862733f, 0.192195433f, 0.192528133f, 0.192860833f,
        0.192997978f, 0.193135122f, 0.193272266f, 0.19340941f, 0.193546554f, 0.193683698f, 0.193820843f, 0.193957987f,
        0.194095131f, 0.194232275f, 0.194187385f, 0.194142495f, 0.194097605f, 0.194052715f, 0.194007825f, 0.193962935f,
        0.193918045f, 0.193873155f, 0.193828265f, 0.193783375f, 0.193637849f, 0.193492323f, 0.193346798f, 0.193201272f,
        0.193055746f, 0.19291022f, 0.192764694f, 0.192619168f, 0.192473643f, 0.192328117f, 0.192201278f, 0.19207444f,
        0.191947602f, 0.191820763f, 0.191693925f, 0.191567087f, 0.191440248f, 0.19131341f, 0.191186572f, 0.191059733f,
        0.191038389f, 0.191017045f, 0.190995701f, 0.190974357f, 0.190953013f, 0.190931668f, 0.190910324f, 0.19088898f,
        0.190867636f, 0.190846292f, 0.190886298f, 0.190926305f, 0.190966312f, 0.191006318f, 0.191046325f, 0.191086332f,
        0.191126338f, 0.191166345f, 0.191206352f, 0.191246358f, 0.191279655f, 0.191312952f, 0.191346248f, 0.191379545f,
        0.191412842f, 0.191446138f, 0.191479435f, 0.191512732f, 0.191546028f, 0.191579325f, 0.191626543f, 0.19167376f,
        0.191720978f, 0.191768195f, 0.191815413f, 0.19186263f, 0.191909848f, 0.191957065f, 0.192004283f, 0.1920515f,
        0.192056469f, 0.192061438f, 0.192066408f, 0.192071377f, 0.192076346f, 0.192081315f, 0.192086284f, 0.192091253f,
        0.192096223f, 0.192101192f, 0.192108111f, 0.19211503f, 0.192121949f, 0.192128868f, 0.192135788f, 0.192142707f,
        0.192149626f, 0.192156545f, 0.192163464f, 0.192170383f, 0.192184616f, 0.192198848f, 0.192213081f, 0.192227313f,
        0.192241546f, 0.192255778f, 0.192270011f, 0.192284243f, 0.192298476f, 0.192312708f, 0.192287254f, 0.1922618f,
        0.192236346f, 0.192210892f, 0.192185438f, 0.192159983f, 0.192134529f, 0.192109075f, 0.192083621f, 0.192058167f,
        0.192115572f, 0.192172977f, 0.192230382f, 0.192287787f, 0.192345192f, 0.192402597f, 0.192460002f, 0.192517407f,
        0.192574812f, 0.192632217f, 0.192607146f, 0.192582075f, 0.192557004f, 0.192531933f, 0.192506863f, 0.192481792f,
        0.192456721f, 0.19243165f, 0.192406579f, 0.192381508f, 0.192336265f, 0.192291022f, 0.192245778f, 0.192200535f,
        0.192155292f, 0.192110048f, 0.192064805f, 0.192019562f, 0.191974318f, 0.191929075f, 0.191795643f, 0.19166221f,
        0.191528778f, 0.191395345f, 0.191261913f, 0.19112848f, 0.190995048f, 0.190861615f, 0.190728183f, 0.19059475f,
        0.190473893f, 0.190353037f, 0.19023218f, 0.190111323f, 0.189990467f, 0.18986961f, 0.189748753f, 0.189627897f,
        0.18950704f, 0.189386183f, 0.189246166f, 0.189106148f, 0.188966131f, 0.188826113f, 0.188686096f, 0.188546078f,
        0.188406061f, 0.188266043f, 0.188126026f, 0.187986008f, 0.187774149f, 0.18756229f, 0.187350431f, 0.187138572f,
        0.186926713f, 0.186714853f, 0.186502994f, 0.186291135f, 0.186079276f, 0.185867417f, 0.185678197f, 0.185488977f,
        0.185299757f, 0.185110537f, 0.184921317f, 0.184732097f, 0.184542877f, 0.184353657f, 0.184164437f, 0.183975217f,
        0.183809435f, 0.183643653f, 0.183477872f, 0.18331209f, 0.183146308f, 0.182980527f, 0.182814745f, 0.182648963f,
        0.182483182f, 0.1823174f, 0.182170823f, 0.182024245f, 0.181877668f, 0.18173109f, 0.181584513f, 0.181437935f,
        0.181291358f, 0.18114478f, 0.180998203f, 0.180851625f, 0.180691691f, 0.180531757f, 0.180371823f, 0.180211888f,
        0.180051954f, 0.17989202f, 0.179732086f, 0.179572152f, 0.179412218f, 0.179252283f, 0.179104241f, 0.178956198f,
        0.178808156f, 0.178660113f, 0.178512071f, 0.178364028f, 0.178215986f, 0.178067943f, 0.177919901f, 0.177771858f,
        0.177599149f, 0.17742644f, 0.177253731f, 0.177081022f, 0.176908313f, 0.176735603f, 0.176562894f, 0.176390185f,
        0.176217476f, 0.176044767f, 0.175874165f, 0.175703563f, 0.175532962f, 0.17536236f, 0.175191758f, 0.175021157f,
        0.174850555f, 0.174679953f, 0.174509352f, 0.17433875f, 0.174241704f, 0.174144658f, 0.174047613f, 0.173950567f,
        0.173853521f, 0.173756475f, 0.173659429f, 0.173562383f, 0.173465338f, 0.173368292f, 0.173250319f, 0.173132347f,
        0.173014374f, 0.172896402f, 0.172778429f, 0.172660457f, 0.172542484f, 0.172424512f, 0.172306539f, 0.172188567f,
        0.172108589f, 0.172028612f, 0.171948634f, 0.171868657f, 0.171788679f, 0.171708702f, 0.171628724f, 0.171548747f,
        0.171468769f, 0.171388792f, 0.154249913f, 0.137111033f, 0.119972154f, 0.102833275f, 0.085694396f, 0.068555517f,
        0.051416638f, 0.034277758f, 0.017138879f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth35({0.067959987f, 0.068836284f, 0.069712581f, 0.070588878f, 0.071465175f,
        0.072341473f, 0.07321777f, 0.074094067f, 0.074970364f, 0.075846661f, 0.076722958f, 0.07743884f, 0.078154722f,
        0.078870604f, 0.079586486f, 0.080302368f, 0.081018249f, 0.081734131f, 0.082450013f, 0.083165895f, 0.083881777f,
        0.084234856f, 0.084587934f, 0.084941013f, 0.085294092f, 0.085647171f, 0.08600025f, 0.086353329f, 0.086706407f,
        0.087059486f, 0.087412565f, 0.087559802f, 0.087707038f, 0.087854275f, 0.088001511f, 0.088148748f, 0.088295984f,
        0.088443221f, 0.088590457f, 0.088737694f, 0.08888493f, 0.089040721f, 0.089196512f, 0.089352304f, 0.089508095f,
        0.089663886f, 0.089819677f, 0.089975468f, 0.090131259f, 0.090287051f, 0.090442842f, 0.090585145f, 0.090727448f,
        0.090869752f, 0.091012055f, 0.091154358f, 0.091296662f, 0.091438965f, 0.091581268f, 0.091723572f, 0.091865875f,
        0.091883666f, 0.091901456f, 0.091919247f, 0.091937038f, 0.091954828f, 0.091972619f, 0.09199041f, 0.0920082f,
        0.092025991f, 0.092043782f, 0.091974486f, 0.09190519f, 0.091835894f, 0.091766598f, 0.091697303f, 0.091628007f,
        0.091558711f, 0.091489415f, 0.091420119f, 0.091350823f, 0.091255008f, 0.091159193f, 0.091063378f, 0.090967563f,
        0.090871748f, 0.090775932f, 0.090680117f, 0.090584302f, 0.090488487f, 0.090392672f, 0.090328126f, 0.090263581f,
        0.090199035f, 0.09013449f, 0.090069944f, 0.090005399f, 0.089940853f, 0.089876308f, 0.089811762f, 0.089747217f,
        0.089737453f, 0.089727689f, 0.089717925f, 0.089708161f, 0.089698397f, 0.089688633f, 0.089678869f, 0.089669105f,
        0.089659341f, 0.089649577f, 0.08966585f, 0.089682123f, 0.089698396f, 0.089714669f, 0.089730943f, 0.089747216f,
        0.089763489f, 0.089779762f, 0.089796035f, 0.089812308f, 0.089818988f, 0.089825667f, 0.089832346f, 0.089839025f,
        0.089845704f, 0.089852383f, 0.089859063f, 0.089865742f, 0.089872421f, 0.0898791f, 0.089888167f, 0.089897234f,
        0.089906301f, 0.089915367f, 0.089924434f, 0.089933501f, 0.089942568f, 0.089951635f, 0.089960702f, 0.089969768f,
        0.089969133f, 0.089968498f, 0.089967863f, 0.089967228f, 0.089966593f, 0.089965958f, 0.089965323f, 0.089964688f,
        0.089964053f, 0.089963418f, 0.089968449f, 0.08997348f, 0.089978511f, 0.089983542f, 0.089988573f, 0.089993604f,
        0.089998635f, 0.090003666f, 0.090008697f, 0.090013728f, 0.090018559f, 0.09002339f, 0.09002822f, 0.090033051f,
        0.090037882f, 0.090042712f, 0.090047543f, 0.090052374f, 0.090057204f, 0.090062035f, 0.090042038f, 0.090022041f,
        0.090002044f, 0.089982047f, 0.08996205f, 0.089942053f, 0.089922056f, 0.089902059f, 0.089882062f, 0.089862065f,
        0.089866237f, 0.089870408f, 0.08987458f, 0.089878752f, 0.089882923f, 0.089887095f, 0.089891267f, 0.089895438f,
        0.08989961f, 0.089903782f, 0.089864664f, 0.089825545f, 0.089786427f, 0.089747309f, 0.089708191f, 0.089669073f,
        0.089629955f, 0.089590836f, 0.089551718f, 0.0895126f, 0.089477832f, 0.089443064f, 0.089408296f, 0.089373528f,
        0.08933876f, 0.089303992f, 0.089269224f, 0.089234456f, 0.089199688f, 0.08916492f, 0.089100586f, 0.089036251f,
        0.088971917f, 0.088907582f, 0.088843248f, 0.088778913f, 0.088714579f, 0.088650244f, 0.08858591f, 0.088521575f,
        0.08847569f, 0.088429804f, 0.088383919f, 0.088338033f, 0.088292148f, 0.088246262f, 0.088200377f, 0.088154491f,
        0.088108606f, 0.08806272f, 0.088005834f, 0.087948947f, 0.087892061f, 0.087835175f, 0.087778288f, 0.087721402f,
        0.087664516f, 0.087607629f, 0.087550743f, 0.087493857f, 0.087389368f, 0.087284879f, 0.08718039f, 0.087075901f,
        0.086971413f, 0.086866924f, 0.086762435f, 0.086657946f, 0.086553457f, 0.086448968f, 0.086355073f, 0.086261178f,
        0.086167283f, 0.086073388f, 0.085979493f, 0.085885597f, 0.085791702f, 0.085697807f, 0.085603912f, 0.085510017f,
        0.085436835f, 0.085363653f, 0.085290471f, 0.085217289f, 0.085144108f, 0.085070926f, 0.084997744f, 0.084924562f,
        0.08485138f, 0.084778198f, 0.084720503f, 0.084662807f, 0.084605112f, 0.084547416f, 0.084489721f, 0.084432025f,
        0.08437433f, 0.084316634f, 0.084258939f, 0.084201243f, 0.084149172f, 0.0840971f, 0.084045029f, 0.083992957f,
        0.083940886f, 0.083888814f, 0.083836743f, 0.083784671f, 0.0837326f, 0.083680528f, 0.083633793f, 0.083587058f,
        0.083540322f, 0.083493587f, 0.083446852f, 0.083400116f, 0.083353381f, 0.083306646f, 0.08325991f, 0.083213175f,
        0.08314572f, 0.083078266f, 0.083010811f, 0.082943356f, 0.082875902f, 0.082808447f, 0.082740992f, 0.082673538f,
        0.082606083f, 0.082538628f, 0.082461047f, 0.082383465f, 0.082305883f, 0.082228302f, 0.08215072f, 0.082073138f,
        0.081995557f, 0.081917975f, 0.081840393f, 0.081762812f, 0.081729338f, 0.081695864f, 0.08166239f, 0.081628916f,
        0.081595443f, 0.081561969f, 0.081528495f, 0.081495021f, 0.081461547f, 0.081428073f, 0.081378168f, 0.081328263f,
        0.081278358f, 0.081228453f, 0.081178548f, 0.081128642f, 0.081078737f, 0.081028832f, 0.080978927f, 0.080929022f,
        0.080903787f, 0.080878553f, 0.080853319f, 0.080828084f, 0.08080285f, 0.080777616f, 0.080752381f, 0.080727147f,
        0.080701913f, 0.080676678f, 0.072609011f, 0.064541343f, 0.056473675f, 0.048406007f, 0.040338339f, 0.032270671f,
        0.024203004f, 0.016135336f, 0.008067668f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth36({0.031023515f, 0.031120073f, 0.031216632f, 0.03131319f, 0.031409748f,
        0.031506307f, 0.031602865f, 0.031699423f, 0.031795982f, 0.03189254f, 0.031989098f, 0.032018177f, 0.032047256f,
        0.032076335f, 0.032105414f, 0.032134493f, 0.032163572f, 0.032192651f, 0.03222173f, 0.032250809f, 0.032279888f,
        0.0323079f, 0.032335911f, 0.032363923f, 0.032391934f, 0.032419946f, 0.032447957f, 0.032475969f, 0.03250398f,
        0.032531992f, 0.032560003f, 0.032579362f, 0.032598721f, 0.03261808f, 0.032637439f, 0.032656798f, 0.032676157f,
        0.032695516f, 0.032714875f, 0.032734234f, 0.032753593f, 0.032760412f, 0.03276723f, 0.032774048f, 0.032780866f,
        0.032787684f, 0.032794502f, 0.032801321f, 0.032808139f, 0.032814957f, 0.032821775f, 0.032821853f, 0.03282193f,
        0.032822008f, 0.032822085f, 0.032822163f, 0.03282224f, 0.032822318f, 0.032822395f, 0.032822473f, 0.03282255f,
        0.032802458f, 0.032782366f, 0.032762274f, 0.032742181f, 0.032722089f, 0.032701997f, 0.032681905f, 0.032661813f,
        0.032641721f, 0.032621628f, 0.03260782f, 0.032594011f, 0.032580202f, 0.032566393f, 0.032552584f, 0.032538775f,
        0.032524967f, 0.032511158f, 0.032497349f, 0.03248354f, 0.032475542f, 0.032467544f, 0.032459547f, 0.032451549f,
        0.032443551f, 0.032435553f, 0.032427555f, 0.032419557f, 0.03241156f, 0.032403562f, 0.032395802f, 0.032388042f,
        0.032380282f, 0.032372522f, 0.032364762f, 0.032357002f, 0.032349242f, 0.032341482f, 0.032333722f, 0.032325962f,
        0.03232411f, 0.032322258f, 0.032320407f, 0.032318555f, 0.032316703f, 0.032314852f, 0.032313f, 0.032311148f,
        0.032309297f, 0.032307445f, 0.032304859f, 0.032302274f, 0.032299688f, 0.032297102f, 0.032294517f, 0.032291931f,
        0.032289345f, 0.03228676f, 0.032284174f, 0.032281588f, 0.032268883f, 0.032256177f, 0.032243471f, 0.032230766f,
        0.03221806f, 0.032205354f, 0.032192649f, 0.032179943f, 0.032167237f, 0.032154532f, 0.032147745f, 0.032140959f,
        0.032134173f, 0.032127386f, 0.0321206f, 0.032113814f, 0.032107027f, 0.032100241f, 0.032093455f, 0.032086668f,
        0.032076223f, 0.032065778f, 0.032055333f, 0.032044888f, 0.032034443f, 0.032023997f, 0.032013552f, 0.032003107f,
        0.031992662f, 0.031982217f, 0.031980094f, 0.03197797f, 0.031975847f, 0.031973724f, 0.031971601f, 0.031969478f,
        0.031967355f, 0.031965231f, 0.031963108f, 0.031960985f, 0.031962523f, 0.03196406f, 0.031965598f, 0.031967136f,
        0.031968673f, 0.031970211f, 0.031971749f, 0.031973286f, 0.031974824f, 0.031976362f, 0.03197069f, 0.031965018f,
        0.031959346f, 0.031953674f, 0.031948003f, 0.031942331f, 0.031936659f, 0.031930987f, 0.031925315f, 0.031919643f,
        0.031922461f, 0.031925278f, 0.031928095f, 0.031930913f, 0.03193373f, 0.031936547f, 0.031939365f, 0.031942182f,
        0.031944999f, 0.031947817f, 0.031934769f, 0.031921721f, 0.031908674f, 0.031895626f, 0.031882578f, 0.031869531f,
        0.031856483f, 0.031843435f, 0.031830388f, 0.03181734f, 0.031814895f, 0.03181245f, 0.031810006f, 0.031807561f,
        0.031805116f, 0.031802671f, 0.031800226f, 0.031797781f, 0.031795337f, 0.031792892f, 0.031786646f, 0.031780399f,
        0.031774153f, 0.031767907f, 0.031761661f, 0.031755415f, 0.031749169f, 0.031742922f, 0.031736676f, 0.03173043f,
        0.031742255f, 0.031754081f, 0.031765906f, 0.031777731f, 0.031789557f, 0.031801382f, 0.031813207f, 0.031825033f,
        0.031836858f, 0.031848683f, 0.031857231f, 0.031865779f, 0.031874327f, 0.031882875f, 0.031891423f, 0.031899971f,
        0.031908519f, 0.031917067f, 0.031925615f, 0.031934163f, 0.031931159f, 0.031928154f, 0.03192515f, 0.031922145f,
        0.031919141f, 0.031916136f, 0.031913132f, 0.031910127f, 0.031907123f, 0.031904118f, 0.03190447f, 0.031904822f,
        0.031905174f, 0.031905526f, 0.031905878f, 0.03190623f, 0.031906582f, 0.031906934f, 0.031907286f, 0.031907638f,
        0.031910956f, 0.031914274f, 0.031917592f, 0.03192091f, 0.031924228f, 0.031927545f, 0.031930863f, 0.031934181f,
        0.031937499f, 0.031940817f, 0.031948264f, 0.031955711f, 0.031963158f, 0.031970605f, 0.031978053f, 0.0319855f,
        0.031992947f, 0.032000394f, 0.032007841f, 0.032015288f, 0.032025932f, 0.032036576f, 0.032047219f, 0.032057863f,
        0.032068507f, 0.03207915f, 0.032089794f, 0.032100438f, 0.032111081f, 0.032121725f, 0.032131217f, 0.032140709f,
        0.032150201f, 0.032159692f, 0.032169184f, 0.032178676f, 0.032188168f, 0.03219766f, 0.032207152f, 0.032216643f,
        0.032217578f, 0.032218512f, 0.032219446f, 0.03222038f, 0.032221314f, 0.032222248f, 0.032223183f, 0.032224117f,
        0.032225051f, 0.032225985f, 0.032227669f, 0.032229353f, 0.032231037f, 0.032232721f, 0.032234405f, 0.032236089f,
        0.032237773f, 0.032239457f, 0.032241141f, 0.032242825f, 0.032251451f, 0.032260077f, 0.032268704f, 0.03227733f,
        0.032285956f, 0.032294582f, 0.032303208f, 0.032311834f, 0.032320461f, 0.032329087f, 0.032333441f, 0.032337796f,
        0.032342151f, 0.032346505f, 0.03235086f, 0.032355215f, 0.032359569f, 0.032363924f, 0.032368279f, 0.032372633f,
        0.032385684f, 0.032398735f, 0.032411786f, 0.032424837f, 0.032437888f, 0.032450939f, 0.03246399f, 0.032477041f,
        0.032490092f, 0.032503143f, 0.029252829f, 0.026002515f, 0.0227522f, 0.019501886f, 0.016251572f, 0.013001257f,
        0.009750943f, 0.006500629f, 0.003250314f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    // Brown skin, light skin, sky, folliage, ...
    static constexpr std::array<Spectrum, 24> Patch{
        Macbeth01,
        Macbeth02,
        Macbeth03,
        Macbeth04,
        Macbeth05,
        Macbeth06,
        Macbeth11,
        Macbeth12,
        Macbeth13,
        Macbeth14,
        Macbeth15,
        Macbeth16,
        Macbeth21,
        Macbeth22,
        Macbeth23,
        Macbeth24,
        Macbeth25,
        Macbeth26,
        Macbeth31,
        Macbeth32,
        Macbeth33,
        Macbeth34,
        Macbeth35,
        Macbeth36,
    };

    static std::vector<Tristimulus> reference(const Spectrum &light, const Observer &obs = CIE1931)
    {
        std::vector<Tristimulus> result(Patch.size());
        for (int i = 0; i < Patch.size(); i++)
        {
            result[i] = obs.fromReflectanceAndLight(Patch[i], light);
        }
        return result;
    }

} // namespace Macbeth



// ---
namespace SOLVER
{
    /*
     * svdcomp - SVD decomposition routine.
     * Takes an mxn matrix a and decomposes it into udv, where u,v are
     * left and right orthogonal transformation matrices, and d is a
     * diagonal matrix of singular values.
     *
     * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
     * code from Numerical Recipes adapted by Luke Tierney and David Betz.
     *
     * Input to dsvd is as follows:
     *   a = mxn matrix to be decomposed, gets overwritten with u
     *   m = row dimension of a
     *   n = column dimension of a
     *   w = returns the vector of singular values of a
     *   v = returns the right orthogonal transformation matrix
     */

    class Vector
    {
      public:
        std::vector<float> v_;

        Vector() { ; }
        Vector(int d) : v_(d) { ; }
        float &operator[](int i) { return v_[i]; }

        const float *data(void) { return v_.data(); }
    };

    class Matrix
    {
      public:
        std::vector<std::vector<float>> v_;
        Matrix() { ; }
        Matrix(int r, int c)
        {
            v_.resize(r);
            for (auto &v : v_)
                v.resize(c);
        }
        const int           rows(void) { return (int)(v_.size()); }
        const int           cols(void) { return (int)(v_.size() ? v_[0].size() : 0); }
        float &             v(int r, int c) { return v_[r][c]; }
        std::vector<float> &operator[](int r) { return v_[r]; }
    };

    static inline float DSIGN(float a, float b) { return (b > 0.f) ? fabsf(a) : -fabsf(a); }
    static inline float PYTHAG(float a, float b)
    {
        float at = fabsf(a);
        float bt = fabsf(b);
        if (at > bt)
        {
            float ct = bt / at;
            return at * sqrtf(1.0f + ct * ct);
        }
        else if (bt > 0.0)
        {
            float ct = at / bt;
            return bt * sqrtf(1.0f + ct * ct);
        }
        return 0.f;
    }

    static int svdcmp(Matrix &a, Vector &w, Matrix &v)
    {
        int   flag, its, j, jj, k, l, nm;
        float c, f, h, s, x, y, z;
        float anorm = 0.0, g = 0.0, scale = 0.0;

        int m = a.rows();
        int n = a.cols();

        if (m < n)
        {
            return (-1);
        }

        Vector rv1(n);

        /* Householder reduction to bidiagonal form */
        for (int i = 0; i < n; i++)
        {
            /* left-hand reduction */
            l      = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < m)
            {
                for (k = i; k < m; k++)
                    scale += fabs(a[k][i]);
                if (scale)
                {
                    for (k = i; k < m; k++)
                    {
                        a[k][i] = (a[k][i] / scale);
                        s += (a[k][i] * a[k][i]);
                    }
                    f       = a[i][i];
                    g       = -DSIGN(sqrtf(s), f);
                    h       = f * g - s;
                    a[i][i] = (f - g);
                    if (i != n - 1)
                    {
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0, k = i; k < m; k++)
                                s += (a[k][i] * a[k][j]);
                            f = s / h;
                            for (k = i; k < m; k++)
                                a[k][j] += (f * a[k][i]);
                        }
                    }
                    for (k = i; k < m; k++)
                        a[k][i] = (a[k][i] * scale);
                }
            }
            w[i] = (scale * g);

            /* right-hand reduction */
            g = s = scale = 0.0;
            if (i < m && i != n - 1)
            {
                for (k = l; k < n; k++)
                    scale += fabs(a[i][k]);
                if (scale)
                {
                    for (k = l; k < n; k++)
                    {
                        a[i][k] = (a[i][k] / scale);
                        s += (a[i][k] * a[i][k]);
                    }
                    f       = a[i][l];
                    g       = -DSIGN(sqrtf(s), f);
                    h       = f * g - s;
                    a[i][l] = (f - g);
                    for (k = l; k < n; k++)
                        rv1[k] = a[i][k] / h;
                    if (i != m - 1)
                    {
                        for (j = l; j < m; j++)
                        {
                            for (s = 0.0, k = l; k < n; k++)
                                s += (a[j][k] * a[i][k]);
                            for (k = l; k < n; k++)
                                a[j][k] += (s * rv1[k]);
                        }
                    }
                    for (k = l; k < n; k++)
                        a[i][k] = (a[i][k] * scale);
                }
            }
            anorm = std::max(anorm, (fabsf(w[i]) + fabsf(rv1[i])));
        }

        /* accumulate the right-hand transformation */
        for (int i = n - 1; i >= 0; i--)
        {
            if (i < n - 1)
            {
                if (g)
                {
                    for (j = l; j < n; j++)
                        v[j][i] = ((a[i][j] / a[i][l]) / g);
                    /* float division to avoid underflow */
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += (a[i][k] * v[k][j]);
                        for (k = l; k < n; k++)
                            v[k][j] += (s * v[k][i]);
                    }
                }
                for (j = l; j < n; j++)
                    v[i][j] = v[j][i] = 0.0;
            }
            v[i][i] = 1.0;
            g       = rv1[i];
            l       = i;
        }

        /* accumulate the left-hand transformation */
        for (int i = n - 1; i >= 0; i--)
        {
            l = i + 1;
            g = w[i];
            if (i < n - 1)
                for (j = l; j < n; j++)
                    a[i][j] = 0.0;
            if (g)
            {
                g = 1.0f / g;
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (a[k][i] * a[k][j]);
                        f = (s / a[i][i]) * g;
                        for (k = i; k < m; k++)
                            a[k][j] += (f * a[k][i]);
                    }
                }
                for (j = i; j < m; j++)
                    a[j][i] = (a[j][i] * g);
            }
            else
            {
                for (j = i; j < m; j++)
                    a[j][i] = 0.0;
            }
            ++a[i][i];
        }

        /* diagonalize the bidiagonal form */
        for (k = n - 1; k >= 0; k--)
        { /* loop over singular values */
            for (its = 0; its < 30; its++)
            { /* loop over allowed iterations */
                flag = 1;
                for (l = k; l >= 0; l--)
                { /* test for splitting */
                    nm = l - 1;
                    if (fabs(rv1[l]) + anorm == anorm)
                    {
                        flag = 0;
                        break;
                    }
                    if (fabs(w[nm]) + anorm == anorm)
                        break;
                }
                if (flag)
                {
                    c = 0.0;
                    s = 1.0;
                    for (int i = l; i <= k; i++)
                    {
                        f = s * rv1[i];
                        if (fabs(f) + anorm != anorm)
                        {
                            g    = w[i];
                            h    = PYTHAG(f, g);
                            w[i] = h;
                            h    = 1.0f / h;
                            c    = g * h;
                            s    = (-f * h);
                            for (j = 0; j < m; j++)
                            {
                                y        = a[j][nm];
                                z        = a[j][i];
                                a[j][nm] = (y * c + z * s);
                                a[j][i]  = (z * c - y * s);
                            }
                        }
                    }
                }
                z = w[k];
                if (l == k)
                { /* convergence */
                    if (z < 0.0)
                    { /* make singular value nonnegative */
                        w[k] = (-z);
                        for (j = 0; j < n; j++)
                            v[j][k] = (-v[j][k]);
                    }
                    break;
                }
                if (its >= 30)
                {
                    return (-1);
                }

                /* shift from bottom 2 x 2 minor */
                x  = w[l];
                nm = k - 1;
                y  = w[nm];
                g  = rv1[nm];
                h  = rv1[k];
                f  = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
                g  = PYTHAG(f, 1.0f);
                f  = ((x - z) * (x + z) + h * ((y / (f + DSIGN(g, f))) - h)) / x;

                /* next QR transformation */
                c = s = 1.0;
                for (j = l; j <= nm; j++)
                {
                    int i  = j + 1;
                    g      = rv1[i];
                    y      = w[i];
                    h      = s * g;
                    g      = c * g;
                    z      = PYTHAG(f, h);
                    rv1[j] = z;
                    c      = f / z;
                    s      = h / z;
                    f      = x * c + g * s;
                    g      = g * c - x * s;
                    h      = y * s;
                    y      = y * c;
                    for (jj = 0; jj < n; jj++)
                    {
                        x        = v[jj][j];
                        z        = v[jj][i];
                        v[jj][j] = (x * c + z * s);
                        v[jj][i] = (z * c - x * s);
                    }
                    z    = PYTHAG(f, h);
                    w[j] = z;
                    if (z)
                    {
                        z = 1.0f / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = (c * g) + (s * y);
                    x = (c * y) - (s * g);
                    for (jj = 0; jj < m; jj++)
                    {
                        y        = a[jj][j];
                        z        = a[jj][i];
                        a[jj][j] = (y * c + z * s);
                        a[jj][i] = (z * c - y * s);
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k]   = x;
            }
        }
        return 1;
    }

    static void svbksb(Matrix &u, Vector &w, Matrix &v, Vector &b, Vector &x)
    {
        int m = u.rows();
        int n = u.cols();

        Vector tmp(n);
        for (int j = 0; j < n; j++)
        {
            float s = 0.0;
            if (w[j])
            {
                for (int i = 0; i < m; i++)
                    s += u[i][j] * b[i];
                s /= w[j];
            }
            tmp[j] = s;
        }
        for (int j = 0; j < n; j++)
        {
            float s = 0.0;
            for (int jj = 0; jj < n; jj++)
                s += v[j][jj] * tmp[jj];
            x[j] = s;
        }
    }

    // solve x for Ax=b.
    static int solve(Matrix A, Vector b, Vector &x, float TOL = 1e-5f)
    {
        // Matrix TMP(rows, cols);
        int    rows = A.rows();
        int    cols = A.cols();
        Matrix V(cols, cols);
        Vector W(cols);
        float  wmax, thresh;

        int r = svdcmp(A, W, V);
        if (r == 0)
        {
            return r;
        }

        wmax = 0.0;
        for (int j = 0; j < cols; j++)
        {
            if (W[j] > wmax)
                wmax = W[j];
        }

        thresh = TOL * wmax;
        for (int j = 0; j < cols; j++)
        {
            if (W[j] < thresh)
                W[j] = 0.0;
        }

        svbksb(A, W, V, b, x);

        return 0;
    }
} // namespace SOLVER

// --------------- Color correction.

class Corrector
{
  public:
    Corrector()          = default;
    virtual ~Corrector() = default;

    static Matrix3 solve(const std::vector<Tristimulus> &patch, const std::vector<Tristimulus> &target)
    {
        // M*A = B, M and B are known. solve A.
        int row_count = (int)(patch.size() + 1) * 3; // forex, (24+1)*3 = 75
        int col_count = 9;

        SOLVER::Matrix matrix(row_count, col_count); // M. (75*9)
        SOLVER::Vector result(row_count);            // B.  75
        SOLVER::Vector x(col_count);                 // A.     9

        for (int i = 0; i < patch.size() * 3; i += 3)
        {
            // each patch
            for (int c = 0; c < col_count; c++)
            {
                matrix[i + 0][c] = 0.f;
                matrix[i + 1][c] = 0.f;
                matrix[i + 2][c] = 0.f;
            }
            matrix[i + 0][3 * 0 + 0] = patch[i / 3][0];
            matrix[i + 0][3 * 0 + 1] = patch[i / 3][1];
            matrix[i + 0][3 * 0 + 2] = patch[i / 3][2];
            matrix[i + 1][3 * 1 + 0] = patch[i / 3][0];
            matrix[i + 1][3 * 1 + 1] = patch[i / 3][1];
            matrix[i + 1][3 * 1 + 2] = patch[i / 3][2];
            matrix[i + 2][3 * 2 + 0] = patch[i / 3][0];
            matrix[i + 2][3 * 2 + 1] = patch[i / 3][1];
            matrix[i + 2][3 * 2 + 2] = patch[i / 3][2];
            result[i + 0]            = target[(i / 3)][0];
            result[i + 1]            = target[(i / 3)][1];
            result[i + 2]            = target[(i / 3)][2];
        }

        // 0 restriction.
        {
            int i = (int)patch.size() * 3;
            for (int c = 0; c < col_count; c++)
            {
                matrix[i + 0][c] = 0.f;
                matrix[i + 1][c] = 0.f;
                matrix[i + 2][c] = 0.f;
            }
            matrix[i + 0][3 * 0 + 0] = 0.f;
            matrix[i + 0][3 * 0 + 1] = 0.f;
            matrix[i + 0][3 * 0 + 2] = 0.f;
            matrix[i + 1][3 * 1 + 0] = 0.f;
            matrix[i + 1][3 * 1 + 1] = 0.f;
            matrix[i + 1][3 * 1 + 2] = 0.f;
            matrix[i + 2][3 * 2 + 0] = 0.f;
            matrix[i + 2][3 * 2 + 1] = 0.f;
            matrix[i + 2][3 * 2 + 2] = 0.f;
            result[i + 0]            = 0.f;
            result[i + 1]            = 0.f;
            result[i + 2]            = 0.f;
        }

        SOLVER::solve(matrix, result, x);
        return Matrix3::fromArray(x.data());
    }

    static Matrix3 solve(
        std::vector<Tristimulus> &patch, const Spectrum &light = CIE_D65, const Observer &observer = CIE1931)
    {
        return solve(patch, Macbeth::reference(light, observer));
    }

    static Gamut makeGamut(const char *name, std::vector<Tristimulus> &patch, const Spectrum &light = CIE_D65,
        const Observer &observer = CIE1931)
    {
        return Gamut(name, solve(patch, light, observer));
    }
};

} // namespace ColorSystem

#endif /* colorsystem_hpp__abf47c16efbc4a80838738ff9b8a0eea */
