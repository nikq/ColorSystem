
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
//
// policy: currently single-headered, simple, low overhead, no ICC support.
//

#include <array>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#pragma once

namespace ColorSystem
{
class Vector3
{
  public:
    typedef std::array<float, 3> vec3;
    vec3 v_;
    constexpr Vector3(const float a, const float b, const float c) : v_({a, b, c}) { ; }
    constexpr float operator[](const int &i) const
    {
        return v_[i];
    }
    constexpr float x() const { return v_[0]; }
    constexpr float y() const { return v_[1]; }
    constexpr float z() const { return v_[2]; }
    constexpr auto  size() const { return v_.size(); }
    auto            begin() const { return v_.begin(); }
    auto            end() const { return v_.end(); }
};

class Matrix3
{
  public:
    typedef std::array<float, 9> matrix;
    matrix m_;

  private:
    static constexpr int M(const int x, const int y) { return x + y * 3; }
    static constexpr int I(const int y, const int x) { return ((x - 1) + (y - 1) * 3); }

  public:
    constexpr Matrix3(
        const float &a00, const float &a01, const float &a02,
        const float &a10, const float &a11, const float &a12,
        const float &a20, const float &a21, const float &a22) : m_({a00, a01, a02, a10, a11, a12, a20, a21, a22}) { ; }
    constexpr Matrix3(void) : m_({1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f}) { ; }
    static constexpr Matrix3 diag(const Vector3 &v)
    {
        return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]);
    }

    constexpr float operator[](const int &i) const
    {
        return m_[i];
    }

    static constexpr Matrix3 mul(const Matrix3 &a, const Matrix3 &b)
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

    static constexpr Vector3 apply(const Matrix3 &m, const Vector3 &v)
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

    static constexpr float det(const Matrix3 &m)
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

    static constexpr Matrix3 mul(const Matrix3 &m, const float &a)
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
    static constexpr Matrix3 div(const Matrix3 &m, const float &a)
    {
        return mul(m, 1.f / a);
    }
    constexpr Matrix3 div(const float a) const
    {
        return div(*this, a);
    }

    static constexpr Matrix3 add(const Matrix3 &a, const Matrix3 &b)
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

    static constexpr Matrix3 invert(const Matrix3 &a)
    {
        return mul(
            Matrix3(
                (a.m_[I(2, 2)] * a.m_[I(3, 3)] - a.m_[I(2, 3)] * a.m_[I(3, 2)]),
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
    constexpr Matrix3 invert(void) const
    {
        return invert(*this);
    }
    constexpr auto size() const { return m_.size(); }
    auto           begin() const { return m_.begin(); }
    auto           end() const { return m_.end(); }
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
    constexpr Tristimulus(const float &v) : v_(v, v, v) { ; }
    constexpr float operator[](const int &i) const
    {
        return v_[i];
    }
    constexpr auto               size() const { return v_.size(); }
    auto                         begin() const { return v_.begin(); }
    auto                         end() const { return v_.end(); }
    constexpr const Vector3 &    vec3(void) const { return v_; }
    constexpr float              a() const { return v_[0]; }
    constexpr float              b() const { return v_[1]; }
    constexpr float              c() const { return v_[2]; }
    static constexpr Tristimulus scale(const Tristimulus &t, const float &s)
    {
        return Tristimulus(t[0] * s, t[1] * s, t[2] * s);
    }
    constexpr Tristimulus scale(const float &s) const
    {
        return scale(*this, s);
    }
    constexpr Tristimulus operator*(const float &s) const
    {
        return scale(s);
    }
    constexpr Tristimulus operator/(const float &s) const
    {
        return scale(1.f / s);
    }

    // apply color transform matrix
    static constexpr Tristimulus apply(const Tristimulus &t, const Matrix3 &m)
    {
        return Tristimulus(m.apply(t.vec3()));
    }
    constexpr Tristimulus apply(const Matrix3 &m) const
    {
        return apply(*this, m);
    }

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
    constexpr const Tristimulus operator*(const Tristimulus &b) const
    {
        return mul(*this, b);
    }
    static constexpr float dot(const Tristimulus &a, const Tristimulus &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    constexpr float dot(const Tristimulus &b) const
    {
        return dot(*this, b);
    }
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

    constexpr float min(void) { return mini(mini(a(), b()), c()); }
    constexpr float max(void) { return maxi(maxi(a(), b()), c()); }

    constexpr Tristimulus clip(const float &l, const float &h) const
    {
        return max(min(*this, Tristimulus(h)), Tristimulus(l));
    }
    static constexpr Tristimulus clip(const Tristimulus &t, const float &l, const float &h) { return t.clip(l, h); }
    constexpr Tristimulus positive() const
    {
        return max(*this, Tristimulus(0.f));
    }
    static constexpr Tristimulus positive(const Tristimulus &t) { return t.positive(); }
    constexpr bool isNegative(const float &a) const
    {
        return (a < 0.f);
    }
    constexpr bool hasNegative() const
    {
        return isNegative(v_[0]) || isNegative(v_[1]) || isNegative(v_[2]);
    }
    static constexpr bool hasNegative(const Tristimulus &t)
    {
        return t.hasNegative();
    }
    static constexpr float z_from_xy(const float &x, const float &y) { return 1 - x - y; }
    static constexpr float X_from_Yxy(const float &Y, const float &x, const float &y) { return (Y < 1e-8f) ? 0.f : x * Y / y; }
    static constexpr float Y_from_Yxy(const float &Y, const float &x, const float &y) { return (Y < 1e-8f) ? 0.f : Y; }
    static constexpr float Z_from_Yxy(const float &Y, const float &x, const float &y) { return (Y < 1e-8f) ? 0.f : z_from_xy(x, y) * Y / y; }
    static constexpr float Y_from_XYZ(const float &x, const float &y, const float &z) { return (y < 1e-8f) ? 0.f : y; }
    static constexpr float x_from_XYZ(const float &x, const float &y, const float &z) { return (y < 1e-8f) ? 0.3127f : x / (x + y + z); }
    static constexpr float y_from_XYZ(const float &x, const float &y, const float &z) { return (y < 1e-8f) ? 0.3290f : y / (x + y + z); }
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
    static constexpr Tristimulus YxyToYuv(const Tristimulus &Yxy) { return Tristimulus(Yxy[0], u_from_xy(Yxy[1], Yxy[2]), v_from_xy(Yxy[1], Yxy[2])); }
    static constexpr Tristimulus YuvToYxy(const Tristimulus &Yuv) { return Tristimulus(Yuv[0], x_from_uv(Yuv[1], Yuv[2]), y_from_uv(Yuv[1], Yuv[2])); }
    static constexpr Tristimulus toYuv(const Tristimulus &XYZ)
    {
        return YxyToYuv(toYxy(XYZ));
    }
    static constexpr Tristimulus toYuv(const float &X, const float &Y, const float &Z)
    {
        return toYuv(Tristimulus(X, Y, Z));
    }
    static constexpr Tristimulus fromYuv(const Tristimulus &Yuv)
    {
        return fromYxy(YuvToYxy(Yuv));
    }
    static constexpr Tristimulus fromYuv(const float &X, const float &Y, const float &Z)
    {
        return fromYuv(Tristimulus(X, Y, Z));
    }
    constexpr Tristimulus toYuv(void) const { return toYuv(*this); }
    constexpr Tristimulus fromYuv(void) const { return fromYuv(*this); }
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
    //HSV
    static constexpr float mod360(const float &r)
    {
        return (r < 0.f) ? mod360(r + 360.f) : ((r > 360.f) ? mod360(r - 360.f) : r);
    }
    static constexpr Tristimulus toHSV(const Tristimulus &t)
    {
        const float max = maxi(maxi(t[0], t[1]), t[2]);
        const float min = mini(mini(t[0], t[1]), t[2]);
        return Tristimulus(
            mod360(((max == min) ? 0.f : ((max == t[0]) ? (60.f * (t[1] - t[2]) / (max - min)) : ((max == t[1]) ? (60.f * (t[2] - t[0]) / (max - min) + 120.f) : (60.f * (t[0] - t[1]) / (max - min) + 240.f))))),
            (max == 0.f) ? 0.f : (max - min) / max,
            max);
    }
    constexpr Tristimulus        toHSV(void) const { return toHSV(*this); }
    static constexpr Tristimulus fromHSV(const Tristimulus &t)
    {
        const float h  = mod360(t[0]);
        const int   d  = (int)(h / 60.f);
        const int   i  = d % 6;
        const float r  = (h / 60.f) - d;
        const float hi = t[2];
        const float lo = t[2] * (1.f - t[1]);
        const float m1 = lo + (hi - lo) * r;
        const float m2 = hi + (hi - lo) * r;
        return (i == 0) ? Tristimulus(hi, m1, lo) : ((i == 1) ? Tristimulus(m2, hi, lo) : ((i == 2) ? Tristimulus(lo, hi, m1) : ((i == 3) ? Tristimulus(lo, m2, hi) : ((i == 4) ? Tristimulus(m1, lo, hi) : ((i == 5) ? Tristimulus(hi, lo, m2) : Tristimulus(0, 0, 0))))));
    }
    constexpr Tristimulus fromHSV(void) const { return fromHSV(*this); }
};

class Delta
{
  public:
    static const float UV(const Tristimulus &a_Yuv, const Tristimulus &b_Yuv) // a, b both are XYZ
    {
        return sqrtf((a_Yuv[1] - b_Yuv[1]) * (a_Yuv[1] - b_Yuv[1]) + (a_Yuv[2] - b_Yuv[2]) * (a_Yuv[2] - b_Yuv[2]));
    }
    static const float E76(const Tristimulus &a_LAB, const Tristimulus &b_LAB)
    {
        return sqrtf(
            (a_LAB[0] - b_LAB[0]) * (a_LAB[0] - b_LAB[0]) +
            (a_LAB[1] - b_LAB[1]) * (a_LAB[1] - b_LAB[1]) +
            (a_LAB[2] - b_LAB[2]) * (a_LAB[2] - b_LAB[2]));
    }
    static const float E00(const Tristimulus &lab1, const Tristimulus &lab2, const float &Kl = 1.f, const float &Kc = 1.f, const float &Kh = 1.f)
    {
        const float PI         = 3.14159265358979323846264338327950288f;
        const float L1         = lab1[0];
        const float a1         = lab1[1];
        const float b1         = lab1[2];
        const float L2         = lab2[0];
        const float a2         = lab2[1];
        const float b2         = lab2[2];
        const float Lbar       = (L1 + L2) / 2.f;
        const float C1         = sqrtf(a1 * a1 + b1 * b1);
        const float C2         = sqrtf(a2 * a2 + b2 * b2);
        const float Cbar       = (C1 + C2) / 2.f;
        const float C7         = powf(Cbar, 7.f);
        const float pow25_7    = 25.f * 25.f * 25.f * 25.f * 25.f * 25.f * 25.f;
        const float G          = (1.f - sqrtf(C7 / (C7 + pow25_7))) / 2.f;
        const float ad1        = a1 * (1.f + G);
        const float ad2        = a2 * (1.f + G);
        const float Cd1        = sqrtf(ad1 * ad1 + b1 * b1);
        const float Cd2        = sqrtf(ad2 * ad2 + b2 * b2);
        const float CdBar      = (Cd1 + Cd2) / 2.f;
        const float h1         = fmod(360.f + atan2f(b1, ad1) * 180.0f / PI, 360.f);
        const float h2         = fmod(360.f + atan2f(b2, ad2) * 180.0f / PI, 360.f);
        const float HdBar      = (fabs(h1 - h2) > 180.f ? (h1 + h2 + 360.f) : (h1 + h2)) / 2.f;
        const float T          = 1.f - 0.17f * cosf(PI * (1.f * HdBar - 30.f) / 180.f) + 0.24f * cosf(PI * (2.f * HdBar) / 180.f) + 0.32f * cosf(PI * (3.f * HdBar + 6.f) / 180.f) - 0.20f * cosf(PI * (4.f * HdBar - 63.f) / 180.f);
        const float deltah     = (fabs(h2 - h1) <= 180.f) ? h2 - h1 : ((h2 <= h1) ? h2 - h1 + 360.f : h2 - h1 - 360.f);
        const float deltaL     = L2 - L1;
        const float deltaC     = Cd2 - Cd1;
        const float deltaH     = 2.f * sqrtf(Cd1 * Cd2) * sinf(PI * deltah / (180.f * 2.f));
        const float Lbar2      = (Lbar - 50.f) * (Lbar - 50.f);
        const float Sl         = 1.f + 0.015f * Lbar2 / sqrtf(20.f + Lbar2);
        const float Sc         = 1.f + 0.045f * CdBar;
        const float Sh         = 1.f + 0.015f * CdBar * T;
        const float HdBar2     = (HdBar - 275.f) * (HdBar - 275.f) / (25.f * 25.f);
        const float deltaTheta = 30.f * expf(-HdBar2);
        const float CdBar7     = powf(CdBar, 7.f);
        const float Rc         = 2.f * sqrtf(CdBar7 / (CdBar7 + pow25_7));
        const float Rt         = -Rc * sinf(2.f * deltaTheta * PI / 180.f);
        const float dl         = deltaL / (Kl * Sl);
        const float dc         = deltaC / (Kc * Sc);
        const float dh         = deltaH / (Kh * Sh);
        return sqrtf(dl * dl + dc * dc + dh * dh + Rt * dc * dh);
    }
};

class Gamut
{
  public:
    Matrix3 toXYZ_;
    Matrix3 fromXYZ_;

    static constexpr Matrix3 primMat(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB)
    {
        return Matrix3(xR, xG, xB, yR, yG, yB, Tristimulus::z_from_xy(xR, yR), Tristimulus::z_from_xy(xG, yG), Tristimulus::z_from_xy(xB, yB));
    }
    static constexpr Matrix3 diag(const Vector3 &v)
    {
        return Matrix3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]);
    }
    static constexpr Matrix3 mulDiag(const Matrix3 &m, const Vector3 &v)
    {
        return m.mul(diag(m.invert().apply(v)));
    }
    static constexpr Matrix3 fromPrimaries(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB, const float &xW, const float &yW)
    {
        return mulDiag(primMat(xR, yR, xG, yG, xB, yB), Tristimulus::fromYxy(1.f, xW, yW).vec3());
    }
    constexpr Gamut(const Matrix3 &fromXYZ) : toXYZ_(fromXYZ.invert()), fromXYZ_(fromXYZ)
    {
        ;
    }
    constexpr Gamut(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB, const float &xW, const float &yW)
        : toXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW)), fromXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW).invert())
    {
        ;
    }
    constexpr Matrix3 toXYZ(void) const { return toXYZ_; }
    constexpr Matrix3 fromXYZ(void) const { return fromXYZ_; }
    constexpr Tristimulus
    toXYZ(const Tristimulus &tri) const
    {
        return Tristimulus(toXYZ_.apply(tri.vec3()));
    }
    constexpr Tristimulus
    fromXYZ(const Tristimulus &tri) const
    {
        return Tristimulus(fromXYZ_.apply(tri.vec3()));
    }

    constexpr Vector3 primaryVector(void) const
    {
        return Vector3(
            toXYZ_[0] + toXYZ_[3] + toXYZ_[6],
            toXYZ_[1] + toXYZ_[4] + toXYZ_[7],
            toXYZ_[2] + toXYZ_[5] + toXYZ_[8]);
    }
    constexpr Matrix3 primaryMatrix(void) const
    {
        const Vector3 t(primaryVector());
        return Matrix3(
            toXYZ_[0] / t[0], toXYZ_[1] / t[1], toXYZ_[2] / t[2],
            toXYZ_[3] / t[0], toXYZ_[4] / t[1], toXYZ_[5] / t[2],
            toXYZ_[6] / t[0], toXYZ_[7] / t[1], toXYZ_[8] / t[2]);
    }
    constexpr Tristimulus primaryWhite() const
    {
        return Tristimulus(primaryMatrix().apply(primaryVector()));
    }
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
    typedef enum {
        LINEAR,
        GAMMA, // simplest gamma
        SRGB,
        BT709,
        ST2084,
        // OTF_HLG // Hybrid-log-gamma
    } TYPE;

    static float gamma(const float &v, const float &g) { return powf(v, 1.f / g); }
    static float degamma(const float &v, const float &g) { return powf(v, g); }
    static const float ST2084_to_Y(const float &pixel) // pixel should be 0-1
    {
        const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
        const float pq_m2 = 78.84375;        // ( 2523.0 / 4096.0 ) * 128.0;
        const float pq_c1 = 0.8359375;       // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
        const float pq_c2 = 18.8515625;      // ( 2413.0 / 4096.0 ) * 32.0;
        const float pq_c3 = 18.6875;         // ( 2392.0 / 4096.0 ) * 32.0;
        const float pq_C  = 10000.0;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this assumes full range (0 - 1)
        float Np = powf(pixel, 1.0f / pq_m2);
        float L  = Np - pq_c1;
        if (L < 0.0)
            L = 0.0;
        L     = L / (pq_c2 - pq_c3 * Np);
        L     = powf(L, 1.0f / pq_m1);
        return L * pq_C; // returns cd/m^2
    }
    static const float Y_to_ST2084(const float &nit) // nit should be 0-10000(cd/m^2)
    {
        const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
        const float pq_m2 = 78.84375;        // ( 2523.0 / 4096.0 ) * 128.0;
        const float pq_c1 = 0.8359375;       // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
        const float pq_c2 = 18.8515625;      // ( 2413.0 / 4096.0 ) * 32.0;
        const float pq_c3 = 18.6875;         // ( 2392.0 / 4096.0 ) * 32.0;
        const float pq_C  = 10000.0;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this returns full range (0 - 1)
        float L  = nit / pq_C;
        float Lm = powf(L, pq_m1);
        float N  = (pq_c1 + pq_c2 * Lm) / (1.0f + pq_c3 * Lm);
        N        = powf(N, pq_m2);
        return N;
    }
    static const float Y_to_sRGB(const float &nits) // returns signal, 0-1, input nits [0-100]
    {
        const float C = nits / 100.f;
        return (C < 0.0031308f) ? C * 12.92f : (1.055f * powf(C, 1.0f / 2.4f) - 0.055f);
    }
    static const float sRGB_to_Y(const float &C) // returns nits, 0-100[cd/m^2]
    {
        return (C < 0.04045f) ? C / 12.92f : powf((C + 0.055f) / 1.055f, 2.4f);
    }
    static const float Y_to_BT709(const float &nits) // returns signal, 0-1, input nits [0-100]
    {
        const float C = nits / 100.f;
        return (C < 0.018f) ? C * 4.50f : (1.099f * powf(C, 0.45f) - 0.099f);
    }
    static const float BT709_to_Y(const float &C) // returns nits, 0-100[cd/m^2]
    {
        return (C < 0.081f) ? C / 4.50f : powf((C + 0.099f) / 1.099f, 1.f / 0.45f);
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
        break;
        case LINEAR:
        default:
            return screen;
        }
    }
};

class Spectrum
{
  public:
    typedef std::array<float, 400> spectrum; // 380-780, 1nm, fixed.
    spectrum s_;

    constexpr Spectrum() : s_{} { ; }
    constexpr Spectrum(const spectrum &s) : s_(s) { ; }
    static constexpr float lerp(const float a, const float b, const float r)
    {
        return a * (1 - r) + b * r;
    }
    static constexpr float fetch(const float *sample, const int samples, const float &lo, const float &hi, const float &t)
    {
        const int   count = samples - 1;
        const int   index = (int)(count * (t - lo) / (hi - lo));
        const float l     = index * (hi - lo) / count;
        const float h     = (index + 1) * (hi - lo) / count;
        return (index >= count) ? sample[count] : (index <= 0) ? sample[0] : lerp(sample[index], sample[index + 1], (t - l - lo) / (h - l));
    }
    constexpr float fetch(const float &lambda) const
    {
        return fetch(s_.data(), 400, 380.f, 779.f, lambda);
    }

    Spectrum(
        const float *sample,
        const int    samples = 400,
        const float &lo      = 380.f,
        const float &hi      = 779.f)
    {
        for (int i = 0; i < 400; i++)
        {
            s_[i] = fetch(sample, samples, lo, hi, 380.f + i);
        }
    }
    constexpr float operator[](const int i) const
    {
        return s_[i];
    }

    const spectrum &s(void) const
    {
        return s_;
    }

    static const double planck(
        const double &T, // temperature (Kelvin)
        const double &l) // wavelength (meter)
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
    const Spectrum operator*(const Spectrum &b) const { return mul(*this, b); }
    const Spectrum operator+(const Spectrum &b) const { return add(*this, b); }
    static constexpr float sumHelper(const Spectrum &a, const int i)
    {
        return (i > 0) ? sumHelper(a, i - 1) + a[i] : a[0];
    }
    static constexpr float sum(const Spectrum &a)
    {
        return sumHelper(a, 399);
    }
    constexpr float sum() const
    {
        return sum(*this);
    }
    static constexpr float dotHelper(const Spectrum &a, const Spectrum &b, const int i)
    {
        return (i > 0) ? a[i] * b[i] + dotHelper(a, b, i - 1) : a[0] * b[0];
    }
    static constexpr float dot(const Spectrum &a, const Spectrum &b)
    {
        return dotHelper(a, b, 399);
    }
    constexpr float dot(const Spectrum &s) const
    {
        return dotHelper(*this, s, 399);
    }
};

class Observer
{
  public:
    Spectrum    X_, Y_, Z_;
    Tristimulus normalize_;
    constexpr Observer(const Spectrum &X, const Spectrum &Y, const Spectrum &Z) : X_(X), Y_(Y), Z_(Z), normalize_(1.f / X.sum(), 1.f / Y.sum(), 1.f / Z.sum()) { ; }
    static constexpr Tristimulus SpectrumIntegrate(const Spectrum &s, const Spectrum &x, const Spectrum &y, const Spectrum &z)
    {
        return Tristimulus(Spectrum::dot(s, x), Spectrum::dot(s, y), Spectrum::dot(s, z));
    }
    constexpr Tristimulus fromSpectrum(const Spectrum &s) const
    {
        return SpectrumIntegrate(s, X_, Y_, Z_) * normalize_;
    }
    constexpr float lumensFromMonochromaticFlux(const float &lambda, const float &watt) const // return lm
    {
        return Y_.fetch(lambda) * watt * 683.0f; // photoptic luminance efficiency
    }
    constexpr Tristimulus candellasFromMonochromaticRadiance(const float &lambda, const float &watt_per_sr_m2) const // return cd/m^2
    {
        return Tristimulus(X_.fetch(lambda), Y_.fetch(lambda), Z_.fetch(lambda)) * 683.0f * watt_per_sr_m2;
    }
};

static constexpr Spectrum CIE1931_X({0.00136800000f, 0.00150205000f, 0.00164232800f, 0.00180238200f, 0.00199575700f, 0.00223600000f, 0.00253538500f, 0.00289260300f, 0.00330082900f, 0.00375323600f,
    0.00424300000f, 0.00476238900f, 0.00533004800f, 0.00597871200f, 0.00674111700f, 0.00765000000f, 0.00875137300f, 0.01002888000f, 0.01142170000f, 0.01286901000f,
    0.01431000000f, 0.01570443000f, 0.01714744000f, 0.01878122000f, 0.02074801000f, 0.02319000000f, 0.02620736000f, 0.02978248000f, 0.03388092000f, 0.03846824000f,
    0.04351000000f, 0.04899560000f, 0.05502260000f, 0.06171880000f, 0.06921200000f, 0.07763000000f, 0.08695811000f, 0.09717672000f, 0.10840630000f, 0.12076720000f,
    0.13438000000f, 0.14935820000f, 0.16539570000f, 0.18198310000f, 0.19861100000f, 0.21477000000f, 0.23018680000f, 0.24487970000f, 0.25877730000f, 0.27180790000f,
    0.28390000000f, 0.29494380000f, 0.30489650000f, 0.31378730000f, 0.32164540000f, 0.32850000000f, 0.33435130000f, 0.33921010000f, 0.34312130000f, 0.34612960000f,
    0.34828000000f, 0.34959990000f, 0.35014740000f, 0.35001300000f, 0.34928700000f, 0.34806000000f, 0.34637330000f, 0.34426240000f, 0.34180880000f, 0.33909410000f,
    0.33620000000f, 0.33319770000f, 0.33004110000f, 0.32663570000f, 0.32288680000f, 0.31870000000f, 0.31402510000f, 0.30888400000f, 0.30329040000f, 0.29725790000f,
    0.29080000000f, 0.28397010000f, 0.27672140000f, 0.26891780000f, 0.26042270000f, 0.25110000000f, 0.24084750000f, 0.22985120000f, 0.21840720000f, 0.20681150000f,
    0.19536000000f, 0.18421360000f, 0.17332730000f, 0.16268810000f, 0.15228330000f, 0.14210000000f, 0.13217860000f, 0.12256960000f, 0.11327520000f, 0.10429790000f,
    0.09564000000f, 0.08729955000f, 0.07930804000f, 0.07171776000f, 0.06458099000f, 0.05795001000f, 0.05186211000f, 0.04628152000f, 0.04115088000f, 0.03641283000f,
    0.03201000000f, 0.02791720000f, 0.02414440000f, 0.02068700000f, 0.01754040000f, 0.01470000000f, 0.01216179000f, 0.00991996000f, 0.00796724000f, 0.00629634600f,
    0.00490000000f, 0.00377717300f, 0.00294532000f, 0.00242488000f, 0.00223629300f, 0.00240000000f, 0.00292552000f, 0.00383656000f, 0.00517484000f, 0.00698208000f,
    0.00930000000f, 0.01214949000f, 0.01553588000f, 0.01947752000f, 0.02399277000f, 0.02910000000f, 0.03481485000f, 0.04112016000f, 0.04798504000f, 0.05537861000f,
    0.06327000000f, 0.07163501000f, 0.08046224000f, 0.08973996000f, 0.09945645000f, 0.10960000000f, 0.12016740000f, 0.13111450000f, 0.14236790000f, 0.15385420000f,
    0.16550000000f, 0.17725710000f, 0.18914000000f, 0.20116940000f, 0.21336580000f, 0.22574990000f, 0.23832090000f, 0.25106680000f, 0.26399220000f, 0.27710170000f,
    0.29040000000f, 0.30389120000f, 0.31757260000f, 0.33143840000f, 0.34548280000f, 0.35970000000f, 0.37408390000f, 0.38863960000f, 0.40337840000f, 0.41831150000f,
    0.43344990000f, 0.44879530000f, 0.46433600000f, 0.48006400000f, 0.49597130000f, 0.51205010000f, 0.52829590000f, 0.54469160000f, 0.56120940000f, 0.57782150000f,
    0.59450000000f, 0.61122090000f, 0.62797580000f, 0.64476020000f, 0.66156970000f, 0.67840000000f, 0.69523920000f, 0.71205860000f, 0.72882840000f, 0.74551880000f,
    0.76210000000f, 0.77854320000f, 0.79482560000f, 0.81092640000f, 0.82682480000f, 0.84250000000f, 0.85793250000f, 0.87308160000f, 0.88789440000f, 0.90231810000f,
    0.91630000000f, 0.92979950000f, 0.94279840000f, 0.95527760000f, 0.96721790000f, 0.97860000000f, 0.98938560000f, 0.99954880000f, 1.00908920000f, 1.01800640000f,
    1.02630000000f, 1.03398270000f, 1.04098600000f, 1.04718800000f, 1.05246670000f, 1.05670000000f, 1.05979440000f, 1.06179920000f, 1.06280680000f, 1.06290960000f,
    1.06220000000f, 1.06073520000f, 1.05844360000f, 1.05522440000f, 1.05097680000f, 1.04560000000f, 1.03903690000f, 1.03136080000f, 1.02266620000f, 1.01304770000f,
    1.00260000000f, 0.99136750000f, 0.97933140000f, 0.96649160000f, 0.95284790000f, 0.93840000000f, 0.92319400000f, 0.90724400000f, 0.89050200000f, 0.87292000000f,
    0.85444990000f, 0.83508400000f, 0.81494600000f, 0.79418600000f, 0.77295400000f, 0.75140000000f, 0.72958360000f, 0.70758880000f, 0.68560220000f, 0.66381040000f,
    0.64240000000f, 0.62151490000f, 0.60111380000f, 0.58110520000f, 0.56139770000f, 0.54190000000f, 0.52259950000f, 0.50354640000f, 0.48474360000f, 0.46619390000f,
    0.44790000000f, 0.42986130000f, 0.41209800000f, 0.39464400000f, 0.37753330000f, 0.36080000000f, 0.34445630000f, 0.32851680000f, 0.31301920000f, 0.29800110000f,
    0.28350000000f, 0.26954480000f, 0.25611840000f, 0.24318960000f, 0.23072720000f, 0.21870000000f, 0.20709710000f, 0.19592320000f, 0.18517080000f, 0.17483230000f,
    0.16490000000f, 0.15536670000f, 0.14623000000f, 0.13749000000f, 0.12914670000f, 0.12120000000f, 0.11363970000f, 0.10646500000f, 0.09969044000f, 0.09333061000f,
    0.08740000000f, 0.08190096000f, 0.07680428000f, 0.07207712000f, 0.06768664000f, 0.06360000000f, 0.05980685000f, 0.05628216000f, 0.05297104000f, 0.04981861000f,
    0.04677000000f, 0.04378405000f, 0.04087536000f, 0.03807264000f, 0.03540461000f, 0.03290000000f, 0.03056419000f, 0.02838056000f, 0.02634484000f, 0.02445275000f,
    0.02270000000f, 0.02108429000f, 0.01959988000f, 0.01823732000f, 0.01698717000f, 0.01584000000f, 0.01479064000f, 0.01383132000f, 0.01294868000f, 0.01212920000f,
    0.01135916000f, 0.01062935000f, 0.00993884600f, 0.00928842200f, 0.00867885400f, 0.00811091600f, 0.00758238800f, 0.00708874600f, 0.00662731300f, 0.00619540800f,
    0.00579034600f, 0.00540982600f, 0.00505258300f, 0.00471751200f, 0.00440350700f, 0.00410945700f, 0.00383391300f, 0.00357574800f, 0.00333434200f, 0.00310907500f,
    0.00289932700f, 0.00270434800f, 0.00252302000f, 0.00235416800f, 0.00219661600f, 0.00204919000f, 0.00191096000f, 0.00178143800f, 0.00166011000f, 0.00154645900f,
    0.00143997100f, 0.00134004200f, 0.00124627500f, 0.00115847100f, 0.00107643000f, 0.00099994930f, 0.00092873580f, 0.00086243320f, 0.00080075030f, 0.00074339600f,
    0.00069007860f, 0.00064051560f, 0.00059450210f, 0.00055186460f, 0.00051242900f, 0.00047602130f, 0.00044245360f, 0.00041151170f, 0.00038298140f, 0.00035664910f,
    0.00033230110f, 0.00030975860f, 0.00028888710f, 0.00026953940f, 0.00025156820f, 0.00023482610f, 0.00021917100f, 0.00020452580f, 0.00019084050f, 0.00017806540f,
    0.00016615050f, 0.00015502360f, 0.00014462190f, 0.00013490980f, 0.00012585200f, 0.00011741300f, 0.00010955150f, 0.00010222450f, 0.00009539445f, 0.00008902390f,
    0.00008307527f, 0.00007751269f, 0.00007231304f, 0.00006745778f, 0.00006292844f, 0.00005870652f, 0.00005477028f, 0.00005109918f, 0.00004767654f, 0.00004448567f});
static constexpr Spectrum CIE1931_Y({0.00003900000f, 0.00004282640f, 0.00004691460f, 0.00005158960f, 0.00005717640f, 0.00006400000f, 0.00007234421f, 0.00008221224f, 0.00009350816f, 0.00010613610f,
    0.00012000000f, 0.00013498400f, 0.00015149200f, 0.00017020800f, 0.00019181600f, 0.00021700000f, 0.00024690670f, 0.00028124000f, 0.00031852000f, 0.00035726670f,
    0.00039600000f, 0.00043371470f, 0.00047302400f, 0.00051787600f, 0.00057221870f, 0.00064000000f, 0.00072456000f, 0.00082550000f, 0.00094116000f, 0.00106988000f,
    0.00121000000f, 0.00136209100f, 0.00153075200f, 0.00172036800f, 0.00193532300f, 0.00218000000f, 0.00245480000f, 0.00276400000f, 0.00311780000f, 0.00352640000f,
    0.00400000000f, 0.00454624000f, 0.00515932000f, 0.00582928000f, 0.00654616000f, 0.00730000000f, 0.00808650700f, 0.00890872000f, 0.00976768000f, 0.01066443000f,
    0.01160000000f, 0.01257317000f, 0.01358272000f, 0.01462968000f, 0.01571509000f, 0.01684000000f, 0.01800736000f, 0.01921448000f, 0.02045392000f, 0.02171824000f,
    0.02300000000f, 0.02429461000f, 0.02561024000f, 0.02695857000f, 0.02835125000f, 0.02980000000f, 0.03131083000f, 0.03288368000f, 0.03452112000f, 0.03622571000f,
    0.03800000000f, 0.03984667000f, 0.04176800000f, 0.04376600000f, 0.04584267000f, 0.04800000000f, 0.05024368000f, 0.05257304000f, 0.05498056000f, 0.05745872000f,
    0.06000000000f, 0.06260197000f, 0.06527752000f, 0.06804208000f, 0.07091109000f, 0.07390000000f, 0.07701600000f, 0.08026640000f, 0.08366680000f, 0.08723280000f,
    0.09098000000f, 0.09491755000f, 0.09904584000f, 0.10336740000f, 0.10788460000f, 0.11260000000f, 0.11753200000f, 0.12267440000f, 0.12799280000f, 0.13345280000f,
    0.13902000000f, 0.14467640000f, 0.15046930000f, 0.15646190000f, 0.16271770000f, 0.16930000000f, 0.17624310000f, 0.18355810000f, 0.19127350000f, 0.19941800000f,
    0.20802000000f, 0.21711990000f, 0.22673450000f, 0.23685710000f, 0.24748120000f, 0.25860000000f, 0.27018490000f, 0.28229390000f, 0.29505050000f, 0.30857800000f,
    0.32300000000f, 0.33840210000f, 0.35468580000f, 0.37169860000f, 0.38928750000f, 0.40730000000f, 0.42562990000f, 0.44430960000f, 0.46339440000f, 0.48293950000f,
    0.50300000000f, 0.52356930000f, 0.54451200000f, 0.56569000000f, 0.58696530000f, 0.60820000000f, 0.62934560000f, 0.65030680000f, 0.67087520000f, 0.69084240000f,
    0.71000000000f, 0.72818520000f, 0.74546360000f, 0.76196940000f, 0.77783680000f, 0.79320000000f, 0.80811040000f, 0.82249620000f, 0.83630680000f, 0.84949160000f,
    0.86200000000f, 0.87381080000f, 0.88496240000f, 0.89549360000f, 0.90544320000f, 0.91485010000f, 0.92373480000f, 0.93209240000f, 0.93992260000f, 0.94722520000f,
    0.95400000000f, 0.96025610000f, 0.96600740000f, 0.97126060000f, 0.97602250000f, 0.98030000000f, 0.98409240000f, 0.98741820000f, 0.99031280000f, 0.99281160000f,
    0.99495010000f, 0.99671080000f, 0.99809830000f, 0.99911200000f, 0.99974820000f, 1.00000000000f, 0.99985670000f, 0.99930460000f, 0.99832550000f, 0.99689870000f,
    0.99500000000f, 0.99260050000f, 0.98974260000f, 0.98644440000f, 0.98272410000f, 0.97860000000f, 0.97408370000f, 0.96917120000f, 0.96385680000f, 0.95813490000f,
    0.95200000000f, 0.94545040000f, 0.93849920000f, 0.93116280000f, 0.92345760000f, 0.91540000000f, 0.90700640000f, 0.89827720000f, 0.88920480000f, 0.87978160000f,
    0.87000000000f, 0.85986130000f, 0.84939200000f, 0.83862200000f, 0.82758130000f, 0.81630000000f, 0.80479470000f, 0.79308200000f, 0.78119200000f, 0.76915470000f,
    0.75700000000f, 0.74475410000f, 0.73242240000f, 0.72000360000f, 0.70749650000f, 0.69490000000f, 0.68221920000f, 0.66947160000f, 0.65667440000f, 0.64384480000f,
    0.63100000000f, 0.61815550000f, 0.60531440000f, 0.59247560000f, 0.57963790000f, 0.56680000000f, 0.55396110000f, 0.54113720000f, 0.52835280000f, 0.51563230000f,
    0.50300000000f, 0.49046880000f, 0.47803040000f, 0.46567760000f, 0.45340320000f, 0.44120000000f, 0.42908000000f, 0.41703600000f, 0.40503200000f, 0.39303200000f,
    0.38100000000f, 0.36891840000f, 0.35682720000f, 0.34477680000f, 0.33281760000f, 0.32100000000f, 0.30933810000f, 0.29785040000f, 0.28659360000f, 0.27562450000f,
    0.26500000000f, 0.25476320000f, 0.24488960000f, 0.23533440000f, 0.22605280000f, 0.21700000000f, 0.20816160000f, 0.19954880000f, 0.19115520000f, 0.18297440000f,
    0.17500000000f, 0.16722350000f, 0.15964640000f, 0.15227760000f, 0.14512590000f, 0.13820000000f, 0.13150030000f, 0.12502480000f, 0.11877920000f, 0.11276910000f,
    0.10700000000f, 0.10147620000f, 0.09618864000f, 0.09112296000f, 0.08626485000f, 0.08160000000f, 0.07712064000f, 0.07282552000f, 0.06871008000f, 0.06476976000f,
    0.06100000000f, 0.05739621000f, 0.05395504000f, 0.05067376000f, 0.04754965000f, 0.04458000000f, 0.04175872000f, 0.03908496000f, 0.03656384000f, 0.03420048000f,
    0.03200000000f, 0.02996261000f, 0.02807664000f, 0.02632936000f, 0.02470805000f, 0.02320000000f, 0.02180077000f, 0.02050112000f, 0.01928108000f, 0.01812069000f,
    0.01700000000f, 0.01590379000f, 0.01483718000f, 0.01381068000f, 0.01283478000f, 0.01192000000f, 0.01106831000f, 0.01027339000f, 0.00953331100f, 0.00884615700f,
    0.00821000000f, 0.00762378100f, 0.00708542400f, 0.00659147600f, 0.00613848500f, 0.00572300000f, 0.00534305900f, 0.00499579600f, 0.00467640400f, 0.00438007500f,
    0.00410200000f, 0.00383845300f, 0.00358909900f, 0.00335421900f, 0.00313409300f, 0.00292900000f, 0.00273813900f, 0.00255987600f, 0.00239324400f, 0.00223727500f,
    0.00209100000f, 0.00195358700f, 0.00182458000f, 0.00170358000f, 0.00159018700f, 0.00148400000f, 0.00138449600f, 0.00129126800f, 0.00120409200f, 0.00112274400f,
    0.00104700000f, 0.00097658960f, 0.00091110880f, 0.00085013320f, 0.00079323840f, 0.00074000000f, 0.00069008270f, 0.00064331000f, 0.00059949600f, 0.00055845470f,
    0.00052000000f, 0.00048391360f, 0.00045005280f, 0.00041834520f, 0.00038871840f, 0.00036110000f, 0.00033538350f, 0.00031144040f, 0.00028916560f, 0.00026845390f,
    0.00024920000f, 0.00023130190f, 0.00021468560f, 0.00019928840f, 0.00018504750f, 0.00017190000f, 0.00015977810f, 0.00014860440f, 0.00013830160f, 0.00012879250f,
    0.00012000000f, 0.00011185950f, 0.00010432240f, 0.00009733560f, 0.00009084587f, 0.00008480000f, 0.00007914667f, 0.00007385800f, 0.00006891600f, 0.00006430267f,
    0.00006000000f, 0.00005598187f, 0.00005222560f, 0.00004871840f, 0.00004544747f, 0.00004240000f, 0.00003956104f, 0.00003691512f, 0.00003444868f, 0.00003214816f,
    0.00003000000f, 0.00002799125f, 0.00002611356f, 0.00002436024f, 0.00002272461f, 0.00002120000f, 0.00001977855f, 0.00001845285f, 0.00001721687f, 0.00001606459f});
static constexpr Spectrum CIE1931_Z({0.00645000100f, 0.00708321600f, 0.00774548800f, 0.00850115200f, 0.00941454400f, 0.01054999000f, 0.01196580000f, 0.01365587000f, 0.01558805000f, 0.01773015000f,
    0.02005001000f, 0.02251136000f, 0.02520288000f, 0.02827972000f, 0.03189704000f, 0.03621000000f, 0.04143771000f, 0.04750372000f, 0.05411988000f, 0.06099803000f,
    0.06785001000f, 0.07448632000f, 0.08136156000f, 0.08915364000f, 0.09854048000f, 0.11020000000f, 0.12461330000f, 0.14170170000f, 0.16130350000f, 0.18325680000f,
    0.20740000000f, 0.23369210000f, 0.26261140000f, 0.29477460000f, 0.33079850000f, 0.37130000000f, 0.41620910000f, 0.46546420000f, 0.51969480000f, 0.57953030000f,
    0.64560000000f, 0.71848380000f, 0.79671330000f, 0.87784590000f, 0.95943900000f, 1.03905010000f, 1.11536730000f, 1.18849710000f, 1.25812330000f, 1.32392960000f,
    1.38560000000f, 1.44263520000f, 1.49480350000f, 1.54219030000f, 1.58488070000f, 1.62296000000f, 1.65640480000f, 1.68529590000f, 1.70987450000f, 1.73038210000f,
    1.74706000000f, 1.76004460000f, 1.76962330000f, 1.77626370000f, 1.78043340000f, 1.78260000000f, 1.78296820000f, 1.78169980000f, 1.77919820000f, 1.77586710000f,
    1.77211000000f, 1.76825890000f, 1.76403900000f, 1.75894380000f, 1.75246630000f, 1.74410000000f, 1.73355950000f, 1.72085810000f, 1.70593690000f, 1.68873720000f,
    1.66920000000f, 1.64752870000f, 1.62341270000f, 1.59602230000f, 1.56452800000f, 1.52810000000f, 1.48611140000f, 1.43952150000f, 1.38987990000f, 1.33873620000f,
    1.28764000000f, 1.23742230000f, 1.18782430000f, 1.13876110000f, 1.09014800000f, 1.04190000000f, 0.99419760000f, 0.94734730000f, 0.90145310000f, 0.85661930000f,
    0.81295010000f, 0.77051730000f, 0.72944480000f, 0.68991360000f, 0.65210490000f, 0.61620000000f, 0.58232860000f, 0.55041620000f, 0.52033760000f, 0.49196730000f,
    0.46518000000f, 0.43992460000f, 0.41618360000f, 0.39388220000f, 0.37294590000f, 0.35330000000f, 0.33485780000f, 0.31755210000f, 0.30133750000f, 0.28616860000f,
    0.27200000000f, 0.25881710000f, 0.24648380000f, 0.23477180000f, 0.22345330000f, 0.21230000000f, 0.20116920000f, 0.19011960000f, 0.17922540000f, 0.16856080000f,
    0.15820000000f, 0.14813830000f, 0.13837580000f, 0.12899420000f, 0.12007510000f, 0.11170000000f, 0.10390480000f, 0.09666748000f, 0.08998272000f, 0.08384531000f,
    0.07824999000f, 0.07320899000f, 0.06867816000f, 0.06456784000f, 0.06078835000f, 0.05725001000f, 0.05390435000f, 0.05074664000f, 0.04775276000f, 0.04489859000f,
    0.04216000000f, 0.03950728000f, 0.03693564000f, 0.03445836000f, 0.03208872000f, 0.02984000000f, 0.02771181000f, 0.02569444000f, 0.02378716000f, 0.02198925000f,
    0.02030000000f, 0.01871805000f, 0.01724036000f, 0.01586364000f, 0.01458461000f, 0.01340000000f, 0.01230723000f, 0.01130188000f, 0.01037792000f, 0.00952930600f,
    0.00874999900f, 0.00803520000f, 0.00738160000f, 0.00678540000f, 0.00624280000f, 0.00574999900f, 0.00530360000f, 0.00489980000f, 0.00453420000f, 0.00420240000f,
    0.00390000000f, 0.00362320000f, 0.00337060000f, 0.00314140000f, 0.00293480000f, 0.00274999900f, 0.00258520000f, 0.00243860000f, 0.00230940000f, 0.00219680000f,
    0.00210000000f, 0.00201773300f, 0.00194820000f, 0.00188980000f, 0.00184093300f, 0.00180000000f, 0.00176626700f, 0.00173780000f, 0.00171120000f, 0.00168306700f,
    0.00165000100f, 0.00161013300f, 0.00156440000f, 0.00151360000f, 0.00145853300f, 0.00140000000f, 0.00133666700f, 0.00127000000f, 0.00120500000f, 0.00114666700f,
    0.00110000000f, 0.00106880000f, 0.00104940000f, 0.00103560000f, 0.00102120000f, 0.00100000000f, 0.00096864000f, 0.00092992000f, 0.00088688000f, 0.00084256000f,
    0.00080000000f, 0.00076096000f, 0.00072368000f, 0.00068592000f, 0.00064544000f, 0.00060000000f, 0.00054786670f, 0.00049160000f, 0.00043540000f, 0.00038346670f,
    0.00034000000f, 0.00030725330f, 0.00028316000f, 0.00026544000f, 0.00025181330f, 0.00024000000f, 0.00022954670f, 0.00022064000f, 0.00021196000f, 0.00020218670f,
    0.00019000000f, 0.00017421330f, 0.00015564000f, 0.00013596000f, 0.00011685330f, 0.00010000000f, 0.00008613333f, 0.00007460000f, 0.00006500000f, 0.00005693333f,
    0.00004999999f, 0.00004416000f, 0.00003948000f, 0.00003572000f, 0.00003264000f, 0.00003000000f, 0.00002765333f, 0.00002556000f, 0.00002364000f, 0.00002181333f,
    0.00002000000f, 0.00001813333f, 0.00001620000f, 0.00001420000f, 0.00001213333f, 0.00001000000f, 0.00000773333f, 0.00000540000f, 0.00000320000f, 0.00000133333f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f,
    0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f, 0.00000000000f});

static constexpr Spectrum CIE2012_X({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0037696470f, 0.0045324160f, 0.0054465530f, 0.0065388680f, 0.0078396990f, 0.0093829670f, 0.0112060800f, 0.0133496500f, 0.0158569000f, 0.0187728600f,
    0.0221430200f, 0.0260128500f, 0.0304303600f, 0.0354432500f, 0.0410964000f, 0.0474298600f, 0.0544739400f, 0.0622361200f, 0.0707004800f, 0.0798251300f,
    0.0895380300f, 0.0997484800f, 0.1104019000f, 0.1214566000f, 0.1328741000f, 0.1446214000f, 0.1566468000f, 0.1687901000f, 0.1808328000f, 0.1925216000f,
    0.2035729000f, 0.2137531000f, 0.2231348000f, 0.2319245000f, 0.2403892000f, 0.2488523000f, 0.2575896000f, 0.2664991000f, 0.2753532000f, 0.2838921000f,
    0.2918246000f, 0.2989200000f, 0.3052993000f, 0.3112031000f, 0.3169047000f, 0.3227087000f, 0.3288194000f, 0.3349242000f, 0.3405452000f, 0.3451688000f,
    0.3482554000f, 0.3494153000f, 0.3489075000f, 0.3471746000f, 0.3446705000f, 0.3418483000f, 0.3390240000f, 0.3359926000f, 0.3324276000f, 0.3280157000f,
    0.3224637000f, 0.3156225000f, 0.3078201000f, 0.2994771000f, 0.2909776000f, 0.2826646000f, 0.2747962000f, 0.2674312000f, 0.2605847000f, 0.2542749000f,
    0.2485254000f, 0.2433039000f, 0.2383414000f, 0.2333253000f, 0.2279619000f, 0.2219781000f, 0.2151735000f, 0.2075619000f, 0.1992183000f, 0.1902290000f,
    0.1806905000f, 0.1707154000f, 0.1604471000f, 0.1500244000f, 0.1395705000f, 0.1291920000f, 0.1189859000f, 0.1090615000f, 0.0995142400f, 0.0904185000f,
    0.0818289500f, 0.0737681700f, 0.0661947700f, 0.0590638000f, 0.0523424200f, 0.0460086500f, 0.0400615400f, 0.0345437300f, 0.0294909100f, 0.0249214000f,
    0.0208398100f, 0.0172359100f, 0.0140792400f, 0.0113451600f, 0.0090196580f, 0.0070977310f, 0.0055711450f, 0.0043945660f, 0.0035163030f, 0.0028876380f,
    0.0024615880f, 0.0022063480f, 0.0021495590f, 0.0023370910f, 0.0028189310f, 0.0036491780f, 0.0048913590f, 0.0066293640f, 0.0089429020f, 0.0119022400f,
    0.0155698900f, 0.0199766800f, 0.0250469800f, 0.0306753000f, 0.0367499900f, 0.0431517100f, 0.0497858400f, 0.0566855400f, 0.0639165100f, 0.0715435200f,
    0.0796291700f, 0.0882147300f, 0.0972697800f, 0.1067504000f, 0.1166192000f, 0.1268468000f, 0.1374060000f, 0.1482471000f, 0.1593076000f, 0.1705181000f,
    0.1818026000f, 0.1931090000f, 0.2045085000f, 0.2161166000f, 0.2280650000f, 0.2405015000f, 0.2535441000f, 0.2671300000f, 0.2811351000f, 0.2954164000f,
    0.3098117000f, 0.3241678000f, 0.3384319000f, 0.3525786000f, 0.3665839000f, 0.3804244000f, 0.3940988000f, 0.4076972000f, 0.4213484000f, 0.4352003000f,
    0.4494206000f, 0.4641616000f, 0.4794395000f, 0.4952180000f, 0.5114395000f, 0.5280233000f, 0.5448696000f, 0.5618898000f, 0.5790137000f, 0.5961882000f,
    0.6133784000f, 0.6305897000f, 0.6479223000f, 0.6654866000f, 0.6833782000f, 0.7016774000f, 0.7204110000f, 0.7394495000f, 0.7586285000f, 0.7777885000f,
    0.7967750000f, 0.8154530000f, 0.8337389000f, 0.8515493000f, 0.8687862000f, 0.8853376000f, 0.9011588000f, 0.9165278000f, 0.9318245000f, 0.9474524000f,
    0.9638388000f, 0.9812596000f, 0.9992953000f, 1.0173430000f, 1.0347900000f, 1.0510110000f, 1.0655220000f, 1.0784210000f, 1.0899440000f, 1.1003200000f,
    1.1097670000f, 1.1184380000f, 1.1262660000f, 1.1331380000f, 1.1389520000f, 1.1436200000f, 1.1470950000f, 1.1494640000f, 1.1508380000f, 1.1513260000f,
    1.1510330000f, 1.1500020000f, 1.1480610000f, 1.1449980000f, 1.1406220000f, 1.1347570000f, 1.1272980000f, 1.1183420000f, 1.1080330000f, 1.0965150000f,
    1.0839280000f, 1.0703870000f, 1.0559340000f, 1.0405920000f, 1.0243850000f, 1.0073440000f, 0.9895268000f, 0.9711213000f, 0.9523257000f, 0.9333248000f,
    0.9142877000f, 0.8952798000f, 0.8760157000f, 0.8561607000f, 0.8354235000f, 0.8135565000f, 0.7904565000f, 0.7664364000f, 0.7418777000f, 0.7171219000f,
    0.6924717000f, 0.6681600000f, 0.6442697000f, 0.6208450000f, 0.5979243000f, 0.5755410000f, 0.5537296000f, 0.5325412000f, 0.5120218000f, 0.4922070000f,
    0.4731224000f, 0.4547417000f, 0.4368719000f, 0.4193121000f, 0.4018980000f, 0.3844986000f, 0.3670592000f, 0.3497167000f, 0.3326305000f, 0.3159341000f,
    0.2997374000f, 0.2841189000f, 0.2691053000f, 0.2547077000f, 0.2409319000f, 0.2277792000f, 0.2152431000f, 0.2033010000f, 0.1919276000f, 0.1810987000f,
    0.1707914000f, 0.1609842000f, 0.1516577000f, 0.1427936000f, 0.1343737000f, 0.1263808000f, 0.1187979000f, 0.1116088000f, 0.1047975000f, 0.0983483500f,
    0.0922459700f, 0.0864750600f, 0.0810198600f, 0.0758651400f, 0.0709963300f, 0.0663996000f, 0.0620622500f, 0.0579740900f, 0.0541253300f, 0.0505060000f,
    0.0471060600f, 0.0439141100f, 0.0409141100f, 0.0380906700f, 0.0354303400f, 0.0329213800f, 0.0305567200f, 0.0283414600f, 0.0262803300f, 0.0243746500f,
    0.0226230600f, 0.0210193500f, 0.0195464700f, 0.0181872700f, 0.0169272700f, 0.0157541700f, 0.0146585400f, 0.0136357100f, 0.0126820500f, 0.0117939400f,
    0.0109677800f, 0.0101996400f, 0.0094843170f, 0.0088168510f, 0.0081929210f, 0.0076087500f, 0.0070613910f, 0.0065495090f, 0.0060719700f, 0.0056274760f,
    0.0052146080f, 0.0048318480f, 0.0044775790f, 0.0041501660f, 0.0038479880f, 0.0035694520f, 0.0033128570f, 0.0030760220f, 0.0028568940f, 0.0026536810f,
    0.0024648210f, 0.0022890600f, 0.0021256940f, 0.0019741210f, 0.0018337230f, 0.0017038760f, 0.0015839040f, 0.0014729390f, 0.0013701510f, 0.0012748030f,
    0.0011862380f, 0.0011038710f, 0.0010271940f, 0.0009557493f, 0.0008891262f, 0.0008269535f, 0.0007689351f, 0.0007149425f, 0.0006648590f, 0.0006185421f,
    0.0005758303f, 0.0005365046f, 0.0005001842f, 0.0004665005f, 0.0004351386f, 0.0004058303f, 0.0003783733f, 0.0003526892f, 0.0003287199f, 0.0003063998f,
    0.0002856577f, 0.0002664108f, 0.0002485462f, 0.0002319529f, 0.0002165300f, 0.0002021853f, 0.0001888338f, 0.0001763935f, 0.0001647895f, 0.0001539542f,
    0.0001438270f, 0.0001343572f, 0.0001255141f, 0.0001172706f, 0.0001095983f, 0.0001024685f, 0.0000958472f, 0.0000896832f, 0.0000839273f, 0.0000785371f,
    0.0000734755f, 0.0000687158f, 0.0000642526f, 0.0000600829f, 0.0000562010f, 0.0000525987f, 0.0000492628f, 0.0000461662f, 0.0000432821f, 0.0000405872f});

static constexpr Spectrum CIE2012_Y({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0004146161f, 0.0005028333f, 0.0006084991f, 0.0007344436f, 0.0008837389f, 0.0010596460f, 0.0012655320f, 0.0015047530f, 0.0017804930f, 0.0020955720f,
    0.0024521940f, 0.0028522160f, 0.0032991150f, 0.0037974660f, 0.0043527680f, 0.0049717170f, 0.0056610140f, 0.0064216150f, 0.0072503120f, 0.0081401730f,
    0.0090798600f, 0.0100560800f, 0.0110645600f, 0.0121052200f, 0.0131801400f, 0.0142937700f, 0.0154500400f, 0.0166409300f, 0.0178530200f, 0.0190701800f,
    0.0202736900f, 0.0214480500f, 0.0226004100f, 0.0237478900f, 0.0249124700f, 0.0261210600f, 0.0273992300f, 0.0287499300f, 0.0301690900f, 0.0316514500f,
    0.0331903800f, 0.0347791200f, 0.0364149500f, 0.0380956900f, 0.0398184300f, 0.0415794000f, 0.0433709800f, 0.0451718000f, 0.0469542000f, 0.0486871800f,
    0.0503365700f, 0.0518761100f, 0.0533221800f, 0.0547060300f, 0.0560633500f, 0.0574339300f, 0.0588510700f, 0.0603080900f, 0.0617864400f, 0.0632657000f,
    0.0647235200f, 0.0661474900f, 0.0675725600f, 0.0690492800f, 0.0706328000f, 0.0723833900f, 0.0743596000f, 0.0765938300f, 0.0791143600f, 0.0819534500f,
    0.0851481600f, 0.0887265700f, 0.0926600800f, 0.0968972300f, 0.1013746000f, 0.1060145000f, 0.1107377000f, 0.1155111000f, 0.1203122000f, 0.1251161000f,
    0.1298957000f, 0.1346299000f, 0.1393309000f, 0.1440235000f, 0.1487372000f, 0.1535066000f, 0.1583644000f, 0.1633199000f, 0.1683761000f, 0.1735365000f,
    0.1788048000f, 0.1841819000f, 0.1896559000f, 0.1952101000f, 0.2008259000f, 0.2064828000f, 0.2121826000f, 0.2180279000f, 0.2241586000f, 0.2307302000f,
    0.2379160000f, 0.2458706000f, 0.2546023000f, 0.2640760000f, 0.2742490000f, 0.2850680000f, 0.2964837000f, 0.3085010000f, 0.3211393000f, 0.3344175000f,
    0.3483536000f, 0.3629601000f, 0.3782275000f, 0.3941359000f, 0.4106582000f, 0.4277595000f, 0.4453993000f, 0.4635396000f, 0.4821376000f, 0.5011430000f,
    0.5204972000f, 0.5401387000f, 0.5600208000f, 0.5800972000f, 0.6003172000f, 0.6206256000f, 0.6409398000f, 0.6610772000f, 0.6808134000f, 0.6999044000f,
    0.7180890000f, 0.7351593000f, 0.7511821000f, 0.7663143000f, 0.7807352000f, 0.7946448000f, 0.8082074000f, 0.8213817000f, 0.8340701000f, 0.8461711000f,
    0.8575799000f, 0.8682408000f, 0.8783061000f, 0.8879907000f, 0.8975211000f, 0.9071347000f, 0.9169947000f, 0.9269295000f, 0.9366731000f, 0.9459482000f,
    0.9544675000f, 0.9619834000f, 0.9684390000f, 0.9738289000f, 0.9781519000f, 0.9814106000f, 0.9836669000f, 0.9852081000f, 0.9863813000f, 0.9875357000f,
    0.9890228000f, 0.9910811000f, 0.9934913000f, 0.9959172000f, 0.9980205000f, 0.9994608000f, 0.9999930000f, 0.9997557000f, 0.9989839000f, 0.9979123000f,
    0.9967737000f, 0.9957356000f, 0.9947115000f, 0.9935534000f, 0.9921156000f, 0.9902549000f, 0.9878596000f, 0.9849324000f, 0.9815036000f, 0.9776035000f,
    0.9732611000f, 0.9684764000f, 0.9631369000f, 0.9571062000f, 0.9502540000f, 0.9424569000f, 0.9336897000f, 0.9242893000f, 0.9146707000f, 0.9052333000f,
    0.8963613000f, 0.8883069000f, 0.8808462000f, 0.8736445000f, 0.8663755000f, 0.8587203000f, 0.8504295000f, 0.8415047000f, 0.8320109000f, 0.8220154000f,
    0.8115868000f, 0.8007874000f, 0.7896515000f, 0.7782053000f, 0.7664733000f, 0.7544785000f, 0.7422473000f, 0.7298229000f, 0.7172525000f, 0.7045818000f,
    0.6918553000f, 0.6791009000f, 0.6662846000f, 0.6533595000f, 0.6402807000f, 0.6270066000f, 0.6135148000f, 0.5998494000f, 0.5860682000f, 0.5722261000f,
    0.5583746000f, 0.5445535000f, 0.5307673000f, 0.5170130000f, 0.5032889000f, 0.4895950000f, 0.4759442000f, 0.4623958000f, 0.4490154000f, 0.4358622000f,
    0.4229897000f, 0.4104152000f, 0.3980356000f, 0.3857300000f, 0.3733907000f, 0.3609245000f, 0.3482860000f, 0.3355702000f, 0.3228963000f, 0.3103704000f,
    0.2980865000f, 0.2861160000f, 0.2744822000f, 0.2631953000f, 0.2522628000f, 0.2416902000f, 0.2314809000f, 0.2216378000f, 0.2121622000f, 0.2030542000f,
    0.1943124000f, 0.1859227000f, 0.1778274000f, 0.1699654000f, 0.1622841000f, 0.1547397000f, 0.1473081000f, 0.1400169000f, 0.1329013000f, 0.1259913000f,
    0.1193120000f, 0.1128820000f, 0.1067113000f, 0.1008052000f, 0.0951665300f, 0.0897959400f, 0.0846904400f, 0.0798400900f, 0.0752337200f, 0.0708606100f,
    0.0667104500f, 0.0627736000f, 0.0590417900f, 0.0555070300f, 0.0521613900f, 0.0489969900f, 0.0460057800f, 0.0431788500f, 0.0405075500f, 0.0379837600f,
    0.0355998200f, 0.0333485600f, 0.0312233200f, 0.0292178000f, 0.0273260100f, 0.0255422300f, 0.0238612100f, 0.0222785900f, 0.0207902000f, 0.0193918500f,
    0.0180793900f, 0.0168481700f, 0.0156918800f, 0.0146044600f, 0.0135806200f, 0.0126157300f, 0.0117069600f, 0.0108560800f, 0.0100647600f, 0.0093333760f,
    0.0086612840f, 0.0080460480f, 0.0074811300f, 0.0069599870f, 0.0064770700f, 0.0060276770f, 0.0056081690f, 0.0052166910f, 0.0048517850f, 0.0045120080f,
    0.0041959410f, 0.0039020570f, 0.0036283710f, 0.0033730050f, 0.0031343150f, 0.0029108640f, 0.0027015280f, 0.0025057960f, 0.0023232310f, 0.0021533330f,
    0.0019955570f, 0.0018493160f, 0.0017139760f, 0.0015888990f, 0.0014734530f, 0.0013670220f, 0.0012689540f, 0.0011784210f, 0.0010946440f, 0.0010169430f,
    0.0009447269f, 0.0008775171f, 0.0008150438f, 0.0007570755f, 0.0007033755f, 0.0006537050f, 0.0006078048f, 0.0005653435f, 0.0005260046f, 0.0004895061f,
    0.0004555970f, 0.0004240548f, 0.0003946860f, 0.0003673178f, 0.0003417941f, 0.0003179738f, 0.0002957441f, 0.0002750558f, 0.0002558640f, 0.0002381142f,
    0.0002217445f, 0.0002066711f, 0.0001927474f, 0.0001798315f, 0.0001678023f, 0.0001565566f, 0.0001460168f, 0.0001361535f, 0.0001269451f, 0.0001183671f,
    0.0001103928f, 0.0001029908f, 0.0000961184f, 0.0000897332f, 0.0000837969f, 0.0000782744f, 0.0000731331f, 0.0000683414f, 0.0000638704f, 0.0000596939f,
    0.0000557886f, 0.0000521351f, 0.0000487218f, 0.0000455385f, 0.0000425744f, 0.0000398188f, 0.0000372588f, 0.0000348747f, 0.0000326477f, 0.0000305614f,
    0.0000286018f, 0.0000267584f, 0.0000250294f, 0.0000234137f, 0.0000219091f, 0.0000205126f, 0.0000192190f, 0.0000180180f, 0.0000168990f, 0.0000158531f});
static constexpr Spectrum CIE2012_Z({0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0184726000f, 0.0222110100f, 0.0266981900f, 0.0320693700f, 0.0384783200f, 0.0460978400f, 0.0551195300f, 0.0657525700f, 0.0782211300f, 0.0927601300f,
    0.1096090000f, 0.1290077000f, 0.1512047000f, 0.1764441000f, 0.2049517000f, 0.2369246000f, 0.2725123000f, 0.3117820000f, 0.3547064000f, 0.4011473000f,
    0.4508369000f, 0.5034164000f, 0.5586361000f, 0.6162734000f, 0.6760982000f, 0.7378822000f, 0.8013019000f, 0.8655573000f, 0.9295791000f, 0.9921293000f,
    1.0518210000f, 1.1075090000f, 1.1595270000f, 1.2088690000f, 1.2568340000f, 1.3050080000f, 1.3547580000f, 1.4055940000f, 1.4564140000f, 1.5059600000f,
    1.5528260000f, 1.5959020000f, 1.6357680000f, 1.6735730000f, 1.7106040000f, 1.7482800000f, 1.7875040000f, 1.8266090000f, 1.8631080000f, 1.8943320000f,
    1.9174790000f, 1.9305290000f, 1.9348190000f, 1.9326500000f, 1.9263950000f, 1.9184370000f, 1.9104300000f, 1.9012240000f, 1.8890000000f, 1.8719960000f,
    1.8485450000f, 1.8177920000f, 1.7816270000f, 1.7425140000f, 1.7027490000f, 1.6644390000f, 1.6292070000f, 1.5973600000f, 1.5688960000f, 1.5438230000f,
    1.5221570000f, 1.5036110000f, 1.4866730000f, 1.4695950000f, 1.4507090000f, 1.4284400000f, 1.4015870000f, 1.3700940000f, 1.3342200000f, 1.2942750000f,
    1.2506100000f, 1.2036960000f, 1.1543160000f, 1.1032840000f, 1.0513470000f, 0.9991789000f, 0.9473958000f, 0.8966222000f, 0.8473981000f, 0.8001576000f,
    0.7552379000f, 0.7127879000f, 0.6725198000f, 0.6340976000f, 0.5972433000f, 0.5617313000f, 0.5274921000f, 0.4948809000f, 0.4642586000f, 0.4358841000f,
    0.4099313000f, 0.3864261000f, 0.3650566000f, 0.3454812000f, 0.3274095000f, 0.3105939000f, 0.2948102000f, 0.2798194000f, 0.2654100000f, 0.2514084000f,
    0.2376753000f, 0.2241211000f, 0.2107484000f, 0.1975839000f, 0.1846574000f, 0.1720018000f, 0.1596918000f, 0.1479415000f, 0.1369428000f, 0.1268279000f,
    0.1176796000f, 0.1094970000f, 0.1020943000f, 0.0952799300f, 0.0889007500f, 0.0828354800f, 0.0770098200f, 0.0714400100f, 0.0661543600f, 0.0611719900f,
    0.0565040700f, 0.0521512100f, 0.0480956600f, 0.0443172000f, 0.0407973400f, 0.0375191200f, 0.0344684600f, 0.0316376400f, 0.0290190100f, 0.0266036400f,
    0.0243816400f, 0.0223409700f, 0.0204641500f, 0.0187345600f, 0.0171378800f, 0.0156617400f, 0.0142964400f, 0.0130370200f, 0.0118789700f, 0.0108172500f,
    0.0098464700f, 0.0089606870f, 0.0081528110f, 0.0074160250f, 0.0067441150f, 0.0061314210f, 0.0055727780f, 0.0050634630f, 0.0045991690f, 0.0041759710f,
    0.0037902910f, 0.0034389520f, 0.0031193410f, 0.0028290380f, 0.0025657220f, 0.0023271860f, 0.0021112800f, 0.0019157660f, 0.0017385890f, 0.0015779200f,
    0.0014321280f, 0.0012997810f, 0.0011796670f, 0.0010706940f, 0.0009718623f, 0.0008822531f, 0.0008010231f, 0.0007273884f, 0.0006606347f, 0.0006001146f,
    0.0005452416f, 0.0004954847f, 0.0004503642f, 0.0004094455f, 0.0003723345f, 0.0003386739f, 0.0003081396f, 0.0002804370f, 0.0002552996f, 0.0002324859f,
    0.0002117772f, 0.0001929758f, 0.0001759024f, 0.0001603947f, 0.0001463059f, 0.0001335031f, 0.0001218660f, 0.0001112857f, 0.0001016634f, 0.0000929100f,
    0.0000849447f, 0.0000776943f, 0.0000710925f, 0.0000650794f, 0.0000596006f, 0.0000546071f, 0.0000500542f, 0.0000459016f, 0.0000421127f, 0.0000386544f,
    0.0000354966f, 0.0000326122f, 0.0000299764f, 0.0000275669f, 0.0000253634f, 0.0000233474f, 0.0000215022f, 0.0000198127f, 0.0000182650f, 0.0000168467f,
    0.0000155463f, 0.0000143536f, 0.0000132592f, 0.0000122544f, 0.0000113317f, 0.0000104839f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f,
    0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f, 0.0000000000f});

// standard illuminant D65( linear interpolated to 1nm )
static constexpr Spectrum CIE_D65({49.97550f, 50.44276f, 50.91002f, 51.37728f, 51.84454f, 52.31180f, 52.77908f, 53.24636f, 53.71364f, 54.18092f,
    54.64820f, 57.45886f, 60.26952f, 63.08018f, 65.89084f, 68.70150f, 71.51218f, 74.32286f, 77.13354f, 79.94422f,
    82.75490f, 83.62800f, 84.50110f, 85.37420f, 86.24730f, 87.12040f, 87.99352f, 88.86664f, 89.73976f, 90.61288f,
    91.48600f, 91.68058f, 91.87516f, 92.06974f, 92.26432f, 92.45890f, 92.65348f, 92.84806f, 93.04264f, 93.23722f,
    93.43180f, 92.75684f, 92.08188f, 91.40692f, 90.73196f, 90.05700f, 89.38206f, 88.70712f, 88.03218f, 87.35724f,
    86.68230f, 88.50056f, 90.31882f, 92.13708f, 93.95534f, 95.77360f, 97.59188f, 99.41016f, 101.22844f, 103.04672f,
    104.86500f, 106.07920f, 107.29340f, 108.50760f, 109.72180f, 110.93600f, 112.15040f, 113.36480f, 114.57920f, 115.79360f,
    117.00800f, 117.08840f, 117.16880f, 117.24920f, 117.32960f, 117.41000f, 117.49040f, 117.57080f, 117.65120f, 117.73160f,
    117.81200f, 117.51680f, 117.22160f, 116.92640f, 116.63120f, 116.33600f, 116.04100f, 115.74600f, 115.45100f, 115.15600f,
    114.86100f, 114.96720f, 115.07340f, 115.17960f, 115.28580f, 115.39200f, 115.49820f, 115.60440f, 115.71060f, 115.81680f,
    115.92300f, 115.21180f, 114.50060f, 113.78940f, 113.07820f, 112.36700f, 111.65580f, 110.94460f, 110.23340f, 109.52220f,
    108.81100f, 108.86520f, 108.91940f, 108.97360f, 109.02780f, 109.08200f, 109.13640f, 109.19080f, 109.24520f, 109.29960f,
    109.35400f, 109.19880f, 109.04360f, 108.88840f, 108.73320f, 108.57800f, 108.42280f, 108.26760f, 108.11240f, 107.95720f,
    107.80200f, 107.50080f, 107.19960f, 106.89840f, 106.59720f, 106.29600f, 105.99480f, 105.69360f, 105.39240f, 105.09120f,
    104.79000f, 105.07980f, 105.36960f, 105.65940f, 105.94920f, 106.23900f, 106.52900f, 106.81900f, 107.10900f, 107.39900f,
    107.68900f, 107.36060f, 107.03220f, 106.70380f, 106.37540f, 106.04700f, 105.71860f, 105.39020f, 105.06180f, 104.73340f,
    104.40500f, 104.36900f, 104.33300f, 104.29700f, 104.26100f, 104.22500f, 104.18920f, 104.15340f, 104.11760f, 104.08180f,
    104.04600f, 103.64140f, 103.23680f, 102.83220f, 102.42760f, 102.02300f, 101.61840f, 101.21380f, 100.80920f, 100.40460f,
    100.00000f, 99.63342f, 99.26684f, 98.90026f, 98.53368f, 98.16710f, 97.80052f, 97.43394f, 97.06736f, 96.70078f,
    96.33420f, 96.27958f, 96.22496f, 96.17034f, 96.11572f, 96.06110f, 96.00648f, 95.95186f, 95.89724f, 95.84262f,
    95.78800f, 95.07776f, 94.36752f, 93.65728f, 92.94704f, 92.23680f, 91.52656f, 90.81632f, 90.10608f, 89.39584f,
    88.68560f, 88.81766f, 88.94972f, 89.08178f, 89.21384f, 89.34590f, 89.47796f, 89.61002f, 89.74208f, 89.87414f,
    90.00620f, 89.96548f, 89.92476f, 89.88404f, 89.84332f, 89.80260f, 89.76190f, 89.72120f, 89.68050f, 89.63980f,
    89.59910f, 89.40906f, 89.21902f, 89.02898f, 88.83894f, 88.64890f, 88.45886f, 88.26882f, 88.07878f, 87.88874f,
    87.69870f, 87.25768f, 86.81666f, 86.37564f, 85.93462f, 85.49360f, 85.05260f, 84.61160f, 84.17060f, 83.72960f,
    83.28860f, 83.32966f, 83.37072f, 83.41178f, 83.45284f, 83.49390f, 83.53496f, 83.57602f, 83.61708f, 83.65814f,
    83.69920f, 83.33196f, 82.96472f, 82.59748f, 82.23024f, 81.86300f, 81.49576f, 81.12852f, 80.76128f, 80.39404f,
    80.02680f, 80.04558f, 80.06436f, 80.08314f, 80.10192f, 80.12070f, 80.13948f, 80.15826f, 80.17704f, 80.19582f,
    80.21460f, 80.42092f, 80.62724f, 80.83356f, 81.03988f, 81.24620f, 81.45252f, 81.65884f, 81.86516f, 82.07148f,
    82.27780f, 81.87844f, 81.47908f, 81.07972f, 80.68036f, 80.28100f, 79.88164f, 79.48228f, 79.08292f, 78.68356f,
    78.28420f, 77.42790f, 76.57160f, 75.71530f, 74.85900f, 74.00270f, 73.14642f, 72.29014f, 71.43386f, 70.57758f,
    69.72130f, 69.91008f, 70.09886f, 70.28764f, 70.47642f, 70.66520f, 70.85398f, 71.04276f, 71.23154f, 71.42032f,
    71.60910f, 71.88308f, 72.15706f, 72.43104f, 72.70502f, 72.97900f, 73.25300f, 73.52700f, 73.80100f, 74.07500f,
    74.34900f, 73.07450f, 71.80000f, 70.52550f, 69.25100f, 67.97650f, 66.70200f, 65.42750f, 64.15300f, 62.87850f,
    61.60400f, 62.43216f, 63.26032f, 64.08848f, 64.91664f, 65.74480f, 66.57296f, 67.40112f, 68.22928f, 69.05744f,
    69.88560f, 70.40574f, 70.92588f, 71.44602f, 71.96616f, 72.48630f, 73.00644f, 73.52658f, 74.04672f, 74.56686f,
    75.08700f, 73.93756f, 72.78812f, 71.63868f, 70.48924f, 69.33980f, 68.19038f, 67.04096f, 65.89154f, 64.74212f,
    63.59270f, 61.87524f, 60.15778f, 58.44032f, 56.72286f, 55.00540f, 53.28796f, 51.57052f, 49.85308f, 48.13564f,
    46.41820f, 48.45692f, 50.49564f, 52.53436f, 54.57308f, 56.61180f, 58.65052f, 60.68924f, 62.72796f, 64.76668f,
    66.80540f, 66.46314f, 66.12088f, 65.77862f, 65.43636f, 65.09410f, 64.75184f, 64.40958f, 64.06732f, 63.72506f});

namespace Macbeth
{
    static constexpr Spectrum Macbeth01({0.054745777f, 0.055103879f, 0.055461981f, 0.055820083f, 0.056178185f, 0.056536288f, 0.05689439f, 0.057252492f, 0.057610594f, 0.057968696f,
        0.058326798f, 0.058610495f, 0.058894191f, 0.059177887f, 0.059461584f, 0.05974528f, 0.060028976f, 0.060312673f, 0.060596369f, 0.060880065f,
        0.061163762f, 0.061285862f, 0.061407963f, 0.061530064f, 0.061652164f, 0.061774265f, 0.061896366f, 0.062018466f, 0.062140567f, 0.062262668f,
        0.062384768f, 0.062377367f, 0.062369965f, 0.062362563f, 0.062355162f, 0.06234776f, 0.062340358f, 0.062332957f, 0.062325555f, 0.062318153f,
        0.062310752f, 0.062286679f, 0.062262607f, 0.062238535f, 0.062214462f, 0.06219039f, 0.062166318f, 0.062142245f, 0.062118173f, 0.062094101f,
        0.062070028f, 0.062045879f, 0.06202173f, 0.06199758f, 0.061973431f, 0.061949282f, 0.061925132f, 0.061900983f, 0.061876834f, 0.061852684f,
        0.061828535f, 0.061804253f, 0.061779972f, 0.06175569f, 0.061731408f, 0.061707127f, 0.061682845f, 0.061658563f, 0.061634282f, 0.06161f,
        0.061585718f, 0.061581288f, 0.061576857f, 0.061572426f, 0.061567995f, 0.061563564f, 0.061559133f, 0.061554703f, 0.061550272f, 0.061545841f,
        0.06154141f, 0.061549581f, 0.061557751f, 0.061565922f, 0.061574092f, 0.061582263f, 0.061590433f, 0.061598604f, 0.061606774f, 0.061614945f,
        0.061623115f, 0.061663561f, 0.061704008f, 0.061744454f, 0.0617849f, 0.061825347f, 0.061865793f, 0.061906239f, 0.061946686f, 0.061987132f,
        0.062027578f, 0.062121243f, 0.062214907f, 0.062308571f, 0.062402235f, 0.062495899f, 0.062589563f, 0.062683228f, 0.062776892f, 0.062870556f,
        0.06296422f, 0.063185879f, 0.063407537f, 0.063629196f, 0.063850855f, 0.064072513f, 0.064294172f, 0.064515831f, 0.064737489f, 0.064959148f,
        0.065180807f, 0.065689872f, 0.066198937f, 0.066708003f, 0.067217068f, 0.067726133f, 0.068235199f, 0.068744264f, 0.069253329f, 0.069762395f,
        0.07027146f, 0.070883902f, 0.071496344f, 0.072108786f, 0.072721227f, 0.073333669f, 0.073946111f, 0.074558553f, 0.075170995f, 0.075783437f,
        0.076395878f, 0.076705358f, 0.077014837f, 0.077324316f, 0.077633795f, 0.077943274f, 0.078252753f, 0.078562233f, 0.078871712f, 0.079181191f,
        0.07949067f, 0.079669186f, 0.079847702f, 0.080026218f, 0.080204734f, 0.08038325f, 0.080561766f, 0.080740282f, 0.080918798f, 0.081097314f,
        0.08127583f, 0.081577204f, 0.081878577f, 0.082179951f, 0.082481324f, 0.082782698f, 0.083084071f, 0.083385445f, 0.083686818f, 0.083988192f,
        0.084289565f, 0.084918979f, 0.085548393f, 0.086177807f, 0.086807221f, 0.087436635f, 0.088066049f, 0.088695463f, 0.089324877f, 0.089954291f,
        0.090583705f, 0.091815509f, 0.093047312f, 0.094279116f, 0.09551092f, 0.096742723f, 0.097974527f, 0.099206331f, 0.100438134f, 0.101669938f,
        0.102901742f, 0.104516701f, 0.10613166f, 0.107746619f, 0.109361578f, 0.110976538f, 0.112591497f, 0.114206456f, 0.115821415f, 0.117436374f,
        0.119051333f, 0.120572658f, 0.122093982f, 0.123615306f, 0.12513663f, 0.126657954f, 0.128179278f, 0.129700603f, 0.131221927f, 0.132743251f,
        0.134264575f, 0.135157913f, 0.13605125f, 0.136944588f, 0.137837925f, 0.138731263f, 0.1396246f, 0.140517938f, 0.141411275f, 0.142304613f,
        0.14319795f, 0.143566098f, 0.143934245f, 0.144302393f, 0.14467054f, 0.145038688f, 0.145406835f, 0.145774983f, 0.14614313f, 0.146511278f,
        0.146879425f, 0.147268983f, 0.14765854f, 0.148048098f, 0.148437655f, 0.148827213f, 0.14921677f, 0.149606328f, 0.149995885f, 0.150385443f,
        0.150775f, 0.151507093f, 0.152239185f, 0.152971278f, 0.15370337f, 0.154435463f, 0.155167555f, 0.155899648f, 0.15663174f, 0.157363833f,
        0.158095925f, 0.159105434f, 0.160114943f, 0.161124453f, 0.162133962f, 0.163143471f, 0.16415298f, 0.165162489f, 0.166171998f, 0.167181508f,
        0.168191017f, 0.169262342f, 0.170333667f, 0.171404992f, 0.172476317f, 0.173547642f, 0.174618967f, 0.175690292f, 0.176761617f, 0.177832942f,
        0.178904267f, 0.179768972f, 0.180633677f, 0.181498382f, 0.182363087f, 0.183227792f, 0.184092497f, 0.184957202f, 0.185821907f, 0.186686612f,
        0.187551317f, 0.187759782f, 0.187968247f, 0.188176712f, 0.188385177f, 0.188593642f, 0.188802107f, 0.189010572f, 0.189219037f, 0.189427502f,
        0.189635967f, 0.189249621f, 0.188863275f, 0.188476929f, 0.188090583f, 0.187704238f, 0.187317892f, 0.186931546f, 0.1865452f, 0.186158854f,
        0.185772508f, 0.185344593f, 0.184916678f, 0.184488763f, 0.184060848f, 0.183632933f, 0.183205018f, 0.182777103f, 0.182349188f, 0.181921273f,
        0.181493358f, 0.181504526f, 0.181515693f, 0.181526861f, 0.181538028f, 0.181549196f, 0.181560363f, 0.181571531f, 0.181582698f, 0.181593866f,
        0.181605033f, 0.18216564f, 0.182726247f, 0.183286853f, 0.18384746f, 0.184408067f, 0.184968673f, 0.18552928f, 0.186089887f, 0.186650493f,
        0.1872111f, 0.188095342f, 0.188979583f, 0.189863825f, 0.190748067f, 0.191632308f, 0.19251655f, 0.193400792f, 0.194285033f, 0.195169275f,
        0.196053517f, 0.197397435f, 0.198741353f, 0.200085272f, 0.20142919f, 0.202773108f, 0.204117027f, 0.205460945f, 0.206804863f, 0.208148782f,
        0.2094927f, 0.18854343f, 0.16759416f, 0.14664489f, 0.12569562f, 0.10474635f, 0.08379708f, 0.06284781f, 0.04189854f, 0.02094927f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth02({0.117132943f, 0.119765069f, 0.122397195f, 0.12502932f, 0.127661446f, 0.130293572f, 0.132925697f, 0.135557823f, 0.138189949f, 0.140822074f,
        0.1434542f, 0.146561686f, 0.149669172f, 0.152776658f, 0.155884143f, 0.158991629f, 0.162099115f, 0.165206601f, 0.168314087f, 0.171421573f,
        0.174529058f, 0.176169633f, 0.177810208f, 0.179450783f, 0.181091358f, 0.182731933f, 0.184372508f, 0.186013083f, 0.187653658f, 0.189294233f,
        0.190934808f, 0.19140154f, 0.191868272f, 0.192335003f, 0.192801735f, 0.193268467f, 0.193735198f, 0.19420193f, 0.194668662f, 0.195135393f,
        0.195602125f, 0.195941709f, 0.196281293f, 0.196620878f, 0.196960462f, 0.197300046f, 0.19763963f, 0.197979214f, 0.198318798f, 0.198658383f,
        0.198997967f, 0.199521663f, 0.20004536f, 0.200569057f, 0.201092753f, 0.20161645f, 0.202140147f, 0.202663843f, 0.20318754f, 0.203711237f,
        0.204234933f, 0.20512993f, 0.206024927f, 0.206919923f, 0.20781492f, 0.208709917f, 0.209604913f, 0.21049991f, 0.211394907f, 0.212289903f,
        0.2131849f, 0.214708742f, 0.216232583f, 0.217756425f, 0.219280267f, 0.220804108f, 0.22232795f, 0.223851792f, 0.225375633f, 0.226899475f,
        0.228423317f, 0.230708312f, 0.232993307f, 0.235278302f, 0.237563297f, 0.239848292f, 0.242133287f, 0.244418282f, 0.246703277f, 0.248988272f,
        0.251273267f, 0.254150625f, 0.257027983f, 0.259905342f, 0.2627827f, 0.265660058f, 0.268537417f, 0.271414775f, 0.274292133f, 0.277169492f,
        0.28004685f, 0.282920274f, 0.285793698f, 0.288667123f, 0.291540547f, 0.294413971f, 0.297287395f, 0.300160819f, 0.303034243f, 0.305907668f,
        0.308781092f, 0.310848142f, 0.312915192f, 0.314982242f, 0.317049292f, 0.319116342f, 0.321183392f, 0.323250442f, 0.325317492f, 0.327384542f,
        0.329451592f, 0.329842878f, 0.330234163f, 0.330625449f, 0.331016735f, 0.331408021f, 0.331799307f, 0.332190593f, 0.332581878f, 0.332973164f,
        0.33336445f, 0.331487665f, 0.32961088f, 0.327734095f, 0.32585731f, 0.323980525f, 0.32210374f, 0.320226955f, 0.31835017f, 0.316473385f,
        0.3145966f, 0.31176484f, 0.30893308f, 0.30610132f, 0.30326956f, 0.3004378f, 0.29760604f, 0.29477428f, 0.29194252f, 0.28911076f, 0.286279f,
        0.284999998f, 0.283720997f, 0.282441995f, 0.281162993f, 0.279883992f, 0.27860499f, 0.277325988f, 0.276046987f, 0.274767985f, 0.273488983f,
        0.27378576f, 0.274082537f, 0.274379313f, 0.27467609f, 0.274972867f, 0.275269643f, 0.27556642f, 0.275863197f, 0.276159973f, 0.27645675f,
        0.276531368f, 0.276605987f, 0.276680605f, 0.276755223f, 0.276829842f, 0.27690446f, 0.276979078f, 0.277053697f, 0.277128315f, 0.277202933f,
        0.278412386f, 0.279621838f, 0.280831291f, 0.282040743f, 0.283250196f, 0.284459648f, 0.285669101f, 0.286878553f, 0.288088006f, 0.289297458f,
        0.294305343f, 0.299313227f, 0.304321111f, 0.309328995f, 0.314336879f, 0.319344763f, 0.324352648f, 0.329360532f, 0.334368416f, 0.3393763f,
        0.347460277f, 0.355544253f, 0.36362823f, 0.371712207f, 0.379796183f, 0.38788016f, 0.395964137f, 0.404048113f, 0.41213209f, 0.420216067f,
        0.426973407f, 0.433730747f, 0.440488087f, 0.447245427f, 0.454002767f, 0.460760107f, 0.467517447f, 0.474274787f, 0.481032127f, 0.487789467f,
        0.491521339f, 0.495253212f, 0.498985084f, 0.502716957f, 0.506448829f, 0.510180702f, 0.513912574f, 0.517644447f, 0.521376319f, 0.525108192f,
        0.527171099f, 0.529234007f, 0.531296914f, 0.533359822f, 0.535422729f, 0.537485637f, 0.539548544f, 0.541611452f, 0.543674359f, 0.545737267f,
        0.547319515f, 0.548901763f, 0.550484012f, 0.55206626f, 0.553648508f, 0.555230757f, 0.556813005f, 0.558395253f, 0.559977502f, 0.56155975f,
        0.563192083f, 0.564824417f, 0.56645675f, 0.568089083f, 0.569721417f, 0.57135375f, 0.572986083f, 0.574618417f, 0.57625075f, 0.577883083f,
        0.579591896f, 0.581300708f, 0.583009521f, 0.584718333f, 0.586427146f, 0.588135958f, 0.589844771f, 0.591553583f, 0.593262396f, 0.594971208f,
        0.596653813f, 0.598336417f, 0.600019021f, 0.601701625f, 0.603384229f, 0.605066833f, 0.606749438f, 0.608432042f, 0.610114646f, 0.61179725f,
        0.613092223f, 0.614387195f, 0.615682168f, 0.61697714f, 0.618272113f, 0.619567085f, 0.620862058f, 0.62215703f, 0.623452003f, 0.624746975f,
        0.626082417f, 0.627417858f, 0.6287533f, 0.630088742f, 0.631424183f, 0.632759625f, 0.634095067f, 0.635430508f, 0.63676595f, 0.638101392f,
        0.639887074f, 0.641672757f, 0.643458439f, 0.645244122f, 0.647029804f, 0.648815487f, 0.650601169f, 0.652386852f, 0.654172534f, 0.655958217f,
        0.65818396f, 0.660409703f, 0.662635447f, 0.66486119f, 0.667086933f, 0.669312677f, 0.67153842f, 0.673764163f, 0.675989907f, 0.67821565f,
        0.680352183f, 0.682488717f, 0.68462525f, 0.686761783f, 0.688898317f, 0.69103485f, 0.693171383f, 0.695307917f, 0.69744445f, 0.699580983f,
        0.701331758f, 0.703082532f, 0.704833306f, 0.70658408f, 0.708334854f, 0.710085628f, 0.711836403f, 0.713587177f, 0.715337951f, 0.717088725f,
        0.718762279f, 0.720435833f, 0.722109388f, 0.723782942f, 0.725456496f, 0.72713005f, 0.728803604f, 0.730477158f, 0.732150713f, 0.733824267f,
        0.66044184f, 0.587059413f, 0.513676987f, 0.44029456f, 0.366912133f, 0.293529707f, 0.22014728f, 0.146764853f, 0.073382427f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth03({0.130362525f, 0.135033678f, 0.13970483f, 0.144375983f, 0.149047135f, 0.153718288f, 0.15838944f, 0.163060593f, 0.167731745f, 0.172402898f,
        0.17707405f, 0.184467636f, 0.191861222f, 0.199254808f, 0.206648393f, 0.214041979f, 0.221435565f, 0.228829151f, 0.236222737f, 0.243616323f,
        0.251009908f, 0.256533426f, 0.262056943f, 0.267580461f, 0.273103978f, 0.278627496f, 0.284151013f, 0.289674531f, 0.295198048f, 0.300721566f,
        0.306245083f, 0.308012125f, 0.309779167f, 0.311546208f, 0.31331325f, 0.315080292f, 0.316847333f, 0.318614375f, 0.320381417f, 0.322148458f,
        0.3239155f, 0.324517105f, 0.32511871f, 0.325720315f, 0.32632192f, 0.326923525f, 0.32752513f, 0.328126735f, 0.32872834f, 0.329329945f,
        0.32993155f, 0.330221306f, 0.330511062f, 0.330800818f, 0.331090573f, 0.331380329f, 0.331670085f, 0.331959841f, 0.332249597f, 0.332539353f,
        0.332829108f, 0.332643691f, 0.332458273f, 0.332272856f, 0.332087438f, 0.331902021f, 0.331716603f, 0.331531186f, 0.331345768f, 0.331160351f,
        0.330974933f, 0.330219834f, 0.329464735f, 0.328709636f, 0.327954537f, 0.327199438f, 0.326444338f, 0.325689239f, 0.32493414f, 0.324179041f,
        0.323423942f, 0.322215453f, 0.321006963f, 0.319798474f, 0.318589985f, 0.317381496f, 0.316173007f, 0.314964518f, 0.313756028f, 0.312547539f,
        0.31133905f, 0.310027743f, 0.308716437f, 0.30740513f, 0.306093823f, 0.304782517f, 0.30347121f, 0.302159903f, 0.300848597f, 0.29953729f,
        0.298225983f, 0.296936445f, 0.295646907f, 0.294357368f, 0.29306783f, 0.291778292f, 0.290488753f, 0.289199215f, 0.287909677f, 0.286620138f,
        0.2853306f, 0.283740298f, 0.282149995f, 0.280559693f, 0.27896939f, 0.277379088f, 0.275788785f, 0.274198483f, 0.27260818f, 0.271017878f,
        0.269427575f, 0.26752225f, 0.265616925f, 0.2637116f, 0.261806275f, 0.25990095f, 0.257995625f, 0.2560903f, 0.254184975f, 0.25227965f,
        0.250374325f, 0.248481174f, 0.246588023f, 0.244694873f, 0.242801722f, 0.240908571f, 0.23901542f, 0.237122269f, 0.235229118f, 0.233335968f,
        0.231442817f, 0.229724124f, 0.228005432f, 0.226286739f, 0.224568047f, 0.222849354f, 0.221130662f, 0.219411969f, 0.217693277f, 0.215974584f,
        0.214255892f, 0.212772542f, 0.211289192f, 0.209805842f, 0.208322492f, 0.206839142f, 0.205355792f, 0.203872442f, 0.202389092f, 0.200905742f,
        0.199422392f, 0.197931317f, 0.196440242f, 0.194949167f, 0.193458092f, 0.191967017f, 0.190475942f, 0.188984867f, 0.187493792f, 0.186002717f,
        0.184511642f, 0.182998908f, 0.181486173f, 0.179973439f, 0.178460705f, 0.176947971f, 0.175435237f, 0.173922503f, 0.172409768f, 0.170897034f,
        0.1693843f, 0.168174663f, 0.166965025f, 0.165755388f, 0.16454575f, 0.163336113f, 0.162126475f, 0.160916838f, 0.1597072f, 0.158497563f,
        0.157287925f, 0.156469936f, 0.155651947f, 0.154833958f, 0.154015968f, 0.153197979f, 0.15237999f, 0.151562001f, 0.150744012f, 0.149926023f,
        0.149108033f, 0.14867956f, 0.148251087f, 0.147822613f, 0.14739414f, 0.146965667f, 0.146537193f, 0.14610872f, 0.145680247f, 0.145251773f,
        0.1448233f, 0.144526793f, 0.144230287f, 0.14393378f, 0.143637273f, 0.143340767f, 0.14304426f, 0.142747753f, 0.142451247f, 0.14215474f,
        0.141858233f, 0.141729483f, 0.141600733f, 0.141471983f, 0.141343233f, 0.141214483f, 0.141085733f, 0.140956983f, 0.140828233f, 0.140699483f,
        0.140570733f, 0.140581125f, 0.140591517f, 0.140601908f, 0.1406123f, 0.140622692f, 0.140633083f, 0.140643475f, 0.140653867f, 0.140664258f,
        0.14067465f, 0.140716511f, 0.140758372f, 0.140800233f, 0.140842093f, 0.140883954f, 0.140925815f, 0.140967676f, 0.141009537f, 0.141051398f,
        0.141093258f, 0.141240748f, 0.141388238f, 0.141535728f, 0.141683218f, 0.141830708f, 0.141978198f, 0.142125688f, 0.142273178f, 0.142420668f,
        0.142568158f, 0.142965302f, 0.143362445f, 0.143759588f, 0.144156732f, 0.144553875f, 0.144951018f, 0.145348162f, 0.145745305f, 0.146142448f,
        0.146539592f, 0.147069393f, 0.147599195f, 0.148128997f, 0.148658798f, 0.1491886f, 0.149718402f, 0.150248203f, 0.150778005f, 0.151307807f,
        0.151837608f, 0.152005322f, 0.152173035f, 0.152340748f, 0.152508462f, 0.152676175f, 0.152843888f, 0.153011602f, 0.153179315f, 0.153347028f,
        0.153514742f, 0.153172203f, 0.152829663f, 0.152487124f, 0.152144585f, 0.151802046f, 0.151459507f, 0.151116968f, 0.150774428f, 0.150431889f,
        0.15008935f, 0.149475664f, 0.148861978f, 0.148248293f, 0.147634607f, 0.147020921f, 0.146407235f, 0.145793549f, 0.145179863f, 0.144566178f,
        0.143952492f, 0.143196713f, 0.142440933f, 0.141685154f, 0.140929375f, 0.140173596f, 0.139417817f, 0.138662038f, 0.137906258f, 0.137150479f,
        0.1363947f, 0.135990375f, 0.135586049f, 0.135181724f, 0.134777398f, 0.134373073f, 0.133968747f, 0.133564422f, 0.133160096f, 0.132755771f,
        0.132351445f, 0.132612045f, 0.132872644f, 0.133133244f, 0.133393844f, 0.133654443f, 0.133915043f, 0.134175643f, 0.134436242f, 0.134696842f,
        0.134957442f, 0.136135052f, 0.137312662f, 0.138490272f, 0.139667882f, 0.140845493f, 0.142023103f, 0.143200713f, 0.144378323f, 0.145555933f,
        0.146733543f, 0.132060189f, 0.117386835f, 0.10271348f, 0.088040126f, 0.073366772f, 0.058693417f, 0.044020063f, 0.029346709f, 0.014673354f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth04({0.05123929f, 0.051538762f, 0.051838234f, 0.052137706f, 0.052437178f, 0.05273665f, 0.053036122f, 0.053335594f, 0.053635066f, 0.053934538f,
        0.05423401f, 0.054409636f, 0.054585261f, 0.054760887f, 0.054936512f, 0.055112138f, 0.055287763f, 0.055463389f, 0.055639014f, 0.05581464f,
        0.055990265f, 0.056095125f, 0.056199985f, 0.056304845f, 0.056409705f, 0.056514565f, 0.056619425f, 0.056724285f, 0.056829145f, 0.056934005f,
        0.057038865f, 0.0571212f, 0.057203536f, 0.057285871f, 0.057368206f, 0.057450542f, 0.057532877f, 0.057615212f, 0.057697548f, 0.057779883f,
        0.057862218f, 0.05797104f, 0.058079861f, 0.058188682f, 0.058297503f, 0.058406324f, 0.058515145f, 0.058623967f, 0.058732788f, 0.058841609f,
        0.05895043f, 0.059085396f, 0.059220361f, 0.059355327f, 0.059490292f, 0.059625258f, 0.059760223f, 0.059895189f, 0.060030154f, 0.06016512f,
        0.060300085f, 0.060401482f, 0.060502878f, 0.060604275f, 0.060705671f, 0.060807068f, 0.060908464f, 0.061009861f, 0.061111257f, 0.061212654f,
        0.06131405f, 0.0614103f, 0.06150655f, 0.0616028f, 0.06169905f, 0.0617953f, 0.06189155f, 0.0619878f, 0.06208405f, 0.0621803f, 0.06227655f,
        0.062373402f, 0.062470254f, 0.062567106f, 0.062663958f, 0.06276081f, 0.062857662f, 0.062954514f, 0.063051366f, 0.063148218f, 0.06324507f,
        0.063398907f, 0.063552743f, 0.06370658f, 0.063860416f, 0.064014253f, 0.064168089f, 0.064321926f, 0.064475762f, 0.064629599f, 0.064783435f,
        0.065042848f, 0.065302261f, 0.065561674f, 0.065821086f, 0.066080499f, 0.066339912f, 0.066599325f, 0.066858738f, 0.067118151f, 0.067377563f,
        0.068171087f, 0.068964611f, 0.069758135f, 0.070551659f, 0.071345183f, 0.072138707f, 0.072932231f, 0.073725755f, 0.074519279f, 0.075312803f,
        0.077901374f, 0.080489944f, 0.083078515f, 0.085667085f, 0.088255656f, 0.090844226f, 0.093432797f, 0.096021367f, 0.098609938f, 0.101198508f,
        0.105614732f, 0.110030955f, 0.114447178f, 0.118863402f, 0.123279625f, 0.127695848f, 0.132112072f, 0.136528295f, 0.140944518f, 0.145360742f,
        0.148651089f, 0.151941437f, 0.155231784f, 0.158522132f, 0.161812479f, 0.165102827f, 0.168393174f, 0.171683522f, 0.174973869f, 0.178264217f,
        0.178832204f, 0.179400192f, 0.179968179f, 0.180536167f, 0.181104154f, 0.181672142f, 0.182240129f, 0.182808117f, 0.183376104f, 0.183944092f,
        0.182560788f, 0.181177485f, 0.179794182f, 0.178410878f, 0.177027575f, 0.175644272f, 0.174260968f, 0.172877665f, 0.171494362f, 0.170111058f,
        0.168038043f, 0.165965027f, 0.163892011f, 0.161818995f, 0.159745979f, 0.157672963f, 0.155599948f, 0.153526932f, 0.151453916f, 0.1493809f,
        0.147717223f, 0.146053547f, 0.14438987f, 0.142726193f, 0.141062517f, 0.13939884f, 0.137735163f, 0.136071487f, 0.13440781f, 0.132744133f,
        0.131655826f, 0.130567518f, 0.129479211f, 0.128390903f, 0.127302596f, 0.126214288f, 0.125125981f, 0.124037673f, 0.122949366f, 0.121861058f,
        0.121192263f, 0.120523468f, 0.119854673f, 0.119185878f, 0.118517083f, 0.117848288f, 0.117179493f, 0.116510698f, 0.115841903f, 0.115173108f,
        0.114603934f, 0.11403476f, 0.113465586f, 0.112896412f, 0.112327238f, 0.111758063f, 0.111188889f, 0.110619715f, 0.110050541f, 0.109481367f,
        0.109069129f, 0.108656892f, 0.108244654f, 0.107832417f, 0.107420179f, 0.107007942f, 0.106595704f, 0.106183467f, 0.105771229f, 0.105358992f,
        0.105257307f, 0.105155623f, 0.105053938f, 0.104952254f, 0.104850569f, 0.104748885f, 0.1046472f, 0.104545516f, 0.104443831f, 0.104342147f,
        0.104506853f, 0.104671559f, 0.104836265f, 0.105000971f, 0.105165678f, 0.105330384f, 0.10549509f, 0.105659796f, 0.105824502f, 0.105989208f,
        0.106280923f, 0.106572637f, 0.106864351f, 0.107156065f, 0.107447779f, 0.107739493f, 0.108031208f, 0.108322922f, 0.108614636f, 0.10890635f,
        0.109205153f, 0.109503956f, 0.10980276f, 0.110101563f, 0.110400366f, 0.110699169f, 0.110997972f, 0.111296775f, 0.111595579f, 0.111894382f,
        0.112111166f, 0.11232795f, 0.112544735f, 0.112761519f, 0.112978303f, 0.113195088f, 0.113411872f, 0.113628656f, 0.113845441f, 0.114062225f,
        0.114051304f, 0.114040383f, 0.114029463f, 0.114018542f, 0.114007621f, 0.1139967f, 0.113985779f, 0.113974858f, 0.113963938f, 0.113953017f,
        0.113797518f, 0.11364202f, 0.113486522f, 0.113331023f, 0.113175525f, 0.113020027f, 0.112864528f, 0.11270903f, 0.112553532f, 0.112398033f,
        0.112373548f, 0.112349063f, 0.112324578f, 0.112300093f, 0.112275608f, 0.112251123f, 0.112226638f, 0.112202153f, 0.112177668f, 0.112153183f,
        0.112419765f, 0.112686347f, 0.112952928f, 0.11321951f, 0.113486092f, 0.113752673f, 0.114019255f, 0.114285837f, 0.114552418f, 0.114819f,
        0.115313643f, 0.115808285f, 0.116302928f, 0.11679757f, 0.117292213f, 0.117786855f, 0.118281498f, 0.11877614f, 0.119270783f, 0.119765425f,
        0.120248358f, 0.120731292f, 0.121214225f, 0.121697158f, 0.122180092f, 0.122663025f, 0.123145958f, 0.123628892f, 0.124111825f, 0.124594758f,
        0.125165738f, 0.125736717f, 0.126307696f, 0.126878675f, 0.127449654f, 0.128020633f, 0.128591613f, 0.129162592f, 0.129733571f, 0.13030455f,
        0.117274095f, 0.10424364f, 0.091213185f, 0.07818273f, 0.065152275f, 0.05212182f, 0.039091365f, 0.02606091f, 0.013030455f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth05({0.144234675f, 0.149638121f, 0.155041567f, 0.160445013f, 0.165848458f, 0.171251904f, 0.17665535f, 0.182058796f, 0.187462242f, 0.192865688f,
        0.198269133f, 0.207885608f, 0.217502082f, 0.227118556f, 0.23673503f, 0.246351504f, 0.255967978f, 0.265584453f, 0.275200927f, 0.284817401f,
        0.294433875f, 0.30253424f, 0.310634605f, 0.31873497f, 0.326835335f, 0.3349357f, 0.343036065f, 0.35113643f, 0.359236795f, 0.36733716f,
        0.375437525f, 0.378731121f, 0.382024717f, 0.385318313f, 0.388611908f, 0.391905504f, 0.3951991f, 0.398492696f, 0.401786292f, 0.405079888f,
        0.408373483f, 0.409631083f, 0.410888683f, 0.412146283f, 0.413403883f, 0.414661483f, 0.415919083f, 0.417176683f, 0.418434283f, 0.419691883f,
        0.420949483f, 0.421472757f, 0.42199603f, 0.422519303f, 0.423042577f, 0.42356585f, 0.424089123f, 0.424612397f, 0.42513567f, 0.425658943f,
        0.426182217f, 0.426172657f, 0.426163097f, 0.426153537f, 0.426143977f, 0.426134417f, 0.426124857f, 0.426115297f, 0.426105737f, 0.426096177f,
        0.426086617f, 0.425410258f, 0.424733898f, 0.424057539f, 0.42338118f, 0.422704821f, 0.422028462f, 0.421352103f, 0.420675743f, 0.419999384f,
        0.419323025f, 0.417734222f, 0.416145418f, 0.414556615f, 0.412967812f, 0.411379008f, 0.409790205f, 0.408201402f, 0.406612598f, 0.405023795f,
        0.403434992f, 0.401018862f, 0.398602732f, 0.396186602f, 0.393770472f, 0.391354342f, 0.388938212f, 0.386522082f, 0.384105952f, 0.381689822f,
        0.379273692f, 0.375982226f, 0.37269076f, 0.369399294f, 0.366107828f, 0.362816363f, 0.359524897f, 0.356233431f, 0.352941965f, 0.349650499f,
        0.346359033f, 0.342835295f, 0.339311557f, 0.335787818f, 0.33226408f, 0.328740342f, 0.325216603f, 0.321692865f, 0.318169127f, 0.314645388f,
        0.31112165f, 0.308132998f, 0.305144345f, 0.302155693f, 0.29916704f, 0.296178388f, 0.293189735f, 0.290201083f, 0.28721243f, 0.284223778f,
        0.281235125f, 0.278499214f, 0.275763303f, 0.273027393f, 0.270291482f, 0.267555571f, 0.26481966f, 0.262083749f, 0.259347838f, 0.256611928f,
        0.253876017f, 0.251377398f, 0.248878778f, 0.246380159f, 0.24388154f, 0.241382921f, 0.238884302f, 0.236385683f, 0.233887063f, 0.231388444f,
        0.228889825f, 0.227421307f, 0.225952788f, 0.22448427f, 0.223015752f, 0.221547233f, 0.220078715f, 0.218610197f, 0.217141678f, 0.21567316f,
        0.214204642f, 0.213619655f, 0.213034668f, 0.212449682f, 0.211864695f, 0.211279708f, 0.210694722f, 0.210109735f, 0.209524748f, 0.208939762f,
        0.208354775f, 0.207681531f, 0.207008287f, 0.206335043f, 0.205661798f, 0.204988554f, 0.20431531f, 0.203642066f, 0.202968822f, 0.202295578f,
        0.201622333f, 0.20089973f, 0.200177127f, 0.199454523f, 0.19873192f, 0.198009317f, 0.197286713f, 0.19656411f, 0.195841507f, 0.195118903f,
        0.1943963f, 0.194213338f, 0.194030375f, 0.193847413f, 0.19366445f, 0.193481488f, 0.193298525f, 0.193115563f, 0.1929326f, 0.192749638f,
        0.192566675f, 0.193328157f, 0.194089638f, 0.19485112f, 0.195612602f, 0.196374083f, 0.197135565f, 0.197897047f, 0.198658528f, 0.19942001f,
        0.200181492f, 0.201604331f, 0.20302717f, 0.204450009f, 0.205872848f, 0.207295688f, 0.208718527f, 0.210141366f, 0.211564205f, 0.212987044f,
        0.214409883f, 0.215920689f, 0.217431495f, 0.218942301f, 0.220453107f, 0.221963913f, 0.223474718f, 0.224985524f, 0.22649633f, 0.228007136f,
        0.229517942f, 0.230623808f, 0.231729675f, 0.232835542f, 0.233941408f, 0.235047275f, 0.236153142f, 0.237259008f, 0.238364875f, 0.239470742f,
        0.240576608f, 0.241915044f, 0.24325348f, 0.244591916f, 0.245930352f, 0.247268788f, 0.248607223f, 0.249945659f, 0.251284095f, 0.252622531f,
        0.253960967f, 0.256415994f, 0.258871022f, 0.261326049f, 0.263781077f, 0.266236104f, 0.268691132f, 0.271146159f, 0.273601187f, 0.276056214f,
        0.278511242f, 0.281982214f, 0.285453187f, 0.288924159f, 0.292395132f, 0.295866104f, 0.299337077f, 0.302808049f, 0.306279022f, 0.309749994f,
        0.313220967f, 0.316678026f, 0.320135085f, 0.323592144f, 0.327049203f, 0.330506263f, 0.333963322f, 0.337420381f, 0.34087744f, 0.344334499f,
        0.347791558f, 0.349599475f, 0.351407392f, 0.353215308f, 0.355023225f, 0.356831142f, 0.358639058f, 0.360446975f, 0.362254892f, 0.364062808f,
        0.365870725f, 0.36586278f, 0.365854835f, 0.36584689f, 0.365838945f, 0.365831f, 0.365823055f, 0.36581511f, 0.365807165f, 0.36579922f,
        0.365791275f, 0.36515456f, 0.364517845f, 0.36388113f, 0.363244415f, 0.3626077f, 0.361970985f, 0.36133427f, 0.360697555f, 0.36006084f,
        0.359424125f, 0.359281038f, 0.35913795f, 0.358994863f, 0.358851775f, 0.358708688f, 0.3585656f, 0.358422513f, 0.358279425f, 0.358136338f,
        0.35799325f, 0.35868661f, 0.35937997f, 0.36007333f, 0.36076669f, 0.36146005f, 0.36215341f, 0.36284677f, 0.36354013f, 0.36423349f,
        0.36492685f, 0.366156758f, 0.367386667f, 0.368616575f, 0.369846483f, 0.371076392f, 0.3723063f, 0.373536208f, 0.374766117f, 0.375996025f,
        0.377225933f, 0.379285967f, 0.381346f, 0.383406033f, 0.385466067f, 0.3875261f, 0.389586133f, 0.391646167f, 0.3937062f, 0.395766233f,
        0.397826267f, 0.35804364f, 0.318261013f, 0.278478387f, 0.23869576f, 0.198913133f, 0.159130507f, 0.11934788f, 0.079565253f, 0.039782627f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth06({0.136268951f, 0.140587868f, 0.144906785f, 0.149225702f, 0.153544619f, 0.157863536f, 0.162182452f, 0.166501369f, 0.170820286f, 0.175139203f,
        0.17945812f, 0.186201515f, 0.19294491f, 0.199688305f, 0.2064317f, 0.213175095f, 0.21991849f, 0.226661885f, 0.23340528f, 0.240148675f,
        0.24689207f, 0.251884762f, 0.256877454f, 0.261870146f, 0.266862839f, 0.271855531f, 0.276848223f, 0.281840915f, 0.286833607f, 0.291826299f,
        0.296818991f, 0.299165242f, 0.301511493f, 0.303857744f, 0.306203995f, 0.308550246f, 0.310896496f, 0.313242747f, 0.315588998f, 0.317935249f,
        0.3202815f, 0.321961795f, 0.323642089f, 0.325322384f, 0.327002679f, 0.328682974f, 0.330363268f, 0.332043563f, 0.333723858f, 0.335404153f,
        0.337084447f, 0.338925824f, 0.3407672f, 0.342608576f, 0.344449953f, 0.346291329f, 0.348132705f, 0.349974082f, 0.351815458f, 0.353656834f,
        0.355498211f, 0.358067014f, 0.360635818f, 0.363204621f, 0.365773425f, 0.368342228f, 0.370911032f, 0.373479835f, 0.376048639f, 0.378617442f,
        0.381186246f, 0.38498025f, 0.388774254f, 0.392568259f, 0.396362263f, 0.400156268f, 0.403950272f, 0.407744276f, 0.411538281f, 0.415332285f,
        0.419126289f, 0.423809552f, 0.428492814f, 0.433176076f, 0.437859339f, 0.442542601f, 0.447225863f, 0.451909125f, 0.456592388f, 0.46127565f,
        0.465958912f, 0.470411019f, 0.474863126f, 0.479315233f, 0.48376734f, 0.488219447f, 0.492671554f, 0.497123661f, 0.501575768f, 0.506027875f,
        0.510479982f, 0.514012641f, 0.5175453f, 0.521077959f, 0.524610618f, 0.528143276f, 0.531675935f, 0.535208594f, 0.538741253f, 0.542273911f,
        0.54580657f, 0.547945156f, 0.550083742f, 0.552222328f, 0.554360914f, 0.5564995f, 0.558638086f, 0.560776672f, 0.562915258f, 0.565053844f,
        0.56719243f, 0.567899658f, 0.568606886f, 0.569314114f, 0.570021342f, 0.57072857f, 0.571435798f, 0.572143026f, 0.572850254f, 0.573557482f,
        0.574264711f, 0.573746582f, 0.573228454f, 0.572710326f, 0.572192198f, 0.57167407f, 0.571155942f, 0.570637814f, 0.570119686f, 0.569601558f,
        0.56908343f, 0.567243478f, 0.565403526f, 0.563563575f, 0.561723623f, 0.559883671f, 0.558043719f, 0.556203768f, 0.554363816f, 0.552523864f,
        0.550683912f, 0.547966075f, 0.545248237f, 0.542530399f, 0.539812561f, 0.537094724f, 0.534376886f, 0.531659048f, 0.528941211f, 0.526223373f,
        0.523505535f, 0.519998335f, 0.516491135f, 0.512983935f, 0.509476735f, 0.505969535f, 0.502462335f, 0.498955135f, 0.495447935f, 0.491940735f,
        0.488433535f, 0.48411129f, 0.479789046f, 0.475466801f, 0.471144556f, 0.466822311f, 0.462500067f, 0.458177822f, 0.453855577f, 0.449533332f,
        0.445211088f, 0.440677259f, 0.43614343f, 0.431609601f, 0.427075772f, 0.422541943f, 0.418008114f, 0.413474285f, 0.408940456f, 0.404406627f,
        0.399872798f, 0.394928727f, 0.389984656f, 0.385040585f, 0.380096514f, 0.375152443f, 0.370208372f, 0.365264301f, 0.36032023f, 0.355376159f,
        0.350432088f, 0.345327729f, 0.34022337f, 0.335119011f, 0.330014653f, 0.324910294f, 0.319805935f, 0.314701576f, 0.309597218f, 0.304492859f,
        0.2993885f, 0.294693059f, 0.289997618f, 0.285302176f, 0.280606735f, 0.275911294f, 0.271215853f, 0.266520411f, 0.26182497f, 0.257129529f,
        0.252434088f, 0.249286877f, 0.246139667f, 0.242992456f, 0.239845246f, 0.236698035f, 0.233550825f, 0.230403614f, 0.227256404f, 0.224109193f,
        0.220961982f, 0.219297194f, 0.217632405f, 0.215967617f, 0.214302828f, 0.212638039f, 0.210973251f, 0.209308462f, 0.207643674f, 0.205978885f,
        0.204314096f, 0.20346194f, 0.202609784f, 0.201757628f, 0.200905472f, 0.200053316f, 0.19920116f, 0.198349004f, 0.197496847f, 0.196644691f,
        0.195792535f, 0.195300989f, 0.194809442f, 0.194317896f, 0.193826349f, 0.193334803f, 0.192843256f, 0.19235171f, 0.191860163f, 0.191368617f,
        0.19087707f, 0.190612275f, 0.190347479f, 0.190082683f, 0.189817888f, 0.189553092f, 0.189288296f, 0.189023501f, 0.188758705f, 0.18849391f,
        0.188229114f, 0.188478167f, 0.188727219f, 0.188976272f, 0.189225325f, 0.189474377f, 0.18972343f, 0.189972482f, 0.190221535f, 0.190470588f,
        0.19071964f, 0.191589889f, 0.192460137f, 0.193330385f, 0.194200633f, 0.195070882f, 0.19594113f, 0.196811378f, 0.197681626f, 0.198551875f,
        0.199422123f, 0.200638499f, 0.201854875f, 0.203071252f, 0.204287628f, 0.205504004f, 0.206720381f, 0.207936757f, 0.209153133f, 0.21036951f,
        0.211585886f, 0.212737445f, 0.213889004f, 0.215040562f, 0.216192121f, 0.21734368f, 0.218495239f, 0.219646797f, 0.220798356f, 0.221949915f,
        0.223101474f, 0.223955268f, 0.224809061f, 0.225662855f, 0.226516649f, 0.227370443f, 0.228224237f, 0.229078031f, 0.229931825f, 0.230785618f,
        0.231639412f, 0.231807176f, 0.23197494f, 0.232142705f, 0.232310469f, 0.232478233f, 0.232645997f, 0.232813761f, 0.232981525f, 0.233149289f,
        0.233317053f, 0.232926626f, 0.232536198f, 0.232145771f, 0.231755343f, 0.231364916f, 0.230974488f, 0.23058406f, 0.230193633f, 0.229803205f,
        0.229412778f, 0.229406829f, 0.22940088f, 0.229394931f, 0.229388982f, 0.229383033f, 0.229377084f, 0.229371136f, 0.229365187f, 0.229359238f,
        0.229353289f, 0.20641796f, 0.183482631f, 0.160547302f, 0.137611973f, 0.114676644f, 0.091741316f, 0.068805987f, 0.045870658f, 0.022935329f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth11({0.053807993f, 0.053796048f, 0.053784102f, 0.053772156f, 0.053760211f, 0.053748265f, 0.053736319f, 0.053724374f, 0.053712428f, 0.053700482f,
        0.053688537f, 0.053645508f, 0.053602479f, 0.053559451f, 0.053516422f, 0.053473393f, 0.053430365f, 0.053387336f, 0.053344307f, 0.053301279f,
        0.05325825f, 0.053302142f, 0.053346034f, 0.053389926f, 0.053433817f, 0.053477709f, 0.053521601f, 0.053565493f, 0.053609385f, 0.053653277f,
        0.053697168f, 0.053729406f, 0.053761644f, 0.053793882f, 0.05382612f, 0.053858358f, 0.053890596f, 0.053922834f, 0.053955072f, 0.05398731f,
        0.054019548f, 0.054069433f, 0.054119318f, 0.054169203f, 0.054219088f, 0.054268973f, 0.054318858f, 0.054368743f, 0.054418628f, 0.054468513f,
        0.054518398f, 0.054561244f, 0.05460409f, 0.054646936f, 0.054689782f, 0.054732628f, 0.054775473f, 0.054818319f, 0.054861165f, 0.054904011f,
        0.054946857f, 0.054967891f, 0.054988926f, 0.055009961f, 0.055030995f, 0.05505203f, 0.055073065f, 0.055094099f, 0.055115134f, 0.055136169f,
        0.055157203f, 0.055209092f, 0.05526098f, 0.055312868f, 0.055364756f, 0.055416644f, 0.055468532f, 0.055520421f, 0.055572309f, 0.055624197f,
        0.055676085f, 0.0557723f, 0.055868515f, 0.055964731f, 0.056060946f, 0.056157161f, 0.056253376f, 0.056349591f, 0.056445806f, 0.056542022f,
        0.056638237f, 0.056814154f, 0.056990071f, 0.057165988f, 0.057341905f, 0.057517822f, 0.057693739f, 0.057869656f, 0.058045573f, 0.05822149f,
        0.058397407f, 0.058679266f, 0.058961125f, 0.059242984f, 0.059524843f, 0.059806703f, 0.060088562f, 0.060370421f, 0.06065228f, 0.060934139f,
        0.061215998f, 0.061916158f, 0.062616317f, 0.063316477f, 0.064016636f, 0.064716796f, 0.065416955f, 0.066117115f, 0.066817274f, 0.067517434f,
        0.068217593f, 0.070337773f, 0.072457953f, 0.074578133f, 0.076698313f, 0.078818493f, 0.080938673f, 0.083058853f, 0.085179033f, 0.087299213f,
        0.089419393f, 0.092938903f, 0.096458412f, 0.099977921f, 0.10349743f, 0.107016939f, 0.110536448f, 0.114055958f, 0.117575467f, 0.121094976f,
        0.124614485f, 0.12750334f, 0.130392195f, 0.13328105f, 0.136169904f, 0.139058759f, 0.141947614f, 0.144836469f, 0.147725324f, 0.150614179f,
        0.153503033f, 0.155532077f, 0.15756112f, 0.159590163f, 0.161619207f, 0.16364825f, 0.165677293f, 0.167706337f, 0.16973538f, 0.171764423f,
        0.173793467f, 0.176358407f, 0.178923347f, 0.181488287f, 0.184053227f, 0.186618167f, 0.189183107f, 0.191748047f, 0.194312987f, 0.196877927f,
        0.199442867f, 0.204325613f, 0.20920836f, 0.214091107f, 0.218973853f, 0.2238566f, 0.228739347f, 0.233622093f, 0.23850484f, 0.243387587f,
        0.248270333f, 0.256985634f, 0.265700935f, 0.274416236f, 0.283131537f, 0.291846838f, 0.300562138f, 0.309277439f, 0.31799274f, 0.326708041f,
        0.335423342f, 0.34628038f, 0.357137418f, 0.367994457f, 0.378851495f, 0.389708533f, 0.400565572f, 0.41142261f, 0.422279648f, 0.433136687f,
        0.443993725f, 0.453441677f, 0.462889628f, 0.47233758f, 0.481785532f, 0.491233483f, 0.500681435f, 0.510129387f, 0.519577338f, 0.52902529f,
        0.538473242f, 0.543293256f, 0.54811327f, 0.552933284f, 0.557753298f, 0.562573313f, 0.567393327f, 0.572213341f, 0.577033355f, 0.581853369f,
        0.586673383f, 0.587489945f, 0.588306507f, 0.589123068f, 0.58993963f, 0.590756192f, 0.591572753f, 0.592389315f, 0.593205877f, 0.594022438f,
        0.594839f, 0.594414013f, 0.593989027f, 0.59356404f, 0.593139053f, 0.592714067f, 0.59228908f, 0.591864093f, 0.591439107f, 0.59101412f,
        0.590589133f, 0.590191936f, 0.589794738f, 0.589397541f, 0.589000343f, 0.588603146f, 0.588205948f, 0.587808751f, 0.587411553f, 0.587014356f,
        0.586617158f, 0.586372875f, 0.586128592f, 0.585884308f, 0.585640025f, 0.585395742f, 0.585151458f, 0.584907175f, 0.584662892f, 0.584418608f,
        0.584174325f, 0.584142693f, 0.58411106f, 0.584079428f, 0.584047795f, 0.584016163f, 0.58398453f, 0.583952898f, 0.583921265f, 0.583889633f,
        0.583858f, 0.584446924f, 0.585035848f, 0.585624773f, 0.586213697f, 0.586802621f, 0.587391545f, 0.587980469f, 0.588569393f, 0.589158318f,
        0.589747242f, 0.591023664f, 0.592300087f, 0.593576509f, 0.594852932f, 0.596129354f, 0.597405777f, 0.598682199f, 0.599958622f, 0.601235044f,
        0.602511467f, 0.604299346f, 0.606087225f, 0.607875104f, 0.609662983f, 0.611450863f, 0.613238742f, 0.615026621f, 0.6168145f, 0.618602379f,
        0.620390258f, 0.622230919f, 0.62407158f, 0.625912241f, 0.627752902f, 0.629593563f, 0.631434223f, 0.633274884f, 0.635115545f, 0.636956206f,
        0.638796867f, 0.640398172f, 0.641999477f, 0.643600782f, 0.645202087f, 0.646803392f, 0.648404697f, 0.650006002f, 0.651607307f, 0.653208612f,
        0.654809917f, 0.655583534f, 0.656357152f, 0.657130769f, 0.657904387f, 0.658678004f, 0.659451622f, 0.660225239f, 0.660998857f, 0.661772474f,
        0.662546092f, 0.662546326f, 0.66254656f, 0.662546794f, 0.662547028f, 0.662547263f, 0.662547497f, 0.662547731f, 0.662547965f, 0.662548199f,
        0.662548433f, 0.662975086f, 0.663401738f, 0.663828391f, 0.664255043f, 0.664681696f, 0.665108348f, 0.665535001f, 0.665961653f, 0.666388306f,
        0.666814958f, 0.600133463f, 0.533451967f, 0.466770471f, 0.400088975f, 0.333407479f, 0.266725983f, 0.200044488f, 0.133362992f, 0.066681496f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth12({0.122364925f, 0.126576671f, 0.130788417f, 0.135000163f, 0.139211908f, 0.143423654f, 0.1476354f, 0.151847146f, 0.156058892f, 0.160270638f,
        0.164482383f, 0.170884575f, 0.177286767f, 0.183688958f, 0.19009115f, 0.196493342f, 0.202895533f, 0.209297725f, 0.215699917f, 0.222102108f,
        0.2285043f, 0.234261801f, 0.240019302f, 0.245776803f, 0.251534303f, 0.257291804f, 0.263049305f, 0.268806806f, 0.274564307f, 0.280321808f,
        0.286079308f, 0.290201841f, 0.294324373f, 0.298446906f, 0.302569438f, 0.306691971f, 0.310814503f, 0.314937036f, 0.319059568f, 0.323182101f,
        0.327304633f, 0.330682438f, 0.334060242f, 0.337438046f, 0.34081585f, 0.344193654f, 0.347571458f, 0.350949263f, 0.354327067f, 0.357704871f,
        0.361082675f, 0.363731076f, 0.366379477f, 0.369027878f, 0.371676278f, 0.374324679f, 0.37697308f, 0.379621481f, 0.382269882f, 0.384918283f,
        0.387566683f, 0.388773437f, 0.38998019f, 0.391186943f, 0.392393697f, 0.39360045f, 0.394807203f, 0.396013957f, 0.39722071f, 0.398427463f,
        0.399634217f, 0.398827703f, 0.398021188f, 0.397214674f, 0.39640816f, 0.395601646f, 0.394795132f, 0.393988618f, 0.393182103f, 0.392375589f,
        0.391569075f, 0.388654937f, 0.385740798f, 0.38282666f, 0.379912522f, 0.376998383f, 0.374084245f, 0.371170107f, 0.368255968f, 0.36534183f,
        0.362427692f, 0.357796759f, 0.353165827f, 0.348534894f, 0.343903962f, 0.339273029f, 0.334642097f, 0.330011164f, 0.325380232f, 0.320749299f,
        0.316118367f, 0.310530894f, 0.304943422f, 0.299355949f, 0.293768477f, 0.288181004f, 0.282593532f, 0.277006059f, 0.271418587f, 0.265831114f,
        0.260243642f, 0.255077422f, 0.249911202f, 0.244744982f, 0.239578762f, 0.234412542f, 0.229246322f, 0.224080102f, 0.218913882f, 0.213747662f,
        0.208581442f, 0.204554178f, 0.200526913f, 0.196499649f, 0.192472385f, 0.188445121f, 0.184417857f, 0.180390593f, 0.176363328f, 0.172336064f,
        0.1683088f, 0.165246354f, 0.162183908f, 0.159121463f, 0.156059017f, 0.152996571f, 0.149934125f, 0.146871679f, 0.143809233f, 0.140746788f,
        0.137684342f, 0.135571868f, 0.133459393f, 0.131346919f, 0.129234445f, 0.127121971f, 0.125009497f, 0.122897023f, 0.120784548f, 0.118672074f,
        0.1165596f, 0.115328577f, 0.114097553f, 0.11286653f, 0.111635507f, 0.110404483f, 0.10917346f, 0.107942437f, 0.106711413f, 0.10548039f,
        0.104249367f, 0.103460951f, 0.102672535f, 0.10188412f, 0.101095704f, 0.100307288f, 0.099518873f, 0.098730457f, 0.097942041f, 0.097153626f,
        0.09636521f, 0.095708316f, 0.095051423f, 0.094394529f, 0.093737635f, 0.093080742f, 0.092423848f, 0.091766954f, 0.091110061f, 0.090453167f,
        0.089796273f, 0.08936739f, 0.088938507f, 0.088509623f, 0.08808074f, 0.087651857f, 0.087222973f, 0.08679409f, 0.086365207f, 0.085936323f,
        0.08550744f, 0.085329032f, 0.085150623f, 0.084972215f, 0.084793807f, 0.084615398f, 0.08443699f, 0.084258582f, 0.084080173f, 0.083901765f,
        0.083723357f, 0.083747411f, 0.083771464f, 0.083795518f, 0.083819572f, 0.083843626f, 0.08386768f, 0.083891734f, 0.083915787f, 0.083939841f,
        0.083963895f, 0.083999012f, 0.084034128f, 0.084069245f, 0.084104362f, 0.084139478f, 0.084174595f, 0.084209712f, 0.084244828f, 0.084279945f,
        0.084315062f, 0.084294631f, 0.084274199f, 0.084253768f, 0.084233337f, 0.084212906f, 0.084192475f, 0.084172044f, 0.084151612f, 0.084131181f,
        0.08411075f, 0.084085745f, 0.08406074f, 0.084035736f, 0.084010731f, 0.083985726f, 0.083960721f, 0.083935716f, 0.083910711f, 0.083885707f,
        0.083860702f, 0.083991644f, 0.084122585f, 0.084253527f, 0.084384469f, 0.084515411f, 0.084646353f, 0.084777295f, 0.084908236f, 0.085039178f,
        0.08517012f, 0.085629799f, 0.086089477f, 0.086549156f, 0.087008834f, 0.087468513f, 0.087928191f, 0.08838787f, 0.088847548f, 0.089307227f,
        0.089766905f, 0.090575043f, 0.091383181f, 0.092191319f, 0.092999456f, 0.093807594f, 0.094615732f, 0.09542387f, 0.096232008f, 0.097040146f,
        0.097848283f, 0.098975833f, 0.100103383f, 0.101230933f, 0.102358483f, 0.103486033f, 0.104613583f, 0.105741133f, 0.106868683f, 0.107996233f,
        0.109123783f, 0.110557826f, 0.111991868f, 0.113425911f, 0.114859953f, 0.116293996f, 0.117728038f, 0.119162081f, 0.120596123f, 0.122030166f,
        0.123464208f, 0.125387027f, 0.127309845f, 0.129232663f, 0.131155482f, 0.1330783f, 0.135001118f, 0.136923937f, 0.138846755f, 0.140769573f,
        0.142692392f, 0.145353644f, 0.148014897f, 0.150676149f, 0.153337402f, 0.155998654f, 0.158659907f, 0.161321159f, 0.163982412f, 0.166643664f,
        0.169304917f, 0.172839598f, 0.176374278f, 0.179908959f, 0.18344364f, 0.186978321f, 0.190513002f, 0.194047683f, 0.197582363f, 0.201117044f,
        0.204651725f, 0.208581847f, 0.212511968f, 0.21644209f, 0.220372212f, 0.224302333f, 0.228232455f, 0.232162577f, 0.236092698f, 0.24002282f,
        0.243952942f, 0.248276378f, 0.252599813f, 0.256923249f, 0.261246685f, 0.265570121f, 0.269893557f, 0.274216993f, 0.278540428f, 0.282863864f,
        0.2871873f, 0.291717082f, 0.296246863f, 0.300776645f, 0.305306427f, 0.309836208f, 0.31436599f, 0.318895772f, 0.323425553f, 0.327955335f,
        0.332485117f, 0.299236605f, 0.265988093f, 0.232739582f, 0.19949107f, 0.166242558f, 0.132994047f, 0.099745535f, 0.066497023f, 0.033248512f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth13({0.096004948f, 0.097870687f, 0.099736425f, 0.101602164f, 0.103467902f, 0.105333641f, 0.107199379f, 0.109065118f, 0.110930856f, 0.112796595f,
        0.114662333f, 0.116254343f, 0.117846353f, 0.119438363f, 0.121030373f, 0.122622383f, 0.124214393f, 0.125806403f, 0.127398413f, 0.128990423f,
        0.130582433f, 0.131032368f, 0.131482303f, 0.131932238f, 0.132382173f, 0.132832108f, 0.133282043f, 0.133731978f, 0.134181913f, 0.134631848f,
        0.135081783f, 0.134918414f, 0.134755045f, 0.134591676f, 0.134428307f, 0.134264938f, 0.134101568f, 0.133938199f, 0.13377483f, 0.133611461f,
        0.133448092f, 0.13326246f, 0.133076828f, 0.132891197f, 0.132705565f, 0.132519933f, 0.132334302f, 0.13214867f, 0.131963038f, 0.131777407f,
        0.131591775f, 0.131453653f, 0.13131553f, 0.131177408f, 0.131039285f, 0.130901163f, 0.13076304f, 0.130624918f, 0.130486795f, 0.130348673f,
        0.13021055f, 0.130000184f, 0.129789818f, 0.129579453f, 0.129369087f, 0.129158721f, 0.128948355f, 0.128737989f, 0.128527623f, 0.128317258f,
        0.128106892f, 0.127801283f, 0.127495673f, 0.127190064f, 0.126884455f, 0.126578846f, 0.126273237f, 0.125967628f, 0.125662018f, 0.125356409f,
        0.1250508f, 0.124593303f, 0.124135807f, 0.12367831f, 0.123220813f, 0.122763317f, 0.12230582f, 0.121848323f, 0.121390827f, 0.12093333f,
        0.120475833f, 0.119940628f, 0.119405423f, 0.118870218f, 0.118335013f, 0.117799808f, 0.117264603f, 0.116729398f, 0.116194193f, 0.115658988f,
        0.115123783f, 0.11459622f, 0.114068657f, 0.113541093f, 0.11301353f, 0.112485967f, 0.111958403f, 0.11143084f, 0.110903277f, 0.110375713f,
        0.10984815f, 0.109357321f, 0.108866492f, 0.108375663f, 0.107884833f, 0.107394004f, 0.106903175f, 0.106412346f, 0.105921517f, 0.105430688f,
        0.104939858f, 0.104427923f, 0.103915987f, 0.103404051f, 0.102892116f, 0.10238018f, 0.101868244f, 0.101356309f, 0.100844373f, 0.100332437f,
        0.099820502f, 0.099354759f, 0.098889016f, 0.098423273f, 0.09795753f, 0.097491787f, 0.097026044f, 0.096560301f, 0.096094558f, 0.095628815f,
        0.095163072f, 0.094912235f, 0.094661398f, 0.094410561f, 0.094159724f, 0.093908888f, 0.093658051f, 0.093407214f, 0.093156377f, 0.09290554f,
        0.092654703f, 0.09263596f, 0.092617217f, 0.092598474f, 0.092579731f, 0.092560988f, 0.092542244f, 0.092523501f, 0.092504758f, 0.092486015f,
        0.092467272f, 0.092539468f, 0.092611665f, 0.092683862f, 0.092756058f, 0.092828255f, 0.092900452f, 0.092972648f, 0.093044845f, 0.093117042f,
        0.093189238f, 0.093490886f, 0.093792533f, 0.094094181f, 0.094395828f, 0.094697476f, 0.094999123f, 0.095300771f, 0.095602418f, 0.095904066f,
        0.096205713f, 0.097396776f, 0.098587839f, 0.099778902f, 0.100969965f, 0.102161028f, 0.10335209f, 0.104543153f, 0.105734216f, 0.106925279f,
        0.108116342f, 0.112861428f, 0.117606515f, 0.122351602f, 0.127096688f, 0.131841775f, 0.136586862f, 0.141331948f, 0.146077035f, 0.150822122f,
        0.155567208f, 0.166549495f, 0.177531782f, 0.188514068f, 0.199496355f, 0.210478642f, 0.221460928f, 0.232443215f, 0.243425502f, 0.254407788f,
        0.265390075f, 0.278722221f, 0.292054367f, 0.305386513f, 0.318718658f, 0.332050804f, 0.34538295f, 0.358715096f, 0.372047242f, 0.385379388f,
        0.398711533f, 0.408848528f, 0.418985522f, 0.429122516f, 0.43925951f, 0.449396504f, 0.459533498f, 0.469670493f, 0.479807487f, 0.489944481f,
        0.500081475f, 0.505705238f, 0.511329f, 0.516952763f, 0.522576525f, 0.528200288f, 0.53382405f, 0.539447813f, 0.545071575f, 0.550695338f,
        0.5563191f, 0.558632255f, 0.56094541f, 0.563258565f, 0.56557172f, 0.567884875f, 0.57019803f, 0.572511185f, 0.57482434f, 0.577137495f,
        0.57945065f, 0.580278085f, 0.58110552f, 0.581932955f, 0.58276039f, 0.583587825f, 0.58441526f, 0.585242695f, 0.58607013f, 0.586897565f,
        0.587725f, 0.588015169f, 0.588305338f, 0.588595508f, 0.588885677f, 0.589175846f, 0.589466015f, 0.589756184f, 0.590046353f, 0.590336523f,
        0.590626692f, 0.590815022f, 0.591003352f, 0.591191682f, 0.591380012f, 0.591568342f, 0.591756672f, 0.591945002f, 0.592133332f, 0.592321662f,
        0.592509992f, 0.592704069f, 0.592898147f, 0.593092224f, 0.593286302f, 0.593480379f, 0.593674457f, 0.593868534f, 0.594062612f, 0.594256689f,
        0.594450767f, 0.594790606f, 0.595130445f, 0.595470284f, 0.595810123f, 0.596149963f, 0.596489802f, 0.596829641f, 0.59716948f, 0.597509319f,
        0.597849158f, 0.598283257f, 0.598717355f, 0.599151453f, 0.599585552f, 0.60001965f, 0.600453748f, 0.600887847f, 0.601321945f, 0.601756043f,
        0.602190142f, 0.602661484f, 0.603132827f, 0.603604169f, 0.604075512f, 0.604546854f, 0.605018197f, 0.605489539f, 0.605960882f, 0.606432224f,
        0.606903567f, 0.607138388f, 0.607373208f, 0.607608029f, 0.60784285f, 0.608077671f, 0.608312492f, 0.608547313f, 0.608782133f, 0.609016954f,
        0.609251775f, 0.609222884f, 0.609193993f, 0.609165103f, 0.609136212f, 0.609107321f, 0.60907843f, 0.609049539f, 0.609020648f, 0.608991758f,
        0.608962867f, 0.609090122f, 0.609217377f, 0.609344632f, 0.609471887f, 0.609599142f, 0.609726397f, 0.609853652f, 0.609980907f, 0.610108162f,
        0.610235417f, 0.549211875f, 0.488188333f, 0.427164792f, 0.36614125f, 0.305117708f, 0.244094167f, 0.183070625f, 0.122047083f, 0.061023542f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    static constexpr Spectrum Macbeth14({0.091993695f, 0.094395723f, 0.096797751f, 0.099199779f, 0.101601807f, 0.104003835f, 0.106405863f, 0.108807891f, 0.111209919f, 0.113611947f,
        0.116013975f, 0.118973427f, 0.121932878f, 0.12489233f, 0.127851782f, 0.130811233f, 0.133770685f, 0.136730137f, 0.139689588f, 0.14264904f,
        0.145608492f, 0.147901087f, 0.150193682f, 0.152486277f, 0.154778872f, 0.157071467f, 0.159364062f, 0.161656657f, 0.163949252f, 0.166241847f,
        0.168534442f, 0.169527616f, 0.17052079f, 0.171513964f, 0.172507138f, 0.173500313f, 0.174493487f, 0.175486661f, 0.176479835f, 0.177473009f,
        0.178466183f, 0.177920935f, 0.177375687f, 0.176830438f, 0.17628519f, 0.175739942f, 0.175194693f, 0.174649445f, 0.174104197f, 0.173558948f,
        0.1730137f, 0.171508923f, 0.170004147f, 0.16849937f, 0.166994593f, 0.165489817f, 0.16398504f, 0.162480263f, 0.160975487f, 0.15947071f,
        0.157965933f, 0.156047798f, 0.154129663f, 0.152211528f, 0.150293393f, 0.148375258f, 0.146457123f, 0.144538988f, 0.142620853f, 0.140702718f,
        0.138784583f, 0.136819271f, 0.134853958f, 0.132888646f, 0.130923333f, 0.128958021f, 0.126992708f, 0.125027396f, 0.123062083f, 0.121096771f,
        0.119131458f, 0.117357856f, 0.115584253f, 0.113810651f, 0.112037048f, 0.110263446f, 0.108489843f, 0.106716241f, 0.104942638f, 0.103169036f,
        0.101395433f, 0.099950526f, 0.098505619f, 0.097060712f, 0.095615805f, 0.094170898f, 0.09272599f, 0.091281083f, 0.089836176f, 0.088391269f,
        0.086946362f, 0.085769857f, 0.084593351f, 0.083416846f, 0.082240341f, 0.081063836f, 0.079887331f, 0.078710826f, 0.07753432f, 0.076357815f,
        0.07518131f, 0.074272676f, 0.073364041f, 0.072455407f, 0.071546773f, 0.070638138f, 0.069729504f, 0.06882087f, 0.067912235f, 0.067003601f,
        0.066094967f, 0.065517321f, 0.064939675f, 0.064362029f, 0.063784383f, 0.063206737f, 0.062629091f, 0.062051445f, 0.061473799f, 0.060896153f,
        0.060318507f, 0.059932656f, 0.059546804f, 0.059160953f, 0.058775102f, 0.058389251f, 0.0580034f, 0.057617549f, 0.057231697f, 0.056845846f,
        0.056459995f, 0.056125823f, 0.055791651f, 0.05545748f, 0.055123308f, 0.054789136f, 0.054454964f, 0.054120792f, 0.05378662f, 0.053452449f,
        0.053118277f, 0.052927354f, 0.052736432f, 0.05254551f, 0.052354587f, 0.052163665f, 0.051972743f, 0.05178182f, 0.051590898f, 0.051399976f,
        0.051209053f, 0.05121224f, 0.051215426f, 0.051218612f, 0.051221798f, 0.051224984f, 0.05122817f, 0.051231357f, 0.051234543f, 0.051237729f,
        0.051240915f, 0.051311842f, 0.051382769f, 0.051453697f, 0.051524624f, 0.051595551f, 0.051666478f, 0.051737405f, 0.051808332f, 0.05187926f,
        0.051950187f, 0.05194242f, 0.051934653f, 0.051926886f, 0.051919119f, 0.051911352f, 0.051903585f, 0.051895818f, 0.051888051f, 0.051880284f,
        0.051872517f, 0.051805078f, 0.051737639f, 0.051670201f, 0.051602762f, 0.051535323f, 0.051467885f, 0.051400446f, 0.051333007f, 0.051265569f,
        0.05119813f, 0.051320318f, 0.051442506f, 0.051564695f, 0.051686883f, 0.051809071f, 0.051931259f, 0.052053447f, 0.052175635f, 0.052297824f,
        0.052420012f, 0.053019264f, 0.053618516f, 0.054217768f, 0.05481702f, 0.055416273f, 0.056015525f, 0.056614777f, 0.057214029f, 0.057813281f,
        0.058412533f, 0.0598892f, 0.061365867f, 0.062842534f, 0.064319201f, 0.065795868f, 0.067272535f, 0.068749202f, 0.070225869f, 0.071702536f,
        0.073179203f, 0.075412842f, 0.077646481f, 0.07988012f, 0.082113759f, 0.084347398f, 0.086581037f, 0.088814676f, 0.091048315f, 0.093281954f,
        0.095515593f, 0.097856567f, 0.100197541f, 0.102538515f, 0.104879489f, 0.107220463f, 0.109561437f, 0.111902411f, 0.114243385f, 0.116584359f,
        0.118925333f, 0.121172248f, 0.123419163f, 0.125666078f, 0.127912993f, 0.130159908f, 0.132406823f, 0.134653738f, 0.136900653f, 0.139147568f,
        0.141394483f, 0.143809122f, 0.14622376f, 0.148638398f, 0.151053037f, 0.153467675f, 0.155882313f, 0.158296952f, 0.16071159f, 0.163126228f,
        0.165540867f, 0.168391362f, 0.171241857f, 0.174092352f, 0.176942847f, 0.179793342f, 0.182643837f, 0.185494332f, 0.188344827f, 0.191195322f,
        0.194045817f, 0.197347446f, 0.200649075f, 0.203950704f, 0.207252333f, 0.210553963f, 0.213855592f, 0.217157221f, 0.22045885f, 0.223760479f,
        0.227062108f, 0.230895335f, 0.234728562f, 0.238561788f, 0.242395015f, 0.246228242f, 0.250061468f, 0.253894695f, 0.257727922f, 0.261561148f,
        0.265394375f, 0.269747375f, 0.274100375f, 0.278453375f, 0.282806375f, 0.287159375f, 0.291512375f, 0.295865375f, 0.300218375f, 0.304571375f,
        0.308924375f, 0.313487192f, 0.318050008f, 0.322612825f, 0.327175642f, 0.331738458f, 0.336301275f, 0.340864092f, 0.345426908f, 0.349989725f,
        0.354552542f, 0.358673831f, 0.36279512f, 0.366916409f, 0.371037698f, 0.375158988f, 0.379280277f, 0.383401566f, 0.387522855f, 0.391644144f,
        0.395765433f, 0.399772599f, 0.403779765f, 0.407786931f, 0.411794097f, 0.415801263f, 0.419808428f, 0.423815594f, 0.42782276f, 0.431829926f,
        0.435837092f, 0.440100286f, 0.44436348f, 0.448626674f, 0.452889868f, 0.457153063f, 0.461416257f, 0.465679451f, 0.469942645f, 0.474205839f,
        0.478469033f, 0.43062213f, 0.382775227f, 0.334928323f, 0.28708142f, 0.239234517f, 0.191387613f, 0.14354071f, 0.095693807f, 0.047846903f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth15({0.061031217f, 0.061052864f, 0.061074512f, 0.061096159f, 0.061117807f, 0.061139454f, 0.061161102f, 0.061182749f, 0.061204397f, 0.061226044f,
        0.061247692f, 0.061315304f, 0.061382916f, 0.061450528f, 0.06151814f, 0.061585752f, 0.061653364f, 0.061720976f, 0.061788588f, 0.0618562f,
        0.061923812f, 0.062022906f, 0.062122f, 0.062221094f, 0.062320188f, 0.062419282f, 0.062518376f, 0.06261747f, 0.062716564f, 0.062815658f,
        0.062914752f, 0.063020304f, 0.063125856f, 0.063231408f, 0.06333696f, 0.063442512f, 0.063548064f, 0.063653616f, 0.063759168f, 0.06386472f,
        0.063970272f, 0.064166072f, 0.064361873f, 0.064557674f, 0.064753474f, 0.064949275f, 0.065145076f, 0.065340876f, 0.065536677f, 0.065732478f,
        0.065928278f, 0.066256456f, 0.066584634f, 0.066912811f, 0.067240989f, 0.067569167f, 0.067897344f, 0.068225522f, 0.0685537f, 0.068881877f,
        0.069210055f, 0.069761582f, 0.07031311f, 0.070864637f, 0.071416164f, 0.071967692f, 0.072519219f, 0.073070746f, 0.073622274f, 0.074173801f,
        0.074725328f, 0.075801767f, 0.076878206f, 0.077954645f, 0.079031084f, 0.080107523f, 0.081183961f, 0.0822604f, 0.083336839f, 0.084413278f,
        0.085489717f, 0.087446972f, 0.089404227f, 0.091361482f, 0.093318737f, 0.095275992f, 0.097233247f, 0.099190502f, 0.101147757f, 0.103105012f,
        0.105062267f, 0.108422853f, 0.111783438f, 0.115144024f, 0.11850461f, 0.121865196f, 0.125225782f, 0.128586368f, 0.131946953f, 0.135307539f,
        0.138668125f, 0.144010023f, 0.14935192f, 0.154693818f, 0.160035715f, 0.165377613f, 0.17071951f, 0.176061408f, 0.181403305f, 0.186745203f,
        0.1920871f, 0.199951233f, 0.207815367f, 0.2156795f, 0.223543633f, 0.231407767f, 0.2392719f, 0.247136033f, 0.255000167f, 0.2628643f,
        0.270728433f, 0.281266263f, 0.291804092f, 0.302341921f, 0.31287975f, 0.323417579f, 0.333955408f, 0.344493238f, 0.355031067f, 0.365568896f,
        0.376106725f, 0.386073688f, 0.39604065f, 0.406007613f, 0.415974575f, 0.425941538f, 0.4359085f, 0.445875463f, 0.455842425f, 0.465809388f,
        0.47577635f, 0.481321102f, 0.486865853f, 0.492410605f, 0.497955357f, 0.503500108f, 0.50904486f, 0.514589612f, 0.520134363f, 0.525679115f,
        0.531223867f, 0.533017229f, 0.534810592f, 0.536603954f, 0.538397317f, 0.540190679f, 0.541984042f, 0.543777404f, 0.545570767f, 0.547364129f,
        0.549157492f, 0.548813237f, 0.548468982f, 0.548124727f, 0.547780472f, 0.547436217f, 0.547091962f, 0.546747707f, 0.546403452f, 0.546059197f,
        0.545714942f, 0.543950689f, 0.542186437f, 0.540422184f, 0.538657932f, 0.536893679f, 0.535129427f, 0.533365174f, 0.531600922f, 0.529836669f,
        0.528072417f, 0.525710778f, 0.523349138f, 0.520987499f, 0.51862586f, 0.516264221f, 0.513902582f, 0.511540943f, 0.509179303f, 0.506817664f,
        0.504456025f, 0.501062431f, 0.497668837f, 0.494275243f, 0.490881648f, 0.487488054f, 0.48409446f, 0.480700866f, 0.477307272f, 0.473913678f,
        0.470520083f, 0.466231866f, 0.461943648f, 0.457655431f, 0.453367213f, 0.449078996f, 0.444790778f, 0.440502561f, 0.436214343f, 0.431926126f,
        0.427637908f, 0.422998625f, 0.418359342f, 0.413720058f, 0.409080775f, 0.404441492f, 0.399802208f, 0.395162925f, 0.390523642f, 0.385884358f,
        0.381245075f, 0.377800438f, 0.374355802f, 0.370911165f, 0.367466528f, 0.364021892f, 0.360577255f, 0.357132618f, 0.353687982f, 0.350243345f,
        0.346798708f, 0.344863173f, 0.342927637f, 0.340992101f, 0.339056565f, 0.337121029f, 0.335185493f, 0.333249958f, 0.331314422f, 0.329378886f,
        0.32744335f, 0.326469735f, 0.32549612f, 0.324522505f, 0.32354889f, 0.322575275f, 0.32160166f, 0.320628045f, 0.31965443f, 0.318680815f,
        0.3177072f, 0.317183179f, 0.316659158f, 0.316135138f, 0.315611117f, 0.315087096f, 0.314563075f, 0.314039054f, 0.313515033f, 0.312991013f,
        0.312466992f, 0.312213849f, 0.311960707f, 0.311707564f, 0.311454422f, 0.311201279f, 0.310948137f, 0.310694994f, 0.310441852f, 0.310188709f,
        0.309935567f, 0.31038308f, 0.310830593f, 0.311278107f, 0.31172562f, 0.312173133f, 0.312620647f, 0.31306816f, 0.313515673f, 0.313963187f,
        0.3144107f, 0.315710298f, 0.317009897f, 0.318309495f, 0.319609093f, 0.320908692f, 0.32220829f, 0.323507888f, 0.324807487f, 0.326107085f,
        0.327406683f, 0.329188613f, 0.330970543f, 0.332752473f, 0.334534403f, 0.336316333f, 0.338098263f, 0.339880193f, 0.341662123f, 0.343444053f,
        0.345225983f, 0.346958421f, 0.348690858f, 0.350423296f, 0.352155733f, 0.353888171f, 0.355620608f, 0.357353046f, 0.359085483f, 0.360817921f,
        0.362550358f, 0.363916925f, 0.365283492f, 0.366650058f, 0.368016625f, 0.369383192f, 0.370749758f, 0.372116325f, 0.373482892f, 0.374849458f,
        0.376216025f, 0.376648712f, 0.377081398f, 0.377514085f, 0.377946772f, 0.378379458f, 0.378812145f, 0.379244832f, 0.379677518f, 0.380110205f,
        0.380542892f, 0.380256029f, 0.379969167f, 0.379682304f, 0.379395442f, 0.379108579f, 0.378821717f, 0.378534854f, 0.378247992f, 0.377961129f,
        0.377674267f, 0.377847467f, 0.378020667f, 0.378193867f, 0.378367067f, 0.378540267f, 0.378713467f, 0.378886667f, 0.379059867f, 0.379233067f,
        0.379406267f, 0.34146564f, 0.303525013f, 0.265584387f, 0.22764376f, 0.189703133f, 0.151762507f, 0.11382188f, 0.075881253f, 0.037940627f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth16({0.06281745f, 0.062819804f, 0.062822157f, 0.062824511f, 0.062826865f, 0.062829218f, 0.062831572f, 0.062833926f, 0.062836279f, 0.062838633f,
        0.062840987f, 0.062890569f, 0.062940151f, 0.062989733f, 0.063039315f, 0.063088897f, 0.063138479f, 0.063188061f, 0.063237643f, 0.063287225f,
        0.063336807f, 0.063357351f, 0.063377895f, 0.063398439f, 0.063418983f, 0.063439528f, 0.063460072f, 0.063480616f, 0.06350116f, 0.063521704f,
        0.063542248f, 0.063558928f, 0.063575608f, 0.063592288f, 0.063608968f, 0.063625648f, 0.063642327f, 0.063659007f, 0.063675687f, 0.063692367f,
        0.063709047f, 0.063779971f, 0.063850896f, 0.063921821f, 0.063992745f, 0.06406367f, 0.064134595f, 0.064205519f, 0.064276444f, 0.064347369f,
        0.064418293f, 0.064512504f, 0.064606714f, 0.064700924f, 0.064795135f, 0.064889345f, 0.064983555f, 0.065077766f, 0.065171976f, 0.065266186f,
        0.065360397f, 0.065423373f, 0.065486349f, 0.065549325f, 0.065612301f, 0.065675277f, 0.065738253f, 0.065801229f, 0.065864205f, 0.065927181f,
        0.065990157f, 0.066084808f, 0.066179459f, 0.066274111f, 0.066368762f, 0.066463413f, 0.066558065f, 0.066652716f, 0.066747367f, 0.066842019f,
        0.06693667f, 0.06708427f, 0.06723187f, 0.06737947f, 0.06752707f, 0.06767467f, 0.06782227f, 0.06796987f, 0.06811747f, 0.06826507f, 0.06841267f,
        0.068699645f, 0.068986619f, 0.069273594f, 0.069560569f, 0.069847543f, 0.070134518f, 0.070421493f, 0.070708467f, 0.070995442f, 0.071282417f,
        0.071724769f, 0.072167121f, 0.072609473f, 0.073051825f, 0.073494177f, 0.073936529f, 0.074378881f, 0.074821233f, 0.075263585f, 0.075705937f,
        0.076857519f, 0.078009102f, 0.079160684f, 0.080312267f, 0.081463849f, 0.082615432f, 0.083767014f, 0.084918597f, 0.086070179f, 0.087221762f,
        0.091030431f, 0.094839101f, 0.098647771f, 0.10245644f, 0.10626511f, 0.11007378f, 0.113882449f, 0.117691119f, 0.121499789f, 0.125308458f,
        0.133360902f, 0.141413345f, 0.149465788f, 0.157518232f, 0.165570675f, 0.173623118f, 0.181675562f, 0.189728005f, 0.197780448f, 0.205832892f,
        0.215775898f, 0.225718903f, 0.235661909f, 0.245604915f, 0.255547921f, 0.265490927f, 0.275433933f, 0.285376938f, 0.295319944f, 0.30526295f,
        0.31305211f, 0.32084127f, 0.32863043f, 0.33641959f, 0.34420875f, 0.35199791f, 0.35978707f, 0.36757623f, 0.37536539f, 0.38315455f, 0.387933383f,
        0.392712215f, 0.397491048f, 0.40226988f, 0.407048713f, 0.411827545f, 0.416606378f, 0.42138521f, 0.426164043f, 0.430942875f, 0.434763645f,
        0.438584415f, 0.442405185f, 0.446225955f, 0.450046725f, 0.453867495f, 0.457688265f, 0.461509035f, 0.465329805f, 0.469150575f, 0.474024545f,
        0.478898515f, 0.483772485f, 0.488646455f, 0.493520425f, 0.498394395f, 0.503268365f, 0.508142335f, 0.513016305f, 0.517890275f, 0.522894058f,
        0.52789784f, 0.532901623f, 0.537905405f, 0.542909188f, 0.54791297f, 0.552916753f, 0.557920535f, 0.562924318f, 0.5679281f, 0.571823276f,
        0.575718452f, 0.579613628f, 0.583508803f, 0.587403979f, 0.591299155f, 0.595194331f, 0.599089507f, 0.602984683f, 0.606879858f, 0.608996518f,
        0.611113177f, 0.613229836f, 0.615346495f, 0.617463154f, 0.619579813f, 0.621696473f, 0.623813132f, 0.625929791f, 0.62804645f, 0.628944887f,
        0.629843323f, 0.63074176f, 0.631640197f, 0.632538633f, 0.63343707f, 0.634335507f, 0.635233943f, 0.63613238f, 0.637030817f, 0.637326863f,
        0.637622908f, 0.637918954f, 0.638215f, 0.638511046f, 0.638807092f, 0.639103138f, 0.639399183f, 0.639695229f, 0.639991275f, 0.640190068f,
        0.640388862f, 0.640587655f, 0.640786448f, 0.640985242f, 0.641184035f, 0.641382828f, 0.641581622f, 0.641780415f, 0.641979208f, 0.642326169f,
        0.64267313f, 0.643020091f, 0.643367052f, 0.643714013f, 0.644060973f, 0.644407934f, 0.644754895f, 0.645101856f, 0.645448817f, 0.645727897f,
        0.646006977f, 0.646286057f, 0.646565137f, 0.646844217f, 0.647123297f, 0.647402377f, 0.647681457f, 0.647960537f, 0.648239617f, 0.648517497f,
        0.648795377f, 0.649073257f, 0.649351137f, 0.649629017f, 0.649906897f, 0.650184777f, 0.650462657f, 0.650740537f, 0.651018417f, 0.651223242f,
        0.651428067f, 0.651632892f, 0.651837717f, 0.652042542f, 0.652247367f, 0.652452192f, 0.652657017f, 0.652861842f, 0.653066667f, 0.65349558f,
        0.653924493f, 0.654353407f, 0.65478232f, 0.655211233f, 0.655640147f, 0.65606906f, 0.656497973f, 0.656926887f, 0.6573558f, 0.658023494f,
        0.658691188f, 0.659358883f, 0.660026577f, 0.660694271f, 0.661361965f, 0.662029659f, 0.662697353f, 0.663365048f, 0.664032742f, 0.66489469f,
        0.665756638f, 0.666618587f, 0.667480535f, 0.668342483f, 0.669204432f, 0.67006638f, 0.670928328f, 0.671790277f, 0.672652225f, 0.67335715f,
        0.674062075f, 0.674767f, 0.675471925f, 0.67617685f, 0.676881775f, 0.6775867f, 0.678291625f, 0.67899655f, 0.679701475f, 0.680107636f,
        0.680513797f, 0.680919958f, 0.681326118f, 0.681732279f, 0.68213844f, 0.682544601f, 0.682950762f, 0.683356923f, 0.683763083f, 0.684215498f,
        0.684667913f, 0.685120328f, 0.685572743f, 0.686025158f, 0.686477573f, 0.686929988f, 0.687382403f, 0.687834818f, 0.688287233f, 0.61945851f,
        0.550629787f, 0.481801063f, 0.41297234f, 0.344143617f, 0.275314893f, 0.20648617f, 0.137657447f, 0.068828723f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth21({0.066244793f, 0.067483929f, 0.068723065f, 0.069962201f, 0.071201337f, 0.072440473f, 0.073679608f, 0.074918744f, 0.07615788f, 0.077397016f,
        0.078636152f, 0.080931439f, 0.083226725f, 0.085522012f, 0.087817299f, 0.090112586f, 0.092407873f, 0.09470316f, 0.096998446f, 0.099293733f,
        0.10158902f, 0.105984264f, 0.110379508f, 0.114774752f, 0.119169995f, 0.123565239f, 0.127960483f, 0.132355727f, 0.136750971f, 0.141146215f,
        0.145541458f, 0.150938738f, 0.156336017f, 0.161733296f, 0.167130575f, 0.172527854f, 0.177925133f, 0.183322413f, 0.188719692f, 0.194116971f,
        0.19951425f, 0.204002765f, 0.20849128f, 0.212979795f, 0.21746831f, 0.221956825f, 0.22644534f, 0.230933855f, 0.23542237f, 0.239910885f,
        0.2443994f, 0.248209173f, 0.252018945f, 0.255828718f, 0.25963849f, 0.263448263f, 0.267258035f, 0.271067808f, 0.27487758f, 0.278687353f,
        0.282497125f, 0.28518319f, 0.287869255f, 0.29055532f, 0.293241385f, 0.29592745f, 0.298613515f, 0.30129958f, 0.303985645f, 0.30667171f,
        0.309357775f, 0.309180931f, 0.309004087f, 0.308827243f, 0.308650398f, 0.308473554f, 0.30829671f, 0.308119866f, 0.307943022f, 0.307766178f,
        0.307589333f, 0.304611171f, 0.301633008f, 0.298654846f, 0.295676683f, 0.292698521f, 0.289720358f, 0.286742196f, 0.283764033f, 0.280785871f,
        0.277807708f, 0.273113663f, 0.268419618f, 0.263725573f, 0.259031528f, 0.254337483f, 0.249643438f, 0.244949393f, 0.240255348f, 0.235561303f,
        0.230867258f, 0.225534513f, 0.220201767f, 0.214869021f, 0.209536275f, 0.204203529f, 0.198870783f, 0.193538038f, 0.188205292f, 0.182872546f,
        0.1775398f, 0.172757168f, 0.167974537f, 0.163191905f, 0.158409273f, 0.153626642f, 0.14884401f, 0.144061378f, 0.139278747f, 0.134496115f,
        0.129713483f, 0.126170058f, 0.122626632f, 0.119083206f, 0.11553978f, 0.111996354f, 0.108452928f, 0.104909503f, 0.101366077f, 0.097822651f,
        0.094279225f, 0.091805435f, 0.089331646f, 0.086857856f, 0.084384066f, 0.081910277f, 0.079436487f, 0.076962697f, 0.074488908f, 0.072015118f,
        0.069541328f, 0.067986135f, 0.066430942f, 0.064875749f, 0.063320556f, 0.061765363f, 0.06021017f, 0.058654977f, 0.057099784f, 0.055544591f,
        0.053989398f, 0.053172617f, 0.052355836f, 0.051539055f, 0.050722274f, 0.049905493f, 0.049088711f, 0.04827193f, 0.047455149f, 0.046638368f,
        0.045821587f, 0.04540613f, 0.044990673f, 0.044575217f, 0.04415976f, 0.043744303f, 0.043328847f, 0.04291339f, 0.042497933f, 0.042082477f,
        0.04166702f, 0.041444174f, 0.041221328f, 0.040998483f, 0.040775637f, 0.040552791f, 0.040329945f, 0.040107099f, 0.039884253f, 0.039661408f,
        0.039438562f, 0.03932571f, 0.039212859f, 0.039100007f, 0.038987156f, 0.038874304f, 0.038761453f, 0.038648601f, 0.03853575f, 0.038422898f,
        0.038310047f, 0.038254032f, 0.038198017f, 0.038142002f, 0.038085987f, 0.038029973f, 0.037973958f, 0.037917943f, 0.037861928f, 0.037805913f,
        0.037749898f, 0.037750194f, 0.03775049f, 0.037750786f, 0.037751082f, 0.037751378f, 0.037751674f, 0.03775197f, 0.037752266f, 0.037752562f,
        0.037752858f, 0.037774517f, 0.037796176f, 0.037817835f, 0.037839494f, 0.037861153f, 0.037882811f, 0.03790447f, 0.037926129f, 0.037947788f,
        0.037969447f, 0.038024552f, 0.038079657f, 0.038134763f, 0.038189868f, 0.038244973f, 0.038300079f, 0.038355184f, 0.038410289f, 0.038465395f,
        0.0385205f, 0.038593144f, 0.038665789f, 0.038738433f, 0.038811077f, 0.038883722f, 0.038956366f, 0.03902901f, 0.039101655f, 0.039174299f,
        0.039246943f, 0.039308816f, 0.039370688f, 0.039432561f, 0.039494433f, 0.039556306f, 0.039618178f, 0.039680051f, 0.039741923f, 0.039803796f,
        0.039865668f, 0.039958673f, 0.040051677f, 0.040144682f, 0.040237686f, 0.040330691f, 0.040423695f, 0.0405167f, 0.040609704f, 0.040702709f,
        0.040795713f, 0.040948515f, 0.041101316f, 0.041254117f, 0.041406919f, 0.04155972f, 0.041712521f, 0.041865323f, 0.042018124f, 0.042170925f,
        0.042323727f, 0.042509255f, 0.042694784f, 0.042880312f, 0.043065841f, 0.043251369f, 0.043436898f, 0.043622426f, 0.043807955f, 0.043993483f,
        0.044179012f, 0.04430846f, 0.044437908f, 0.044567356f, 0.044696804f, 0.044826252f, 0.0449557f, 0.045085148f, 0.045214596f, 0.045344044f,
        0.045473492f, 0.045510572f, 0.045547653f, 0.045584734f, 0.045621814f, 0.045658895f, 0.045695976f, 0.045733056f, 0.045770137f, 0.045807218f,
        0.045844298f, 0.045901646f, 0.045958994f, 0.046016342f, 0.04607369f, 0.046131038f, 0.046188385f, 0.046245733f, 0.046303081f, 0.046360429f,
        0.046417777f, 0.046612773f, 0.04680777f, 0.047002766f, 0.047197763f, 0.047392759f, 0.047587756f, 0.047782752f, 0.047977749f, 0.048172745f,
        0.048367742f, 0.048748155f, 0.049128568f, 0.049508982f, 0.049889395f, 0.050269808f, 0.050650222f, 0.051030635f, 0.051411048f, 0.051791462f,
        0.052171875f, 0.05268579f, 0.053199704f, 0.053713619f, 0.054227534f, 0.054741448f, 0.055255363f, 0.055769278f, 0.056283192f, 0.056797107f,
        0.057311022f, 0.058078237f, 0.058845453f, 0.059612669f, 0.060379884f, 0.0611471f, 0.061914316f, 0.062681531f, 0.063448747f, 0.064215963f,
        0.064983178f, 0.058484861f, 0.051986543f, 0.045488225f, 0.038989907f, 0.032491589f, 0.025993271f, 0.019494954f, 0.012996636f, 0.006498318f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth22({0.051950585f, 0.052061574f, 0.052172564f, 0.052283553f, 0.052394542f, 0.052505532f, 0.052616521f, 0.05272751f, 0.0528385f, 0.052949489f,
        0.053060478f, 0.053174412f, 0.053288345f, 0.053402279f, 0.053516212f, 0.053630146f, 0.053744079f, 0.053858013f, 0.053971946f, 0.05408588f,
        0.054199813f, 0.054324554f, 0.054449295f, 0.054574036f, 0.054698777f, 0.054823518f, 0.054948259f, 0.055073f, 0.055197741f, 0.055322482f,
        0.055447223f, 0.055579124f, 0.055711025f, 0.055842925f, 0.055974826f, 0.056106727f, 0.056238627f, 0.056370528f, 0.056502429f, 0.056634329f,
        0.05676623f, 0.056945853f, 0.057125476f, 0.057305099f, 0.057484722f, 0.057664345f, 0.057843968f, 0.058023591f, 0.058203214f, 0.058382837f,
        0.05856246f, 0.058843538f, 0.059124616f, 0.059405694f, 0.059686771f, 0.059967849f, 0.060248927f, 0.060530005f, 0.060811083f, 0.061092161f,
        0.061373238f, 0.061811859f, 0.062250479f, 0.062689099f, 0.06312772f, 0.06356634f, 0.06400496f, 0.064443581f, 0.064882201f, 0.065320821f,
        0.065759442f, 0.066666911f, 0.06757438f, 0.068481849f, 0.069389318f, 0.070296788f, 0.071204257f, 0.072111726f, 0.073019195f, 0.073926664f,
        0.074834133f, 0.076619936f, 0.078405738f, 0.08019154f, 0.081977342f, 0.083763144f, 0.085548946f, 0.087334749f, 0.089120551f, 0.090906353f,
        0.092692155f, 0.095911162f, 0.099130169f, 0.102349176f, 0.105568183f, 0.10878719f, 0.112006197f, 0.115225204f, 0.118444211f, 0.121663218f,
        0.124882225f, 0.130183128f, 0.135484032f, 0.140784935f, 0.146085838f, 0.151386742f, 0.156687645f, 0.161988548f, 0.167289452f, 0.172590355f,
        0.177891258f, 0.184680648f, 0.191470038f, 0.198259428f, 0.205048818f, 0.211838208f, 0.218627598f, 0.225416988f, 0.232206378f, 0.238995768f,
        0.245785158f, 0.25193172f, 0.258078282f, 0.264224843f, 0.270371405f, 0.276517967f, 0.282664528f, 0.28881109f, 0.294957652f, 0.301104213f,
        0.307250775f, 0.310241468f, 0.31323216f, 0.316222853f, 0.319213545f, 0.322204238f, 0.32519493f, 0.328185623f, 0.331176315f, 0.334167008f,
        0.3371577f, 0.336795513f, 0.336433325f, 0.336071138f, 0.33570895f, 0.335346763f, 0.334984575f, 0.334622388f, 0.3342602f, 0.333898013f,
        0.333535825f, 0.331835113f, 0.3301344f, 0.328433688f, 0.326732975f, 0.325032263f, 0.32333155f, 0.321630838f, 0.319930125f, 0.318229413f,
        0.3165287f, 0.314175215f, 0.31182173f, 0.309468245f, 0.30711476f, 0.304761275f, 0.30240779f, 0.300054305f, 0.29770082f, 0.295347335f,
        0.29299385f, 0.289880905f, 0.28676796f, 0.283655015f, 0.28054207f, 0.277429125f, 0.27431618f, 0.271203235f, 0.26809029f, 0.264977345f,
        0.2618644f, 0.258677209f, 0.255490018f, 0.252302828f, 0.249115637f, 0.245928446f, 0.242741255f, 0.239554064f, 0.236366873f, 0.233179683f,
        0.229992492f, 0.22675861f, 0.223524728f, 0.220290847f, 0.217056965f, 0.213823083f, 0.210589202f, 0.20735532f, 0.204121438f, 0.200887557f,
        0.197653675f, 0.194392049f, 0.191130423f, 0.187868798f, 0.184607172f, 0.181345546f, 0.17808392f, 0.174822294f, 0.171560668f, 0.168299043f,
        0.165037417f, 0.162034949f, 0.159032482f, 0.156030014f, 0.153027547f, 0.150025079f, 0.147022612f, 0.144020144f, 0.141017677f, 0.138015209f,
        0.135012742f, 0.133001506f, 0.13099027f, 0.128979034f, 0.126967798f, 0.124956563f, 0.122945327f, 0.120934091f, 0.118922855f, 0.116911619f,
        0.114900383f, 0.113807521f, 0.112714659f, 0.111621797f, 0.110528935f, 0.109436073f, 0.108343211f, 0.107250349f, 0.106157487f, 0.105064625f,
        0.103971763f, 0.103366052f, 0.10276034f, 0.102154628f, 0.101548916f, 0.100943204f, 0.100337492f, 0.099731781f, 0.099126069f, 0.098520357f,
        0.097914645f, 0.097562141f, 0.097209637f, 0.096857133f, 0.096504629f, 0.096152125f, 0.095799621f, 0.095447117f, 0.095094613f, 0.094742109f,
        0.094389605f, 0.094185508f, 0.093981411f, 0.093777314f, 0.093573216f, 0.093369119f, 0.093165022f, 0.092960925f, 0.092756828f, 0.092552731f,
        0.092348633f, 0.092390764f, 0.092432894f, 0.092475024f, 0.092517155f, 0.092559285f, 0.092601415f, 0.092643546f, 0.092685676f, 0.092727806f,
        0.092769937f, 0.09314561f, 0.093521283f, 0.093896957f, 0.09427263f, 0.094648303f, 0.095023977f, 0.09539965f, 0.095775323f, 0.096150997f,
        0.09652667f, 0.097114159f, 0.097701649f, 0.098289138f, 0.098876627f, 0.099464117f, 0.100051606f, 0.100639095f, 0.101226585f, 0.101814074f,
        0.102401563f, 0.103003696f, 0.103605828f, 0.10420796f, 0.104810093f, 0.105412225f, 0.106014357f, 0.10661649f, 0.107218622f, 0.107820754f,
        0.108422887f, 0.108925907f, 0.109428927f, 0.109931948f, 0.110434968f, 0.110937988f, 0.111441009f, 0.111944029f, 0.112447049f, 0.11295007f,
        0.11345309f, 0.11364037f, 0.11382765f, 0.11401493f, 0.11420221f, 0.11438949f, 0.11457677f, 0.11476405f, 0.11495133f, 0.11513861f, 0.11532589f,
        0.115185591f, 0.115045292f, 0.114904993f, 0.114764693f, 0.114624394f, 0.114484095f, 0.114343796f, 0.114203497f, 0.114063198f, 0.113922898f,
        0.113957618f, 0.113992338f, 0.114027058f, 0.114061778f, 0.114096498f, 0.114131218f, 0.114165938f, 0.114200658f, 0.114235378f, 0.114270098f,
        0.102843089f, 0.091416079f, 0.079989069f, 0.068562059f, 0.057135049f, 0.045708039f, 0.03428103f, 0.02285402f, 0.01142701f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth23({0.049923225f, 0.049808127f, 0.049693028f, 0.04957793f, 0.049462832f, 0.049347733f, 0.049232635f, 0.049117537f, 0.049002438f, 0.04888734f,
        0.048772242f, 0.048654367f, 0.048536492f, 0.048418618f, 0.048300743f, 0.048182868f, 0.048064994f, 0.047947119f, 0.047829244f, 0.04771137f, 0.047593495f, 0.047557731f,
        0.047521968f, 0.047486204f, 0.04745044f, 0.047414677f, 0.047378913f, 0.047343149f, 0.047307386f, 0.047271622f, 0.047235858f, 0.047228584f, 0.047221309f, 0.047214034f,
        0.047206759f, 0.047199484f, 0.047192209f, 0.047184935f, 0.04717766f, 0.047170385f, 0.04716311f, 0.047182243f, 0.047201376f, 0.047220509f, 0.047239641f, 0.047258774f,
        0.047277907f, 0.04729704f, 0.047316173f, 0.047335306f, 0.047354438f, 0.047360874f, 0.047367309f, 0.047373745f, 0.04738018f, 0.047386616f, 0.047393051f, 0.047399487f,
        0.047405922f, 0.047412358f, 0.047418793f, 0.047368927f, 0.047319061f, 0.047269195f, 0.047219329f, 0.047169463f, 0.047119597f, 0.047069731f, 0.047019865f, 0.046969999f,
        0.046920133f, 0.046835367f, 0.046750601f, 0.046665835f, 0.046581069f, 0.046496303f, 0.046411537f, 0.046326771f, 0.046242005f, 0.046157239f, 0.046072473f, 0.045982751f,
        0.045893029f, 0.045803307f, 0.045713585f, 0.045623863f, 0.04553414f, 0.045444418f, 0.045354696f, 0.045264974f, 0.045175252f, 0.045101658f, 0.045028063f, 0.044954469f,
        0.044880875f, 0.044807281f, 0.044733687f, 0.044660093f, 0.044586498f, 0.044512904f, 0.04443931f, 0.044424251f, 0.044409192f, 0.044394133f, 0.044379074f, 0.044364015f,
        0.044348956f, 0.044333897f, 0.044318838f, 0.044303779f, 0.04428872f, 0.04432793f, 0.04436714f, 0.04440635f, 0.04444556f, 0.04448477f, 0.04452398f, 0.04456319f,
        0.0446024f, 0.04464161f, 0.04468082f, 0.044772309f, 0.044863799f, 0.044955288f, 0.045046777f, 0.045138267f, 0.045229756f, 0.045321245f, 0.045412735f, 0.045504224f,
        0.045595713f, 0.045714163f, 0.045832613f, 0.045951063f, 0.046069513f, 0.046187963f, 0.046306413f, 0.046424863f, 0.046543313f, 0.046661763f, 0.046780213f, 0.046866198f,
        0.046952183f, 0.047038167f, 0.047124152f, 0.047210137f, 0.047296121f, 0.047382106f, 0.047468091f, 0.047554075f, 0.04764006f, 0.047734641f, 0.047829222f, 0.047923803f,
        0.048018383f, 0.048112964f, 0.048207545f, 0.048302126f, 0.048396707f, 0.048491288f, 0.048585868f, 0.048763807f, 0.048941745f, 0.049119684f, 0.049297622f, 0.049475561f,
        0.049653499f, 0.049831438f, 0.050009376f, 0.050187315f, 0.050365253f, 0.050713413f, 0.051061572f, 0.051409731f, 0.051757891f, 0.05210605f, 0.052454209f, 0.052802369f,
        0.053150528f, 0.053498687f, 0.053846847f, 0.054448226f, 0.055049604f, 0.055650983f, 0.056252362f, 0.056853741f, 0.05745512f, 0.058056499f, 0.058657877f, 0.059259256f,
        0.059860635f, 0.061086313f, 0.062311991f, 0.063537669f, 0.064763347f, 0.065989025f, 0.067214703f, 0.068440381f, 0.069666059f, 0.070891737f, 0.072117415f, 0.07526156f,
        0.078405705f, 0.08154985f, 0.084693995f, 0.08783814f, 0.090982285f, 0.09412643f, 0.097270575f, 0.10041472f, 0.103558865f, 0.110955209f, 0.118351554f, 0.125747898f,
        0.133144242f, 0.140540587f, 0.147936931f, 0.155333275f, 0.16272962f, 0.170125964f, 0.177522308f, 0.190977247f, 0.204432185f, 0.217887123f, 0.231342062f, 0.244797f,
        0.258251938f, 0.271706877f, 0.285161815f, 0.298616753f, 0.312071692f, 0.32754781f, 0.343023928f, 0.358500047f, 0.373976165f, 0.389452283f, 0.404928402f, 0.42040452f,
        0.435880638f, 0.451356757f, 0.466832875f, 0.478232743f, 0.489632612f, 0.50103248f, 0.512432348f, 0.523832217f, 0.535232085f, 0.546631953f, 0.558031822f, 0.56943169f,
        0.580831558f, 0.587191888f, 0.593552218f, 0.599912548f, 0.606272878f, 0.612633208f, 0.618993538f, 0.625353868f, 0.631714198f, 0.638074528f, 0.644434858f, 0.647475183f,
        0.650515507f, 0.653555831f, 0.656596155f, 0.659636479f, 0.662676803f, 0.665717128f, 0.668757452f, 0.671797776f, 0.6748381f, 0.676371843f, 0.677905586f, 0.679439329f,
        0.680973071f, 0.682506814f, 0.684040557f, 0.6855743f, 0.687108043f, 0.688641786f, 0.690175528f, 0.690982426f, 0.691789324f, 0.692596222f, 0.69340312f, 0.694210018f,
        0.695016916f, 0.695823814f, 0.696630712f, 0.69743761f, 0.698244508f, 0.699011744f, 0.69977898f, 0.700546216f, 0.701313452f, 0.702080688f, 0.702847923f, 0.703615159f,
        0.704382395f, 0.705149631f, 0.705916867f, 0.706820624f, 0.707724382f, 0.708628139f, 0.709531897f, 0.710435654f, 0.711339412f, 0.712243169f, 0.713146927f, 0.714050684f,
        0.714954442f, 0.715829284f, 0.716704127f, 0.717578969f, 0.718453812f, 0.719328654f, 0.720203497f, 0.721078339f, 0.721953182f, 0.722828024f, 0.723702867f, 0.724343013f,
        0.724983158f, 0.725623304f, 0.72626345f, 0.726903596f, 0.727543742f, 0.728183888f, 0.728824033f, 0.729464179f, 0.730104325f, 0.730465293f, 0.730826262f, 0.73118723f,
        0.731548198f, 0.731909167f, 0.732270135f, 0.732631103f, 0.732992072f, 0.73335304f, 0.733714008f, 0.734183318f, 0.734652627f, 0.735121936f, 0.735591245f, 0.736060554f,
        0.736529863f, 0.736999173f, 0.737468482f, 0.737937791f, 0.7384071f, 0.66456639f, 0.59072568f, 0.51688497f, 0.44304426f, 0.36920355f, 0.29536284f, 0.22152213f, 0.14768142f,
        0.07384071f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth24({0.057984715f, 0.057628052f, 0.057271388f, 0.056914725f, 0.056558062f, 0.056201398f, 0.055844735f, 0.055488072f, 0.055131408f, 0.054774745f,
        0.054418082f, 0.054192412f, 0.053966742f, 0.053741073f, 0.053515403f, 0.053289733f, 0.053064064f, 0.052838394f, 0.052612724f, 0.052387055f, 0.052161385f, 0.052143699f,
        0.052126012f, 0.052108326f, 0.052090639f, 0.052072953f, 0.052055266f, 0.05203758f, 0.052019893f, 0.052002207f, 0.05198452f, 0.052048576f, 0.052112633f, 0.052176689f,
        0.052240745f, 0.052304802f, 0.052368858f, 0.052432914f, 0.052496971f, 0.052561027f, 0.052625083f, 0.052760931f, 0.052896779f, 0.053032627f, 0.053168475f, 0.053304323f,
        0.05344017f, 0.053576018f, 0.053711866f, 0.053847714f, 0.053983562f, 0.054192932f, 0.054402303f, 0.054611674f, 0.054821044f, 0.055030415f, 0.055239786f, 0.055449156f,
        0.055658527f, 0.055867898f, 0.056077268f, 0.056411799f, 0.05674633f, 0.057080861f, 0.057415392f, 0.057749923f, 0.058084454f, 0.058418985f, 0.058753516f, 0.059088047f,
        0.059422578f, 0.060139201f, 0.060855824f, 0.061572446f, 0.062289069f, 0.063005692f, 0.063722314f, 0.064438937f, 0.06515556f, 0.065872182f, 0.066588805f, 0.067997994f,
        0.069407183f, 0.070816372f, 0.07222556f, 0.073634749f, 0.075043938f, 0.076453127f, 0.077862316f, 0.079271505f, 0.080680693f, 0.083301016f, 0.085921339f, 0.088541662f,
        0.091161985f, 0.093782308f, 0.09640263f, 0.099022953f, 0.101643276f, 0.104263599f, 0.106883922f, 0.111399713f, 0.115915504f, 0.120431295f, 0.124947086f, 0.129462878f,
        0.133978669f, 0.13849446f, 0.143010251f, 0.147526042f, 0.152041833f, 0.159344925f, 0.166648017f, 0.173951108f, 0.1812542f, 0.188557292f, 0.195860383f, 0.203163475f,
        0.210466567f, 0.217769658f, 0.22507275f, 0.236118538f, 0.247164327f, 0.258210115f, 0.269255903f, 0.280301692f, 0.29134748f, 0.302393268f, 0.313439057f, 0.324484845f,
        0.335530633f, 0.348216821f, 0.360903008f, 0.373589196f, 0.386275383f, 0.398961571f, 0.411647758f, 0.424333946f, 0.437020133f, 0.449706321f, 0.462392508f, 0.472026513f,
        0.481660517f, 0.491294521f, 0.500928525f, 0.510562529f, 0.520196533f, 0.529830538f, 0.539464542f, 0.549098546f, 0.55873255f, 0.564431976f, 0.570131402f, 0.575830828f,
        0.581530253f, 0.587229679f, 0.592929105f, 0.598628531f, 0.604327957f, 0.610027383f, 0.615726808f, 0.619126633f, 0.622526458f, 0.625926283f, 0.629326108f, 0.632725933f,
        0.636125758f, 0.639525583f, 0.642925408f, 0.646325233f, 0.649725058f, 0.651974476f, 0.654223893f, 0.656473311f, 0.658722728f, 0.660972146f, 0.663221563f, 0.665470981f,
        0.667720398f, 0.669969816f, 0.672219233f, 0.674384135f, 0.676549037f, 0.678713938f, 0.68087884f, 0.683043742f, 0.685208643f, 0.687373545f, 0.689538447f, 0.691703348f,
        0.69386825f, 0.695476616f, 0.697084982f, 0.698693348f, 0.700301713f, 0.701910079f, 0.703518445f, 0.705126811f, 0.706735177f, 0.708343543f, 0.709951908f, 0.711275701f,
        0.712599493f, 0.713923286f, 0.715247078f, 0.716570871f, 0.717894663f, 0.719218456f, 0.720542248f, 0.721866041f, 0.723189833f, 0.724014358f, 0.724838883f, 0.725663408f,
        0.726487933f, 0.727312458f, 0.728136983f, 0.728961508f, 0.729786033f, 0.730610558f, 0.731435083f, 0.732195703f, 0.732956322f, 0.733716941f, 0.73447756f, 0.735238179f,
        0.735998798f, 0.736759418f, 0.737520037f, 0.738280656f, 0.739041275f, 0.739757163f, 0.74047305f, 0.741188938f, 0.741904825f, 0.742620713f, 0.7433366f, 0.744052488f,
        0.744768375f, 0.745484263f, 0.74620015f, 0.746760517f, 0.747320883f, 0.74788125f, 0.748441617f, 0.749001983f, 0.74956235f, 0.750122717f, 0.750683083f, 0.75124345f,
        0.751803817f, 0.752439286f, 0.753074755f, 0.753710224f, 0.754345693f, 0.754981163f, 0.755616632f, 0.756252101f, 0.75688757f, 0.757523039f, 0.758158508f, 0.758736744f,
        0.75931498f, 0.759893216f, 0.760471452f, 0.761049688f, 0.761627923f, 0.762206159f, 0.762784395f, 0.763362631f, 0.763940867f, 0.764415688f, 0.76489051f, 0.765365332f,
        0.765840153f, 0.766314975f, 0.766789797f, 0.767264618f, 0.76773944f, 0.768214262f, 0.768689083f, 0.768918358f, 0.769147632f, 0.769376906f, 0.76960618f, 0.769835454f,
        0.770064728f, 0.770294003f, 0.770523277f, 0.770752551f, 0.770981825f, 0.771434485f, 0.771887145f, 0.772339805f, 0.772792465f, 0.773245125f, 0.773697785f, 0.774150445f,
        0.774603105f, 0.775055765f, 0.775508425f, 0.776197791f, 0.776887157f, 0.777576523f, 0.778265888f, 0.778955254f, 0.77964462f, 0.780333986f, 0.781023352f, 0.781712718f,
        0.782402083f, 0.783179574f, 0.783957065f, 0.784734556f, 0.785512047f, 0.786289538f, 0.787067028f, 0.787844519f, 0.78862201f, 0.789399501f, 0.790176992f, 0.79077827f,
        0.791379548f, 0.791980827f, 0.792582105f, 0.793183383f, 0.793784662f, 0.79438594f, 0.794987218f, 0.795588497f, 0.796189775f, 0.796500418f, 0.796811062f, 0.797121705f,
        0.797432348f, 0.797742992f, 0.798053635f, 0.798364278f, 0.798674922f, 0.798985565f, 0.799296208f, 0.799732681f, 0.800169153f, 0.800605626f, 0.801042098f, 0.801478571f,
        0.801915043f, 0.802351516f, 0.802787988f, 0.803224461f, 0.803660933f, 0.72329484f, 0.642928747f, 0.562562653f, 0.48219656f, 0.401830467f, 0.321464373f, 0.24109828f,
        0.160732187f, 0.080366093f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth25({0.144547292f, 0.149603995f, 0.154660698f, 0.159717402f, 0.164774105f, 0.169830808f, 0.174887512f, 0.179944215f, 0.185000918f, 0.190057622f,
        0.195114325f, 0.20386221f, 0.212610095f, 0.22135798f, 0.230105865f, 0.23885375f, 0.247601635f, 0.25634952f, 0.265097405f, 0.27384529f, 0.282593175f, 0.288910735f,
        0.295228295f, 0.301545855f, 0.307863415f, 0.314180975f, 0.320498535f, 0.326816095f, 0.333133655f, 0.339451215f, 0.345768775f, 0.347374304f, 0.348979833f, 0.350585363f,
        0.352190892f, 0.353796421f, 0.35540195f, 0.357007479f, 0.358613008f, 0.360218538f, 0.361824067f, 0.361073996f, 0.360323925f, 0.359573854f, 0.358823783f, 0.358073713f,
        0.357323642f, 0.356573571f, 0.3558235f, 0.355073429f, 0.354323358f, 0.352252065f, 0.350180772f, 0.348109478f, 0.346038185f, 0.343966892f, 0.341895598f, 0.339824305f,
        0.337753012f, 0.335681718f, 0.333610425f, 0.330820404f, 0.328030383f, 0.325240363f, 0.322450342f, 0.319660321f, 0.3168703f, 0.314080279f, 0.311290258f, 0.308500238f,
        0.305710217f, 0.302761862f, 0.299813507f, 0.296865152f, 0.293916797f, 0.290968442f, 0.288020087f, 0.285071732f, 0.282123377f, 0.279175022f, 0.276226667f, 0.273359616f,
        0.270492565f, 0.267625514f, 0.264758463f, 0.261891413f, 0.259024362f, 0.256157311f, 0.25329026f, 0.250423209f, 0.247556158f, 0.244605828f, 0.241655498f, 0.238705168f,
        0.235754838f, 0.232804508f, 0.229854178f, 0.226903848f, 0.223953518f, 0.221003188f, 0.218052858f, 0.215235745f, 0.212418632f, 0.209601518f, 0.206784405f, 0.203967292f,
        0.201150178f, 0.198333065f, 0.195515952f, 0.192698838f, 0.189881725f, 0.187679393f, 0.185477062f, 0.18327473f, 0.181072398f, 0.178870067f, 0.176667735f, 0.174465403f,
        0.172263072f, 0.17006074f, 0.167858408f, 0.165968076f, 0.164077743f, 0.162187411f, 0.160297078f, 0.158406746f, 0.156516413f, 0.154626081f, 0.152735748f, 0.150845416f,
        0.148955083f, 0.146756609f, 0.144558135f, 0.142359661f, 0.140161187f, 0.137962713f, 0.135764238f, 0.133565764f, 0.13136729f, 0.129168816f, 0.126970342f, 0.12499661f,
        0.123022878f, 0.121049147f, 0.119075415f, 0.117101683f, 0.115127952f, 0.11315422f, 0.111180488f, 0.109206757f, 0.107233025f, 0.106471288f, 0.10570955f, 0.104947813f,
        0.104186075f, 0.103424338f, 0.1026626f, 0.101900863f, 0.101139125f, 0.100377388f, 0.09961565f, 0.099842944f, 0.100070238f, 0.100297533f, 0.100524827f, 0.100752121f,
        0.100979415f, 0.101206709f, 0.101434003f, 0.101661298f, 0.101888592f, 0.102055773f, 0.102222955f, 0.102390137f, 0.102557318f, 0.1027245f, 0.102891682f, 0.103058863f,
        0.103226045f, 0.103393227f, 0.103560408f, 0.104111032f, 0.104661655f, 0.105212278f, 0.105762902f, 0.106313525f, 0.106864148f, 0.107414772f, 0.107965395f, 0.108516018f,
        0.109066642f, 0.111840288f, 0.114613935f, 0.117387582f, 0.120161228f, 0.122934875f, 0.125708522f, 0.128482168f, 0.131255815f, 0.134029462f, 0.136803108f, 0.14308545f,
        0.149367792f, 0.155650133f, 0.161932475f, 0.168214817f, 0.174497158f, 0.1807795f, 0.187061842f, 0.193344183f, 0.199626525f, 0.208676675f, 0.217726825f, 0.226776975f,
        0.235827125f, 0.244877275f, 0.253927425f, 0.262977575f, 0.272027725f, 0.281077875f, 0.290128025f, 0.301121203f, 0.312114382f, 0.32310756f, 0.334100738f, 0.345093917f,
        0.356087095f, 0.367080273f, 0.378073452f, 0.38906663f, 0.400059808f, 0.411633658f, 0.423207507f, 0.434781356f, 0.446355205f, 0.457929054f, 0.469502903f, 0.481076753f,
        0.492650602f, 0.504224451f, 0.5157983f, 0.525704643f, 0.535610985f, 0.545517328f, 0.55542367f, 0.565330013f, 0.575236355f, 0.585142698f, 0.59504904f, 0.604955383f,
        0.614861725f, 0.622030803f, 0.629199882f, 0.63636896f, 0.643538038f, 0.650707117f, 0.657876195f, 0.665045273f, 0.672214352f, 0.67938343f, 0.686552508f, 0.69107417f,
        0.695595832f, 0.700117493f, 0.704639155f, 0.709160817f, 0.713682478f, 0.71820414f, 0.722725802f, 0.727247463f, 0.731769125f, 0.734567121f, 0.737365117f, 0.740163113f,
        0.742961108f, 0.745759104f, 0.7485571f, 0.751355096f, 0.754153092f, 0.756951088f, 0.759749083f, 0.761206923f, 0.762664762f, 0.764122601f, 0.76558044f, 0.767038279f,
        0.768496118f, 0.769953958f, 0.771411797f, 0.772869636f, 0.774327475f, 0.775209151f, 0.776090827f, 0.776972503f, 0.777854178f, 0.778735854f, 0.77961753f, 0.780499206f,
        0.781380882f, 0.782262558f, 0.783144233f, 0.784085772f, 0.78502731f, 0.785968848f, 0.786910387f, 0.787851925f, 0.788793463f, 0.789735002f, 0.79067654f, 0.791618078f,
        0.792559617f, 0.793640925f, 0.794722233f, 0.795803542f, 0.79688485f, 0.797966158f, 0.799047467f, 0.800128775f, 0.801210083f, 0.802291392f, 0.8033727f, 0.804190628f,
        0.805008555f, 0.805826483f, 0.80664441f, 0.807462338f, 0.808280265f, 0.809098193f, 0.80991612f, 0.810734048f, 0.811551975f, 0.812114972f, 0.812677968f, 0.813240965f,
        0.813803962f, 0.814366958f, 0.814929955f, 0.815492952f, 0.816055948f, 0.816618945f, 0.817181942f, 0.818004698f, 0.818827453f, 0.819650209f, 0.820472965f, 0.821295721f,
        0.822118477f, 0.822941233f, 0.823763988f, 0.824586744f, 0.8254095f, 0.74286855f, 0.6603276f, 0.57778665f, 0.4952457f, 0.41270475f, 0.3301638f, 0.24762285f, 0.1650819f,
        0.08254095f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth26({0.107731417f, 0.111077022f, 0.114422627f, 0.117768232f, 0.121113837f, 0.124459442f, 0.127805047f, 0.131150652f, 0.134496257f, 0.137841862f,
        0.141187467f, 0.146315534f, 0.151443602f, 0.156571669f, 0.161699737f, 0.166827804f, 0.171955872f, 0.177083939f, 0.182212007f, 0.187340074f, 0.192468142f, 0.196862589f,
        0.201257037f, 0.205651484f, 0.210045932f, 0.214440379f, 0.218834827f, 0.223229274f, 0.227623722f, 0.232018169f, 0.236412617f, 0.238856745f, 0.241300873f, 0.243745002f,
        0.24618913f, 0.248633258f, 0.251077387f, 0.253521515f, 0.255965643f, 0.258409772f, 0.2608539f, 0.263318711f, 0.265783522f, 0.268248333f, 0.270713143f, 0.273177954f,
        0.275642765f, 0.278107576f, 0.280572387f, 0.283037198f, 0.285502008f, 0.2886915f, 0.291880992f, 0.295070483f, 0.298259975f, 0.301449467f, 0.304638958f, 0.30782845f,
        0.311017942f, 0.314207433f, 0.317396925f, 0.320970598f, 0.32454427f, 0.328117943f, 0.331691615f, 0.335265288f, 0.33883896f, 0.342412633f, 0.345986305f, 0.349559978f,
        0.35313365f, 0.356843817f, 0.360553983f, 0.36426415f, 0.367974317f, 0.371684483f, 0.37539465f, 0.379104817f, 0.382814983f, 0.38652515f, 0.390235317f, 0.393808434f,
        0.397381552f, 0.400954669f, 0.404527787f, 0.408100904f, 0.411674022f, 0.415247139f, 0.418820257f, 0.422393374f, 0.425966492f, 0.427930518f, 0.429894545f, 0.431858572f,
        0.433822598f, 0.435786625f, 0.437750652f, 0.439714678f, 0.441678705f, 0.443642732f, 0.445606758f, 0.445469133f, 0.445331508f, 0.445193883f, 0.445056258f, 0.444918633f,
        0.444781008f, 0.444643383f, 0.444505758f, 0.444368133f, 0.444230508f, 0.442128823f, 0.440027137f, 0.437925451f, 0.435823765f, 0.433722079f, 0.431620393f, 0.429518708f,
        0.427417022f, 0.425315336f, 0.42321365f, 0.419441042f, 0.415668433f, 0.411895825f, 0.408123217f, 0.404350608f, 0.400578f, 0.396805392f, 0.393032783f, 0.389260175f,
        0.385487567f, 0.380610718f, 0.375733868f, 0.370857019f, 0.36598017f, 0.361103321f, 0.356226472f, 0.351349623f, 0.346472773f, 0.341595924f, 0.336719075f, 0.331320268f,
        0.325921462f, 0.320522655f, 0.315123848f, 0.309725042f, 0.304326235f, 0.298927428f, 0.293528622f, 0.288129815f, 0.282731008f, 0.277586156f, 0.272441303f, 0.267296451f,
        0.262151598f, 0.257006746f, 0.251861893f, 0.246717041f, 0.241572188f, 0.236427336f, 0.231282483f, 0.226659771f, 0.222037058f, 0.217414346f, 0.212791633f, 0.208168921f,
        0.203546208f, 0.198923496f, 0.194300783f, 0.189678071f, 0.185055358f, 0.181103803f, 0.177152248f, 0.173200693f, 0.169249138f, 0.165297583f, 0.161346028f, 0.157394473f,
        0.153442918f, 0.149491363f, 0.145539808f, 0.142792583f, 0.140045357f, 0.137298131f, 0.134550905f, 0.131803679f, 0.129056453f, 0.126309228f, 0.123562002f, 0.120814776f,
        0.11806755f, 0.116313407f, 0.114559264f, 0.112805121f, 0.111050978f, 0.109296835f, 0.107542692f, 0.105788549f, 0.104034406f, 0.102280263f, 0.10052612f, 0.099431749f,
        0.098337378f, 0.097243008f, 0.096148637f, 0.095054266f, 0.093959895f, 0.092865524f, 0.091771153f, 0.090676783f, 0.089582412f, 0.088780049f, 0.087977686f, 0.087175323f,
        0.08637296f, 0.085570597f, 0.084768234f, 0.083965871f, 0.083163508f, 0.082361145f, 0.081558782f, 0.08104252f, 0.080526259f, 0.080009997f, 0.079493736f, 0.078977474f,
        0.078461213f, 0.077944951f, 0.07742869f, 0.076912428f, 0.076396167f, 0.076162454f, 0.075928742f, 0.075695029f, 0.075461317f, 0.075227604f, 0.074993892f, 0.074760179f,
        0.074526467f, 0.074292754f, 0.074059042f, 0.07395845f, 0.073857859f, 0.073757268f, 0.073656676f, 0.073556085f, 0.073455494f, 0.073354902f, 0.073254311f, 0.07315372f,
        0.073053128f, 0.073042195f, 0.073031262f, 0.073020328f, 0.073009395f, 0.072998462f, 0.072987528f, 0.072976595f, 0.072965662f, 0.072954728f, 0.072943795f, 0.073030288f,
        0.07311678f, 0.073203273f, 0.073289766f, 0.073376258f, 0.073462751f, 0.073549244f, 0.073635736f, 0.073722229f, 0.073808722f, 0.073986753f, 0.074164784f, 0.074342815f,
        0.074520846f, 0.074698878f, 0.074876909f, 0.07505494f, 0.075232971f, 0.075411002f, 0.075589033f, 0.075705172f, 0.075821311f, 0.07593745f, 0.076053589f, 0.076169728f,
        0.076285867f, 0.076402006f, 0.076518145f, 0.076634284f, 0.076750423f, 0.07672324f, 0.076696057f, 0.076668874f, 0.076641691f, 0.076614508f, 0.076587324f, 0.076560141f,
        0.076532958f, 0.076505775f, 0.076478592f, 0.076329891f, 0.07618119f, 0.076032489f, 0.075883788f, 0.075735087f, 0.075586386f, 0.075437685f, 0.075288984f, 0.075140283f,
        0.074991582f, 0.074769379f, 0.074547177f, 0.074324975f, 0.074102772f, 0.07388057f, 0.073658368f, 0.073436165f, 0.073213963f, 0.072991761f, 0.072769558f, 0.072692466f,
        0.072615373f, 0.07253828f, 0.072461188f, 0.072384095f, 0.072307002f, 0.07222991f, 0.072152817f, 0.072075724f, 0.071998632f, 0.072173265f, 0.072347897f, 0.07252253f,
        0.072697163f, 0.072871796f, 0.073046429f, 0.073221062f, 0.073395694f, 0.073570327f, 0.07374496f, 0.074305955f, 0.074866949f, 0.075427944f, 0.075988939f, 0.076549933f,
        0.077110928f, 0.077671923f, 0.078232917f, 0.078793912f, 0.079354907f, 0.071419416f, 0.063483925f, 0.055548435f, 0.047612944f, 0.039677453f, 0.031741963f, 0.023806472f,
        0.015870981f, 0.007935491f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth31({0.189362217f, 0.195889904f, 0.202417592f, 0.208945279f, 0.215472967f, 0.222000654f, 0.228528342f, 0.235056029f, 0.241583717f, 0.248111404f,
        0.254639092f, 0.271435503f, 0.288231913f, 0.305028324f, 0.321824735f, 0.338621146f, 0.355417557f, 0.372213968f, 0.389010378f, 0.405806789f, 0.4226032f, 0.446363808f,
        0.470124417f, 0.493885025f, 0.517645633f, 0.541406242f, 0.56516685f, 0.588927458f, 0.612688067f, 0.636448675f, 0.660209283f, 0.675285888f, 0.690362493f, 0.705439098f,
        0.720515703f, 0.735592308f, 0.750668913f, 0.765745518f, 0.780822123f, 0.795898728f, 0.810975333f, 0.816089584f, 0.821203835f, 0.826318086f, 0.831432337f, 0.836546588f,
        0.841660838f, 0.846775089f, 0.85188934f, 0.857003591f, 0.862117842f, 0.863563629f, 0.865009417f, 0.866455204f, 0.867900992f, 0.869346779f, 0.870792567f, 0.872238354f,
        0.873684142f, 0.875129929f, 0.876575717f, 0.877335507f, 0.878095297f, 0.878855087f, 0.879614877f, 0.880374667f, 0.881134457f, 0.881894247f, 0.882654037f, 0.883413827f,
        0.884173617f, 0.884859832f, 0.885546047f, 0.886232262f, 0.886918477f, 0.887604692f, 0.888290907f, 0.888977122f, 0.889663337f, 0.890349552f, 0.891035767f, 0.891498596f,
        0.891961425f, 0.892424254f, 0.892887083f, 0.893349913f, 0.893812742f, 0.894275571f, 0.8947384f, 0.895201229f, 0.895664058f, 0.896029293f, 0.896394528f, 0.896759763f,
        0.897124998f, 0.897490233f, 0.897855468f, 0.898220703f, 0.898585938f, 0.898951173f, 0.899316408f, 0.899754399f, 0.90019239f, 0.900630381f, 0.901068372f, 0.901506363f,
        0.901944353f, 0.902382344f, 0.902820335f, 0.903258326f, 0.903696317f, 0.904045051f, 0.904393785f, 0.904742519f, 0.905091253f, 0.905439988f, 0.905788722f, 0.906137456f,
        0.90648619f, 0.906834924f, 0.907183658f, 0.907373738f, 0.907563817f, 0.907753896f, 0.907943975f, 0.908134054f, 0.908324133f, 0.908514213f, 0.908704292f, 0.908894371f,
        0.90908445f, 0.90926711f, 0.90944977f, 0.90963243f, 0.90981509f, 0.90999775f, 0.91018041f, 0.91036307f, 0.91054573f, 0.91072839f, 0.91091105f, 0.910825224f, 0.910739398f,
        0.910653573f, 0.910567747f, 0.910481921f, 0.910396095f, 0.910310269f, 0.910224443f, 0.910138618f, 0.910052792f, 0.910169041f, 0.91028529f, 0.910401539f, 0.910517788f,
        0.910634038f, 0.910750287f, 0.910866536f, 0.910982785f, 0.911099034f, 0.911215283f, 0.911495718f, 0.911776152f, 0.912056586f, 0.91233702f, 0.912617454f, 0.912897888f,
        0.913178323f, 0.913458757f, 0.913739191f, 0.914019625f, 0.913960303f, 0.913900982f, 0.91384166f, 0.913782338f, 0.913723017f, 0.913663695f, 0.913604373f, 0.913545052f,
        0.91348573f, 0.913426408f, 0.913685639f, 0.91394487f, 0.914204101f, 0.914463332f, 0.914722563f, 0.914981793f, 0.915241024f, 0.915500255f, 0.915759486f, 0.916018717f,
        0.915964848f, 0.915910978f, 0.915857109f, 0.91580324f, 0.915749371f, 0.915695502f, 0.915641633f, 0.915587763f, 0.915533894f, 0.915480025f, 0.915516173f, 0.915552322f,
        0.91558847f, 0.915624618f, 0.915660767f, 0.915696915f, 0.915733063f, 0.915769212f, 0.91580536f, 0.915841508f, 0.91569081f, 0.915540112f, 0.915389413f, 0.915238715f,
        0.915088017f, 0.914937318f, 0.91478662f, 0.914635922f, 0.914485223f, 0.914334525f, 0.914447788f, 0.91456105f, 0.914674313f, 0.914787575f, 0.914900838f, 0.9150141f,
        0.915127363f, 0.915240625f, 0.915353888f, 0.91546715f, 0.915684095f, 0.91590104f, 0.916117985f, 0.91633493f, 0.916551875f, 0.91676882f, 0.916985765f, 0.91720271f,
        0.917419655f, 0.9176366f, 0.917735608f, 0.917834617f, 0.917933625f, 0.918032633f, 0.918131642f, 0.91823065f, 0.918329658f, 0.918428667f, 0.918527675f, 0.918626683f,
        0.918865424f, 0.919104165f, 0.919342906f, 0.919581647f, 0.919820388f, 0.920059128f, 0.920297869f, 0.92053661f, 0.920775351f, 0.921014092f, 0.921203662f, 0.921393232f,
        0.921582802f, 0.921772372f, 0.921961942f, 0.922151512f, 0.922341082f, 0.922530652f, 0.922720222f, 0.922909792f, 0.923004354f, 0.923098917f, 0.923193479f, 0.923288042f,
        0.923382604f, 0.923477167f, 0.923571729f, 0.923666292f, 0.923760854f, 0.923855417f, 0.923668726f, 0.923482036f, 0.923295345f, 0.923108655f, 0.922921964f, 0.922735274f,
        0.922548583f, 0.922361893f, 0.922175202f, 0.921988512f, 0.922031691f, 0.922074871f, 0.922118051f, 0.92216123f, 0.92220441f, 0.92224759f, 0.922290769f, 0.922333949f,
        0.922377129f, 0.922420308f, 0.922655382f, 0.922890455f, 0.923125528f, 0.923360602f, 0.923595675f, 0.923830748f, 0.924065822f, 0.924300895f, 0.924535968f, 0.924771042f,
        0.925042622f, 0.925314202f, 0.925585782f, 0.925857362f, 0.926128942f, 0.926400522f, 0.926672102f, 0.926943682f, 0.927215262f, 0.927486842f, 0.927715248f, 0.927943653f,
        0.928172059f, 0.928400465f, 0.928628871f, 0.928857277f, 0.929085683f, 0.929314088f, 0.929542494f, 0.9297709f, 0.929834968f, 0.929899037f, 0.929963105f, 0.930027173f,
        0.930091242f, 0.93015531f, 0.930219378f, 0.930283447f, 0.930347515f, 0.930411583f, 0.930699356f, 0.930987128f, 0.931274901f, 0.931562673f, 0.931850446f, 0.932138218f,
        0.932425991f, 0.932713763f, 0.933001536f, 0.933289308f, 0.839960378f, 0.746631447f, 0.653302516f, 0.559973585f, 0.466644654f, 0.373315723f, 0.279986793f, 0.186657862f,
        0.093328931f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth32({0.17084875f, 0.176969716f, 0.183090682f, 0.189211648f, 0.195332613f, 0.201453579f, 0.207574545f, 0.213695511f, 0.219816477f, 0.225937443f,
        0.232058408f, 0.245359243f, 0.258660078f, 0.271960913f, 0.285261748f, 0.298562583f, 0.311863418f, 0.325164253f, 0.338465088f, 0.351765923f, 0.365066758f, 0.379215755f,
        0.393364752f, 0.407513748f, 0.421662745f, 0.435811742f, 0.449960738f, 0.464109735f, 0.478258732f, 0.492407728f, 0.506556725f, 0.512649794f, 0.518742863f, 0.524835933f,
        0.530929002f, 0.537022071f, 0.54311514f, 0.549208209f, 0.555301278f, 0.561394348f, 0.567487417f, 0.569008233f, 0.57052905f, 0.572049867f, 0.573570683f, 0.5750915f,
        0.576612317f, 0.578133133f, 0.57965395f, 0.581174767f, 0.582695583f, 0.583196129f, 0.583696675f, 0.584197221f, 0.584697767f, 0.585198313f, 0.585698858f, 0.586199404f,
        0.58669995f, 0.587200496f, 0.587701042f, 0.587940117f, 0.588179192f, 0.588418267f, 0.588657342f, 0.588896417f, 0.589135492f, 0.589374567f, 0.589613642f, 0.589852717f,
        0.590091792f, 0.590181816f, 0.59027184f, 0.590361864f, 0.590451888f, 0.590541913f, 0.590631937f, 0.590721961f, 0.590811985f, 0.590902009f, 0.590992033f, 0.590869947f,
        0.59074786f, 0.590625773f, 0.590503687f, 0.5903816f, 0.590259513f, 0.590137427f, 0.59001534f, 0.589893253f, 0.589771167f, 0.589634701f, 0.589498235f, 0.589361769f,
        0.589225303f, 0.589088838f, 0.588952372f, 0.588815906f, 0.58867944f, 0.588542974f, 0.588406508f, 0.588408395f, 0.588410282f, 0.588412168f, 0.588414055f, 0.588415942f,
        0.588417828f, 0.588419715f, 0.588421602f, 0.588423488f, 0.588425375f, 0.588480734f, 0.588536093f, 0.588591453f, 0.588646812f, 0.588702171f, 0.58875753f, 0.588812889f,
        0.588868248f, 0.588923608f, 0.588978967f, 0.589029458f, 0.58907995f, 0.589130442f, 0.589180933f, 0.589231425f, 0.589281917f, 0.589332408f, 0.5893829f, 0.589433392f,
        0.589483883f, 0.589594587f, 0.58970529f, 0.589815993f, 0.589926697f, 0.5900374f, 0.590148103f, 0.590258807f, 0.59036951f, 0.590480213f, 0.590590917f, 0.590533679f,
        0.590476442f, 0.590419204f, 0.590361967f, 0.590304729f, 0.590247492f, 0.590190254f, 0.590133017f, 0.590075779f, 0.590018542f, 0.590006389f, 0.589994237f, 0.589982084f,
        0.589969932f, 0.589957779f, 0.589945627f, 0.589933474f, 0.589921322f, 0.589909169f, 0.589897017f, 0.58993754f, 0.589978063f, 0.590018587f, 0.59005911f, 0.590099633f,
        0.590140157f, 0.59018068f, 0.590221203f, 0.590261727f, 0.59030225f, 0.590201258f, 0.590100265f, 0.589999273f, 0.58989828f, 0.589797288f, 0.589696295f, 0.589595303f,
        0.58949431f, 0.589393318f, 0.589292325f, 0.589457419f, 0.589622513f, 0.589787608f, 0.589952702f, 0.590117796f, 0.59028289f, 0.590447984f, 0.590613078f, 0.590778173f,
        0.590943267f, 0.590879753f, 0.59081624f, 0.590752727f, 0.590689213f, 0.5906257f, 0.590562187f, 0.590498673f, 0.59043516f, 0.590371647f, 0.590308133f, 0.590248664f,
        0.590189195f, 0.590129726f, 0.590070257f, 0.590010788f, 0.589951318f, 0.589891849f, 0.58983238f, 0.589772911f, 0.589713442f, 0.589455258f, 0.589197073f, 0.588938889f,
        0.588680705f, 0.588422521f, 0.588164337f, 0.587906153f, 0.587647968f, 0.587389784f, 0.5871316f, 0.586933215f, 0.58673483f, 0.586536445f, 0.58633806f, 0.586139675f,
        0.58594129f, 0.585742905f, 0.58554452f, 0.585346135f, 0.58514775f, 0.584937078f, 0.584726407f, 0.584515735f, 0.584305063f, 0.584094392f, 0.58388372f, 0.583673048f,
        0.583462377f, 0.583251705f, 0.583041033f, 0.582733138f, 0.582425242f, 0.582117346f, 0.58180945f, 0.581501554f, 0.581193658f, 0.580885763f, 0.580577867f, 0.580269971f,
        0.579962075f, 0.579744815f, 0.579527555f, 0.579310295f, 0.579093035f, 0.578875775f, 0.578658515f, 0.578441255f, 0.578223995f, 0.578006735f, 0.577789475f, 0.577605334f,
        0.577421193f, 0.577237053f, 0.577052912f, 0.576868771f, 0.57668463f, 0.576500489f, 0.576316348f, 0.576132208f, 0.575948067f, 0.575793163f, 0.57563826f, 0.575483357f,
        0.575328453f, 0.57517355f, 0.575018647f, 0.574863743f, 0.57470884f, 0.574553937f, 0.574399033f, 0.574180169f, 0.573961305f, 0.573742441f, 0.573523577f, 0.573304713f,
        0.573085848f, 0.572866984f, 0.57264812f, 0.572429256f, 0.572210392f, 0.572050543f, 0.571890693f, 0.571730844f, 0.571570995f, 0.571411146f, 0.571251297f, 0.571091448f,
        0.570931598f, 0.570771749f, 0.5706119f, 0.570473113f, 0.570334325f, 0.570195538f, 0.57005675f, 0.569917963f, 0.569779175f, 0.569640388f, 0.5695016f, 0.569362813f,
        0.569224025f, 0.569129288f, 0.569034552f, 0.568939815f, 0.568845078f, 0.568750342f, 0.568655605f, 0.568560868f, 0.568466132f, 0.568371395f, 0.568276658f, 0.568245926f,
        0.568215193f, 0.568184461f, 0.568153728f, 0.568122996f, 0.568092263f, 0.568061531f, 0.568030798f, 0.568000066f, 0.567969333f, 0.567820551f, 0.567671768f, 0.567522986f,
        0.567374203f, 0.567225421f, 0.567076638f, 0.566927856f, 0.566779073f, 0.566630291f, 0.566481508f, 0.566464265f, 0.566447022f, 0.566429778f, 0.566412535f, 0.566395292f,
        0.566378048f, 0.566360805f, 0.566343562f, 0.566326318f, 0.566309075f, 0.509678168f, 0.45304726f, 0.396416353f, 0.339785445f, 0.283154538f, 0.22652363f, 0.169892723f,
        0.113261815f, 0.056630908f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth33({0.1442082f, 0.149032892f, 0.153857583f, 0.158682275f, 0.163506967f, 0.168331658f, 0.17315635f, 0.177981042f, 0.182805733f, 0.187630425f,
        0.192455117f, 0.200393591f, 0.208332065f, 0.216270539f, 0.224209013f, 0.232147488f, 0.240085962f, 0.248024436f, 0.25596291f, 0.263901384f, 0.271839858f, 0.27773725f,
        0.283634642f, 0.289532033f, 0.295429425f, 0.301326817f, 0.307224208f, 0.3131216f, 0.319018992f, 0.324916383f, 0.330813775f, 0.33277435f, 0.334734925f, 0.3366955f,
        0.338656075f, 0.34061665f, 0.342577225f, 0.3445378f, 0.346498375f, 0.34845895f, 0.350419525f, 0.351069437f, 0.351719348f, 0.35236926f, 0.353019172f, 0.353669083f,
        0.354318995f, 0.354968907f, 0.355618818f, 0.35626873f, 0.356918642f, 0.357349918f, 0.357781193f, 0.358212469f, 0.358643745f, 0.359075021f, 0.359506297f, 0.359937573f,
        0.360368848f, 0.360800124f, 0.3612314f, 0.361434438f, 0.361637477f, 0.361840515f, 0.362043553f, 0.362246592f, 0.36244963f, 0.362652668f, 0.362855707f, 0.363058745f,
        0.363261783f, 0.363232583f, 0.363203382f, 0.363174181f, 0.36314498f, 0.363115779f, 0.363086578f, 0.363057378f, 0.363028177f, 0.362998976f, 0.362969775f, 0.3627538f,
        0.362537825f, 0.36232185f, 0.362105875f, 0.3618899f, 0.361673925f, 0.36145795f, 0.361241975f, 0.361026f, 0.360810025f, 0.360602629f, 0.360395233f, 0.360187838f,
        0.359980442f, 0.359773046f, 0.35956565f, 0.359358254f, 0.359150858f, 0.358943463f, 0.358736067f, 0.35867327f, 0.358610473f, 0.358547677f, 0.35848488f, 0.358422083f,
        0.358359287f, 0.35829649f, 0.358233693f, 0.358170897f, 0.3581081f, 0.358142819f, 0.358177538f, 0.358212258f, 0.358246977f, 0.358281696f, 0.358316415f, 0.358351134f,
        0.358385853f, 0.358420573f, 0.358455292f, 0.358528363f, 0.358601435f, 0.358674507f, 0.358747578f, 0.35882065f, 0.358893722f, 0.358966793f, 0.359039865f, 0.359112937f,
        0.359186008f, 0.359308604f, 0.3594312f, 0.359553796f, 0.359676392f, 0.359798988f, 0.359921583f, 0.360044179f, 0.360166775f, 0.360289371f, 0.360411967f, 0.360416353f,
        0.360420738f, 0.360425124f, 0.36042951f, 0.360433896f, 0.360438282f, 0.360442668f, 0.360447053f, 0.360451439f, 0.360455825f, 0.360466266f, 0.360476707f, 0.360487148f,
        0.360497588f, 0.360508029f, 0.36051847f, 0.360528911f, 0.360539352f, 0.360549793f, 0.360560233f, 0.360587145f, 0.360614057f, 0.360640968f, 0.36066788f, 0.360694792f,
        0.360721703f, 0.360748615f, 0.360775527f, 0.360802438f, 0.36082935f, 0.360789041f, 0.360748732f, 0.360708423f, 0.360668113f, 0.360627804f, 0.360587495f, 0.360547186f,
        0.360506877f, 0.360466568f, 0.360426258f, 0.360568448f, 0.360710637f, 0.360852826f, 0.360995015f, 0.361137204f, 0.361279393f, 0.361421583f, 0.361563772f, 0.361705961f,
        0.36184815f, 0.361838397f, 0.361828643f, 0.36181889f, 0.361809137f, 0.361799383f, 0.36178963f, 0.361779877f, 0.361770123f, 0.36176037f, 0.361750617f, 0.361708134f,
        0.361665652f, 0.361623169f, 0.361580687f, 0.361538204f, 0.361495722f, 0.361453239f, 0.361410757f, 0.361368274f, 0.361325792f, 0.361125494f, 0.360925197f, 0.360724899f,
        0.360524602f, 0.360324304f, 0.360124007f, 0.359923709f, 0.359723412f, 0.359523114f, 0.359322817f, 0.359143894f, 0.358964972f, 0.358786049f, 0.358607127f, 0.358428204f,
        0.358249282f, 0.358070359f, 0.357891437f, 0.357712514f, 0.357533592f, 0.357323623f, 0.357113653f, 0.356903684f, 0.356693715f, 0.356483746f, 0.356273777f, 0.356063808f,
        0.355853838f, 0.355643869f, 0.3554339f, 0.355132009f, 0.354830118f, 0.354528228f, 0.354226337f, 0.353924446f, 0.353622555f, 0.353320664f, 0.353018773f, 0.352716883f,
        0.352414992f, 0.352163242f, 0.351911492f, 0.351659742f, 0.351407992f, 0.351156242f, 0.350904492f, 0.350652742f, 0.350400992f, 0.350149242f, 0.349897492f, 0.349671781f,
        0.34944607f, 0.349220359f, 0.348994648f, 0.348768938f, 0.348543227f, 0.348317516f, 0.348091805f, 0.347866094f, 0.347640383f, 0.347423158f, 0.347205932f, 0.346988706f,
        0.34677148f, 0.346554254f, 0.346337028f, 0.346119803f, 0.345902577f, 0.345685351f, 0.345468125f, 0.345193043f, 0.344917962f, 0.34464288f, 0.344367798f, 0.344092717f,
        0.343817635f, 0.343542553f, 0.343267472f, 0.34299239f, 0.342717308f, 0.342462254f, 0.3422072f, 0.341952146f, 0.341697092f, 0.341442038f, 0.341186983f, 0.340931929f,
        0.340676875f, 0.340421821f, 0.340166767f, 0.33991013f, 0.339653493f, 0.339396857f, 0.33914022f, 0.338883583f, 0.338626947f, 0.33837031f, 0.338113673f, 0.337857037f,
        0.3376004f, 0.337371173f, 0.337141945f, 0.336912718f, 0.33668349f, 0.336454263f, 0.336225035f, 0.335995808f, 0.33576658f, 0.335537353f, 0.335308125f, 0.335160589f,
        0.335013053f, 0.334865518f, 0.334717982f, 0.334570446f, 0.33442291f, 0.334275374f, 0.334127838f, 0.333980303f, 0.333832767f, 0.333629798f, 0.33342683f, 0.333223862f,
        0.333020893f, 0.332817925f, 0.332614957f, 0.332411988f, 0.33220902f, 0.332006052f, 0.331803083f, 0.331677119f, 0.331551155f, 0.331425191f, 0.331299227f, 0.331173263f,
        0.331047298f, 0.330921334f, 0.33079537f, 0.330669406f, 0.330543442f, 0.297489098f, 0.264434753f, 0.231380409f, 0.198326065f, 0.165271721f, 0.132217377f, 0.099163033f,
        0.066108688f, 0.033054344f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth34({0.105192122f, 0.107805515f, 0.110418907f, 0.1130323f, 0.115645693f, 0.118259086f, 0.120872479f, 0.123485872f, 0.126099264f, 0.128712657f,
        0.13132605f, 0.134453238f, 0.137580427f, 0.140707615f, 0.143834803f, 0.146961992f, 0.15008918f, 0.153216368f, 0.156343557f, 0.159470745f, 0.162597933f, 0.164354708f,
        0.166111483f, 0.167868258f, 0.169625033f, 0.171381808f, 0.173138583f, 0.174895358f, 0.176652133f, 0.178408908f, 0.180165683f, 0.180740801f, 0.181315918f, 0.181891036f,
        0.182466153f, 0.183041271f, 0.183616388f, 0.184191506f, 0.184766623f, 0.185341741f, 0.185916858f, 0.186278556f, 0.186640253f, 0.187001951f, 0.187363648f, 0.187725346f,
        0.188087043f, 0.188448741f, 0.188810438f, 0.189172136f, 0.189533833f, 0.189866533f, 0.190199233f, 0.190531933f, 0.190864633f, 0.191197333f, 0.191530033f, 0.191862733f,
        0.192195433f, 0.192528133f, 0.192860833f, 0.192997978f, 0.193135122f, 0.193272266f, 0.19340941f, 0.193546554f, 0.193683698f, 0.193820843f, 0.193957987f, 0.194095131f,
        0.194232275f, 0.194187385f, 0.194142495f, 0.194097605f, 0.194052715f, 0.194007825f, 0.193962935f, 0.193918045f, 0.193873155f, 0.193828265f, 0.193783375f, 0.193637849f,
        0.193492323f, 0.193346798f, 0.193201272f, 0.193055746f, 0.19291022f, 0.192764694f, 0.192619168f, 0.192473643f, 0.192328117f, 0.192201278f, 0.19207444f, 0.191947602f,
        0.191820763f, 0.191693925f, 0.191567087f, 0.191440248f, 0.19131341f, 0.191186572f, 0.191059733f, 0.191038389f, 0.191017045f, 0.190995701f, 0.190974357f, 0.190953013f,
        0.190931668f, 0.190910324f, 0.19088898f, 0.190867636f, 0.190846292f, 0.190886298f, 0.190926305f, 0.190966312f, 0.191006318f, 0.191046325f, 0.191086332f, 0.191126338f,
        0.191166345f, 0.191206352f, 0.191246358f, 0.191279655f, 0.191312952f, 0.191346248f, 0.191379545f, 0.191412842f, 0.191446138f, 0.191479435f, 0.191512732f, 0.191546028f,
        0.191579325f, 0.191626543f, 0.19167376f, 0.191720978f, 0.191768195f, 0.191815413f, 0.19186263f, 0.191909848f, 0.191957065f, 0.192004283f, 0.1920515f, 0.192056469f,
        0.192061438f, 0.192066408f, 0.192071377f, 0.192076346f, 0.192081315f, 0.192086284f, 0.192091253f, 0.192096223f, 0.192101192f, 0.192108111f, 0.19211503f, 0.192121949f,
        0.192128868f, 0.192135788f, 0.192142707f, 0.192149626f, 0.192156545f, 0.192163464f, 0.192170383f, 0.192184616f, 0.192198848f, 0.192213081f, 0.192227313f, 0.192241546f,
        0.192255778f, 0.192270011f, 0.192284243f, 0.192298476f, 0.192312708f, 0.192287254f, 0.1922618f, 0.192236346f, 0.192210892f, 0.192185438f, 0.192159983f, 0.192134529f,
        0.192109075f, 0.192083621f, 0.192058167f, 0.192115572f, 0.192172977f, 0.192230382f, 0.192287787f, 0.192345192f, 0.192402597f, 0.192460002f, 0.192517407f, 0.192574812f,
        0.192632217f, 0.192607146f, 0.192582075f, 0.192557004f, 0.192531933f, 0.192506863f, 0.192481792f, 0.192456721f, 0.19243165f, 0.192406579f, 0.192381508f, 0.192336265f,
        0.192291022f, 0.192245778f, 0.192200535f, 0.192155292f, 0.192110048f, 0.192064805f, 0.192019562f, 0.191974318f, 0.191929075f, 0.191795643f, 0.19166221f, 0.191528778f,
        0.191395345f, 0.191261913f, 0.19112848f, 0.190995048f, 0.190861615f, 0.190728183f, 0.19059475f, 0.190473893f, 0.190353037f, 0.19023218f, 0.190111323f, 0.189990467f,
        0.18986961f, 0.189748753f, 0.189627897f, 0.18950704f, 0.189386183f, 0.189246166f, 0.189106148f, 0.188966131f, 0.188826113f, 0.188686096f, 0.188546078f, 0.188406061f,
        0.188266043f, 0.188126026f, 0.187986008f, 0.187774149f, 0.18756229f, 0.187350431f, 0.187138572f, 0.186926713f, 0.186714853f, 0.186502994f, 0.186291135f, 0.186079276f,
        0.185867417f, 0.185678197f, 0.185488977f, 0.185299757f, 0.185110537f, 0.184921317f, 0.184732097f, 0.184542877f, 0.184353657f, 0.184164437f, 0.183975217f, 0.183809435f,
        0.183643653f, 0.183477872f, 0.18331209f, 0.183146308f, 0.182980527f, 0.182814745f, 0.182648963f, 0.182483182f, 0.1823174f, 0.182170823f, 0.182024245f, 0.181877668f,
        0.18173109f, 0.181584513f, 0.181437935f, 0.181291358f, 0.18114478f, 0.180998203f, 0.180851625f, 0.180691691f, 0.180531757f, 0.180371823f, 0.180211888f, 0.180051954f,
        0.17989202f, 0.179732086f, 0.179572152f, 0.179412218f, 0.179252283f, 0.179104241f, 0.178956198f, 0.178808156f, 0.178660113f, 0.178512071f, 0.178364028f, 0.178215986f,
        0.178067943f, 0.177919901f, 0.177771858f, 0.177599149f, 0.17742644f, 0.177253731f, 0.177081022f, 0.176908313f, 0.176735603f, 0.176562894f, 0.176390185f, 0.176217476f,
        0.176044767f, 0.175874165f, 0.175703563f, 0.175532962f, 0.17536236f, 0.175191758f, 0.175021157f, 0.174850555f, 0.174679953f, 0.174509352f, 0.17433875f, 0.174241704f,
        0.174144658f, 0.174047613f, 0.173950567f, 0.173853521f, 0.173756475f, 0.173659429f, 0.173562383f, 0.173465338f, 0.173368292f, 0.173250319f, 0.173132347f, 0.173014374f,
        0.172896402f, 0.172778429f, 0.172660457f, 0.172542484f, 0.172424512f, 0.172306539f, 0.172188567f, 0.172108589f, 0.172028612f, 0.171948634f, 0.171868657f, 0.171788679f,
        0.171708702f, 0.171628724f, 0.171548747f, 0.171468769f, 0.171388792f, 0.154249913f, 0.137111033f, 0.119972154f, 0.102833275f, 0.085694396f, 0.068555517f, 0.051416638f,
        0.034277758f, 0.017138879f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth35({0.067959987f, 0.068836284f, 0.069712581f, 0.070588878f, 0.071465175f, 0.072341473f, 0.07321777f, 0.074094067f, 0.074970364f, 0.075846661f,
        0.076722958f, 0.07743884f, 0.078154722f, 0.078870604f, 0.079586486f, 0.080302368f, 0.081018249f, 0.081734131f, 0.082450013f, 0.083165895f, 0.083881777f, 0.084234856f,
        0.084587934f, 0.084941013f, 0.085294092f, 0.085647171f, 0.08600025f, 0.086353329f, 0.086706407f, 0.087059486f, 0.087412565f, 0.087559802f, 0.087707038f, 0.087854275f,
        0.088001511f, 0.088148748f, 0.088295984f, 0.088443221f, 0.088590457f, 0.088737694f, 0.08888493f, 0.089040721f, 0.089196512f, 0.089352304f, 0.089508095f, 0.089663886f,
        0.089819677f, 0.089975468f, 0.090131259f, 0.090287051f, 0.090442842f, 0.090585145f, 0.090727448f, 0.090869752f, 0.091012055f, 0.091154358f, 0.091296662f, 0.091438965f,
        0.091581268f, 0.091723572f, 0.091865875f, 0.091883666f, 0.091901456f, 0.091919247f, 0.091937038f, 0.091954828f, 0.091972619f, 0.09199041f, 0.0920082f, 0.092025991f,
        0.092043782f, 0.091974486f, 0.09190519f, 0.091835894f, 0.091766598f, 0.091697303f, 0.091628007f, 0.091558711f, 0.091489415f, 0.091420119f, 0.091350823f, 0.091255008f,
        0.091159193f, 0.091063378f, 0.090967563f, 0.090871748f, 0.090775932f, 0.090680117f, 0.090584302f, 0.090488487f, 0.090392672f, 0.090328126f, 0.090263581f, 0.090199035f,
        0.09013449f, 0.090069944f, 0.090005399f, 0.089940853f, 0.089876308f, 0.089811762f, 0.089747217f, 0.089737453f, 0.089727689f, 0.089717925f, 0.089708161f, 0.089698397f,
        0.089688633f, 0.089678869f, 0.089669105f, 0.089659341f, 0.089649577f, 0.08966585f, 0.089682123f, 0.089698396f, 0.089714669f, 0.089730943f, 0.089747216f, 0.089763489f,
        0.089779762f, 0.089796035f, 0.089812308f, 0.089818988f, 0.089825667f, 0.089832346f, 0.089839025f, 0.089845704f, 0.089852383f, 0.089859063f, 0.089865742f, 0.089872421f,
        0.0898791f, 0.089888167f, 0.089897234f, 0.089906301f, 0.089915367f, 0.089924434f, 0.089933501f, 0.089942568f, 0.089951635f, 0.089960702f, 0.089969768f, 0.089969133f,
        0.089968498f, 0.089967863f, 0.089967228f, 0.089966593f, 0.089965958f, 0.089965323f, 0.089964688f, 0.089964053f, 0.089963418f, 0.089968449f, 0.08997348f, 0.089978511f,
        0.089983542f, 0.089988573f, 0.089993604f, 0.089998635f, 0.090003666f, 0.090008697f, 0.090013728f, 0.090018559f, 0.09002339f, 0.09002822f, 0.090033051f, 0.090037882f,
        0.090042712f, 0.090047543f, 0.090052374f, 0.090057204f, 0.090062035f, 0.090042038f, 0.090022041f, 0.090002044f, 0.089982047f, 0.08996205f, 0.089942053f, 0.089922056f,
        0.089902059f, 0.089882062f, 0.089862065f, 0.089866237f, 0.089870408f, 0.08987458f, 0.089878752f, 0.089882923f, 0.089887095f, 0.089891267f, 0.089895438f, 0.08989961f,
        0.089903782f, 0.089864664f, 0.089825545f, 0.089786427f, 0.089747309f, 0.089708191f, 0.089669073f, 0.089629955f, 0.089590836f, 0.089551718f, 0.0895126f, 0.089477832f,
        0.089443064f, 0.089408296f, 0.089373528f, 0.08933876f, 0.089303992f, 0.089269224f, 0.089234456f, 0.089199688f, 0.08916492f, 0.089100586f, 0.089036251f, 0.088971917f,
        0.088907582f, 0.088843248f, 0.088778913f, 0.088714579f, 0.088650244f, 0.08858591f, 0.088521575f, 0.08847569f, 0.088429804f, 0.088383919f, 0.088338033f, 0.088292148f,
        0.088246262f, 0.088200377f, 0.088154491f, 0.088108606f, 0.08806272f, 0.088005834f, 0.087948947f, 0.087892061f, 0.087835175f, 0.087778288f, 0.087721402f, 0.087664516f,
        0.087607629f, 0.087550743f, 0.087493857f, 0.087389368f, 0.087284879f, 0.08718039f, 0.087075901f, 0.086971413f, 0.086866924f, 0.086762435f, 0.086657946f, 0.086553457f,
        0.086448968f, 0.086355073f, 0.086261178f, 0.086167283f, 0.086073388f, 0.085979493f, 0.085885597f, 0.085791702f, 0.085697807f, 0.085603912f, 0.085510017f, 0.085436835f,
        0.085363653f, 0.085290471f, 0.085217289f, 0.085144108f, 0.085070926f, 0.084997744f, 0.084924562f, 0.08485138f, 0.084778198f, 0.084720503f, 0.084662807f, 0.084605112f,
        0.084547416f, 0.084489721f, 0.084432025f, 0.08437433f, 0.084316634f, 0.084258939f, 0.084201243f, 0.084149172f, 0.0840971f, 0.084045029f, 0.083992957f, 0.083940886f,
        0.083888814f, 0.083836743f, 0.083784671f, 0.0837326f, 0.083680528f, 0.083633793f, 0.083587058f, 0.083540322f, 0.083493587f, 0.083446852f, 0.083400116f, 0.083353381f,
        0.083306646f, 0.08325991f, 0.083213175f, 0.08314572f, 0.083078266f, 0.083010811f, 0.082943356f, 0.082875902f, 0.082808447f, 0.082740992f, 0.082673538f, 0.082606083f,
        0.082538628f, 0.082461047f, 0.082383465f, 0.082305883f, 0.082228302f, 0.08215072f, 0.082073138f, 0.081995557f, 0.081917975f, 0.081840393f, 0.081762812f, 0.081729338f,
        0.081695864f, 0.08166239f, 0.081628916f, 0.081595443f, 0.081561969f, 0.081528495f, 0.081495021f, 0.081461547f, 0.081428073f, 0.081378168f, 0.081328263f, 0.081278358f,
        0.081228453f, 0.081178548f, 0.081128642f, 0.081078737f, 0.081028832f, 0.080978927f, 0.080929022f, 0.080903787f, 0.080878553f, 0.080853319f, 0.080828084f, 0.08080285f,
        0.080777616f, 0.080752381f, 0.080727147f, 0.080701913f, 0.080676678f, 0.072609011f, 0.064541343f, 0.056473675f, 0.048406007f, 0.040338339f, 0.032270671f, 0.024203004f,
        0.016135336f, 0.008067668f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    static constexpr Spectrum Macbeth36({0.031023515f, 0.031120073f, 0.031216632f, 0.03131319f, 0.031409748f, 0.031506307f, 0.031602865f, 0.031699423f, 0.031795982f, 0.03189254f,
        0.031989098f, 0.032018177f, 0.032047256f, 0.032076335f, 0.032105414f, 0.032134493f, 0.032163572f, 0.032192651f, 0.03222173f, 0.032250809f,
        0.032279888f, 0.0323079f, 0.032335911f, 0.032363923f, 0.032391934f, 0.032419946f, 0.032447957f, 0.032475969f, 0.03250398f, 0.032531992f,
        0.032560003f, 0.032579362f, 0.032598721f, 0.03261808f, 0.032637439f, 0.032656798f, 0.032676157f, 0.032695516f, 0.032714875f, 0.032734234f,
        0.032753593f, 0.032760412f, 0.03276723f, 0.032774048f, 0.032780866f, 0.032787684f, 0.032794502f, 0.032801321f, 0.032808139f, 0.032814957f,
        0.032821775f, 0.032821853f, 0.03282193f, 0.032822008f, 0.032822085f, 0.032822163f, 0.03282224f, 0.032822318f, 0.032822395f, 0.032822473f,
        0.03282255f, 0.032802458f, 0.032782366f, 0.032762274f, 0.032742181f, 0.032722089f, 0.032701997f, 0.032681905f, 0.032661813f, 0.032641721f,
        0.032621628f, 0.03260782f, 0.032594011f, 0.032580202f, 0.032566393f, 0.032552584f, 0.032538775f, 0.032524967f, 0.032511158f, 0.032497349f,
        0.03248354f, 0.032475542f, 0.032467544f, 0.032459547f, 0.032451549f, 0.032443551f, 0.032435553f, 0.032427555f, 0.032419557f, 0.03241156f,
        0.032403562f, 0.032395802f, 0.032388042f, 0.032380282f, 0.032372522f, 0.032364762f, 0.032357002f, 0.032349242f, 0.032341482f, 0.032333722f,
        0.032325962f, 0.03232411f, 0.032322258f, 0.032320407f, 0.032318555f, 0.032316703f, 0.032314852f, 0.032313f, 0.032311148f, 0.032309297f,
        0.032307445f, 0.032304859f, 0.032302274f, 0.032299688f, 0.032297102f, 0.032294517f, 0.032291931f, 0.032289345f, 0.03228676f, 0.032284174f,
        0.032281588f, 0.032268883f, 0.032256177f, 0.032243471f, 0.032230766f, 0.03221806f, 0.032205354f, 0.032192649f, 0.032179943f, 0.032167237f,
        0.032154532f, 0.032147745f, 0.032140959f, 0.032134173f, 0.032127386f, 0.0321206f, 0.032113814f, 0.032107027f, 0.032100241f, 0.032093455f,
        0.032086668f, 0.032076223f, 0.032065778f, 0.032055333f, 0.032044888f, 0.032034443f, 0.032023997f, 0.032013552f, 0.032003107f, 0.031992662f,
        0.031982217f, 0.031980094f, 0.03197797f, 0.031975847f, 0.031973724f, 0.031971601f, 0.031969478f, 0.031967355f, 0.031965231f, 0.031963108f,
        0.031960985f, 0.031962523f, 0.03196406f, 0.031965598f, 0.031967136f, 0.031968673f, 0.031970211f, 0.031971749f, 0.031973286f, 0.031974824f,
        0.031976362f, 0.03197069f, 0.031965018f, 0.031959346f, 0.031953674f, 0.031948003f, 0.031942331f, 0.031936659f, 0.031930987f, 0.031925315f,
        0.031919643f, 0.031922461f, 0.031925278f, 0.031928095f, 0.031930913f, 0.03193373f, 0.031936547f, 0.031939365f, 0.031942182f, 0.031944999f,
        0.031947817f, 0.031934769f, 0.031921721f, 0.031908674f, 0.031895626f, 0.031882578f, 0.031869531f, 0.031856483f, 0.031843435f, 0.031830388f,
        0.03181734f, 0.031814895f, 0.03181245f, 0.031810006f, 0.031807561f, 0.031805116f, 0.031802671f, 0.031800226f, 0.031797781f, 0.031795337f,
        0.031792892f, 0.031786646f, 0.031780399f, 0.031774153f, 0.031767907f, 0.031761661f, 0.031755415f, 0.031749169f, 0.031742922f, 0.031736676f,
        0.03173043f, 0.031742255f, 0.031754081f, 0.031765906f, 0.031777731f, 0.031789557f, 0.031801382f, 0.031813207f, 0.031825033f, 0.031836858f,
        0.031848683f, 0.031857231f, 0.031865779f, 0.031874327f, 0.031882875f, 0.031891423f, 0.031899971f, 0.031908519f, 0.031917067f, 0.031925615f,
        0.031934163f, 0.031931159f, 0.031928154f, 0.03192515f, 0.031922145f, 0.031919141f, 0.031916136f, 0.031913132f, 0.031910127f, 0.031907123f,
        0.031904118f, 0.03190447f, 0.031904822f, 0.031905174f, 0.031905526f, 0.031905878f, 0.03190623f, 0.031906582f, 0.031906934f, 0.031907286f,
        0.031907638f, 0.031910956f, 0.031914274f, 0.031917592f, 0.03192091f, 0.031924228f, 0.031927545f, 0.031930863f, 0.031934181f, 0.031937499f,
        0.031940817f, 0.031948264f, 0.031955711f, 0.031963158f, 0.031970605f, 0.031978053f, 0.0319855f, 0.031992947f, 0.032000394f, 0.032007841f,
        0.032015288f, 0.032025932f, 0.032036576f, 0.032047219f, 0.032057863f, 0.032068507f, 0.03207915f, 0.032089794f, 0.032100438f, 0.032111081f,
        0.032121725f, 0.032131217f, 0.032140709f, 0.032150201f, 0.032159692f, 0.032169184f, 0.032178676f, 0.032188168f, 0.03219766f, 0.032207152f,
        0.032216643f, 0.032217578f, 0.032218512f, 0.032219446f, 0.03222038f, 0.032221314f, 0.032222248f, 0.032223183f, 0.032224117f, 0.032225051f,
        0.032225985f, 0.032227669f, 0.032229353f, 0.032231037f, 0.032232721f, 0.032234405f, 0.032236089f, 0.032237773f, 0.032239457f, 0.032241141f,
        0.032242825f, 0.032251451f, 0.032260077f, 0.032268704f, 0.03227733f, 0.032285956f, 0.032294582f, 0.032303208f, 0.032311834f, 0.032320461f,
        0.032329087f, 0.032333441f, 0.032337796f, 0.032342151f, 0.032346505f, 0.03235086f, 0.032355215f, 0.032359569f, 0.032363924f, 0.032368279f,
        0.032372633f, 0.032385684f, 0.032398735f, 0.032411786f, 0.032424837f, 0.032437888f, 0.032450939f, 0.03246399f, 0.032477041f, 0.032490092f,
        0.032503143f, 0.029252829f, 0.026002515f, 0.0227522f, 0.019501886f, 0.016251572f, 0.013001257f, 0.009750943f, 0.006500629f, 0.003250314f,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    // Brown skin, light skin, sky, folliage, ...
    static constexpr std::array<Spectrum, 24> Patch{
        Macbeth01, Macbeth02, Macbeth03, Macbeth04, Macbeth05, Macbeth06,
        Macbeth11, Macbeth12, Macbeth13, Macbeth14, Macbeth15, Macbeth16,
        Macbeth21, Macbeth22, Macbeth23, Macbeth24, Macbeth25, Macbeth26,
        Macbeth31, Macbeth32, Macbeth33, Macbeth34, Macbeth35, Macbeth36,
    };

} // namespace Macbeth

// Standard observers
static constexpr Observer CIE1931(CIE1931_X, CIE1931_Y, CIE1931_Z);
static constexpr Observer CIE2012(CIE2012_X, CIE2012_Y, CIE2012_Z);

//standard illuminants
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

//xR,yR,xG,yG,xB,yB,xW,yW
static constexpr Gamut AdobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
static constexpr Gamut Rec709(0.64f, 0.33f, 0.30f, 0.60f, 0.15f, 0.06f, 0.3127f, 0.3290f);
static constexpr Gamut Rec2020(0.708f, 0.292f, 0.17f, 0.797f, 0.131f, 0.046f, 0.3127f, 0.3290f);
static constexpr Gamut DCI_P3(0.68f, 0.32f, 0.265f, 0.69f, 0.15f, 0.06f, 0.314f, 0.351f);
static constexpr Gamut S_Gamut(0.73f, 0.28f, 0.14f, 0.855f, 0.10f, -0.05f, 0.3127f, 0.3290f);
static constexpr Gamut S_Gamut3_Cine(0.766f, 0.275f, 0.225f, 0.800f, 0.089f, -0.087f, 0.3127f, 0.3290f);
static constexpr Gamut ACEScg(0.713f, 0.293f, 0.165f, 0.830f, 0.128f, 0.044f, 0.32168f, 0.33767f);
static constexpr Gamut ACES2065(0.73470f, 0.26530f, 0.f, 1.f, 0.0001f, -0.077f, 0.32168f, 0.33767f);
static constexpr Gamut LMS(Matrix3(0.8951f, 0.2664f, -0.1614f, -0.7502f, 1.7135f, 0.0367f, 0.0389f, -0.0685f, 1.0296f)); // fromXYZ matrix.
static constexpr Gamut XYZ(Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1));

// returns Gamut convert matrix
static constexpr Matrix3 GamutConvert(const Gamut &src, const Gamut &dst)
{
    return dst.fromXYZ().mul(src.toXYZ());
}

// returns Bradford adaptation matrix
static constexpr Matrix3 Bradford(const Tristimulus &white_src, const Tristimulus &white_dst)
{
    const Tristimulus &lms_src(LMS.fromXYZ(white_src));
    const Tristimulus &lms_dst(LMS.fromXYZ(white_dst));
    const Matrix3      scale(
        Matrix3::diag(
            Vector3(
                lms_dst[0] / lms_src[0],
                lms_dst[1] / lms_src[1],
                lms_dst[2] / lms_src[2])));
    return LMS.toXYZ().mul(scale).mul(LMS.fromXYZ());
}

// cubic spline approx https://en.wikipedia.org/wiki/Planckian_locus#Approximation
static constexpr Tristimulus PlanckianLocus(const float &T, const float &Y = 1.f)
{
    const float x = (T < 4000.f) ? ((-0.2661239e9f / (T * T * T)) - (0.2343580e6f / (T * T)) + (0.8776956e3f / T) + 0.179910f)
                                 : ((-3.0258469e9f / (T * T * T)) + (2.1070379e6f / (T * T)) + (0.2226347e3f / T) + 0.240390f);
    const float y = (T < 2222.f) ? ((-1.1063814f * x * x * x) - (1.34811020f * x * x) + (2.18555832f * x) - 0.20219683f)
                                 : ((T < 4000.f) ? ((-1.1063814f * x * x * x - 1.34811020f * x * x + 2.18555832f * x - 0.20219683f))
                                                 : ((3.0817580f * x * x * x - 5.87338670f * x * x + 3.75112997f * x - 0.37001483f)));
    return Tristimulus(Y, x, y).fromYxy();
}

} // namespace ColorSystem
