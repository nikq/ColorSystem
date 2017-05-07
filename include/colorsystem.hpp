
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
//
// policy: currently single-headered, simple, low overhead, no ICC support.
//

#include <array>
#include <math.h>
#include <stdio.h>
#include <assert.h>

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
    constexpr auto                        size() const { return v_.size(); }
    auto                                  begin() const { return v_.begin(); }
    auto                                  end() const { return v_.end(); }
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

    static constexpr Tristimulus toHSV(const Tristimulus &t)
    {
        const float max = maxi(maxi(t[0], t[1]), t[2]);
        const float min = mini(mini(t[0], t[1]), t[2]);
        return Tristimulus(
            ((max == min) ? 0.f : ((max == t[0]) ? (60.f * (t[1] - t[2]) / (max - min)) : ((max == t[1]) ? (60.f * (t[2] - t[0]) / (max - min) + 120.f) : (60.f * (t[0] - t[1]) / (max - min) + 240.f)))),
            (max == 0.f) ? 0.f : (max - min) / max, max);
    }
    constexpr Tristimulus  toHSV(void) const { return toHSV(*this); }
    static constexpr float mod360(const float &r)
    {
        return (r < 0.f) ? mod360(r + 360.f) : ((r > 360.f) ? mod360(r - 360.f) : r);
    }
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

    constexpr Spectrum() { ; }
    constexpr Spectrum(spectrum &s) : s_(s) { ; }
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
        return (index >= count) ? sample[count] : (index <= 0) ? sample[0] : lerp(sample[index], sample[index + 1], (t - l) / (h - l));
    }

    Spectrum(
        const float *sample,
        const int    samples = 400,
        const float &lo      = 380.f,
        const float &hi      = 780.f)
    {
        for (int i = 0; i < 400; i++)
            s_[i]  = fetch(sample, samples, lo, hi, 380.f + i);
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
    static Spectrum blackbody(const float temp)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = (float)planck(temp, (double)(380.f+i) * 1e-9);
        }
        return Spectrum(s);
    }

    static Spectrum E(const float e = 1.f)
    {
        spectrum s;
        for (int i = 0; i < 400; i++)
        {
            s[i] = e;
        }
        return Spectrum(s);
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

inline void dump(const Matrix3 &m)
{
    printf("-------\n");
    printf("%f,%f,%f\n", m[0], m[1], m[2]);
    printf("%f,%f,%f\n", m[3], m[4], m[5]);
    printf("%f,%f,%f\n", m[6], m[7], m[8]);
}

inline void dump(const Vector3 &v)
{
    printf("v-------\n");
    printf("%f,%f,%f\n", v[0], v[1], v[2]);
}

inline void dump(const Spectrum &s)
{
    for (int i = 0; i < 400; i++)
    {
        printf("%3d:%f\n", i, s[i]);
    }
}

static constexpr Spectrum CIE1931_X(Spectrum::spectrum{
    0.00136800000f, 0.00150205000f, 0.00164232800f, 0.00180238200f, 0.00199575700f, 0.00223600000f, 0.00253538500f, 0.00289260300f, 0.00330082900f, 0.00375323600f,
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
static constexpr Spectrum CIE1931_Y(Spectrum::spectrum{
    0.00003900000f, 0.00004282640f, 0.00004691460f, 0.00005158960f, 0.00005717640f, 0.00006400000f, 0.00007234421f, 0.00008221224f, 0.00009350816f, 0.00010613610f,
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
static constexpr Spectrum CIE1931_Z(Spectrum::spectrum{
    0.00645000100f, 0.00708321600f, 0.00774548800f, 0.00850115200f, 0.00941454400f, 0.01054999000f, 0.01196580000f, 0.01365587000f, 0.01558805000f, 0.01773015000f,
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

static constexpr Tristimulus CIEXYZ1931(const Spectrum &s)
{
    return Tristimulus(
        Spectrum::dot(CIE1931_X, s),
        Spectrum::dot(CIE1931_Y, s),
        Spectrum::dot(CIE1931_Z, s));
}

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
static constexpr Gamut LMS(Matrix3(0.8951f, 0.2664f, -0.1614f, -0.7502f, 1.7135f, 0.0367f, 0.0389f, -0.0685f, 1.0296f));

static constexpr Matrix3 GamutConvert(const Gamut &src, const Gamut &dst)
{
    return dst.fromXYZ().mul(src.toXYZ());
}

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

} // namespace ColorSystem
