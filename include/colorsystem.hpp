
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
//
// policy: currently single-headered, simple, low overhead, no ICC support.
//

#include <array>
#include <stdio.h>
#include <math.h>
#include <tuple>

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

    constexpr auto size() const { return v_.size(); }
    auto begin() const { return v_.begin(); }
    auto end() const { return v_.end(); }
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
    constexpr Matrix3 diag(const Vector3 &v) const
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
    auto begin() const { return m_.begin(); }
    auto end() const { return m_.end(); }
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
    constexpr auto size() const { return v_.size(); }
    auto begin() const { return v_.begin(); }
    auto end() const { return v_.end(); }

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
        const float threshold = 0.008856f;
        const float K         = 903.3f;
        return (f > 1.f) ? 1.f : ((f > threshold) ? powf(f, 1.f / 3.f) : ((K * f + 16.f) / 116.f));
    }
    static constexpr float CIELAB_decurve(const float &f)
    {
        const float K = (3.f / 29.f) * (3.f / 29.f) * (3.f / 29.f);
        return (f > 1.f) ? 1.f : ((f > 6.f / 29.f) ? pow(f, 3.f) : (116.f * f - 16.f) * K);
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
    constexpr Gamut(const Matrix3 &toXYZ) : toXYZ_(toXYZ), fromXYZ_(toXYZ.invert())
    {
        ;
    }
    constexpr Gamut(const float &xR, const float &yR, const float &xG, const float &yG, const float &xB, const float &yB, const float &xW, const float &yW)
        : toXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW)), fromXYZ_(fromPrimaries(xR, yR, xG, yG, xB, yB, xW, yW).invert())
    {
        ;
    }
    const Matrix3 &toXYZ(void) const { return toXYZ_; }
    const Matrix3 &fromXYZ(void) const { return fromXYZ_; }
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

    Spectrum()
    {
        ;
    }

    Spectrum(
        const float*sample,
        const int smaples,
        const float& lo,
        const float& hi)
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

//standard illuminants
static constexpr Tristimulus Illuminant_E(1.f, 1.f, 1.f);
static constexpr Tristimulus Illuminant_D50(0.9642f, 1.0f, 0.8249f);
static constexpr Tristimulus Illuminant_D65(0.95046f, 1.0f, 1.08906f);
static constexpr Tristimulus Illuminant_C(0.98071f, 1.0f, 1.1825f);

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

} // namespace ColorSystem
