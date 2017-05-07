
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

    constexpr const Vector3 &vec3(void) const { return v_; }
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

    constexpr bool isNegative(const float &a) const { return (a < 0.f); }
    constexpr bool hasNegative() const
    {
        return isNegative(v_[0]) || isNegative(v_[1]) || isNegative(v_[2]);
    }
    static constexpr bool hasNegative(const Tristimulus &t) { return t.hasNegative(); }

    static constexpr float z_from_xy(const float &x, const float &y) { return 1 - x - y; }
    static constexpr float X_from_Yxy(const float &Y, const float &x, const float &y) { return x * Y / y; }
    static constexpr float Y_from_Yxy(const float &Y, const float &x, const float &y) { return Y; }
    static constexpr float Z_from_Yxy(const float &Y, const float &x, const float &y) { return z_from_xy(x, y) * Y / y; }
    static constexpr Tristimulus fromYxy(const float &Y, const float &x, const float &y)
    {
        return Tristimulus(
            X_from_Yxy(Y, x, y),
            Y_from_Yxy(Y, x, y),
            Z_from_Yxy(Y, x, y));
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
        const float pq_C = 10000.0;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this assumes full range (0 - 1)
        float Np = powf(pixel, 1.0f / pq_m2);
        float L = Np - pq_c1;
        if (L < 0.0)
            L = 0.0;
        L = L / (pq_c2 - pq_c3 * Np);
        L = powf(L, 1.0f / pq_m1);
        return L * pq_C; // returns cd/m^2
    }
    static const float Y_to_ST2084(const float &nit) // nit should be 0-10000(cd/m^2)
    {
        const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
        const float pq_m2 = 78.84375;        // ( 2523.0 / 4096.0 ) * 128.0;
        const float pq_c1 = 0.8359375;       // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
        const float pq_c2 = 18.8515625;      // ( 2413.0 / 4096.0 ) * 32.0;
        const float pq_c3 = 18.6875;         // ( 2392.0 / 4096.0 ) * 32.0;
        const float pq_C = 10000.0;

        // Note that this does NOT handle any of the signal range
        // considerations from 2084 - this returns full range (0 - 1)
        float L = nit / pq_C;
        float Lm = powf(L, pq_m1);
        float N = (pq_c1 + pq_c2 * Lm) / (1.0f + pq_c3 * Lm);
        N = powf(N, pq_m2);
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
            return Tristimulus(
                Y_to_sRGB(scene[0]),
                Y_to_sRGB(scene[1]),
                Y_to_sRGB(scene[2]));
        }
        break;
        case BT709:
        {
            return Tristimulus(
                Y_to_BT709(scene[0]),
                Y_to_BT709(scene[1]),
                Y_to_BT709(scene[2]));
        }
        break;
        case ST2084:
        {
            return Tristimulus(
                Y_to_ST2084(scene[0]),
                Y_to_ST2084(scene[1]),
                Y_to_ST2084(scene[2]));
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
            return Tristimulus(
                sRGB_to_Y(screen[0]),
                sRGB_to_Y(screen[1]),
                sRGB_to_Y(screen[2]));
        }
        break;
        case BT709:
        {
            return Tristimulus(
                BT709_to_Y(screen[0]),
                BT709_to_Y(screen[1]),
                BT709_to_Y(screen[2]));
        }
        break;
        case ST2084:
        {
            return Tristimulus(
                ST2084_to_Y(screen[0]),
                ST2084_to_Y(screen[1]),
                ST2084_to_Y(screen[2]));
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
//
} // namespace ColorSystem
