
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
#include "common.hpp"
#include "colorsystem.hpp"

namespace
{
const float epsilon = 0.000001f;
}

TEST_CASE("toXYZ", "[XYZ]")
{
    const ColorSystem::Gamut adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);

    SECTION("adobe.toXYZ")
    {
        const ColorSystem::Matrix3 expected{0.576669f, 0.185558f, 0.188229f, 0.297345f, 0.627363f, 0.075291f, 0.027031f, 0.070689f, 0.991337f};
        auto const &xyz = adobeRGB.toXYZ();
        for (int_fast32_t i = 0; i < 9; ++i)
        {
            CAPTURE(i);
            REQUIRE(xyz[i] == Approx(expected[i]).epsilon(epsilon));
        }
    }
    SECTION("adobe.toXYZ (white)")
    {
        const ColorSystem::Tristimulus white(1, 1, 1);
        auto const &v3 = adobeRGB.toXYZ(white).vec3();
        REQUIRE(v3[0] == Approx(0.950456).epsilon(epsilon));
        REQUIRE(v3[1] == Approx(1.000000).epsilon(epsilon));
        REQUIRE(v3[2] == Approx(1.089058).epsilon(epsilon));
        //dump(v3) ;
    }
}

TEST_CASE("fromXYZ", "[XYZ]")
{
    ColorSystem::Gamut adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    SECTION("adobe.fromXYZ")
    {
        const ColorSystem::Matrix3 expected{2.041588f, -0.565007f, -0.344731f, -0.969244f, 1.875968f, 0.041555f, 0.013444f, -0.118362f, 1.015175f};
        auto const &m = adobeRGB.fromXYZ();
        for (int_fast32_t i = 0; i < 9; ++i)
        {
            CAPTURE(i);
            REQUIRE(m[i] == Approx(expected[i]).epsilon(epsilon));
        }
    }
    SECTION("adobe.fromXYZ (white)")
    {
        ColorSystem::Tristimulus white(1, 1, 1);
        auto const &v = adobeRGB.fromXYZ(white).vec3();
        REQUIRE(v[0] == Approx(1.131850).epsilon(epsilon));
        REQUIRE(v[1] == Approx(0.948279).epsilon(epsilon));
        REQUIRE(v[2] == Approx(0.910257).epsilon(epsilon));
    }
}

TEST_CASE("BT709")
{
    ColorSystem::Gamut adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    ColorSystem::Tristimulus white(1, 1, 1);
    auto const &monitor = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, adobeRGB.fromXYZ(white) * 100.f);
    // 1.062991,0.974048,0.954468
    REQUIRE(monitor[0] == Approx(1.062991).epsilon(epsilon));
    REQUIRE(monitor[1] == Approx(0.974048).epsilon(epsilon));
    REQUIRE(monitor[2] == Approx(0.954468).epsilon(epsilon));
}

TEST_CASE("toYxy")
{
    ColorSystem::Tristimulus X(1.f, 0.f, 0.f);
    ColorSystem::Tristimulus Y(0.f, 1.f, 0.f);
    ColorSystem::Tristimulus Z(0.f, 0.f, 1.f);
    ColorSystem::Tristimulus R(19.58f, 11.39f, 4.90f); // 5R

    auto const &X_Yxy = X.toYxy();
    auto const &Y_Yxy = Y.toYxy();
    auto const &Z_Yxy = Z.toYxy();
    auto const &R_Yxy = R.toYxy();
    REQUIRE(X_Yxy[0] == Approx(0.0000f).epsilon(epsilon));
    REQUIRE(X_Yxy[1] == Approx(0.3127f).epsilon(epsilon));
    REQUIRE(X_Yxy[2] == Approx(0.3290f).epsilon(epsilon));

    REQUIRE(Y_Yxy[0] == Approx(1.0000f).epsilon(epsilon));
    REQUIRE(Y_Yxy[1] == Approx(0.0000f).epsilon(epsilon));
    REQUIRE(Y_Yxy[2] == Approx(1.0000f).epsilon(epsilon));

    REQUIRE(Z_Yxy[0] == Approx(0.0000f).epsilon(epsilon));
    REQUIRE(Z_Yxy[1] == Approx(0.3127f).epsilon(epsilon));
    REQUIRE(Z_Yxy[2] == Approx(0.3290f).epsilon(epsilon));

    REQUIRE(R_Yxy[0] == Approx(11.39f).epsilon(epsilon));
    REQUIRE(R_Yxy[1] == Approx(0.54586f).epsilon(1e-5f));
    REQUIRE(R_Yxy[2] == Approx(0.31754f).epsilon(1e-5f));
}

TEST_CASE("fromYxy")
{
    ColorSystem::Tristimulus W(1.f, 0.3127f, 0.3290f);      //D65 w
    ColorSystem::Tristimulus R(11.39f, 0.54586f, 0.31754f); // 5R

    auto const &W_XYZ = W.fromYxy();
    auto const &R_XYZ = R.fromYxy();
    REQUIRE(W_XYZ[0] == Approx(0.95046f).epsilon(1e-5f));
    REQUIRE(W_XYZ[1] == Approx(1.00000f).epsilon(1e-5f));
    REQUIRE(W_XYZ[2] == Approx(1.08906f).epsilon(1e-5f));
    REQUIRE(R_XYZ[0] == Approx(19.58).epsilon(1e-3f));
    REQUIRE(R_XYZ[1] == Approx(11.39f).epsilon(1e-3f));
    REQUIRE(R_XYZ[2] == Approx(4.90f).epsilon(1e-3f));
}

TEST_CASE("ACES2065")
{
    const ColorSystem::Matrix3 &ACES2065_to_XYZ = ColorSystem::ACES2065.toXYZ();
    const ColorSystem::Matrix3 &ACES2065_from_XYZ = ColorSystem::ACES2065.fromXYZ();
    const float EPS=1e-5f;
    REQUIRE( ACES2065_to_XYZ[0] == Approx(0.9525523959f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[1] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[2] == Approx(0.0000936786f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[3] == Approx(0.3439664498f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[4] == Approx(0.7281660966f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[5] == Approx(-0.0721325464f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[6] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[7] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_to_XYZ[8] == Approx(1.0088251844f).epsilon(EPS));
    
    REQUIRE( ACES2065_from_XYZ[0] == Approx(1.0498110175f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[1] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[2] == Approx(-0.0000974845f).epsilon(EPS));

    REQUIRE( ACES2065_from_XYZ[3] == Approx(-0.4959030231f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[4] == Approx(1.3733130458f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[5] == Approx(0.0982400361f).epsilon(EPS));
    
    REQUIRE( ACES2065_from_XYZ[6] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[7] == Approx(0.0000000000f).epsilon(EPS));
    REQUIRE( ACES2065_from_XYZ[8] == Approx(0.9912520182f).epsilon(EPS));
}