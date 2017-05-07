
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
#include "colorsystem.hpp"
#include "TestUtilities.hpp"
#include "common.hpp"

namespace
{
const float epsilon = 0.000001f;
}

TEST_CASE("illuminants")
{
    const float                    EPS     = 1e-4f;
    const ColorSystem::Tristimulus d65_Yxy = ColorSystem::Illuminant_D65.toYxy();
    REQUIRE(d65_Yxy[1] == Approx(0.3127f).epsilon(EPS));
    REQUIRE(d65_Yxy[2] == Approx(0.3290f).epsilon(EPS));

    const ColorSystem::Tristimulus E_Yxy = ColorSystem::Illuminant_E.toYxy();
    REQUIRE(E_Yxy[1] == Approx(1.f / 3.f).epsilon(EPS));
    REQUIRE(E_Yxy[2] == Approx(1.f / 3.f).epsilon(EPS));

    const ColorSystem::Tristimulus d50_Yxy = ColorSystem::Illuminant_D50.toYxy();
    REQUIRE(d50_Yxy[1] == Approx(0.345703f).epsilon(EPS));
    REQUIRE(d50_Yxy[2] == Approx(0.358539f).epsilon(EPS));

    const ColorSystem::Tristimulus d60_Yxy = ColorSystem::Illuminant_D60.toYxy();
    REQUIRE(d60_Yxy[1] == Approx(0.32168f).epsilon(EPS));
    REQUIRE(d60_Yxy[2] == Approx(0.33767f).epsilon(EPS));
}

TEST_CASE("toXYZ", "[XYZ]")
{
    const ColorSystem::Gamut adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);

    SECTION("adobe.toXYZ")
    {
        const ColorSystem::Matrix3 expected{0.576669f, 0.185558f, 0.188229f, 0.297345f, 0.627363f, 0.075291f, 0.027031f, 0.070689f, 0.991337f};
        auto const &               xyz = adobeRGB.toXYZ();
        REQUIRE_THAT(xyz, IsApproxEquals(expected, epsilon));
    }
    SECTION("adobe.toXYZ (white)")
    {
        const ColorSystem::Tristimulus white(1, 1, 1);
        auto const &                   v3 = adobeRGB.toXYZ(white).vec3();
        REQUIRE_THAT(v3, IsApproxEquals(ColorSystem::Vector3{0.950456f, 1.000000f, 1.089058f}, epsilon));
    }
}

TEST_CASE("fromXYZ", "[XYZ]")
{
    ColorSystem::Gamut adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    SECTION("adobe.fromXYZ")
    {
        const ColorSystem::Matrix3 expected{2.041588f, -0.565007f, -0.344731f, -0.969244f, 1.875968f, 0.041555f, 0.013444f, -0.118362f, 1.015175f};
        auto const &               m = adobeRGB.fromXYZ();
        REQUIRE_THAT(m, IsApproxEquals(expected, epsilon));
    }
    SECTION("adobe.fromXYZ (white)")
    {
        ColorSystem::Tristimulus white(1, 1, 1);
        auto const &             v = adobeRGB.fromXYZ(white).vec3();
        REQUIRE_THAT(v, IsApproxEquals(ColorSystem::Vector3{1.131850f, 0.948279f, 0.910257f}, epsilon));
    }
}

TEST_CASE("BT709")
{
    ColorSystem::Gamut       adobeRGB(0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    ColorSystem::Tristimulus white(1, 1, 1);
    auto const &             monitor = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, adobeRGB.fromXYZ(white) * 100.f);
    // 1.062991,0.974048,0.954468
    REQUIRE_THAT(monitor, IsApproxEquals(ColorSystem::Tristimulus{1.062991f, 0.974048f, 0.954468f}, epsilon));
}

TEST_CASE("toYxy")
{
    SECTION("Convert from X")
    {
        ColorSystem::Tristimulus X(1.f, 0.f, 0.f);
        auto const &             X_Yxy = X.toYxy();
        REQUIRE_THAT(X_Yxy, IsApproxEquals(ColorSystem::Tristimulus{0.0000f, 0.3127f, 0.3290f}, epsilon));
    }
    SECTION("Convert from Y")
    {
        ColorSystem::Tristimulus Y(0.f, 1.f, 0.f);
        auto const &             Y_Yxy = Y.toYxy();
        REQUIRE_THAT(Y_Yxy, IsApproxEquals(ColorSystem::Tristimulus{1.0000f, 0.0000f, 1.0000f}, epsilon));
    }
    SECTION("Convert from Z")
    {
        ColorSystem::Tristimulus Z(0.f, 0.f, 1.f);
        auto const &             Z_Yxy = Z.toYxy();
        REQUIRE_THAT(Z_Yxy, IsApproxEquals(ColorSystem::Tristimulus{0.0000f, 0.3127f, 0.3290f}, epsilon));
    }
    SECTION("Convert from 5R")
    {
        ColorSystem::Tristimulus R(19.58f, 11.39f, 4.90f); // 5R
        auto const &             R_Yxy = R.toYxy();
        REQUIRE_THAT(R_Yxy, IsApproxEquals(ColorSystem::Tristimulus{11.39f, 0.54586f, 0.31754f}, 1e-5f));
    }
}

TEST_CASE("fromYxy")
{
    SECTION("Convert from D65 w")
    {
        ColorSystem::Tristimulus W(1.f, 0.3127f, 0.3290f); //D65 w
        auto const &             W_XYZ = W.fromYxy();

        REQUIRE_THAT(W_XYZ, IsApproxEquals(ColorSystem::Tristimulus{0.95046f, 1.00000f, 1.08906f}, 1e-5f));
    }
    SECTION("Convert from 5R")
    {
        ColorSystem::Tristimulus R(11.39f, 0.54586f, 0.31754f); // 5R
        auto const &             R_XYZ = R.fromYxy();

        REQUIRE_THAT(R_XYZ, IsApproxEquals(ColorSystem::Tristimulus{19.58f, 11.39f, 4.90f}, 1e-3f));
    }
}

TEST_CASE("ACES2065")
{
    // http://www.oscars.org/science-technology/aces/aces-documentation
    const float EPS = 1e-5f;
    SECTION("toXYZ")
    {
        const ColorSystem::Matrix3 &ACES2065_to_XYZ = ColorSystem::ACES2065.toXYZ();
        const ColorSystem::Matrix3  expected{0.9525523959f, 0.0000000000f, 0.0000936786f, 0.3439664498f, 0.7281660966f, -0.0721325464f, 0.0000000000f, 0.0000000000f, 1.0088251844f};
        REQUIRE_THAT(ACES2065_to_XYZ, IsApproxEquals(expected, EPS));
    }
    SECTION("fromXYZ")
    {
        const ColorSystem::Matrix3 &ACES2065_from_XYZ = ColorSystem::ACES2065.fromXYZ();
        const ColorSystem::Matrix3  expected{1.0498110175f, 0.0000000000f, -0.0000974845f, -0.4959030231f, 1.3733130458f, 0.0982400361f, 0.0000000000f, 0.0000000000f, 0.9912520182f};
        REQUIRE_THAT(ACES2065_from_XYZ, IsApproxEquals(expected, EPS));
    }
}

TEST_CASE("CIELAB")
{
    SECTION("toLAB")
    {
        const ColorSystem::Tristimulus XYZ1(0.4f, 0.5f, 0.6f);
        const ColorSystem::Tristimulus XYZ2(0.7f, 0.6f, 0.5f);
        const float                    EPS  = 1e-1f; // float is low precision... maybe?
        auto const &                   LAB1 = XYZ1.toCIELAB();
        auto const &                   LAB2 = XYZ2.toCIELAB();
        REQUIRE_THAT(LAB1, IsApproxEquals(ColorSystem::Tristimulus{76.0693f, -23.9455f, -21.1024f}, EPS));
        REQUIRE_THAT(LAB2, IsApproxEquals(ColorSystem::Tristimulus{81.8382f, 27.6605f, -0.5517f}, EPS));
    }
    SECTION("delta")
    {
        //http://qiita.com/tibigame/items/40ab217c863a20cdb264
        const float                    EPS = 1e-6f;
        const ColorSystem::Tristimulus LAB1(50.f, 2.6772f, -79.7751f);
        const ColorSystem::Tristimulus LAB2(50.f, 0.0000f, -82.7485f);
        const float                    de00 = ColorSystem::Delta::E00(LAB1, LAB2);
        REQUIRE(de00 == Approx(2.0424596801565738f).epsilon(EPS));
    }
    SECTION("delta2")
    {
        const float EPS  = 1e-6f;
        const float de00 = ColorSystem::Delta::E00(
            ColorSystem::Tristimulus(50.f, 3.1571f, -77.2803f),
            ColorSystem::Tristimulus(50.f, 0.0000f, -82.7485f));
        REQUIRE(de00 == Approx(2.8615f).epsilon(EPS));
    }
}

TEST_CASE("Matrices")
{
    SECTION("gamut_convert")
    {
        const float                 EPS = 1e-6f;
        const ColorSystem::Matrix3 &ACES2065_to_ACEScg(ColorSystem::GamutConvert(ColorSystem::ACES2065, ColorSystem::ACEScg));
        const ColorSystem::Matrix3  expected{
            1.4514393161f, -0.2365107469f, -0.2149285693f,
            -0.0765537734f, 1.1762296998f, -0.0996759264f,
            0.0083161484f, -0.0060324498f, 0.9977163014f};
        REQUIRE_THAT(ACES2065_to_ACEScg, IsApproxEquals(expected, EPS));
    }

    SECTION("bradford")
    {
        //http://w3.kcua.ac.jp/~fujiwara/infosci/colorspace/bradford.html
        const float                 EPS = 1e-2f; // because some constants are differ from reference.
        const ColorSystem::Matrix3 &adopt_D65_to_D50(ColorSystem::Bradford(ColorSystem::Illuminant_D65, ColorSystem::Illuminant_D50));
        const ColorSystem::Matrix3  expected{
            1.047886f, 0.022919f, -0.050216f,
            0.029582f, 0.990484f, -0.017079f,
            -0.009252f, 0.015073f, 0.751678f};
        REQUIRE_THAT(adopt_D65_to_D50, IsApproxEquals(expected, EPS));
    }
}

TEST_CASE("Spectrum")
{
    SECTION("blackbody")
    {
        const double planck6000 = ColorSystem::Spectrum::planck(6000., 380 * 1e-9 );
        REQUIRE(planck6000 == Approx(27366).epsilon(1e2));
    }
    SECTION("toXYZ")
    {
        const ColorSystem::Spectrum &  D65(ColorSystem::Spectrum::blackbody(6504.f));
        const ColorSystem::Spectrum &  E(ColorSystem::Spectrum::E(1.f));
        const ColorSystem::Tristimulus D65_Yxy = ColorSystem::CIEXYZ1931(D65).toYxy();
        const ColorSystem::Tristimulus E_Yxy   = ColorSystem::CIEXYZ1931(E).toYxy();
        REQUIRE(D65_Yxy[1] == Approx(ColorSystem::Illuminant_D65.toYxy()[1]).epsilon(1e-2f)); // precision problem...
        REQUIRE(D65_Yxy[2] == Approx(ColorSystem::Illuminant_D65.toYxy()[2]).epsilon(1e-2f));
        REQUIRE(E_Yxy[1] == Approx(ColorSystem::Illuminant_E.toYxy()[1]).epsilon(1e-5f));
        REQUIRE(E_Yxy[2] == Approx(ColorSystem::Illuminant_E.toYxy()[2]).epsilon(1e-5f));
    }
}