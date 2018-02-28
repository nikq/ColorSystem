
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

namespace
{
const float epsilon = 0.000001f;
}

TEST_CASE("toXYZ", "[XYZ]")
{
    SECTION("adobe.toXYZ")
    {
        const ColorSystem::Matrix3 expected{0.576669f, 0.185558f, 0.188229f, 0.297345f, 0.627363f, 0.075291f, 0.027031f, 0.070689f, 0.991337f};
        auto const &               xyz = ColorSystem::AdobeRGB.toXYZ();
        REQUIRE_THAT(xyz, IsApproxEquals(expected, epsilon));
    }
    SECTION("adobe.toXYZ (white)")
    {
        const ColorSystem::Tristimulus white(1, 1, 1);
        auto const &                   v3 = ColorSystem::AdobeRGB.toXYZ(white).vec3();
        REQUIRE_THAT(v3, IsApproxEquals(ColorSystem::Vector3{0.950456f, 1.000000f, 1.089058f}, epsilon));
    }
}

TEST_CASE("fromXYZ", "[XYZ]")
{
    SECTION("adobe.fromXYZ")
    {
        const ColorSystem::Matrix3 expected{2.041588f, -0.565007f, -0.344731f, -0.969244f, 1.875968f, 0.041555f, 0.013444f, -0.118362f, 1.015175f};
        auto const &               m = ColorSystem::AdobeRGB.fromXYZ();
        REQUIRE_THAT(m, IsApproxEquals(expected, epsilon));
    }
    SECTION("adobe.fromXYZ (white)")
    {
        ColorSystem::Tristimulus white(1, 1, 1);
        auto const &             v = ColorSystem::AdobeRGB.fromXYZ(white).vec3();
        REQUIRE_THAT(v, IsApproxEquals(ColorSystem::Vector3{1.131850f, 0.948279f, 0.910257f}, epsilon));
    }
}

TEST_CASE("BT709")
{
    ColorSystem::Tristimulus white(1, 1, 1);

    auto const &monitor = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, ColorSystem::AdobeRGB.fromXYZ(white));
    
    // 1.062991,0.974048,0.954468
    REQUIRE_THAT(monitor, IsApproxEquals(ColorSystem::Tristimulus{1.0f, 0.974048f, 0.954468f}, epsilon));
}

TEST_CASE("toYxy")
{
    SECTION("Convert from X")
    {
        ColorSystem::Tristimulus X(1.f, 0.f, 0.f);
        auto const &             X_Yxy = X.toYxy();
        REQUIRE_THAT(X_Yxy, IsApproxEquals(ColorSystem::Tristimulus{0.0000f, 1.0000f, 0.0000f}, epsilon));
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
        REQUIRE_THAT(Z_Yxy, IsApproxEquals(ColorSystem::Tristimulus{0.0000f, 0.0000f, 0.0000f}, epsilon));
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

TEST_CASE("fromCT")
{
    SECTION("1931 6500K = 0.3127 0.3290, WP approx")
    {
        ColorSystem::Tristimulus W     = ColorSystem::Tristimulus::fromPlanckianLocus(6504.f);
        auto const &             W_Yxy = W.toYxy();
        printf("%f,%f,%f\n", W_Yxy[0], W_Yxy[1], W_Yxy[2]);
        REQUIRE_THAT(W_Yxy, IsApproxEquals(ColorSystem::Tristimulus{1.00f, 0.3127f, 0.3290f}, 1e-2f)); // this approximation is not precise.
    }
    SECTION("1931 6500K = 0.3127 0.3290, our approx")
    {
        ColorSystem::Tristimulus W     = ColorSystem::Tristimulus::fromCT(6504.f);
        auto const &             W_Yxy = W.toYxy();
        printf("%f,%f,%f\n", W_Yxy[0], W_Yxy[1], W_Yxy[2]);
        REQUIRE_THAT(W_Yxy, IsApproxEquals(ColorSystem::Tristimulus{1.00f, 0.3127f, 0.3290f}, 1e-2f)); // this approximation is not precise.
    }
}