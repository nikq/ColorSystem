
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

TEST_CASE("oetf", "")
{
    SECTION("sRGB")
    {
        const ColorSystem::Tristimulus v0(0.f, 0.f, 0.f);
        const ColorSystem::Tristimulus v1(0.25f, 0.5f, 0.75f);
        const ColorSystem::Tristimulus v2(1.0f, 1.0f, 1.0f);
        const ColorSystem::Tristimulus s0 = ColorSystem::OTF::toScreen(ColorSystem::OTF::SRGB, v0);
        const ColorSystem::Tristimulus s1 = ColorSystem::OTF::toScreen(ColorSystem::OTF::SRGB, v1);
        const ColorSystem::Tristimulus s2 = ColorSystem::OTF::toScreen(ColorSystem::OTF::SRGB, v2);
        REQUIRE_THAT(s0, IsApproxEquals(v0, epsilon));
        REQUIRE_THAT(s2, IsApproxEquals(v2, epsilon));
        const ColorSystem::Tristimulus c0 = ColorSystem::OTF::toScene(ColorSystem::OTF::SRGB, s0);
        const ColorSystem::Tristimulus c1 = ColorSystem::OTF::toScene(ColorSystem::OTF::SRGB, s1);
        const ColorSystem::Tristimulus c2 = ColorSystem::OTF::toScene(ColorSystem::OTF::SRGB, s2);
        REQUIRE_THAT(v0, IsApproxEquals(c0, epsilon));
        REQUIRE_THAT(v1, IsApproxEquals(c1, epsilon));
        REQUIRE_THAT(v2, IsApproxEquals(c2, epsilon));
    }
    SECTION("BT709")
    {
        const ColorSystem::Tristimulus v0(0.f, 0.f, 0.f);
        const ColorSystem::Tristimulus v1(0.25f, 0.5f, 0.75f);
        const ColorSystem::Tristimulus v2(1.0f, 1.0f, 1.0f);
        const ColorSystem::Tristimulus s0 = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, v0);
        const ColorSystem::Tristimulus s1 = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, v1);
        const ColorSystem::Tristimulus s2 = ColorSystem::OTF::toScreen(ColorSystem::OTF::BT709, v2);
        const ColorSystem::Tristimulus c0 = ColorSystem::OTF::toScene(ColorSystem::OTF::BT709, s0);
        const ColorSystem::Tristimulus c1 = ColorSystem::OTF::toScene(ColorSystem::OTF::BT709, s1);
        const ColorSystem::Tristimulus c2 = ColorSystem::OTF::toScene(ColorSystem::OTF::BT709, s2);
        REQUIRE_THAT(v0, IsApproxEquals(c0, epsilon));
        REQUIRE_THAT(v1, IsApproxEquals(c1, epsilon));
        REQUIRE_THAT(v2, IsApproxEquals(c2, epsilon));
    }
    SECTION("HLG")
    {
        const ColorSystem::Tristimulus v0(0.f, 0.f, 0.f);
        const ColorSystem::Tristimulus v1(0.25f, 0.5f, 0.75f);
        const ColorSystem::Tristimulus v2(1.0f, 2.2f, 4.4f);
        const ColorSystem::Tristimulus s0 = ColorSystem::OTF::toScreen(ColorSystem::OTF::HLG, v0);
        const ColorSystem::Tristimulus s1 = ColorSystem::OTF::toScreen(ColorSystem::OTF::HLG, v1);
        const ColorSystem::Tristimulus s2 = ColorSystem::OTF::toScreen(ColorSystem::OTF::HLG, v2);
        REQUIRE_THAT(s0, IsApproxEquals(v0, epsilon));
        const ColorSystem::Tristimulus c0 = ColorSystem::OTF::toScene(ColorSystem::OTF::HLG, s0);
        const ColorSystem::Tristimulus c1 = ColorSystem::OTF::toScene(ColorSystem::OTF::HLG, s1);
        const ColorSystem::Tristimulus c2 = ColorSystem::OTF::toScene(ColorSystem::OTF::HLG, s2);
        REQUIRE_THAT(v0, IsApproxEquals(c0, epsilon));
        REQUIRE_THAT(v1, IsApproxEquals(c1, epsilon));
        REQUIRE_THAT(v2, IsApproxEquals(c2, epsilon));
    }
}