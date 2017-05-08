#include "common.hpp"
#include "TestUtilities.hpp"

#include <colorsystem.hpp>

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

