
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

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

