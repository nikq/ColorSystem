
#include "common.hpp"
#include "TestUtilities.hpp"

#include <colorsystem.hpp>

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
