
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

TEST_CASE("Units")
{
    SECTION("lumens")
    {
        const float lumens555 = ColorSystem::CIE1931.lumensFromMonochromaticFlux( 555.0f, 1.0f );//
        REQUIRE( lumens555 == Approx(683.0f).epsilon(1e-3f) );
    }
    SECTION("candellas")
    {
        const ColorSystem::Tristimulus cd_m2 = ColorSystem::CIE1931.candellasFromMonochromaticRadiance( 555.0f, 1.0f );
        REQUIRE_THAT(cd_m2, IsApproxEquals(ColorSystem::Tristimulus{349.730225f,683.0f,3.927249f}, 1e-3f));
    }
}
