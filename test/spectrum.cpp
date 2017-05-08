
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

TEST_CASE("Spectrum")
{
    SECTION("blackbody")
    {
        const double planck6000 = ColorSystem::Spectrum::planck(6000., 380 * 1e-9);
        REQUIRE(planck6000 == Approx(27366).epsilon(1e2));
    }
    SECTION("toXYZ")
    {
        const ColorSystem::Spectrum &  BB65(ColorSystem::Spectrum::blackbody(6504.f));
        const ColorSystem::Spectrum &  E(ColorSystem::Spectrum::E(1.f));
        const ColorSystem::Tristimulus D65_Yxy  = ColorSystem::CIE1931.fromSpectrum(ColorSystem::CIE_D65).toYxy();
        const ColorSystem::Tristimulus BB65_Yxy = ColorSystem::CIE1931.fromSpectrum(BB65).toYxy();
        const ColorSystem::Tristimulus E_Yxy    = ColorSystem::CIE1931.fromSpectrum(E).toYxy();
        REQUIRE(D65_Yxy[1] == Approx(ColorSystem::Illuminant_D65.toYxy()[1]).epsilon(1e-4f));
        REQUIRE(D65_Yxy[2] == Approx(ColorSystem::Illuminant_D65.toYxy()[2]).epsilon(1e-4f));
        REQUIRE(BB65_Yxy[1] == Approx(ColorSystem::Illuminant_D65.toYxy()[1]).epsilon(1e-2f)); // precision problem...
        REQUIRE(BB65_Yxy[2] == Approx(ColorSystem::Illuminant_D65.toYxy()[2]).epsilon(1e-2f));
        REQUIRE(E_Yxy[1] == Approx(ColorSystem::Illuminant_E.toYxy()[1]).epsilon(1e-5f));
        REQUIRE(E_Yxy[2] == Approx(ColorSystem::Illuminant_E.toYxy()[2]).epsilon(1e-5f));
    }
}
