
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

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

