
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
#include "common.hpp"
#include "colorsystem.hpp"

namespace {
    const float epsilon = 0.000001 ;
}

TEST_CASE ("toXYZ", "[XYZ]") {
    const ColorSystem::Gamut adobeRGB (0.64, 0.33, 0.21, 0.71, 0.15, 0.06, 0.3127, 0.3290);

    SECTION ("adobe.toXYZ") {
        const ColorSystem::Matrix3    expected { 0.576669, 0.185558, 0.188229
                                               , 0.297345, 0.627363, 0.075291
                                               , 0.027031, 0.070689, 0.991337
                                               };
        auto const &xyz = adobeRGB.toXYZ ();
        for (int_fast32_t i = 0 ; i < 9 ; ++i) {
            CAPTURE (i) ;
            REQUIRE (xyz [i] == Approx (expected [i]).epsilon (epsilon)) ;
        }
    }
    SECTION ("adobe.toXYZ (white)") {
        const ColorSystem::Tristimulus white (1, 1, 1);
        auto const &v3 = adobeRGB.toXYZ (white).vec3 ();
        REQUIRE (v3[0] == Approx (0.950456).epsilon (epsilon));
        REQUIRE (v3[1] == Approx (1.000000).epsilon (epsilon));
        REQUIRE (v3[2] == Approx (1.089058).epsilon (epsilon));
        //dump(v3) ;
    }
}

TEST_CASE ("fromXYZ", "[XYZ]") {
    ColorSystem::Gamut       adobeRGB (0.64, 0.33, 0.21, 0.71, 0.15, 0.06, 0.3127, 0.3290);
    SECTION ("adobe.fromXYZ") {
        const ColorSystem::Matrix3  expected {  2.041588, -0.565007, -0.344731
                                             , -0.969244,  1.875968,  0.041555
                                             ,  0.013444, -0.118362,  1.015175
                                             } ;
        auto const &    m = adobeRGB.fromXYZ () ;
        for (int_fast32_t i = 0 ; i < 9 ; ++i) {
            CAPTURE (i) ;
            REQUIRE (m [i] == Approx (expected [i]).epsilon (epsilon)) ;
        }
    }
    SECTION ("adobe.fromXYZ (white)") {
        ColorSystem::Tristimulus white (1, 1, 1);
        auto const &    v = adobeRGB.fromXYZ (white).vec3 () ;
        REQUIRE (v [0] == Approx (1.131850).epsilon (epsilon)) ;
        REQUIRE (v [1] == Approx (0.948279).epsilon (epsilon)) ;
        REQUIRE (v [2] == Approx (0.910257).epsilon (epsilon)) ;
    }
}
