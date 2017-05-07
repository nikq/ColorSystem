
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved
#include "common.hpp"
#include "colorsystem.hpp"

namespace {
    const float epsilon = 0.000001f ;
}

TEST_CASE ("toXYZ", "[XYZ]") {
    const ColorSystem::Gamut adobeRGB (0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);

    SECTION ("adobe.toXYZ") {
        const ColorSystem::Matrix3    expected { 0.576669f, 0.185558f, 0.188229f
                                               , 0.297345f, 0.627363f, 0.075291f
                                               , 0.027031f, 0.070689f, 0.991337f
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
    ColorSystem::Gamut       adobeRGB (0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    SECTION ("adobe.fromXYZ") {
        const ColorSystem::Matrix3  expected {  2.041588f, -0.565007f, -0.344731f
                                             , -0.969244f,  1.875968f,  0.041555f
                                             ,  0.013444f, -0.118362f,  1.015175f
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

TEST_CASE ("BT709") {
    ColorSystem::Gamut       adobeRGB (0.64f, 0.33f, 0.21f, 0.71f, 0.15f, 0.06f, 0.3127f, 0.3290f);
    ColorSystem::Tristimulus white (1, 1, 1);
    auto const &monitor = ColorSystem::OTF::toScreen( ColorSystem::OTF::BT709, adobeRGB.fromXYZ( white )*100.f) ;
    // 1.062991,0.974048,0.954468
    REQUIRE (monitor [0] == Approx (1.062991).epsilon (epsilon)) ;
    REQUIRE (monitor [1] == Approx (0.974048).epsilon (epsilon)) ;
    REQUIRE (monitor [2] == Approx (0.954468).epsilon (epsilon)) ;
}
