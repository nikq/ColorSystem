#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

void createPairs(const float x, const float y, const float z, const float l, const float a, const float b,
    ColorSystem::Tristimulus &xyz, ColorSystem::Tristimulus &lab, ColorSystem::Tristimulus &xyz_lab,
    ColorSystem::Tristimulus &lab_xyz)
{
    xyz     = ColorSystem::Tristimulus(x, y, z);
    lab     = ColorSystem::Tristimulus(l, a, b);
    xyz_lab = ColorSystem::XYZ_to_OKLAB(xyz);
    lab_xyz = ColorSystem::OKLAB_to_XYZ(lab);
};

TEST_CASE("OKLAB")
{
    const float EPS = 1e-1f;
    SECTION("test1")
    {
        ColorSystem::Tristimulus xyz, lab, xyz_lab, lab_xyz;
        createPairs(0.950, 1.000, 1.089, 1.000, 0.000, 0.000, xyz, lab, xyz_lab, lab_xyz);
        REQUIRE_THAT(xyz, IsApproxEquals(lab_xyz, EPS));
        REQUIRE_THAT(lab, IsApproxEquals(xyz_lab, EPS));
    }
    SECTION("test2")
    {
        ColorSystem::Tristimulus xyz, lab, xyz_lab, lab_xyz;
        createPairs(1.000, 0.000, 0.000, 0.450, 1.236, -0.019, xyz, lab, xyz_lab, lab_xyz);
        REQUIRE_THAT(xyz, IsApproxEquals(lab_xyz, EPS));
        REQUIRE_THAT(lab, IsApproxEquals(xyz_lab, EPS));
    }
    SECTION("test3")
    {
        ColorSystem::Tristimulus xyz, lab, xyz_lab, lab_xyz;
        createPairs(0.000, 1.000, 0.000, 0.922, -0.671, 0.263, xyz, lab, xyz_lab, lab_xyz);
        REQUIRE_THAT(xyz, IsApproxEquals(lab_xyz, EPS));
        REQUIRE_THAT(lab, IsApproxEquals(xyz_lab, EPS));
    }
    SECTION("test4")
    {
        ColorSystem::Tristimulus xyz, lab, xyz_lab, lab_xyz;
        createPairs(0.000, 0.000, 1.000, 0.153, -1.415, -0.449, xyz, lab, xyz_lab, lab_xyz);
        REQUIRE_THAT(xyz, IsApproxEquals(lab_xyz, EPS));
        REQUIRE_THAT(lab, IsApproxEquals(xyz_lab, EPS));
    }
}