
#include "common.hpp"

#include "TestUtilities.hpp"

#include <colorsystem.hpp>

TEST_CASE("Macbeth")
{
    SECTION("ChartToXYZ")
    {
        for( const ColorSystem::Spectrum& s : ColorSystem::Macbeth::Patch )
        {
            const ColorSystem::Tristimulus& xyz = ColorSystem::CIE1931.fromSpectrum( s * ColorSystem::CIE_D65 );
            const ColorSystem::Tristimulus& Yxy = xyz.toYxy();
            printf("%f,%f,%f\n",Yxy[0],Yxy[1],Yxy[2]); // TODO: check with D50 values.
        }

        std::vector<ColorSystem::Tristimulus> v = ColorSystem::Macbeth::reference( ColorSystem::CIE_D65, ColorSystem::CIE1931 );
        for( const auto &p:v){
            const ColorSystem::Tristimulus& Yxy = p.toYxy();
            printf("%f,%f,%f / %f,%f,%f\n",p[0],p[1],p[2],Yxy[0],Yxy[1],Yxy[2]);
        }
    }
}
