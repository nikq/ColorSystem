
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
    }
}
