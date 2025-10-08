#include <stdio.h>

#include "..\include\colorsystem.hpp"

int main(void)
{
    for (int temperature = 1500; temperature < 13000; temperature += 100)
    {
        const ColorSystem::Spectrum    blackbody = ColorSystem::Spectrum::blackbody(temperature);
        const ColorSystem::Tristimulus xyzt      = ColorSystem::CIE31JV.fromSpectrum(blackbody);
        const ColorSystem::Tristimulus xyz       = xyzt / xyzt[1];

        const ColorSystem::Tristimulus rgb2020 = ColorSystem::Rec2020.fromXYZ(xyz);
        const ColorSystem::Tristimulus rgb709  = ColorSystem::Rec709.fromXYZ(xyz);
        printf("temp %d : rgb (Rec2020) %f, %f, %f  rgb(Rec709) %f, %f, %f\n", temperature, rgb2020[0], rgb2020[1],
            rgb2020[2], rgb709[0], rgb709[1], rgb709[2]);
    }
    return 0;
}