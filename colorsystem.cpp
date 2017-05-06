
// "colorsystem"
// copyright 2017 (c) Hajime UCHIMURA / @nikq
// all rights reserved

#include <stdio.h>

#include "colorsystem.hpp"

//sample usage.
int main(int argc, char *argv[])
{
    ColorSystem::Gamut adobeRGB(0.64, 0.33, 0.21, 0.71, 0.15, 0.06, 0.3127, 0.3290);
    ColorSystem::Tristimulus white(1, 1, 1);
    printf("adobe.toXYZ\n");
    dump(adobeRGB.toXYZ());
    dump(adobeRGB.toXYZ(white).vec3());
    printf("adobe.fromXYZ\n");
    dump(adobeRGB.fromXYZ());
    dump(adobeRGB.fromXYZ(white).vec3());
    return 0;
}
