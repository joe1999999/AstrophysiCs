#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CelestialMechanics.h"



int main(){
    printf("Neptune is at a distance of %.5f AU from the sun at perihelion.\n", DistanceFrom(NEPTUNE_SMA,0.01671123, Deg2Rad(0)) );
    printf("Neptune is at a distance of %.5f AU from the sun at ephelion.\n", DistanceFrom(NEPTUNE_SMA,0.008678, Deg2Rad(180)) );
    printf("Neptune takes %.5f years to orbit the sun\n", KeplersThirdLawP(NEPTUNE_SMA*AU, SunMass, NeptuneMass) / 3600/24/365);

    return 0;
}