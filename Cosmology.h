#include <math.h>
#define c 299792458             // speed of light in vacuum (m/s)
#define HUBBLE_CONSTANT 70.0    // Hubble's constant in km/s/Mpc


/* Hubble's Law 
takes H0 : the hubble constant
and r : the distance to the galaxy (in megaparsecs, Mpc)
returns the Recessional velocity of a galaxy (in km/s)

*/
double HubblesLaw(double H, double r){
    return H * r;
}

/* 
For a set of observed “standard candles”, i.e.
galaxies whose absolute magnitudes are close
to some mean M0 a linear relationship between the apparent mag-
nitude m and the logarithm of the redshift, lg z.
This is because a galaxy at distance r has an ap-
parent magnitude m = M0 + 5 lg(r/10 pc), and hence Hubble’s law yields
5 lg z + C.
C: Intercept, related to the intrinsic brightness and distance scale.

*/
double ApparentMagnitude(double redshift, double C){
    return C + ( 5 *  log10(redshift));
}