#include <math.h>

// Constants
#define G 6.67430e-11  // Gravitational constant (m^3 kg^-1 s^-2)
#define c 299792458  // Speed of light (m/s)
#define H0 70.0        // Hubble constant (km/s/Mpc)
#define Mpc_to_m 3.086e22  // Conversion factor: 1 Mpc to meters
#define km_to_m 1000.0     // Conversion factor: 1 km to meters

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
C: Intercept, related to the intrinsic brightness M0 and distance scale.

*/
double ApparentMagnitude(double redshift, double C){
    return C + ( 5 *  log10(redshift));
}



// Function to calculate the Hubble parameter at a given redshift z
double hubble_parameter(double z, double omega_m, double omega_lambda) {
    return H0 * sqrt(omega_m * pow(1 + z, 3) + omega_lambda);
}

// Function to calculate the critical density of the universe
double critical_density(double H) {
    double H_si = H * km_to_m / Mpc_to_m;  // Convert H0 to SI units (s^-1)
    return (3 * pow(H_si, 2)) / (8 * M_PI * G);
}

// Function to calculate the scale factor a(t) from redshift z
double scale_factor(double z) {
    return 1.0 / (1.0 + z);
}

// Function to calculate the Friedmann equation (first equation)
double friedmann_equation(double a, double omega_m, double omega_lambda, double omega_k) {
    return H0 * sqrt(omega_m / pow(a, 3) + omega_k / pow(a, 2) + omega_lambda);
}
