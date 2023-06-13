#ifndef GENERALRELATIVITY_H

#define G 6.67430e-11           // gravitational constant (m^3 kg^-1 s^-2)
#define c 299792458             // speed of light in vacuum (m/s)
#define h 6.62607015e-34        // Planck constant (J s)
#define kB 1.380649e-23         // Boltzmann constant (J K^-1)
#define e 1.602176634e-19       // elementary charge (C)
#define mu0 1.25663706212e-6    // vacuum permeability (N A^-2)
#define eps0 8.8541878128e-12   // vacuum permittivity (F m^-1)
#define alpha 7.2973525693e-3   // fine-structure constant
#define hbar h/(2*M_PI)         // reduced Planck constant (J s)
#define lambda 1.0e-52          // the cosmological constant, in units of m^-2
#define Egc 8.0*M_PI*G/pow(c, 4)  // Einstein gravitational constant (N^-1)

#endif /* GENERALRELATIVITY_H */



/* 
    The Einstein field equations describe the relationship between spacetime curvature and the matter-energy content of the universe.
    G_{mu,nu} = 8 * pi * G / c^4 * T_{mu,nu}
    where G_{mu,nu} is the Einstein tensor, T_{mu,nu} is the stress-energy tensor, and G and c are the gravitational constant and speed of light, respectively. 
*/
double G_mu_nu[4][4];  // Einstein tensor

// Calculates the Einstein Tensor (4x4 matrix) given the Stress-energy tensor (4x4 matrix)
double EinsteinTensor(double T_mu_nu[4][4]){
    for (int mu = 0; mu < 4 ; mu++){
        for (int nu = 0 ; nu < 4 ; nu++){
            G_mu_nu[mu][nu] = Egc * T_mu_nu[mu][nu];
        }
    }
}

// Calculates the Christoffel Symbol 
// 'g' is the metric tensor and 'dg' is its derivatives with respect to the spacetime coordinates x, 'dim' is the number of dimensions duuuuh
double ***calculate_christoffel(double ***g, double **dg, int dim) {
    int i, j, k;
    double ***gamma;
    gamma = (double ***) malloc(dim * sizeof(double **));
    for (i = 0; i < dim; i++) {
        gamma[i] = (double **) malloc(dim * sizeof(double *));
        for (j = 0; j < dim; j++) {
            gamma[i][j] = (double *) malloc(dim * sizeof(double));
            for (k = 0; k < dim; k++) {
                gamma[i][j][k] = 0.5 * (
                    g[i][j][k] * dg[i][i] +
                    g[i][k][j] * dg[i][k] +
                    g[k][j][i] * dg[k][i] -
                    g[k][i][j] * dg[k][i]
                );
            }
        }
    }
    return gamma;
}


// The geodesic equation describes the path of a free-falling particle (i.e., one not subject to any external forces) in curved spacetime.
// d^2x^mu / d lambda^2 + Gamma^mu_{nu,rho} * dx^nu / d lambda * dx^rho / d lambda = 0
// where x^mu is the coordinate of the particle in spacetime, lambda is an affine parameter along the particle's path, and Gamma^mu_{nu,rho} is the Christoffel symbol.

double x[4];             // Particle's spacetime coordinates
double dx[4];            // First derivative of particle's coordinates 
double ddx[4];           // Second derivative of particle's coordinates 
double Gamma[4][4][4];   // Christoffel symbol (calculate with function above)


// Calculates the parameters that describe the motion of a particle in curved spacetime
double GeodesicEquation(double x[4], double dx[4], double ddx[4]){
    for (int mu = 0; mu < 4; mu++) {
        ddx[mu] = 0;
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                ddx[mu] -= Gamma[mu][nu][rho] * dx[nu] * dx[rho];
            }
        }
    }
}

// Calculates Time-dialation factor t/t0 given v; the relative velocity (m/s). t0; the time measured by the stationary observer
double time_dilation(double v, double t0) {
    double gamma = 1 / sqrt(1 - pow(v,2)/pow(c,2)); // Lorentz factor
    double t = gamma * t0; // time measured by moving observer
    return t-t0;
}

/* Calculates the relative frequency change due to gravitational redshift 
Given  : 
    f is the frequency of the light emitted from a source in a weak gravitational field
    f_prime is the frequency of the same light observed at a distance from the source in a strong gravitational field
    GFObjectMass is the mass of the object producing the gravitational field
    r is the distance between the object and the point where the light is emitted
    v is the velocity of the emitting object relative to the observer
*/
double GravitationalRedshift(double f, double GFObjectMass, double r){
    double f_prime = f * sqrt( 1 - ( 2 * G * GFObjectMass / (r * c * c) ) );
    return f_prime;
}