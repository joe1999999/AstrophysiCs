#include "Cosmology.h"

int main() {
    // Cosmological parameters
    double omega_m = 0.3;       // Matter density parameter
    double omega_lambda = 0.7;  // Dark energy density parameter
    double omega_k = 0.0;       // Curvature density parameter (flat universe)
    double z = 1.0;             // Redshift

    // Calculate Hubble parameter at redshift z
    double H_z = hubble_parameter(z, omega_m, omega_lambda);
    printf("Hubble parameter at z = %.2f: %.2f km/s/Mpc\n", z, H_z);

    // Calculate critical density at redshift z
    double rho_crit = critical_density(H_z);
    printf("Critical density at z = %.2f: %.2e kg/m^3\n", z, rho_crit);

    // Calculate scale factor at redshift z
    double a = scale_factor(z);
    printf("Scale factor at z = %.2f: %.2f\n", z, a);

    // Calculate Friedmann equation at scale factor a
    double H_a = friedmann_equation(a, omega_m, omega_lambda, omega_k);
    printf("Friedmann equation at a = %.2f: %.2f km/s/Mpc\n", a, H_a);

    return 0;
}