#ifndef CELESTIALMECHANICS_H

 
#define G 6.6743e-11            // The universal gravitational constant (N.m^2/kg^2)

#define AU 1.496e11             // Astronomical unit (m)

#define SunRadius 6.957e8       // radius of the Sun (m)
#define MercuryRadius 2.44e6    // radius of Mercury (m)
#define VenusRadius 6.052e6     // radius of Venus (m)
#define EarthRadius 6.371e6     // radius of Earth (m)
#define MoonRadius 1.737e6      // radius of the Moon (m)
#define MarsRadius 3.389e6      // radius of Mars (m)
#define JupiterRadius 6.9911e7  // radius of Jupiter (m)
#define SaturnRadius 5.8232e7   // radius of Saturn (m)
#define UranusRadius 2.5362e7   // radius of Uranus (m)
#define NeptuneRadius 2.4622e7  // radius of Neptune (m)
#define PlutoRadius 1.1883e6    // radius of Pluto (m)

#define SunMass 1.989e30        // Mass of the sun (kg)
#define MercuryMass 3.3011e23   // Mass of mercury (kg)
#define VenusMass 4.8675e24     // Mass of Venus (kg)
#define EarthMass 5.9722e24     // Mass of the earth (kg)
#define MoonMass 7.342e22       // Mass of the Moon (kg)
#define MarsMass 6.4171e23      // Mass of Mars (kg)
#define JupiterMass 1.8982e27   // Mass of Jupiter (kg)
#define SaturnMass 5.6834e26    // Mass of Saturn (kg)
#define UranusMass 8.6810e25    // Mass of Uranus (kg)
#define NeptuneMass 1.0243e26   // Mass of Neptune (kg)
#define PultoMass 1.303e22      // Mass of Pluto (kg)

#define SUN_SMA 0 // distance from Sun to Mercury (AU)
#define MERCURY_SMA 0.38709893 // distance from Sun to Mercury (AU)
#define VENUS_SMA 0.72333199 // distance from Sun to Venus (AU)
#define EARTH_SMA 1 // distance from Sun to Earth (AU)
#define MARS_SMA 1.52366231 // distance from Sun to Mars (AU)
#define JUPITER_SMA 5.20336301 // distance from Sun to Jupiter (AU)
#define SATURN_SMA 9.53707032 // distance from Sun to Saturn (AU)
#define URANUS_SMA 19.19126393 // distance from Sun to Uranus (AU)
#define NEPTUNE_SMA 30.06896348 // distance from Sun to Neptune (AU)
#define PLUTO_SMA 39.48168677 // distance from Sun to Pluto (dwarf planet) (AU)

#endif /* CELESTIALMECHANICS_H */

// Converts Degrees to Radians
double Deg2Rad(double Deg){
    return Deg * M_PI / 180;
}
// Calculates the reduced mass µ
double ReducedMass(double Object1Mass, double Object2Mass){
    return (Object1Mass * Object2Mass) / (Object1Mass + Object2Mass);
}

// Returns the gravitational force of attraction in (N) between 2 objects given their masses in (kg) and the distance between them in (m).
double GravitationalForce(double Object1Mass, double Object2Mass, double Distance){
    return G * Object1Mass * Object2Mass / pow(Distance, 2.0);
}

// Newton's second law, returns acceleration (m.s^-2) given the force applied to an object (N) and its mass (Kg)
double NewtonsSecondLaw(double Mass, double Force){
    return Force / Mass;
}

// Returns escape velocity from the surface of a primary body in (m/s) given its mass (kg) and raduis (m)
double EscapeVelocity(double PlanetMass, double PlanetRaduis){
	return sqrt(2 * G * PlanetMass / PlanetRaduis);
}

// Returns the period T (seconds) of a secondary object of mass Object2Mass (kg) orbiting a primary object of mass Object1Mass (kg) in an elliptical orbit given; Object1Mass, Object2Mass and the semi-major-axis
double KeplersThirdLawP(double SemiMajorAxis, double Object1Mass, double Object2Mass){
    return 2 * M_PI * sqrt(pow(SemiMajorAxis,3) / (G*(Object1Mass+Object2Mass)));
}

// Returns the total mass of a binary orbital system of given the semi manjor axis (m) of the orbit and its period (s)
double KeplersThirdLawM(double SemiMajorAxis, double Period){
    return ( 4 * pow(M_PI, 2) / G ) * ( pow( SemiMajorAxis, 3) / pow(Period, 2));
}

// Returns an object's momentum (kg.m.s^-1) given its mass(kg) and velocity (m/s)
double Momentum(double Mass, double Velocity){
    return Mass * Velocity;
}

// Returns the centripetal acceleration of an object moving in a perfect circlar orbit
double CentripetalCircularAcceleration(double OrbitalVelocity, double Distance){
    return pow(OrbitalVelocity, 2 ) / Distance;
} 

// Returns the orbital angular momentum (L) of a secondary object orbiting a primary object given their masses (Object2Mass and Object1Mass respectively), the semi major axis (m), and the orbit's eccentricity (e) 
double OrbitalAngularMomentum(double Object1Mass, double Object2Mass, double SemiMajorAxis, double Eccentricity){
    return ( ReducedMass(Object1Mass, Object2Mass) * sqrt(G * (Object1Mass + Object2Mass) * SemiMajorAxis * (1 - pow(Eccentricity,2))));
}

// Returns the distance (meters) between a secondary object orbiting a primary object in an elleptical orbit (0<=e<1) given the semi-major-axis, the eccentricity of the orbit (e) and the polar coordinates of the secondary object as Theta (radians)
double DistanceFrom(double SemiMajorAxis, double Eccentricity, double Theta){
    /*
        A longer (less efficient) alternative formula 
        return ( pow(OrbitalAngularMomentum(Object1Mass, Object2Mass, SemiMajorAxis, Eccentricity),2) / pow(ReducedMass(Object1Mass,Object2Mass),2) ) / ( G * (Object1Mass+Object2Mass) * (1 + ( Eccentricity * cos(Theta) ) ) );
    */
    return (SemiMajorAxis * ( 1 - pow(Eccentricity, 2)) / (1 + (Eccentricity * cos(Theta))));
}

// Returns the orbital velocity of a secondary object (m/s) given its mass Object2Mass (kg), semi-major-axis (m), the eccentricity of its orbit (e), its distance (m) and the mass of the primary object (Object1Mass)
double OrbitalVelocity(double Object1Mass, double Object2Mass, double SemiMajorAxis, double Eccentricity, double Distance){
    return sqrt ( G * ( Object1Mass + Object2Mass ) * SemiMajorAxis * (1 - pow(Eccentricity, 2)) )  / Distance;
    
}

// Returns the total orbital energy of a binary system given the mass of the primary and secondary object and the semi-major-axis
double TotalOrbitalEnergy(double Object1Mass, double Object2Mass, double SemiMajorAxis){
    return ( -1 * G * Object1Mass * Object2Mass ) / (2 * SemiMajorAxis);
}


