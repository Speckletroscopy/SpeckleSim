#include <iostream>

/**
 * Returns the energy of a photon emitted from hydrogen in eV 
 * \param ni (integer) initial energy level
 * \param nf (integer) final energy level
 * \param kryd (double) Rydberg constant of your preferred flavor (default 13.606 eV)
 */
double rydberg(int ni, int nf, double kryd = 13.606){

    double out;
    
    out = kryd*(1./(nf*nf)-1./(ni*ni));

    return out;
}

/**
 * main function for this script
 */
void hydrogenspec(){
    for(int i = 3; i < 10; i++){
        std::cout << 1.e9*TMath::C()/rydberg(i,2,1.602e-19*13.606/TMath::H()) << "\n";
    }
    double v_ps, i_ps, pow_lamp;
    double eff_lamp, lum_lamp;
    double d_lamp;

    v_ps = 5000.0; // power supply voltage (volts)
    i_ps = 0.01; // power supply current (amperes)
    eff_lamp = 0.25; // lamp efficiency
    d_lamp = 0.2; // lamp distance from aperture in meters

    pow_lamp = i_ps*v_ps; // lamp power (watts)
    lum_lamp = pow_lamp*eff_lamp; // lamp luminosity (watts) 


}
