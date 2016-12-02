// Below gives grid sizes for the interpolation arrays
#ifndef _x_int_VARIABLES_DEFINED
#define x_int_NXHII  14
#define x_int_NENERGY  258
#define _x_int_VARIABLES_DEFINED
#endif

// Initialization; must be called once to 
void initialize_interp_arrays();

// Primary functions to compute heating fractions and number of Lya photons or ionization produced,
// Note that En is the energy of the electron, so the energy in the initial ionization is NOT
// included in all these.
// All energies are in eV.
// xHII_call is the desired ionized fraction.
float interp_fheat(float En, float xHII_call);
float interp_n_Lya(float En, float xHII_call);
float interp_nion_HI(float En, float xHII_call);
float interp_nion_HeI(float En, float xHII_call);
float interp_nion_HeII(float En, float xHII_call);

int locate_energy_index(float En);
int locate_xHII_index(float xHII_call);


