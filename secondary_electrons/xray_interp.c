#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "dnumrecipes.h"
#include "elec_interp.h"

// Interpolates over power-law spectrum to get net energy deposition fractions for secondary 
// electrons in the IGM.  Uses elec_interp plus the tables.  Corresponds to eq. 16, Fig. 6 in 
// Furlanetto & Stoever

// IGM ionized fraction
float xHII_call = 8.0e-3;
// Power law index of spectrum
float alpha = 1.5;

// Analytic Function for initial Photoionization
//Parameters for Function
//Photoionization cross sections, E is in eV. Returns units of cm^-2.
//Species choices are: 0 = HI, 1 = HeI, 2 = HeII.
//Function taken from Verner et al. 1996, ApJ, 465, 487
const double E0[3] = {0.4298, 13.61, 1.72};
//const double sig0[3] = {5.475e2, 9.492, 1.369e2};
const double sig0[3] = {5.475e-14, 9.492e-16, 1.369e-14};
const double P[3] = {2.963, 3.188, 2.963};
const double ya[3] = {32.88, 1.469, 32.88};
const double yw[3] = {0.0, 2.039, 0.0};
const double y0_pi[3] = {0.0, 0.4434, 0.0};
const double y1_pi[3] = {0.0, 2.136, 0.0};
const double Eion[3] = {13.58, 24.586, 54.398}; /* in eV */
const double nuIon[3] = {3.28e15, 5.93e15, 1.31e16};

double photoion_CrossSection(double E, int species);
void fInt(double E, double y[], double deriv[]);

int main( int argc, char *argv[])
{
  double *ans;
  double step;
  int goodSteps,badSteps;
  double Emin(13.61),Emax;
  int i;
  double INT_ACCURACY(1.0e-4);

  double fheat,nlya,nHI;

  while ((argc > 1) && (argv[1][0] == '-')) {
    switch (argv[1][1]) {
    case 'x':
      xHII_call = atof(&argv[1][2]);
      break;
    case 'a':
      alpha = atof(&argv[1][2]);
      break;
    }
    ++argv;
    --argc;
  }

  initialize_interp_arrays();

  Emax = 1.0e6;

  ans = dvector(1,4);
  printf("Emin\tEmax\tfheat\tflya\tfHI-ion\n");
  for (Emin=13.61;Emin<1.0e4;Emin*=1.1) {
    
    for (i=1;i<=4;i++)
      ans[i] = 0.0;

    step = Emin/1.0e4;

    if (Emin < 23.579) {
      odeint(ans,4,Emin,23.579,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
      odeint(ans,4,23.65,54.3,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
      odeint(ans,4,54.5,Emax,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
    } else if (Emin < 54.3) {
      if (Emin < 23.65)
	Emin = 23.65;
      odeint(ans,4,Emin,54.3,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
      odeint(ans,4,54.5,Emax,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
    } else {
      if (Emin < 54.5)
	Emin = 54.5;
      odeint(ans,4,Emin,Emax,INT_ACCURACY,step,0.0,&goodSteps,
	     &badSteps,fInt,bsstep);
    }

    fheat = ans[1]/ans[4];
    nlya = ans[2]/ans[4];
    nHI = ans[3]/ans[4];
    //    printf("%g\t%g\t%g\t%g\t%g\n",Emin,ans[1],ans[2],ans[3],ans[4]);

    printf("%g\t%g\t%g\t%g\t%g\n",Emin,Emax,fheat,nlya,nHI);
  }

  free_dvector(ans,1,4);

  return 0;
}

void fInt(double En, double y[], double deriv[])
{
  double HI_den(0.92),He_den(0.08);
  double ans1;

  // Heating first
  // Photoionization of HI
  deriv[1] = interp_fheat(En-Eion[0],xHII_call);
  //  printf("deriv[1]=%g\tfheat=%g\n",deriv[1],interp_fheat(En-Eion[0],xHII_call));
  deriv[1] *= photoion_CrossSection(En,0);
  //  printf("deriv[1]=%g\n",deriv[1]);
  deriv[1] *= pow(En,-alpha-1.0);
  //  printf("deriv[1]=%g\n",deriv[1]);
  deriv[1] *= (En-Eion[0])*HI_den;  
  //  printf("deriv[1]=%g\n",deriv[1]);
  // Photoionization of HeI
  if (En > Eion[1]) {
    ans1 = interp_fheat(En-Eion[1],xHII_call);
    deriv[1] += (photoion_CrossSection(En,1)*ans1*pow(En,-alpha-1.0)*
		 (En-Eion[1]))*He_den;
  }

  // Now Lyman-alpha photons
  // Photoionization of HI
  ans1 = interp_n_Lya(En-Eion[0],xHII_call);
  deriv[2] = (photoion_CrossSection(En,0)*ans1*pow(En,-alpha-1.0)*
	      10.2)*HI_den;
  // Photoionization of HeI
  if (En > Eion[1]) {
    ans1 = interp_n_Lya(En-Eion[1],xHII_call);
    deriv[2] += (photoion_CrossSection(En,1)*ans1*pow(En,-alpha-1.0)*
		 10.2)*He_den;
  }

  // Now HI ionization
  // Photoionization of HI
  ans1 = (interp_nion_HI(En-Eion[0],xHII_call));
  deriv[3] = (photoion_CrossSection(En,0)*ans1
	      *pow(En,-alpha-1.0)*Eion[0])*HI_den;
  // Photoionization of HeI
  if (En > Eion[1]) {
    ans1 = (interp_nion_HI(En-Eion[1],xHII_call));
    deriv[3] += (photoion_CrossSection(En,1)*ans1
		 *pow(En,-alpha-1.0)*Eion[0])*He_den;
  }

  // Now normalization factor
  deriv[4] = pow(En,-alpha)*photoion_CrossSection(En,0)*HI_den;
  deriv[4] += pow(En,-alpha)*photoion_CrossSection(En,1)*He_den;
  //  printf("%g\t%g\t%g\t%g\t%g\n",En,deriv[1],deriv[2],deriv[3],deriv[4]);
}

double photoion_CrossSection(double E, int species)
{
  double x,y,ans1;

  if ( species == 2 && E < 54.4 ) {
    return 0;
  };
  if (species == 1 && E < 24.586) {
    return 0;
  };
  if (species == 0 && E < 13.6) {
    return 0;
  };
  
  x = E/E0[species] - y0_pi[species];
  y = sqrt(x*x + y1_pi[species]*y1_pi[species]);

  ans1 = sig0[species]*(pow(x-1.0,2.0) + pow(yw[species],2.0));
  ans1 *= pow(y,0.5*P[species]-5.5)/pow(1.0+sqrt(y/ya[species]),P[species]);
  return ans1;
}
