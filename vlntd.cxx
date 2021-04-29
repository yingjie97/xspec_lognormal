#include <stlToCArrays.h>
#include <functionMap.h>
#include <FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <xsTypes.h>
#include <cmath>
#include "funcWrappers.h"

// function from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType,
                        const IntegerArray& Zarray, const RealArray& abun,
                        const Real dens, const Real z, const RealArray& Tarr,
                        const RealArray& DEMarr, const int ifl, const bool qtherm,
                        const Real velocity, RealArray& fluxArray,
                        RealArray& fluxErrArray);

// XSPEC model subroutine to calculate thermal plasma with log-normal temperature distribution:
//    Q(T) = 1/((2*pi)**0.5*sigma)*exp(-(ln(T)-ln(Tmean))**2/(2*sigma**2))
// Parameters:
//    param(1) = Temperature mean
//    param(2) = sigma
//    param(3) = nt (number of grid bins)
//    param(4) = nH (cm^-3)  Fixed at 1 for most applications
//    param(5) = He abundance
//    param(6) = C   "
//    param(7) = N   "
//    param(8) = O   "
//    param(9) = Ne  "
//    param(10) = Na  "
//    param(11)= Mg  "
//    param(12)= Al  "
//    param(13)= Si  "
//    param(14)= S   "
//    param(15)= Ar  "
//    param(16)= Ca  "
//    param(17)= Fe  "
//    param(18)= Ni  " 
//    param(19)= redshift
//    param(20) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=interpolate APEC model)

extern "C" void vlntd(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString)
{

   using namespace XSutility;

   const Real Tmean = params[0]; //mean temperature  
   //define x to be the natural log of mean temperature over sigma
   const Real sigma = params[1]; //standard deviation of x
   const Real xmean = log(Tmean) / sigma; //mean value of x
   const int nt = params[2]; //the number of grid bins

   // *******************************************************************
   // run nx steps from -nSig to +nSig
   int nx = (int) nt;
   int nSig = 3.0;
    
   // set up arrays of temperature and DEM values
   RealArray Tarray(nx);
   RealArray demarray(nx);

   Real TMin = Tmean*exp(-nSig*sigma);
   if ( TMin <= 0.0 ) TMin = 0.001;  // set up the minimum temperature
   const Real xmin = log(TMin) / sigma; // minimum x
   Real dx = 2. * nSig / nx;  // step size

   for (int i=0; i<nx; i++) {
     Real x = xmin + dx * (i + 0.5); // x for each bin
     Tarray[i] = exp(x * sigma); // temperature
     demarray[i] = (erf((x+dx-xmean)/sqrt(2.))-erf((x-xmean)/sqrt(2.)))/2.;
     // the Gaussian weight
   }
   // end of set up arrays of temperature and DEM values
   // *******************************************************************


   // set up all the variables to pass to calcMultiTempPlasma
   int swtch = static_cast<int>(params[19]);
   int plasmaType;
   if ( swtch == 0 ) {
     plasmaType = 3;
   } else if ( swtch == 1 ) {
     plasmaType = 4;
   } else if ( swtch == 2 ) {
     plasmaType = 6;
   } else {
    FunctionUtility::xsWrite("\n VCLUSCOOL: Invalid switch parameter value",2);
    FunctionUtility::xsWrite("            Must be 0, 1, or 2",2);
    return;
  }

   const Real density = params[3];
   const Real redshift = params[18];

   RealArray abun(14);
   for (size_t i=0; i<abun.size(); i++) abun[i] = params[i+4];
   const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
//                        He, C, N，O, Ne, Na, Mg, Al, Si,  S, Ar, Ca, Fe, Ni
   IntegerArray Zarray(14);
   for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

   const bool qtherm = false;
   const Real velocity = 0.0;

   int status=0;
   status = calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, density,
                                redshift, Tarray, demarray, spectrumNumber,
				qtherm, velocity, flux, fluxErr);

   if (status != 0) {
     std::ostringstream msg;
     msg << "vadaf: error status " << status << " returned from calcMultiTempPlasma";
     FunctionUtility::xsWrite(msg.str(),5);
   }
}
