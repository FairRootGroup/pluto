// Author: R. Holzmann & M.A. Kagarlis
// Written: 16.06.00
// Revised: 05.09.05
// Revised: 28.03.06
// Revised: 20.02.07
// Thermal Source Class Header

#ifndef _PThermal_H_
#define _PThermal_H_

#define THERMAL_WIDTH 0.001

#include "TObject.h"
#include "PDynamicData.h"

class PThermal : public TObject {

 public:

  static double dNdE(double *, double *);
  // thermal source + blast

  static double dNdE1(double *, double *);
  // thermal source + blast (E*sqrt(p))

  static double dNdE2(double *, double *);
  // thermal source + blast (E*p**3)

  static double dNdE3(double *, double *);
  // thermal source (p**n)

  static double d2NdEdM(double *, double *);
  // thermal source + blast

  static double d2NdEdM1(double *, double *);
  // thermal source + blast (E*sqrt(p))

  static double d2NdEdM2(double *, double *);
  // thermal source + blast (E*p**3)

  static double d2NdEdM3(double *, double *);
  // thermal source (p**n)

  static double d2NdEdTheta(double *, double *);
  // thermal source (p**n) with T=T(theta_cm) and n = n(theta_cm)

  static double dNdTheta(double*, double *);
  // polar angular distribution


  static double dNdy(double*, double *);
  //  rapidity distribution

  static double dNdMt(double*, double *);
  //  transverse-mass distribution

  static double dNdPt(double*, double *);
  //  transverse-momentum distribution

  static double IntThermal(double, double, int);
  //  normalization integral of exp(-E/T)*E*p**n  where p=sqrt(E**2-m**2)


  static Double_t mtScaleFactor(Int_t id, const Float_t T);
  // Compute the yield enhancement factor due to thermal weighting of a
  // resonance with id at temperature T.

  ClassDef(PThermal,1) //Pluto Thermal Source Class
};
#endif // _PThermal_H_

