///////////////////////////////////////////////////////////////////////////////
//  Pluto Thermal Source Class
//
//  This is a static class containing thermal-source functions
//
//                             Author:  R. Holzmann & M.A. Kagarlis
//                             Written: 16.06.00
//                             Revised: 02.02.07 R.H.
//                             Revised: 20.02.07 R.H.
//
// Ref 1:   Siemens and Rasmussen, PRL 42 (1979) 880 (Note the typos in Eq. 1,
//          corrected e.g. in Reisdorf at al., NP A612 (1997) 512)
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <TMath.h>
#include "PThermal.h"
#include "PData.h"
ClassImp(PThermal)



double PThermal:: dNdE(double* x, double* par) {
  // thermal source + blast (Ref 1)

  double E = x[0];
  double M = par[0];
  double T1 = par[1];
  double T2 = par[2];
  double f = par[3];
  double beta = par[4];
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  if (beta < 0.01 || beta > 0.99) {   // simple thermal source
    val = f*exp(-E/T1);
    if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  } else {                            // thermal + blast
    
    double gamma = 1./sqrt(1.-beta*beta);
    double alpha = beta*gamma*p/T1;
    val = f*exp(-gamma*E/T1)
      * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
    if (f < 1. && T2 > 0.) {
      alpha = beta*gamma*p/T2;
      val += (1.-f)*exp(-gamma*E/T2)
	* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
    }
  }

  // Another definition of f is given if dNdE is normalized to 1:
  // val -> val/[M*M*T*K2(M/T)]  (See e.g. NP A612 (1997) 512) 

  return 100000.*p*E*val;
}

double PThermal:: dNdE1(double* x, double* par) {
  // thermal source + blast (Ref 1) with E*sqrt(p) weight

  double E = x[0];
  double M = par[0];
  double T1 = par[1];
  double T2 = par[2];
  double f = par[3];
  double beta = par[4];
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  if (beta < 0.01 || beta > 0.99) {   // simple thermal source
    val = f*exp(-E/T1);
    if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  } else {                            // thermal + blast
    
    double gamma = 1./sqrt(1.-beta*beta);
    double alpha = beta*gamma*p/T1;
    val = f*exp(-gamma*E/T1)
      * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
    if (f < 1. && T2 > 0.) {
      alpha = beta*gamma*p/T2;
      val += (1.-f)*exp(-gamma*E/T2)
	* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
    }
  }
  
  return 100000.*sqrt(p)*E*val;
}

double PThermal:: dNdE2(double* x, double* par) {
  // thermal source + blast (Ref 1) with E*p**3 weight

  double E = x[0];
  double M = par[0];
  double T1 = par[1];
  double T2 = par[2];
  double f = par[3];
  double beta = par[4];
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  if (beta < 0.01 || beta > 0.99) {   // simple thermal source
    val = f*exp(-E/T1);
    if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  } else {                            // thermal + blast
    
    double gamma = 1./sqrt(1.-beta*beta);
    double alpha = beta*gamma*p/T1;
    val = f*exp(-gamma*E/T1)
      * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
    if (f < 1. && T2 > 0.) {
      alpha = beta*gamma*p/T2;
      val += (1.-f)*exp(-gamma*E/T2)
	* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
    }
  }
  
  return 100000.*p*p*p*E*val;
}

double PThermal:: dNdE3(double* x, double* par) {
  // thermal source with E*p**n weight

  double E = x[0];
  double M = par[0];
  double T1 = par[1];
  double T2 = par[2];
  double f = par[3];
  double n = par[4];
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  val = f*exp(-E/T1);
  if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  return 1000000.*pow(p,n)*val;
}

double PThermal:: d2NdEdM(double* x, double* par) {
  // thermal source + blast (Ref 1)  x  Breit-Wigner

    //cout << "a";

    double T1 = par[0];
    double T2 = par[1];
    double f = par[2];
    double beta = par[3];
    int id = int(par[4]+0.001);
    int idx = int(par[5]+0.001); 

    double val;
    double E = x[0];
    double M = x[1];

    if (int(par[6]+0.001) == 1) {
	//option=rotate
	//cout  << "rot" << endl;
	E = (x[0]*0.5 + 2.*x[1])*0.5;
	M = (2.*x[1] - x[0]*0.5)*0.5;
    }

    double p2 = E*E - M*M;

    if (p2 <= 0.) return 0.;
    
    double p = sqrt(p2);
  
    if (beta < 0.01 || beta > 0.99) {   // simple thermal source
	val = f*exp(-E/T1);
	if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
	
    } else {                            // thermal + blast
	
	double gamma = 1./sqrt(1.-beta*beta);
	double alpha = beta*gamma*p/T1;
	val = f*exp(-gamma*E/T1)
	    * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
	
	if (f < 1. && T2 > 0.) {
	    alpha = beta*gamma*p/T2;
	    val += (1.-f)*exp(-gamma*E/T2)
		* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
	}
    }
    
    if (makeStaticData()->GetParticleTotalWidth(id) > THERMAL_WIDTH) 
	val *= makeDynamicData()->GetParticleTotalWeight(M,id,idx);
  
    // Another definition of f is given if d2NdEdM is normalized to 1
    
    return 100000.*p*E*val;
}



double PThermal:: d2NdEdM1(double* x, double* par) {
  // thermal source + blast (Ref 1) with E*sqrt(p) weight  x  Breit-Wigner

  double E = x[0];
  double M = x[1];
  double T1 = par[0];
  double T2 = par[1];
  double f = par[2];
  double beta = par[3];
  int id = int(par[4]+0.001);
  int idx = int(par[5]+0.001);
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  if (beta < 0.01 || beta > 0.99) {   // simple thermal source
    val = f*exp(-E/T1);
    if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  } else {                            // thermal + blast
    
    double gamma = 1./sqrt(1.-beta*beta);
    double alpha = beta*gamma*p/T1;
    val = f*exp(-gamma*E/T1)
      * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
    if (f < 1. && T2 > 0.) {
      alpha = beta*gamma*p/T2;
      val += (1.-f)*exp(-gamma*E/T2)
	* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
    }
  }
  
  if (makeStaticData()->GetParticleTotalWidth(id) > THERMAL_WIDTH) 
      val *= makeDynamicData()->GetParticleTotalWeight(M,id,idx);
 // else
//       val *= makeDynamicData()->GetParticleTotalWeight(M,id);
  return 100000.*sqrt(p)*E*val;
}

double PThermal:: d2NdEdM2(double* x, double* par) {
  // thermal source + blast (Ref 1) with E*p**3 weight  x  Breit-Wigner

  double E = x[0];
  double M = x[1];
  double T1 = par[0];
  double T2 = par[1];
  double f = par[2];
  double beta = par[3];
  int id = int(par[4]+0.001);
  int idx = int(par[5]+0.001);
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  if (beta < 0.01 || beta > 0.99) {   // simple thermal source
    val = f*exp(-E/T1);
    if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  } else {                            // thermal + blast
    
    double gamma = 1./sqrt(1.-beta*beta);
    double alpha = beta*gamma*p/T1;
    val = f*exp(-gamma*E/T1)
      * ((gamma+T1/E)*sinh(alpha)/alpha -T1/E*cosh(alpha));
    if (f < 1. && T2 > 0.) {
      alpha = beta*gamma*p/T2;
      val += (1.-f)*exp(-gamma*E/T2)
	* ((gamma+T2/E)*sinh(alpha)/alpha -T2/E*cosh(alpha));
    }
  }
  
  if (makeStaticData()->GetParticleTotalWidth(id) > THERMAL_WIDTH) 
      val *= makeDynamicData()->GetParticleTotalWeight(M,id,idx);
//  else
//       val *= makeDynamicData()->GetParticleTotalWeight(M,id);
  return 100000.*p*p*p*E*val;
}

double PThermal:: d2NdEdM3(double* x, double* par) {
  // thermal source with E*p**n weight  x  Breit-Wigner

  double E = x[0];
  double M = x[1];
  double T1 = par[0];
  double T2 = par[1];
  double f = par[2];
  double n = par[3];
  int id = int(par[4]+0.001);
  int idx = int(par[5]+0.001);
  double val;
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  
  val = f*exp(-E/T1);
  if (f < 1. && T2 > 0.) val += (1.-f)*exp(-E/T2);
    
  if (makeStaticData()->GetParticleTotalWidth(id) > THERMAL_WIDTH) 
      val *= makeDynamicData()->GetParticleTotalWeight(M,id,idx);
//  else
//       val *= makeDynamicData()->GetParticleTotalWeight(M,id);
  return 100000.*E*pow(p,n)*val;
}


double PThermal:: d2NdEdTheta(double* x, double* par) {
  // thermal source with E*p**n weight with T=T(theta_cm) and n=n(theta_cm)

  double E = x[0];
  double T1 = par[0];
  double T2 = par[1];
  double n1 = par[2];
  double n2 = par[3];
  double M = par[4];
  double p2 = E*E - M*M;
  if (p2 <= 0.) return 0.;
  double p = sqrt(p2);
  double A2 = par[5];
  double A4 = par[6];
  double cost = cos(x[1]);
  double cost2 = cost*cost;
  double sint = sin(x[1]);
  double val = sint*(1. + (A2 + A4*cost2)*cost2);
//double Theta = 57.29577951*x[1];
//double Thetar = Theta<90. ? Theta : 180.-Theta;
  double T = T1 + (T2-T1)*sint;  // interpolate T and n with angle
  double n = n1 + (n2-n1)*sint;
//double T = T1 + (T2-T1)*Thetar/90.;
//double n = n1 + (n2-n1)*Thetar/90.;

  double norm1 = IntThermal(M,T,int(TMath::Floor(n)));
  double norm2 = IntThermal(M,T,int(TMath::Ceil(n)));
  double frac = n-TMath::Floor(n);  // fractional exponent for interpolation
  double norm = norm1*pow(norm2/norm1,frac);

  // Must be normalized to keep an instrinsically flat angular distribution
  // if summed over all momenta.
  //
  return val*E*pow(p,n)*exp(-E/T)/norm;
}


double PThermal:: dNdTheta(double* x, double* par) {  
  // polar angular distribution

  double theta = x[0];
  double A2 = par[0];
  double A4 = par[1];
  double cost = cos(theta);
  double cost2 = cost*cost;
  return sin(theta)*(1. + A2*cost2 + A4*cost2*cost2);
}

double PThermal:: dNdy(double* x, double* par) {  
  // rapidity distribution

  double y = x[0];
  double m = par[0];
  double T = par[1];
  if (T<=0. || m<=0.) return 0.;
  double chi = T/(m*cosh(y));
  return m*m*T*(1.+2.*chi*(1.+chi))*exp(-1./chi);
}

double PThermal:: dNdMt(double* x, double* par) {  
  // transverse-mass distribution

  double mt = x[0];
  double m = par[0];
  double T1 = par[1];
  double T2 = par[2];
  double f = par[3];
  if (T1<=0.0 || mt<=m) return 0.;
  Double_t val = f*TMath::BesselK1(mt/T1)/T1;
  if (f < 1. && T2 > 0.) val += (1.-f)*TMath::BesselK1(mt/T2)/T2;
  return mt*mt*val;
}

double PThermal:: dNdPt(double* x, double* par) {  
  // transverse-momentum distribution

  double pt = x[0];
  double m = par[0];
  double mt = sqrt(pt*pt + m*m);
  return dNdMt(x,par)*pt/mt;
}


double PThermal:: IntThermal(double m, double T, int n) {
  //
  // return Integral|m,inf| E*p**n exp(-E/T)  n=0,1,2,3,4,5,6,7 
  //
  double sum = 1;

  switch (n) {
    case 0:
      sum = exp(-m/T)*T*(m+T);
      break;
    case 1:
      sum = m*m*T*TMath::BesselK(2,m/T);
      break;
    case 2:
      sum = 2.*exp(-m/T)*pow(T,2.)*(m*m+3.*T*(m+T));
      break;
    case 3:
      sum = 3.*m*pow(m*T,2.)*TMath::BesselK(3,m/T);
      break;
    case 4:
      sum = 8.*exp(-m/T)*pow(T,3.)*(m*m*(m+6.*T)+15.*T*T*(m+T));
      break;
    case 5:
      sum = 15.*m*pow(m*T,3.)*TMath::BesselK(4,m/T);
      break;
    case 6:
      sum = 48.*exp(-m/T)*pow(T,4.)*(m*m*m*(m+10.*T)+45.*m*m*T*T+105.*T*T*T*(m+T));
      break;
    case 7:
      sum = 105.*m*pow(m*T,4.)*TMath::BesselK(5,m/T);
      break;
    default:
      sum = 1;
  }
  return sum;
}



Double_t PThermal:: mtScaleFactor(Int_t id, const Float_t T) {
//
// Compute the yield enhancement factor due to thermal weighting of a
// resonance with id at temperature T.
// 
// Factor = 1/B(M0) x Integral[M=0,Inf](B(M)*BW(M))/Integral[M=0,Inf](BW(M)) 
//
// where BW(M) is the bare resonance shape and B(M) is the integral of the
// Boltzmann distribution: B(M) = Integral[E=M,Inf](exp(-E/T)*E*sqrt(E^2-M^2))
//                              = M^2*T*BesselK(2,M/T)
// 
// Tool function copied from PData (IF)

    Double_t val = 0.;
    if (T<=0.) return val;
    Double_t M0 = makeStaticData()->GetParticleMass(id);  // pole mass
    Double_t BM0 = M0*M0*T*TMath::BesselK(2,M0/T); // normalization
    Double_t I1 = 0.;
    Double_t I2 = 0.;
    Double_t M;
    for (Int_t i=1;i<1000;i++) {  // do a rudimentary integration from 0-5 GeV
	M = float(i)/200;
//	I1 += BreitWigner(id,M)*M*M*T*TMath::BesselK(2,M/T);
//	I2 += BreitWigner(id,M);
	I1 += makeDynamicData()->GetParticleTotalWeight(M,id)*M*M*T*TMath::BesselK(2,M/T);
	I2 += makeDynamicData()->GetParticleTotalWeight(M,id);
    }
    val = I1/(I2*BM0);
    printf("Factor = %f\n\n",val);
    return val;
}







