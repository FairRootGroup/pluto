// Author: M.A. Kagarlis
// Written: 15.12.98
// Revised: 17.10.00 by R. Holzmann
// Utilities Class Header

#ifndef _PUTILS_H_
#define _PUTILS_H_

#include "TRandom3.h"
#include "TMath.h"
using namespace std;
#include <iostream>
#include <cstdlib>
#include "TH1.h"
#include "PKinematics.h"
#include "PStaticData.h"


//#define SEED 65539   
// set to 0 for random initialization of seed at startup
// (this initializes TRandom3 with systime (but only a
// granularity of 1 second)

//This is the new default (IF 20.09.2009)
#define SEED 0


#define PLUTOVERSION_FOR_FAIR 1



class PUtilsREngine;
PUtilsREngine& fPUtilsREngine();
PUtilsREngine* makePUtilsREngine();


class PUtilsREngine : public TObject {
    
 public:
    PUtilsREngine ();
  
    Double_t sampleFlat() {
	// Samples uniformly between 0 and 1
	return rnd->Rndm();
    }
    
    Double_t sampleGaus(Double_t c, Double_t s) {
	// Samples from a Gaussian distribution
	return rnd->Gaus(c,s);
    }
    
    int samplePoisson(Double_t mean) {
	// Samples from a Poisson distribution
	return rnd->Poisson(mean);
    }
    
    Int_t sampleBinomial(Int_t ntot, Float_t prob) {
	// Samples from a binomial distribution
	return rnd->Binomial(ntot,prob);
    }
    
    Double_t sampleBW(Double_t c, Double_t g) {
	// Samples from a Breit Wigner distribution
	if (g==0.) return c;
	return (c+0.5*g*TMath::Tan((2.0*rnd->Rndm()-1.0)*TMath::Pi()*0.5));
    }
    
    void SetSeed(UInt_t s);

    double lambda(double M, double m1, double m2) {
	return PKinematics::lambda(M, m1, m2);
    }
    
    double pcms2(double M, double m1, double m2) {
	return PKinematics::pcms2(M, m1, m2);
    }
    
    double pcms(double M, double m1, double m2) {
	// cm momentum for the decay of M to m1 and m2
	return PKinematics::pcms(M, m1, m2);
    }


 private:
    
    TRandom3 * rnd;
    
    
    ClassDef(PUtilsREngine,0) //Pluto Utilities Class (random wrapper)

};




class PUtils : public TObject {

 public:
  PUtils () { 
      cout << "seed: " << SEED << endl;
      SetSeed(SEED); 

  }
  // create former behaviour

  static void dsort(Double_t *, int);
  // Sort in ascending order the first (int) entries of the array (Double_t *).
  // Adapted from Ref 1

  static void isort(int * i, int n) {
    // Sort in ascending order the first (int) entries of the array (int *).

      //BUGBUG: Quickersort is unstable
      //see example:
      //Int_t a[3]={8,9,9};
      //PUtils::isort(a,3);

      //workaround: add very small number

    Double_t x[n];

    for (int j=0;j<n;++j) x[j]=((Double_t)i[j]) + ((Double_t)j)*0.00001;
    dsort(x,n);
    for (int j=0;j<n;++j) i[j]=(int)x[j];
  }

  static Double_t sampleFlat() {
    // Samples uniformly between 0 and 1
    return makePUtilsREngine()->sampleFlat();
  }

  static Double_t sampleGaus(Double_t c, Double_t s) {
    // Samples from a Gaussian distribution
    return makePUtilsREngine()->sampleGaus(c,s);
  }

  static int samplePoisson(Double_t mean) {
    // Samples from a Poisson distribution

    return makePUtilsREngine()->samplePoisson(mean);
  }

  static Int_t sampleBinomial(Int_t ntot, Float_t prob) {
    // Samples from a binomial distribution

    return makePUtilsREngine()->sampleBinomial(ntot,prob);
  }

  static Double_t sampleBW(Double_t c, Double_t g) {
    // Samples from a Breit Wigner distribution
      return makePUtilsREngine()->sampleBW(c, g);
  }

  static Double_t cgc(const int &, const int &, const int &,
		    const int &, const int &);
  // Clebsch-Gordan coefficients (arguments are 2 x j or m)

  static Double_t s3j(const Double_t &, const Double_t &, const Double_t &,
		    const Double_t &, const Double_t &, Double_t m3=0.);
  // 3j-symbol, related to Clebsch-Gordan coefficient

  static Double_t racah(const int &, const int &, const int &,
		      const int &, const int &, const int &);
  // Racah coefficients (arguments are 2 x j or m)

  static Double_t s6j(const Double_t &, const Double_t &, const Double_t &,
		    const Double_t &, const Double_t &, const Double_t &);
  // 6j-symbol, related to Racah coefficient



  static Int_t FindIndex(Int_t n, Double_t *a, Double_t r);

  static void SetSeed(UInt_t s) { makePUtilsREngine()->SetSeed(s); }

 static Bool_t Tokenize(const char * options, const char *delimiter, char ** array, int *size);
 static void remove_spaces(char ** partc);
 static Int_t remove_brackets(char ** partc, char a, char b);
 static Bool_t ValidVariableName(const char * name);
 static Bool_t IsInt(const char * name);

 static void correct_histo(TH1 * histo);
 static void correct(TH1 * histo);

 private:

  static Double_t phasef(const int & n) { return 1.- 2.*(abs(n)%2); }
  // (-1)**n

  static Double_t j123(const int &, const int &, const int &);
  // used by racah


ClassDef(PUtils,0) //Pluto Utilities Class
};



#endif // _PUTILS_H_

