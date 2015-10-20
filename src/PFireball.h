// Author: Romain Holzmann
// Written: 28.04.00
// Revised: 14.08.00 MK
// Revised: 05.09.05 RH
// Revised: 28.03.06 RH
// Revised: 02.02.07 RH
// Revised: 21.02.07 RH
// PFireball Class Header

#ifndef _PFireball_H_
#define _PFireball_H_

#include "PF2.h"
#include "TH2F.h"
#include "PParticle.h"
#include "PChannel.h"
#include "PChannelModel.h"
#include "PFunction.h"
#include "PAdaptiveMeshN.h"

class PFireball: public PParticle {

 protected:

  int prodId;    // product id
  float T1;      // temperature of component 1 (in GeV)
  float T2;      // temperature of component 2 (in GeV)
  float frac;    // fraction of component 1  (0 < frac < 1) 
  float blast;   // radial blast velocity  (0 < blast < 1)
  float power1;  // power of p**n term in p**n exp(-E/T) sampling at 0 deg
  float power2;  // power of p**n term at 90 deg
  float A2;      // polar anisotropy (dN/dOmega = 1 + A2*cost**2 + A4*cost**4)
  float A4;      // 
  float v1;      // directed flow parameter at mid-rapidity
  float v2;      // elliptic flow (squeeze out) at mid-rapidity
  float Ap;      // projectile mass
  float At;      // target mass
  float prob;    // particle production probability per participant nucleon
  float bmin;    // min impact parameter 
  float bmax;    // max impact parameter
  float meanN;   // mean multiplicity 
  int nProd;     // number of product particles
  bool flag;     //! if set, sample impact parameter first, then Poisson Npart
  float mt_fac;  //! mt scaling of mass distribution
  double sig;    //! sigma of gaussian dNdy distribution
  int spect;     //! if set to 0 participant, if set to 1 target-like spectator,
                 //  if set to 2 projectile-like spectator

  PF2* fE;       //! energy sampling function (true Boltzmann) x BW(M)
  PF2* fE1;      //! energy sampling function with E*sqrt(p) x BW(M)
  PF2* fE2;      //! energy sampling function with E*p**3 x BW(M)
  PF2* fE3;      //! energy sampling function with p**n x W(Theta)
  TF1* fA;       //! polar angle sampling function
  TF1* fMt;      //! mt sampling function
  TH2F* fHisto;  //! dN2/dEdTheta sampling histogram
  bool trueThermal; //! thermal sampling flag
  Int_t npx, npy; //! Granularity for all 2-dim histograms

  PParticle** part; //! array of particle pointers
  bool quasistable; //! long-lived or stable particle flag

  PChannelModel * model; //! Primary model for emitted particle prodId
  int didx_old;          //! cache for didx
  int sample_option;     //! option for a rotation in the E-M-plane

  int mesh_option;     //! option to use the adaptive mesh
  PFunction *pfE;      //!Envelope for thermal channel model E
  PAdaptiveMeshN *afE; //!Mesh for  thermal channel model E

 public:
  PFireball(char* name, float AGeV, float t1, float t2=0.0, float f=1.0,
  float b=0.0, float a2=0.0, float a4=0.0, float w1=0.0, float w2=0.0, int sp=0);

  void setTemperature(float t1, float t2, float f, int id=0) {
    T1 = t1; T2 = t2; frac = f;
    updateFunctions(id);
    if (mt_fac>0.) mt_fac = mtIntegral(makeStaticData()->GetParticleMass(prodId),T1);
  }
  void setSigma(double s) {sig=s;} 
  void setBlast(float b, int id=0) {
    blast = b;
    updateFunctions(id);
  }
  void setAnisotropy(float a2, float a4, int id=0) {
    A2 = a2;
    A4 = a4; 
    updateFunctions(id);
  }
  void setFlow(float a, float b) {v1 = a; v2 = b;}

  Double_t AvApart(Double_t ap, Double_t at, Double_t bn, Double_t bx);
  Double_t AvApart(Double_t ap, Double_t at) {
    return (ap*pow(at,0.667) + at*pow(ap,0.667))
           /pow((pow(ap,0.333)+pow(at,0.333)),2);  // average Apart
  }
  void setRandomB(float ap, float at, float p=0., float bn=0., float bx=20.) {
    Ap = ap; At = at; prob = p; flag = 1; bmin = bn; bmax = bx;
    if (prob<=0.) prob = meanN/AvApart(Ap,At);
    meanN = prob*AvApart(Ap,At,bmin,bmax);
    if (bmax > 1.14*(pow((double)Ap,(double)0.333) + pow((double)At,(double)0.333))+2. )
      bmax = 1.14*(pow((double)Ap,(double)0.333) + pow((double)At,(double)0.333)+2.);
    if (bmin > bmax) {float temp=bmin; bmin=bmax; bmax=temp;}
  }
  void setMeanN(float mN) {if (mN >= 0.) meanN = mN;}
  void setSpectator(int sp){spect=sp;}

  bool IsRandomB() {return flag;}
  bool IsRandomN() {return meanN>0. ? 1 : 0;}

  int getParticleId() {return prodId;}
  float getT1() {return T1;}
  float getT2() {return T1;}
  float getFrac() {return frac;}
  float getBlast() {return blast;}
  float getA2() {return A2;}
  float getA4() {return A4;}
  float getV1() {return v1;}
  float getV2() {return v2;}
  int getSpectator() {return spect;}
  PF2* getFuncE() {return fE;}
  PF2* getFuncE1() {return fE1;}
  PF2* getFuncE2() {return fE2;}
  PF2* getFuncE3() {return fE3;}
  TF1* getFuncA() {return fA;}
  TF1* getFuncMt() {return fMt;}
  void setHisto(TH2F* pH) {
    fHisto = pH;
    A2=A4=0.;
  }
  TH2F* getHisto() {return fHisto;}

  void SetNpx(Int_t my_npx);
  void SetNpy(Int_t my_npy);
  void SetEpsilon(Double_t e);

  void RotateSamplingPlane(void) {
      sample_option =1;
      updateFunctions();
  }
  void UseMesh(void) {
      mesh_option =1;
  }


  void sampleECM(Double_t &E, Double_t &M, int didx);

  void sampleECM1(Double_t &E, Double_t &M, int didx) {
      if (didx != didx_old) {
	  didx_old = didx; //prevent random initialization each time
	  if (didx<0) 
	      fE1->SetParameter(5,-1.1);
	  else
	      fE1->SetParameter(5,(double) didx);
      }
      fE1->GetRandom2(E,M);}
  void sampleECM2(Double_t &E, Double_t &M, int didx) {
      if (didx != didx_old) {
	  didx_old = didx; //prevent random initialization each time
	  if (didx<0) 
	      fE2->SetParameter(5,-1.1);
	  else
	      fE2->SetParameter(5,(double) didx);
      }
      fE2->GetRandom2(E,M);}
  void sampleECM3(Double_t &E, Double_t &Theta) {fE3->GetRandom2(E,Theta);}
  double sampleThetaCM() {return fA->GetRandom();}
  double sampleMt() {return fMt->GetRandom();}
  void samplePartCM(double &px, double &py, double &pz, double &E, int didx);
  float sampleB() const ;
  int sampleNProd();
  int sampleNProd(float b);
  int getLastNProd() {return nProd;}
  int IsFireball() { return 1; }
  PChannel* makeChannel(Int_t nMax, Float_t nAverage=0.);
  //  virtual double KE() {return Rapidity();}
  void setTrueThermal(Bool_t flag=kTRUE) {trueThermal = flag;} 
  void setMtScaling() {mt_fac = mtIntegral(makeStaticData()->GetParticleMass(prodId),T1);}
  float mtScale(double m) {return mt_fac>0. ? mtIntegral(m,T1)/mt_fac : 1.;}
  virtual void Print(const Option_t* delme="") const;
  void printAverages() const ;
  virtual ~PFireball() {
    delete fE;
    delete fE1;
    delete fE2;
    delete fA;
    delete fMt;
  }
  float mtIntegral(double mass, float temperature);

  void SetW(double w=1.) {
      // Warning("SetW","the use of SetW() is depriciated, use SetMultiplicity() and change you macro");
      
      if (*(makeStaticData()->GetBatchValue("_system_weight_version"))) {
	  Info("SetW","Use old weighting method (pure chain based)");
      }

      *(makeStaticData()->GetBatchValue("_system_weight_version")) = 0.;
      
      PParticle::SetW(w);
  };

 protected:
  void updateFunctions(int id=0);
  void setToMidrapidity(float AGeV);
  float Npar(float ap, float at, float b);
  float NparSmeared(float ap, float at, float b) const ;

ClassDef(PFireball,1) // Pluto Fireball Class
};
#endif // _PFireball_H_


