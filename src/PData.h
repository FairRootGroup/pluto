// Author: M.A. Kagarlis
// Written: 31.01.99
// Revised: 23.7.2007  IF
// PData Class Header

#ifndef _PDATA_H_
#define _PDATA_H_



#include "TROOT.h"
#include "TMath.h"
#include <iostream>
#include "PStaticData.h"

using namespace std;

class PData : public TObject {
  
    
 public:

#if 0
  static float eventPlane;  // event plane angle
  static float impactParam; // impact parameter
  static float gVertexX,gVertexY,gVertexZ; // global Vertex x,y,z
  static int   seqNr;       // Seq-Nr used in PEmbeddedParticle
  static int   WeightVersion; // Use old or new weigthing option
#endif

  static void  SetWeightVersion(Int_t i) {
      //PData::WeightVersion = i;
      *(makeStaticData()->GetBatchValue("_system_weight_version")) 
	  = (Double_t) i;
  }

  void SetUnstableWidth(Double_t m) {
      *(makeStaticData()->GetBatchValue("_system_unstable_width")) = m;
  };

#if 0
  static double unstable_width; // threshold when hadrons seems to be stable
  static int dbg;
  static void setQED(Int_t val=1) { 
      cout << "This setQED option is obsolete!" << endl;useQED = val; }
  static void setFlatMD(Int_t val=1) { 
      cout << "This setFlatMD option is obsolete!" << endl;flatMD = val; }
  static void setMtScalingTemp(Float_t T=0.) { // set temperature in GeV
    cout << "This mt-scaling option is obsolete!" << endl;
    tempMtScaling = T;
  } // obsolete since version 4.10
#endif


  static double LMass(const int & id) {
      // lower mass limit, by particle pid      

      Double_t *mass;

      if (makeDataBase()->GetParamDouble ((char*)"pid", id, (char*)"lmass", &mass)) {
	  //cout << "found lmass for: " << id << " with " << *mass << endl;
	  return *mass;
      }

      if ( makeStaticData()->GetParticleTotalWidth(id) < 1.e-3 ) {
	  // (quasi)stable
	  return 0.999*makeStaticData()->GetParticleMass(id);      
      }
//  else if ( id==41 ) return 2.*PMass[8];                 // rho0
      else if ( id==41 )                                   // rho0
	  return 2.*makeStaticData()->GetParticleMass(2);                 
      else if ( id==42 || id==43)                          // rho+,rho-
	  return makeStaticData()->GetParticleMass(8)+
	      makeStaticData()->GetParticleMass(7);  
      else if ( id==52 )                                   // omega
	  return 2.*makeStaticData()->GetParticleMass(2);                 
      else if ( id==55 )                                   // phi
	  return 3.*makeStaticData()->GetParticleMass(2);                 
      else if ( makeStaticData()->GetParticleBaryon(id) )  // resonance
	  return makeStaticData()->GetParticleMass(14);               
      else if ( makeStaticData()->GetParticleTotalWidth(id) > 0.01) 
	  return TMath::Max(0.,0.999*(makeStaticData()->GetParticleMass(id)-
				      5.*makeStaticData()->GetParticleTotalWidth(id)));
      else 
	  return TMath::Max(0.,0.999*(makeStaticData()->GetParticleMass(id)-
				      15.*makeStaticData()->GetParticleTotalWidth(id)));
  }

  
  static double LMass(const char * id) { 
      // lower mass limit, by name
      return LMass(makeStaticData()->IsParticleValid(id)); 
  }

  
  static double UMass(const int & id) {
      // upper mass limit, by particle pid

      Double_t *mass;
      if (makeDataBase()->GetParamDouble ((char*)"pid", id, (char*)"umass", &mass)) return *mass;

      if (id==50 || makeStaticData()->IsParticle(id,"dilepton")) // dilepton
	  return 1.e10;        
      else if ( makeStaticData()->GetParticleTotalWidth(id) < 1.e-3 )   // (quasi)stable
	  return 1.001*makeStaticData()->GetParticleMass(id); 
//      else return makeStaticData()->GetParticleMass(id)+1.;      // just go 1 GeV above pole
      else return makeStaticData()->GetParticleMass(id)+10.*makeStaticData()->GetParticleTotalWidth(id);  //better in units of sigma
  }

  
  static double UMass(const char * id) { 
      // upper mass limit, by name
      return UMass(makeStaticData()->IsParticleValid(id)); 
  }

#if 0
  static Int_t useQED;
  static Int_t flatMD;
  static Float_t tempMtScaling;
#endif



public:


  static int IsDelta(const int & id) {
    // check if id corresponds to a Delta
    return (id==makeStaticData()->GetParticleID("D++")||
	    id==makeStaticData()->GetParticleID("D+") ||
	    id==makeStaticData()->GetParticleID("D0") ||
	    id==makeStaticData()->GetParticleID("D-"));
  }

  static int IsN1535(const int & id) {
    // check if id corresponds to a S11(1535) resonance
    return (id==makeStaticData()->GetParticleID("NS11+")||
	    id==makeStaticData()->GetParticleID("NS110"));
  }
  
  static int IsPi(const int & id) {
    // check if id corresponds to a pion
    return (id==makeStaticData()->GetParticleID("pi0")||
	    id==makeStaticData()->GetParticleID("pi+")||
	    id==makeStaticData()->GetParticleID("pi-"));
  }
  
  static int IsN(const int & id) { 
      // check if id corresponds to a nucleon
      return (id==makeStaticData()->GetParticleID("p")||
	      id==makeStaticData()->GetParticleID("n")); }
  
  static int IsRho(const int & id) {
      // check if id corresponds to a rho meson
      return (id==makeStaticData()->GetParticleID("rho0")||
	      id==makeStaticData()->GetParticleID("rho+")||
	      id==makeStaticData()->GetParticleID("rho-"));
  }
  
  static int IsNPiPi(const int & i1, const int & i2, const int & i3) {
      // checks if i1+i2+i3 are N+pi+pi
      return ((IsN(i1)&&IsPi(i2)&&IsPi(i3))||
	      (IsN(i2)&&IsPi(i1)&&IsPi(i3))|| 
	      (IsN(i3)&&IsPi(i1)&&IsPi(i2)));
  }
  
  static int IsDalitz(const int & id, const int & i1, const int & i2) {
      // checks if the decay id->i1+i2 is a known Dalitz decay
      int eeg = 
	  (makeStaticData()->IsParticle(i1,"dilepton") && 
	   makeStaticData()->IsParticle(i2,"g")) || 
	  (makeStaticData()->IsParticle(i1,"g") && 
	   makeStaticData()->IsParticle(i2,"dilepton")); // e+e-gamma?
      int mumug = 
	  (makeStaticData()->IsParticle(i1,"dimuon") && 
	   makeStaticData()->IsParticle(i2,"g")) || 
	  (makeStaticData()->IsParticle(i1,"g") && 
	   makeStaticData()->IsParticle(i2,"dimuon"));   // mu+mu-gamma?
      int eepi = 
	  (makeStaticData()->IsParticle(i1,"dilepton") && 
	   makeStaticData()->IsParticle(i2,"pi0")) || 
	  (makeStaticData()->IsParticle(i1,"pi0") && 
	   makeStaticData()->IsParticle(i2,"dilepton")); // e+e-pi0?
      int mumupi = 
	  (makeStaticData()->IsParticle(i1,"dimuon") && 
	   makeStaticData()->IsParticle(i2,"pi0")) || 
	  (makeStaticData()->IsParticle(i1,"pi0") && 
	   makeStaticData()->IsParticle(i2,"dimuon"));   // mu+mu-pi0?
//  int eeeta = (Is(i1,"dilepton")&&Is(i2,"eta")) || (Is(i1,"eta")&&Is(i2,"dilepton")); // e+e-eta?
      
      int pseudo = 
	  (makeStaticData()->IsParticle(id,"eta") || 
	   makeStaticData()->IsParticle(id,"eta'") || 
	   makeStaticData()->IsParticle(id,"pi0")); // pseudo-scalar meson?
      int vector = 
	  (makeStaticData()->IsParticle(id,"w") || 
	   makeStaticData()->IsParticle(id,"phi") || 
	   makeStaticData()->IsParticle(id,"J/Psi"));  // vector meson?
      int D0 = 
	  makeStaticData()->IsParticle(id,"D0") &&  // Delta0 Dalitz decay?
	  ((makeStaticData()->IsParticle(i1,"dilepton") && 
	    makeStaticData()->IsParticle(i2,"n")) || 
	   (makeStaticData()->IsParticle(i1,"n") && 
	    makeStaticData()->IsParticle(i2,"dilepton")));
      int Dp = 
	  makeStaticData()->IsParticle(id,"D+") &&  // Delta+ Dalitz decay?
	  ((makeStaticData()->IsParticle(i1,"dilepton") && 
	    makeStaticData()->IsParticle(i2,"p")) || 
	   (makeStaticData()->IsParticle(i1,"p") && 
	    makeStaticData()->IsParticle(i2,"dilepton")));
      int pn = 
	  makeStaticData()->IsParticle(id,"pn") &&   // pn bremsstrahlung?
	  ((makeStaticData()->IsParticle(i1,"dilepton") && 
	    makeStaticData()->IsParticle(i2,"p")) || 
	   (makeStaticData()->IsParticle(i1,"p") && 
	    makeStaticData()->IsParticle(i2,"dilepton")));
      int NS0 = 
	  makeStaticData()->IsParticle(id,"NS110") && // N(1535)0 Dalitz decay?
	  ((makeStaticData()->IsParticle(i1,"dilepton") && 
	    makeStaticData()->IsParticle(i2,"n")) || 
	   (makeStaticData()->IsParticle(i1,"n") && 
	    makeStaticData()->IsParticle(i2,"dilepton")));
      int NSp = 
	  makeStaticData()->IsParticle(id,"NS11+") && // N(1535)+ Dalitz decay?
	  ((makeStaticData()->IsParticle(i1,"dilepton") && 
	    makeStaticData()->IsParticle(i2,"p")) || 
	   (makeStaticData()->IsParticle(i1,"p") && 
	    makeStaticData()->IsParticle(i2,"dilepton")));
      
      return (pseudo&&eeg) || (pseudo&&mumug) || 
	  (vector&&eepi) || (vector&&mumupi)
	  || D0 || Dp || pn || NS0 || NSp;
  }
  
  static int IsDalitz(int * i) {
    // as above, for pid array

    if (!i) return 0;
    return IsDalitz(i[0],i[1],i[2]);
  }

  static int IsDalitz(const int & idx) {
    // checks if idx is the index of a Dalitz decay mode

      int i[10];
      
      i[0]=2; //max 2 particles
      makeStaticData()->GetDecayMode(idx,i);
      if (*i!=2) return 0; // other than two products in this decay mode
      return IsDalitz(makeStaticData()->GetDecayParent(idx),i[1],i[2]);
  }

  static int IsPseudoscalarDalitz(const int & id, const int & i1, const int & i2) {
    // checks if the decay id->i1+i2 is a known pseudoscalar Dalitz decay

    int eeg = 
	(makeStaticData()->IsParticle(i1,"dilepton") && 
	 makeStaticData()->IsParticle(i2,"g")) || 
	(makeStaticData()->IsParticle(i1,"g") && 
	 makeStaticData()->IsParticle(i2,"dilepton")); // e+e-gamma?
    int mumug = 
	(makeStaticData()->IsParticle(i1,"dimuon") && 
	 makeStaticData()->IsParticle(i2,"g")) || 
	(makeStaticData()->IsParticle(i1,"g") && 
	 makeStaticData()->IsParticle(i2,"dimuon"));   // mu+mu-gamma?

    int pseudo = 
	(makeStaticData()->IsParticle(id,"eta") || 
	 makeStaticData()->IsParticle(id,"eta'") || 
	 makeStaticData()->IsParticle(id,"pi0")); // pseudo-scalar meson?
    
    return (pseudo&&eeg) || (pseudo&&mumug);
  }

  static int IsPseudoscalarDalitz(int * i) {
    // as above, for pid array

    if (!i) return 0;
    return IsPseudoscalarDalitz(i[0],i[1],i[2]);
  }

  static int IsPseudoscalarDalitz(const int & idx) {
    // checks if idx is the index of a Dalitz decay mode
        int i[10];

    i[0]=2; //max 2 particles
    makeStaticData()->GetDecayMode(idx,i);
    if (*i!=2) 
	return 0; // other than two products in this decay mode
    return 
	IsPseudoscalarDalitz(makeStaticData()->GetDecayParent(idx),i[1],i[2]);
  }

  static int IsDirectEE(const int & id, const int & i1, const int & i2) {
    // checks if the decay id->i1+i2 is a direct vector-meson ee decay
    
    return 
	(makeStaticData()->IsParticleMeson(id) && 
	 ((makeStaticData()->IsParticle(i1,"e+") && 
	   makeStaticData()->IsParticle(i2,"e-"))||
	  (makeStaticData()->IsParticle(i1,"e-") && 
	   makeStaticData()->IsParticle(i2,"e+"))));
  }

  static int IsDirectEE(const int & idx) {
    // checks if idx is the index of a direct vector-meson ee decay
      int i[10];
    i[0]=2; //max 2 particles
    makeStaticData()->GetDecayMode(idx,i);
    if (*i!=2) 
	return 0; // other than two products in this decay mode
    return 
	IsDirectEE(makeStaticData()->GetDecayParent(idx),i[1],i[2]);
  }

  static int IsDirectMuMu(const int & id, const int & i1, const int & i2) {
    // checks if the decay id->i1+i2 is a direct vector-meson mumu decay
    
    return 
	(makeStaticData()->IsParticleMeson(id) && 
	 ((makeStaticData()->IsParticle(i1,"mu+") && 
	   makeStaticData()->IsParticle(i2,"mu-"))||
	  (makeStaticData()->IsParticle(i1,"mu-") && 
	   makeStaticData()->IsParticle(i2,"mu+"))));
  }

  static int IsDirectMuMu(const int & idx) {
    // checks if idx is the index of a direct vector-meson mumu decay
      int i[10];
 
    i[0]=2; //max 2 particles
    makeStaticData()->GetDecayMode(idx,i);
    if (*i!=2) return 0; // other than two products in this decay mode
    return IsDirectMuMu(makeStaticData()->GetDecayParent(idx),i[1],i[2]);
  }

  // static int isWide(const int &);
  // returns total width index

  static int LPW(const int &, const int &, const int &);
  // lowest allowed transition partial wave, for parent id --> pid = i1 & i2

  static int IsMDalitz(const int &);
  // checks if decay mode with index idx is a mesonic Dalitz decay
  

  static Float_t mtIntegral(Double_t m, Float_t T) { // integral of thermal
                                                     // mt distribution
      if (m<=0. || T<=0.) return 0.;
      return T*exp(-m/T)*sqrt(m)*(1.5*T+m)
	  + 1.329340*pow((double)T,(double)2.5)*(1.-TMath::Erf(sqrt(m/T)));
  }
  

/*   static Double_t mtScaleFactor(const Int_t id, const Float_t T); */
/*   static Double_t mtScaleFactor(const char* id, const Float_t T) { */
/*       return mtScaleFactor(makeStaticData()->IsParticleValid(id),T); */
/*   } */
  //moved to PTHermal
  // compute the yield enhancement due to thermal weighting of a resonance 

  ClassDef(PData,1) //Pluto Particle Data Tool Class
};


class PSplash;
R__EXTERN PSplash *gSplash;

class PSplash : public TObject {

 public:

    PSplash();
    
 protected:
    ClassDef(PSplash,0) //Pluto welcome message

};

#endif // _PDATA_H_







