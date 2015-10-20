// Author: Marios A. Kagarlis
// Written: 15.11.98
// Revised: 20.07.99
// Revised: 30.08.00  R. Holzmann
// Revised: 07.09.00  J. Ritman (added parent-,siblings-,daughter- index
//                    and decayTime)
// Revised: 15.12.00  R. Holzmann  (creation time added)
// Revised: 22.03.05  R. Holzmann  (get/setParent() added)
// Revised: 23.07.07  IF (new framework) 

// PParticle Class Header
 
#ifndef _PPARTICLE_H_
#define _PPARTICLE_H_

#include "TLorentzVector.h"
#include "PData.h"
#include "PValues.h"
#include "PUtils.h"

#include "TMath.h"



class PParticle: public TLorentzVector {
  
 public:
  PParticle(Int_t id=0, Double_t T=0., Double_t w=1.);
  // id, lab kinetic energy (GeV), weight

  PParticle(const char *, Double_t T=0., Double_t w=1.);
  // name, lab kinetic energy (GeV), weight

  PParticle(Int_t, Double_t, Double_t, Double_t, Double_t m=0., Double_t w=1.);
  // id, Px, Py, Pz (GeV/c), mass (GeV/c**2) overrides default, weight

  PParticle(const char *, Double_t, Double_t, Double_t,Double_t m=0., Double_t w=1.);
  // name, Px, Py, Pz (GeV/c), mass (GeV/c**2) overrides default, weight

  PParticle(Int_t, const TVector3 &, Double_t m=0., Double_t w=1.);
  // id, 3-momentum vector (GeV/c), mass (GeV/c**2) overrides default, weight

  PParticle(Int_t, Double_t *, Double_t w=1.);
  // id, pointer to 4-dim array (Px, Py, Pz, E) (GeV/c, GeV), weight
  
  PParticle(Int_t, float *, Double_t w=1.);
  // id, pointer to 4-dim array (Px, Py, Pz, E) (GeV/c, GeV), weight

  PParticle(const PParticle &);
  // copy constructor

  PParticle(const PParticle *);
  // copy constructor

  virtual ~PParticle() {
      if (qParticle1) delete qParticle1;
      if (qParticle2) delete qParticle2;
      if (values)  delete values;
  }

  Int_t Is(const char * id) { 
      return (pid==makeStaticData()->GetParticleID(id)); 
  }
  Int_t is(const char * id) { 
      //kept for compatibility only
      return Is(id); 
  }
  Int_t size() { return (pid<1000)?1:2; }
  Int_t IsNucleon() { return (Is("p")||Is("n")); }
  Int_t IsDelta() { return (Is("D0")||Is("D-")||Is("D+")||Is("D++")); }
  Int_t IsPi() { return (Is("pi0")||Is("pi+")||Is("pi-")); }
  Int_t IsRho() { return (Is("rho0")||Is("rho+")||Is("rho-")); }
  Int_t ID() const { return pid; }

  const char * Name() { 
      return makeStaticData()->GetParticleName(pid); 
  }

  void SetSpectator(Int_t s) {spectator=s;};
  //=0 don't care; =-1 never; =1 force

  Int_t  IsSpectator() {return spectator;};

  void SetID(const Int_t id) { pid=id; }
  Double_t W() const { return wt; }
  void SetW(Double_t w=1.) { wt=w; }
  void SetMultiplicity(Double_t w=1.) { mult=w; }
  Double_t GetMultiplicity() {return mult;};
  Double_t GenW() const { return genwt; }
  void SetGenW(Double_t w=1.) { genwt=w; }
  Double_t InvGenW() const { return invgenwt; }
  void SetInvGenW(Double_t w=1.) { invgenwt=w; }


  virtual Double_t KE() const { return E()-M(); } // kinetic energy
  void SetKE(Double_t T=0.);    // reset by kinetic energy
  void SetM(Double_t m=0.);     // reset by mass
  void SetMom(Double_t mom=0.); // reset by momentum

  void ResetE() {             // reset E to be consistent with mass(ID)
    if (makeStaticData()->GetParticleTotalWidth(pid)>1.e-3) 
	return;      // broad particle
    Double_t m = makeStaticData()->GetParticleMass(pid);            // get tabulated mass
    Double_t px = Px(), py = Py(), pz = Pz();
    Double_t lE = sqrt(m*m + px*px + py*py + pz*pz);
    SetPxPyPzE(px, py, pz, lE);
  }

  Double_t Life(Double_t m=0., Int_t idx=-1); // lifetime in lab frame
  TLorentzVector Vect4() const { return TLorentzVector(Vect(),E()); }
  void SetVect4(const TLorentzVector & v) {SetPxPyPzE(v[0],v[1],v[2],v[3]);}
  Int_t IsMeson() const { return makeStaticData()->IsParticleMeson(pid); }
  Int_t IsHadron() const { return makeStaticData()->IsParticleHadron(pid); }
  Int_t BaryonN() const { return makeStaticData()->GetParticleBaryon(pid); }
  Int_t LeptonN() const { return makeStaticData()->GetParticleLepton(pid); }
  Int_t Charge() const { return makeStaticData()->GetParticleCharge(pid); }
  Int_t Key() const { return makeStaticData()->GetParticleKey(pid); }
  void Reset(const Int_t id, const TLorentzVector & v, const Double_t w=1.) { 
    pid=id;
    wt=w;
    SetVect4(v);
  }
  virtual Int_t IsFireball() { return 0; }
  virtual Int_t IsDilepton() { return 0; }
  virtual Int_t IsFileInput() { return 0; }
  void Reset(const Int_t id=0, const Double_t px=0., const Double_t py=0.,
	     const Double_t pz=0., const Double_t e=0., const Double_t w=1.) {
    SetPxPyPzE(px,py,pz,e);
    pid=id;
    wt=w;
  }
  void Reset(const PParticle & p) {
    *this=p;
  }

  PParticle* Clone(const char*delme=NULL) const;

  bool SetValue(Int_t id , Double_t val) {
      if (!values) values=new PValues();
      return values->SetValue(id , val);
  };

  bool SetValue(Int_t id)  {
      return SetValue(id , 1.);
  }

  bool GetValue(Int_t id , Double_t * val) {
      if (!values) return kFALSE;
      return values->GetValue(id , val);
  };
   

  void SetStatus(Int_t st){
      status=st;
  };
  
  Int_t GetStatus(){
      return status;
  };

  Double_t InvariantT(Double_t m3, Double_t m4, Double_t cos_th_cm);
  // The Mandelstam invariant t in 1 + 2 --> 3 + 4 scattering
  
  PParticle operator + (const PParticle & b) const {
    // "addition" for composite quasi-particles
    PParticle v( *this );
    return v += b;
  }
  void Scatter(PParticle *p1, PParticle *p2);

  PParticle operator - (const PParticle & b) const {
    // "addition" for composite quasi-particles
    PParticle v( *this );
    return v -= b;
  }
  void Reconstruct(void);

  PParticle & operator += ( const PParticle & );
  PParticle & operator -= ( const PParticle & );
  PParticle & operator = ( const PParticle &);
  PParticle & operator *= (const TRotation &);
  PParticle & operator *= (const TLorentzRotation &);
  PParticle & Transform(const TRotation &);
  PParticle & Transform(const TLorentzRotation &);

  PParticle & AddTmp( const PParticle & p);

  void Print(const Option_t* delme= NULL) const;

  inline void SetParentId(Int_t pId) {parentId = pId;}
  inline void SetSourceId(Int_t sId) {sourceId = sId;}
  inline Int_t GetParentId() const {return parentId;}
  inline Int_t GetSourceId() const {return sourceId;}

  inline Int_t getSourceId() const {return sourceId;} //kept for backw. comp.

  inline void  SetParentIndex(Int_t pInd) {parentIndex = pInd;}
  inline Int_t GetParentIndex() const {return parentIndex;}
  inline void  SetDecayModeIndex(Int_t pInd, Int_t i=0) {
      if (i) {
	  decayModeIndex=-1;
	  destroyDecayModeIndex=pInd;	  
      } else {
	  decayModeIndex = pInd;
	  destroyDecayModeIndex = -1;
      }
  }
  
  //Internal models should use opt=1, they will not get the "wrong" getDecayModeIndex
  inline Int_t GetDecayModeIndex(Int_t opt=0) const {
      if ((opt ==0) && (decayModeIndex<0) && (destroyDecayModeIndex>0))
	  return destroyDecayModeIndex;      
      return decayModeIndex;
  }
  inline void  SetDaughterIndex(Int_t dInd) {daughterIndex = dInd;}
  inline Int_t GetDaughterIndex() const {return daughterIndex;}
  inline void  SetSiblingIndex(Int_t sInd) {siblingIndex = sInd;}
  inline Int_t GetSiblingIndex() const {return siblingIndex;}
  inline void SetIndex(Int_t Ind) {index = Ind;}
  inline Int_t GetIndex() const {return index;}

  inline Bool_t IsActive() const {return active;}
  inline void SetActive() {active = kTRUE;}
  inline void SetInActive() {active = kFALSE;}

  void SetProperTime();     // sample time until decay in proper time (sec)
  void SetProperTime(Double_t t) {decayTime = t;}
  inline Double_t GetProperTime() const {return decayTime;}  

  void SetVertex(Double_t x, Double_t y, Double_t z, Double_t t)
                                                 {fV.SetXYZ(x,y,z); fT=t;}
  inline void SetVertex(TVector3& v, Double_t t) {fV = v; fT=t;}
  inline TVector3& GetVertex() {return fV;}
  inline TVector3& getVertex() {return fV;}
  inline Double_t X() const {return fV.X();}
  inline Double_t Y() const {return fV.Y();}
  inline Double_t Z() const {return fV.Z();}
  inline Double_t R() const {return fV.Mag();}
  inline Double_t T() const {return fT;}
  inline void SetT(Double_t time) {fT = time;}
  inline void SetInput(PParticle *p) {pParticle = p;}
  inline PParticle* GetInput() {return pParticle;}
  inline void SetParent(PParticle *p) {pParticle = p;}
  inline PParticle* GetParent() {return pParticle;}
  inline void SetSibling(PParticle *p) {sParticle = p;}
  inline PParticle* GetSibling() {if (sParticle) return sParticle;
  else return this;}

  inline void ResetDaughters(void) {
    for (int i=0; i<(MAX_DAUGHTERS+1); i++)
      daughters[i]=NULL;
  };

  inline void SetDaughter(PParticle *d) {
    for (int i=0; i<MAX_DAUGHTERS; i++) {
      if (daughters[i]==NULL) {
	daughters[i]=d;
	return;
      }
    }
    Error("SetDaughter","MAX_DAUGHTERS reached");
  };

  inline PParticle *GetDaughter(int i) { 
    return daughters[i];
  }


  inline PParticle *GetScattering(Int_t i){
      if (i==0) return qParticle1; else return qParticle2;
  };

  void addDebugString(char * s) {
      debug.Append(s);  
      debug.Append(":");  
  }
  void clearDebugString() {
       debug="";  //BUGBUG memory leak?
  }
  TString * GetDebugString(){
      return &debug;
  }

  Int_t  GetScatterClone(void) {return make_new_qParticle;};
  void SetScatterClone(Int_t t) {make_new_qParticle=t;};

  Double_t weight_divisor; //!Helper for PReaction

  //Something for PBatch:
  Double_t Sample() {return PUtils::sampleFlat();};

 protected:

  void defaults(void);

  Int_t pid;           //  particle code  (partially conforms with Geant 3)
  Int_t sourceId;      //  Source ID
  Int_t parentId;      //  parent ID
  Int_t parentIndex;   //  parent index in particle array
  Int_t decayModeIndex;//  decay Mode index (for decayAll option)
  Int_t destroyDecayModeIndex;// save only for data file
  Int_t daughterIndex; //  daughter index
  Int_t siblingIndex;  //  sibling index
  Int_t spectator;     //! flag that forces particle to be treated as spectator

  Double_t decayTime;  //  proper time until decay  (in mm/c)
  Double_t wt;         //  weight
  Double_t genwt;      //! generator weight
  Double_t invgenwt;   //! inverted generator weight
  Double_t mult;       //! multiplicity

  TVector3 fV;         //  creation vertex (in mm)
  Double_t fT;         //  creation time (in mm/c)
  Bool_t active;       //! internal activity flag
  Int_t index;         //! index in particle array
  PParticle* pParticle;//! pointer to particle object
  PParticle* qParticle1;//!
  PParticle* qParticle2;//!
  PParticle* sParticle; //! pointer to particle object
  PParticle* daughters[MAX_DAUGHTERS+1]; //!pointer to daughter array 
  TString debug;        //! debug string

  PValues * values;   //!pointer to value container
  Int_t status;       //! status of parent particle in PChannel::Decay

  Bool_t make_new_qParticle; //! Workaround

  ClassDef(PParticle,4)  // Pluto Particle Class
};
#endif // _PPARTICLE_H_













