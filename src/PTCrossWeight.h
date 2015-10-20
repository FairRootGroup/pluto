// Author: I. Froehlich
// Written: 31.1.2008
// Revised: 

#ifndef _PTCROSSWEIGHT_H_
#define _PTCROSSWEIGHT_H_

#define TCROSS_MAX_DAUGHTERS 10

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TSpline.h"

class PTCrossWeight : public PChannelModel  {
  
 public:
    PTCrossWeight();
    PTCrossWeight(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    
    using PChannelModel::GetWeight;
	
    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
	
    void SetCrossSection(TF1 * c) {
	//Sets cross section as a function of excess energy
	tcross=c;
    };

    
    void SetCrossSection(TGraph * f, Bool_t useSpline=kFALSE) {
	//Sets cross section as a function of excess energy
	tcrossg=f;
	spline=useSpline;
	if (spline) g_spline=new TSpline3("",f);
    };

    void SetCrossSection(Double_t f) {
	//Sets cross section as a constant -> overrides everything
	tcross=NULL;
	tcrossg=NULL;
	tcrossc=f;
    };

    void SetScaling(Double_t scal) {
	scaling=scal;
    };

    void SetExcessEnergy(Bool_t a,Bool_t b) {
	isExcessEnergy=a;isMevEnergy=b;
    };

 private:

    TF1 * tcross;
    TGraph * tcrossg;
    Double_t tcrossc;
    TSpline3 * g_spline;
    Bool_t spline;
    Double_t scaling;

    Bool_t isExcessEnergy,isMevEnergy;
    
    PParticle *daughters[TCROSS_MAX_DAUGHTERS], *beam, *target, *parent; 
    int daughter_pos;

    ClassDef(PTCrossWeight,0)  //Total cross section weights

};

#endif


