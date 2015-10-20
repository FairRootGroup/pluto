// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PDALITZDECAY_H_
#define _PDALITZDECAY_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PDalitzDecay : public PChannelModel  {
  
 public:
    PDalitzDecay();
    PDalitzDecay(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t SampleMass(void);
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    using PDistribution::GetWeight;
    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    int GetDepth(int i=0);

    Bool_t FreezeOut(void);

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    void SetUseQED(int t){useQED=t;};
    virtual void Print(const Option_t* delme=NULL) const ;  //Debug info
    Double_t draw_parent_mass;

 protected:
  
    PParticle *parent, *dilepton, *other;

    long double alpha;
    double mass_pi0;
    double mass_n;
    double mass_eta;
    double mass_p;
    double p1;

    double mass_e, mass_ee;
    
    double p2;
    double p3;
    double mass_x;   // mass of second (non dilepton) product
    double mass_parent;
    int flag;    // case flag
    int sw;      //0 for ee, 1 for mumu
    double ml;   // dilepton mass threshold

    int rejection_flag; //Case flag for the mass sampling

    int pi0, eta, eta_prime, w, Delta_0, Delta_plus, phi, S11_0, S11_plus; //PIDs
    int parent_id,dilepton_pid,other_pid;
    int dilepton_position;  //is dilepton particle 1 or 2?
    int others_position;
    double photon_br;  //PhotonBR

    int useQED; //use QED FF and not VMD
    int flatMD; // for test purposes sample a flat mass distribution

    PMesh *integral;                 //parent-mass dependent integral
    PChannelModel *formfactor_model; //form factor object


    //The following private methods are copied from PData
    //they are used here by the official wrapper funcions
    virtual double dGdM(const int& id, const double& m, const double& ecm);
    double FDalitz(const int & id, const double & m, double ecm);
    void sampleMD(const double & ecm, const int & id, 
		  double & m, const double & m1);
    double PhotonBR(const int & id);
 
    ClassDef(PDalitzDecay,0)  // Dalitz decays "a -> dilepton/dimuon + b"
};

#endif


