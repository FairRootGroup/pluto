// Author: I. Froehlich
// Written: 22.01.2010
// Revised: 

#ifndef _PHADRONDECAYM1N_H_
#define _PHADRONDECAYM1N_H_

#include "TF1.h"
#include "TF2.h"
#include "TGenPhaseSpace.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"
#include "PMesh.h"


#define M1N_MAX_DAUGHTERS 10
#define M1N_PARENT_GRID 0.1

class PHadronDecayM1N : public PChannelModel  {
  
 public:
    PHadronDecayM1N();
    PHadronDecayM1N(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;


    Bool_t Prepare(void);
    Bool_t Init(void);
    Bool_t SampleMass(void);
    using  PChannelModel::SampleMass;


    int      GetDepth(int i=0);

    Double_t Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    Double_t EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    void SubPrint(Int_t opt) const ;  //Print sub-models

 protected:

    //int didx1,didx2,didx3;
    PParticle * parent, *daughters[M1N_MAX_DAUGHTERS], *unstable_particle; 
    Int_t daughter_pos,unstable_daughter;
    PChannelModel *unstable_model;

    Int_t unstable_pos,daughter_id[M1N_MAX_DAUGHTERS],daughter_didx[M1N_MAX_DAUGHTERS],unstable_id,parent_id,unstable_didx;
    Double_t unstable_width,daughter_masses[M1N_MAX_DAUGHTERS],parent_mass,old_parent_mass;
    
    PMesh *mesh;
    TGenPhaseSpace *event;

    void UpdateMesh(void);
    Int_t maxmesh;

    ClassDef(PHadronDecayM1N,0)  // Decay of Hadron -> 1 unstable + N stable
};

#endif


