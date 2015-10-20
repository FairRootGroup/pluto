// Author: IF
// Written: 25.05.2007
// Revised: 
// PDynamicData Class Header

#ifndef _PDYNAMICDATA_H_
#define _PDYNAMICDATA_H_

#define MAX_PCHANNEL_OBJECTS 10000

#include "TROOT.h"
#include "TF1.h"
#include <iostream>

#include "PDataBase.h"
#include "PDistribution.h"
#include "PMesh.h"
#include "PChannelModel.h"
#include "PParticle.h"

using namespace std;

class PDynamicData;
PDynamicData* makeDynamicData();
PDynamicData& fDynamicData();

class PDynamicData : public TObject {
      
 public:

    PDynamicData();

    //**** Decays
    Double_t GetDecayPartialWidth(Double_t mass,Int_t didx);
   
    PChannelModel *GetDecayModel(Int_t didx);
    PChannelModel *GetDecayModelByKey(Int_t key);
    PChannelModel *GetDecayModelByKey(Int_t key, Int_t defkey);
    bool SetDecayModel(Int_t didx, PChannelModel * model);
    bool SetDecayModelByKey(Int_t didx, PChannelModel * model);
    Double_t GetDecaySCFactor(Int_t didx);
    void SetDecaySCFactor(Int_t didx, Double_t factor);
    bool CheckSCAbort(Int_t didx);  //to stop endless loops
    double GetDecayBR(int , double); 
    // branching ratio by mode index & mass (GeV/c^2) 
    
    void ListDecayBR(const Char_t* pid, double mass) { 
	ListDecayBR(makeStaticData()->IsParticleValid(pid), mass); 
    } 
    void ListDecayBR(int pid, double mass); 
    // list branching ratios by particle id & mass (GeV/c^2) 

  
    int PickDecayChannel(const int &, const double &);
    // decay index for particle pid of given mass (GeV/c**2),
    // randomly selected consistent with the branching ratios

    int PickDecayChannel(const char * id, const double & m) {
	// as above, by particle name
	return PickDecayChannel(makeStaticData()->IsParticleValid(id),m); 
    }

    int PickDecayChannel(const int &, const double &, int *);
    // number and pids of decay products for the decay of particle 
    // pid of given mass(GeV/c**2), via a randomly selected mode 
    
    int PickDecayChannel(const char * id, const double & m, int * array) {
	// as above, by particle name
	return PickDecayChannel(makeStaticData()->IsParticleValid(id),m,array);
    }

    //**** Particles
    PChannelModel *GetParticleModel(Int_t pid);
    PChannelModel *GetParticleSecondaryModel(const char * name, const char * modelname);

    Double_t GetParticleTotalWidth(Double_t mass,Int_t pid);
    // mass-dependent total resonance decay width by id (units GeV/c**2)

    Double_t GetParticleTotalWidthSum(Double_t mass,Int_t id, Int_t flag=0);
    // same as above but by resumming the partial width (slower but more precise)
    // flag=1: take into account only hadronic decays

    Double_t GetParticleTotalWeight(Double_t mass,Int_t pid,Int_t didx=-1);
    Double_t GetParticleScalingFactor(Int_t didx);
    void SetParticleScalingFactor(Int_t didx, Double_t factor);

    int GetParticleDepth(const int & id, int flag=0);
    
    int GetParticleDepth(const char * id, int i=0) { 
	return GetParticleDepth(makeStaticData()->IsParticleValid(id),i);
    }
    
    
    
    double GetParticleLife(const int & id, double m=0, int idx=-1);
    // mean life
    // Arguments: 1. id=particle id 
    //            2. m=mass (GeV/c**2)
    //            3. idx=decay-mode index
    
    double GetParticleLife(const char * id, double m=0, int idx=-1) {
	// as above, by particle name
	return GetParticleLife(makeStaticData()->IsParticleValid(id),m,idx);
    }

    //Debug info... similar to PStaticData
    void PrintParticle(int pid);
    void PrintParticleByKey(int key);
    void PrintDecayByKey(int key);

    PParticle * GetBatchParticle(const char * name, Int_t make_val=1);

    void PushPChannel(TObject * obj) {
        if (num_pchannels==MAX_PCHANNEL_OBJECTS) {
            Warning("Push","MAX_PCHANNEL_OBJECTS reached");
            return;
        }
        pchannel_list[num_pchannels++]=obj;
    };

    TObject ** GetPChannels(Int_t * num) {
        *num=num_pchannels;
        if (num_pchannels) return pchannel_list;
        return NULL;
    }

    void DumpPChannels(void) {
        for (int i=0;i<num_pchannels;i++)
            pchannel_list[i]->Print();

    };


 private:
     Int_t *i_result;
     char  *c_result;
     Double_t *d_result;
     TObject *t_result;
     Int_t  pid_param;
     Int_t  name_param,model_param,didx_param,scfactor_param,sccount_param;
     Int_t  pnmodes_param, link_param;


     TObject * pchannel_list[MAX_PCHANNEL_OBJECTS]; //list of all created PChannels for bookkeeping
     Int_t num_pchannels;

  ClassDef(PDynamicData,0) //Pluto Dynamic Data Interface (for mass-dependent values)
};

#endif // _PDATAUTIL_H_







