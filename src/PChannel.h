// Author: M.A. Kagarlis
// Written: 21.01.99
// Revised: 05.12.05  by R. Holzmann
// Revised: 29.08.06  by R. Holzmann
// PChannel Class Header

#ifndef _PCHANNEL_H_
#define _PCHANNEL_H_

#include "TArrayI.h"
#include "TF1.h"
#include "TF2.h"
#include "PSaid.h"
#include "PParticle.h"
#include "PDistribution.h"
#include "PStdData.h"

#define MAX_DISTRIBUTIONS 10
#define MAX_NUM_NOT_FINALIZED 5

#define STATUS_OK          0
#define STATUS_EMPTY_EVENT 80
#define STATUS_PSEUDO_EMPTY_EVENT 81
#define STATUS_REDO_CHAIN  82
#define STATUS_NOT_DECAYED 99



class PChannel : public TObject {
  
 public:

    PChannel(PParticle **, int nt=2, int mf=0, int af=0, int bf=0);
    // Channel constructor by particle array

    PChannel(int idx, PParticle ** particles=NULL, int mf=0, int af=0, int bf=0);
    // Channel constructor by decay-mode index

    PChannel(int idx, PParticle & parent, int mf=0, int af=0, int bf=0);
    // Channel constructor by decay-mode index and parent reference

    ~PChannel() {
	// Channel destructor

	if (ipid.GetArray()) ipid.~TArrayI();
	// necessary, because of constr 2+3, for constr 1 also better
	if (ptcls && ptcls_extern>=0) {
	    for(int i=ptcls_extern; i<=n; i++)
		if (ptcls[i]) delete ptcls[i];
	    if (ptcls) delete [] ptcls;
	}
    }

    Bool_t Init();
    // sets on the first call all parent, siblings etc IDs. This
    // has to be done only once per PReaction

    Bool_t Reset();
    // cleans the channel, so an init will be done in the next decay()
  

    Bool_t SetDaughters();

    int Decay();
    // N-body phase-space decay function based on the CERNLIB routine 
    // GENBOD. Resonance and dilepton masses (for Dalitz decays)
    // are sampled in PData. Scattering angles are sampled here, for
    // selected channels (see idChannel).

    int decay() {return Decay();}; //kept for comp.

    PParticle ** GetParticles() {return ptcls;}
    // address of array of pointers to the channel particles

    int * GetPids() { return ipid.GetArray(); }
    // address of array of channel-particle ids

    int GetNumPar() { return n; }
    // number of decay products

    
    //******for distribution handling:
    int SetDistribution(PDistribution * distribution);
    // return -1 if distribution does not match the channel
    int CheckSiblings(PParticle * p, PDistribution * dist, int flag); //check for siblings


    double GetBT() {
	PParticle *beam=ptcls[0]->GetScattering(0);
	PParticle *target=ptcls[0]->GetScattering(1);
	
	if (beam) {
	    return (beam->KE() > target->KE()? beam->KE(): target->KE());
	}
	else return ptcls[0]->Rapidity();
    }

    //******for quasi-elastic handling:
    PChannel * GetQuasi(void) {return quasi_pchannel;};
    void ClearQuasi(void) {quasi_pchannel=NULL;};

    int GetParentSize() { return 1+(ipid[0]>=1000); }
    // 1 or 2 for elementary or quasi particles respectively

    static void SetGlobalWeight(double w) { 
	globalWeight = w; 
    } 
    static double GetGlobalWeight() { 
	return globalWeight; 
    }

    int GetDMIndex() { 
	//returns the decay mode index of this channel
	return DMIndex; 
    }

    void GetMessage();
    // status message

    void Print(const Option_t* delme="") const;
    void PrintReaction(Int_t check_key = 1) const;
    void PrintNew();
    

    void PrintReport() const; //Print a final report after sampling
    void CheckDecayKey() const;

    char const *GetName(void) const ;

//  double ds_dt(double);
    // ...> moved to PParticle (IF)
    // ds/dt(cos_th_cm) in the cm for N+N->N+Delta. With the mass of Delta
    // sampled independently in PData, the simulated events are distributed
    // as ds/dOmega (normalized). Ref: NPA459 (1986) 503

    void SetPrintTentative(bool t){ 
	// Print the tentative messages on the 1st call
	print_tentative=t;
    };

    void DisableHelicityAngle() {
	cout << "DisableHelicityAngle() is obsolete -> use DistributionManager" << endl;
    };

    //this is for envelope distributions:
    int GetNumNotFinalized() {return num_not_finalized;};
    PDistribution *GetDistributionNotFinalized(int i) {return distribution_not_finalized[i];};


 private:

    int n, status,DMIndex;
    // number of products, error status
    // decay-mode index matching PData modes

    int thSrc, tcSrc, dlSrc;

    Double_t weight_sum;   //Sum of all weights for debugging

    int decay_key;                       //! Key to data base entry
    TString decay_string; //! For parsing the decay string
    TString decay_string2; 

    // For generic distributions
    PDistribution * dist[MAX_DISTRIBUTIONS];
    Double_t dist_weight[MAX_DISTRIBUTIONS];
    Double_t dist_weight_sum[MAX_DISTRIBUTIONS];
    Double_t dist_counter[MAX_DISTRIBUTIONS];

    Int_t distribution_position;
    bool init_done, thermal_disable_didx;
    bool print_tentative;
    Long_t fEnablePattern;

    int bid,tid;
    PParticle *spectator,*participant,*quasi_composite;
    PChannel *quasi_pchannel;
    // general quasifree particle production on the deuteron  (V. Hejny 13/10/00)

    TArrayI ipid;
    // particle id array

    static double globalWeight;
    // global weight for the event (see PDecayManager, BUGBUG: never used)

    double *event_impact_param , *event_plane;  //from batch system
    double *weight_version;  //from batch system

    PParticle ** ptcls, *parent, *orig_parent, *grandparent, *grandgrandparent;
    // array of channel particles

    int ptcls_extern;
    // flag whether ptcls are coming from external source (do not delete those!)

    double ecm, w, emin; // invariant mass, weight of decay, energy threshold
    double e_cm;  // used to save parameters between calls of ds_dt()


    void IsReaction();
    // if a reaction channel identify beam, target, and decay mode

    void IdChannel();
    // identify reaction channel of known angular distribution
    // these are: p+p->p+p (elastic), N+N->N+Delta, pi+N->N+w, pi+ + P -> pi+ + P + w

    void ThermalSampling();
    // parent is midrapidity thermal source

    Int_t ReadFileInput();
    // parent is file input interface

    void MakeDilepton();
    // parent is a sampled dilepton


    int Genbod(int);
    // N-body phase-space decay of a particle in its rest frame.
    // This follows the CERNLIB routine GENBOD, in C++

    //this is for envelope distributions:
    int num_not_finalized;
    PDistribution *distribution_not_finalized[MAX_NUM_NOT_FINALIZED];

    
    ClassDef(PChannel,0) //Pluto Channel Class

};

#endif // _PCHANNEL_H_
