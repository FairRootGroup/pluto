// Author: I. Froehlich
// Written: 3.7.2006
// Revised: 
// Base class for common distributions

#ifndef _PDISTRIBUTION_H_
#define _PDISTRIBUTION_H_

#define MAX_PARTICLE_LIST      10
#define PARTICLE_LIST_DAUGHTER 1
#define PARTICLE_LIST_PARENT   2
#define PARTICLE_LIST_GRANDPARENT 4
#define PARTICLE_LIST_GRANDGRANDPARENT 8
#define PARTICLE_LIST_GRANDDAUGHTER 16
#define PARTICLE_LIST_SIBLING   128

#define DISTRIBUTION_QUASI_PID 1000
#define DISTRIBUTION_SOMETHING_PID 2000
#define DISTRIBUTION_NUCLEON_PID 3000

//Flags for sampling, generator, etc.
#define VERSION_SAMPLING          1
#define VERSION_WEIGHTING         2
#define VERSION_INVERT_WEIGHTING  4
#define VERSION_MASS_SAMPLING     128
#define VERSION_IS_PRIMARY        256
#define VERSION_GENERATOR_MC      512
#define VERSION_NO_PHASE_SPACE    1024

#include "TF1.h"
#include "TF2.h"
#include "PParticle.h"
#include "PData.h"
#include "PUtils.h"


class PDistribution : public TF1 {
  
 public:

    PDistribution();
    PDistribution(const Char_t *id, const Char_t *de);
    virtual ~PDistribution();
  
    virtual PDistribution* Clone(const char*delme=NULL) const;

    Int_t Add(const Char_t * opt);
    Int_t Add(const Char_t * name, int flag, 
	      const Char_t * pflag);  //Adds a particle placehoulder to the list
  
    Int_t Add(const Char_t * name, const Char_t *  flag1, 
	      const Char_t *  flag2 = NULL, const Char_t *  flag3 = NULL);
    //Helper function
  
    Int_t Set(const Char_t * opt);
    Int_t Set(const Char_t * name, const Char_t * pflag); //Like Add(), but rather replace it

    void Reset(void);
  
    Int_t SetParticle(PParticle *part, int pid, int flag); 
    //Replaces the pointer for the given particle
    //This function should be called from PChannel

    void  ResetRelatives(Int_t flag = 0);
    const char * GetGroupID(void) { //Some wish where I want to go to...
	return group;
    };
    void  SetGroupID(const char * gr) {
	group=gr;
    }

    PParticle *GetParticle(const Char_t *pflag=NULL); //get the particly with the private flag 
  
    Int_t GetStatus(void); //returns 0 if all particles are set
    Int_t CheckDaughters(void) ; //returns 0 if all daughters are set
    void  NoDaughters(void) {no_daughters=1;};
    void ResetStatus(void);  //Clean setted particles
  
    virtual void Print(const Option_t* delme=NULL) const ;  //Debug info
    void BasePrint(void) const ;  //Debug info
    const Char_t * OptString(void) {return opt_string;};
    virtual void SubPrint(Int_t opt) const ;  //Print sub-models
    const Char_t  *Path(void) const; //print path for models

    const Char_t  *GetIdentifier(void) {return identifier;};
    const Char_t  *GetDescription(void) {return description;};

    void  SetEnable(Int_t en){enable=en;};
    Int_t GetEnable(void) {return enable;};
    
    void  SetVersionFlag(Int_t f){version_flag |= f;};
    void  ClearVersionFlag(Int_t f){version_flag &= (~f);};
    UInt_t GetVersionFlag(UInt_t f = 0xffffffff)  const {
	return version_flag & f;
    };
    

    /* the following tool functions set different versions */
    void DisableSampling(void) {
	ClearVersionFlag(VERSION_SAMPLING);
    }; //No Sampling in the GenBod
    
    void EnableGenerator(void) {
	ClearVersionFlag(VERSION_WEIGHTING);
	SetVersionFlag(VERSION_INVERT_WEIGHTING);
    }; //Use distribution as a generator

    void EnableWeighting(void) {
	ClearVersionFlag(VERSION_SAMPLING);
	SetVersionFlag(VERSION_WEIGHTING);
    }; //Use distribution as a generator
    
    virtual Bool_t FreezeOut(void);

    virtual Bool_t Init(void);

    virtual Bool_t Prepare(void);
 
    virtual Bool_t SampleMass(void);

    virtual Bool_t SampleMomentum(void);
  
    virtual Bool_t SampleAngle(void);

    virtual Bool_t IsValid(void);

    virtual Bool_t CheckAbort(void);

    virtual Bool_t Finalize(void);

    virtual Bool_t EndOfChain(void);
 
    virtual Double_t GetWeight(void);

    virtual Bool_t WriteDebugInfo(PParticle *par);

    Int_t GetKey(void) {
	return primary_key;
    };
    void SetDecayIdx(Int_t opt) {didx_option=opt;};

    virtual int GetDepth(int i=0);

    int LinkDBdone(void) {return linkdbdone;};
    void LinkDBdone(int i)  {linkdbdone=i;};

    void IncrementWeightSum(Double_t w, Double_t sc=1.) {
//	if (w_num && (sc>(w_num  *1.1))) return;

	w_sum+=w*sc;
	w_num+=sc;
    };

    Double_t GetExpectedWeightMean() {
	//For PChannelModels: ExpectedWeightMean is set to the
	//static branching ratio

	return exp_w_mean;
    };

    void SetExpectedWeightMean(Double_t w) {
	exp_w_mean=w;
    };

    Double_t CalcWeightMean() {
	return w_sum/w_num;
    }

    void SetWeightMax(Double_t w) {
	w_max = w;
    }

    Double_t GetWeightMax(void) {
	return w_max;
    }

    Double_t GetDynamicRange() {
	//Used to get the Delta_x for the Monte-Carlo Integral
	return dynamic_range;
    };

    Bool_t Exec(const char * command); 
    virtual Bool_t ExecCommand(const char * command, Double_t value); 

    Int_t  debug_flag;  

    void PreParticles(void) {preliminary_particles=1;};

 protected:
  
    const Char_t *identifier, *description;
    Int_t enable;
    char *parse[4];  //Option strings
    Int_t parse_s; //max 4 options

    PParticle * particle[MAX_PARTICLE_LIST];
    const Char_t * names[MAX_PARTICLE_LIST];
    const Char_t * private_flag[MAX_PARTICLE_LIST];
    const Char_t * private_flag_int[MAX_PARTICLE_LIST];
    const Char_t * opt_string;  
    Int_t particle_flag[MAX_PARTICLE_LIST], pid[MAX_PARTICLE_LIST];
    Int_t position, current_flag;

    Int_t primary_key;   //Model coupled to data base key
    Int_t model_def_key; //Model type definition key
    Int_t sec_key;       //Key if secondary model (tmp)

    Int_t didx_option;  //Model has target didx

    Int_t preliminary_particles; //Used for auto-Adds in PChannelModel

    Int_t GetFlag(const Char_t *  flag);
    const Char_t * group;     //Group ID
    Int_t linkdbdone;
    Int_t no_daughters;
    UInt_t version_flag,relative_warning;

    Double_t w_sum,exp_w_mean,w_num,dynamic_range; //For weighting method
    Double_t w_max;


    ClassDef(PDistribution,0) //Base class for all distributions

};

#endif


	
