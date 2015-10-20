// Author: I. Froehlich
// Written: 23.5.2007
// Revised: 

#ifndef _PCHANNELMODEL_H_
#define _PCHANNELMODEL_H_

#include "TComplex.h"

#include "TF1.h"
#include "TF2.h"
#include "PDistribution.h"


class PChannelModel : public PDistribution {
  
 public:
    PChannelModel();
    PChannelModel(const Char_t *id, const Char_t *de, Int_t key=-1);
    PDistribution* Clone(const char*delme=NULL) const;

    using PDistribution::SampleMass;   
    virtual Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    Double_t SampleMass(Int_t didx) {
	Double_t mass;
	if (SampleMass(&mass,&didx)) return mass;
	return 0.;
    }

    //weight/amplitude calculations: To avoid endless loops,
    using PDistribution::GetWeight;   
    virtual Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    
    Double_t GetWeight(Double_t mass, Int_t *didx=NULL) {
	//wrapper for models with single variable only
	return GetWeight(&mass,didx);
    };
    Double_t GetWeight(Double_t mass, Double_t mass2, Int_t *didx=NULL) {
	//wrapper for models with 2 variables only
	Double_t m[2];
	m[0]=mass;m[1]=mass2;
	return GetWeight(m,didx);
    };
    Double_t GetWeight(Double_t mass, Double_t mass2, Double_t mass3, Int_t *didx=NULL) {
	//wrapper for models with 3 variables
	Double_t m[3];
	m[0]=mass;m[1]=mass2;m[2]=mass3;
	return GetWeight(m,didx);
    };
    Double_t GetWeight(Double_t mass, Double_t mass2, Double_t mass3, Double_t mass4, 
		       Int_t *didx=NULL) {
	//wrapper for models with 4 variables
	Double_t m[4];
	m[0]=mass;m[1]=mass2;m[2]=mass3;m[3]=mass4;
	return GetWeight(m,didx);
    };
    Double_t GetWeight(Double_t mass, Double_t mass2, Double_t mass3, Double_t mass4, Double_t mass5, 
		       Int_t *didx=NULL) {
	//wrapper for models with 5 variables
	Double_t m[5];
	m[0]=mass;m[1]=mass2;m[2]=mass3;m[3]=mass4;m[4]=mass5;
	return GetWeight(m,didx);
    };

    
    virtual TComplex GetAmplitude(Double_t *mass, Int_t *didx=NULL);


    virtual Bool_t GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    Double_t GetWidth(Double_t mass) {
	Double_t lWidth;
	GetWidth(mass, &lWidth, -1);
	return lWidth;
    };

    virtual Bool_t GetBR(Double_t mass, Double_t *br, Double_t totalwidth=-1.);

    virtual int GetDepth(int i=0);
    // Returns the total number of embedded recursive decays 
    // (TDepth if flag=0) or the number of embedded recursive
    // hadronic decays (HDepth if flag=1).
    // The default value 0 means the index has not been accessed yet.
    // Zero depth corresponds to index -1.

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    void SetDraw(Int_t opt) {draw_option=opt;};
    void SetDidx(Int_t opt) {didx_option=opt;};

    void SetDynamicRange(Double_t my_mmin,Double_t my_mmax) {
	mmin=my_mmin;
	mmax=my_mmax;
	fXmin      = mmin;
	fXmax      = mmax;
    };

    Int_t GetDef() {return model_def_key;};

    void SetWidthInit (Int_t w){width_init=w;};

 protected:
  
    Double_t width;     //Mass-dependent width container
    PMesh * mesh;       //Linear mesh container for fast lookup
    Int_t width_init;   //flag for the 1st initialization in getWidth
    Int_t is_channel,is_pid;
    Double_t mmin,mmax; //Dynamic range of the model
    Int_t maxmesh;
    Int_t draw_option;
    int mc_max,mc;      //Monte-Carlo intergration params
    Int_t didx_param, scfactor_param;
    Double_t *unstable_width;


    Int_t loop_flag; //To avoid endless loops of GetWidth/GetAmplitude calls

    PChannelModel * GetSecondaryModel(const char *);

    ClassDef(PChannelModel,0)  // Base class for coupled-channel distributions
};

#endif


