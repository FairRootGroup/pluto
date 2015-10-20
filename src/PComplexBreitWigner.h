// Author: I. Froehlich
// Written: 28.6.2007
// Revised: 

#ifndef _PCOMPLEXBREITWIGNER_H_
#define _PCOMPLEXBREITWIGNER_H_

#define COMPLEX_MAX_DECAYCHANNELS 30
#define COMPLEX_MAX_TERMS 10

#include "TF1.h"
#include "TF2.h"
#include "PBreitWigner.h"


class PComplexBreitWigner : public PBreitWigner  {
  
 public:
    PComplexBreitWigner();
    PComplexBreitWigner(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    using PDistribution::GetWeight;   
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    TComplex GetAmplitude(Double_t *mass, Int_t *didx=NULL);

    using PDistribution::SampleMass;
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    void SetUpdateAmplitude(int i){updateAmplitude=i;};

    void AddInterference(int idx,int key,int didx,double ampl,double phase);
    void AddAmplitude(int idx,double ampl,double phase);

    virtual void Print(const Option_t* delme=NULL) const ;  //Debug info

 private:
    
    void ReadModes(void);
    void ReadModels(void);
    int readModesDone,readModelsDone,updateAmplitude;
    TComplex GetAmplitudeLocal(Double_t *mass, Int_t num);

    Int_t num_decaychannels;

    PChannelModel *p[COMPLEX_MAX_DECAYCHANNELS][COMPLEX_MAX_TERMS+1]; //Amplitude models for coherent additions
    Int_t    int_index[COMPLEX_MAX_DECAYCHANNELS]; //Index numbers of decaychannels
    Int_t    num_terms[COMPLEX_MAX_DECAYCHANNELS];
    Double_t int_phase[COMPLEX_MAX_DECAYCHANNELS][COMPLEX_MAX_TERMS+1]; //Phase for each term
    Double_t int_ampl[COMPLEX_MAX_DECAYCHANNELS][COMPLEX_MAX_TERMS+1];  //Ampl for each term
    Int_t    int_key[COMPLEX_MAX_DECAYCHANNELS][COMPLEX_MAX_TERMS+1];  //Key pointing to the term
    Int_t    int_didx[COMPLEX_MAX_DECAYCHANNELS][COMPLEX_MAX_TERMS+1];  //Target didx
    


    ClassDef(PComplexBreitWigner,0)  // Breit Wigner (complex version) with mass-dependent width
};

#endif


