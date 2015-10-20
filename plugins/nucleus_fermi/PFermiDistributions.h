// Author: L. Witthauer &
//         M. Dieterle
// Written 24.03.2009

#ifndef _PFERMIDISTRIBUTIONS_H_
#define _PFERMIDISTRIBUTIONS_H_

#include "PUtils.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>

#include "PChannelModel.h"

#define NUCLEAR_FERMI_D    1
#define NUCLEAR_FERMI_3HE  2
#define NUCLEAR_FERMI_4HE  3
#define NUCLEAR_FERMI_7LI  4
#define NUCLEAR_FERMI_12C  5
#define NUCLEAR_FERMI_40CA 6



class PFermiDistributions : public PChannelModel {
    
 public:
    
    PFermiDistributions();
    PFermiDistributions(const Char_t *id, const Char_t *de, Int_t key);
    
    //PFermiDistributions(TString Target);
    ~PFermiDistributions();

    PDistribution* Clone(const char*delme=NULL) const;

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    
    
    
 private:
    int    version_flag;
    double FermiDisD( double *x, double *par ) const;
    double FermiDis3He( double *x, double *par ) const;
    double FermiDis4He( double *x, double *par ) const;
    double FermiDis7Li( double *x, double *par ) const;
    double FermiDis12C( double *x, double *par ) const;
    double FermiDis40Ca( double *x, double *par ) const;
    void SetDeuteronValues();
    
    //    TF1 *FM;
    double MeV2GeV;

    double *DeutC; 
    double *DeutD; 
    double *DeutM2; 
    
 public:
    
    ClassDef(PFermiDistributions,0) //Helper class for A-Fermi

};

#endif
