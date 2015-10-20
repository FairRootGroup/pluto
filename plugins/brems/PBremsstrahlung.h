// Author: Wuestenfeld/Dohrmann
// Written: 2.2.2008
// Revised: 


#ifndef PBREMSSTRAHLUNG__H
#define PBREMSSTRAHLUNG__H
//Class definition

#include "PChannelModel.h"
#include "PAdaptiveMesh.h"
#include "TGraph2D.h"
#include "TGraph.h"

#define BREMS_KAPTARI_KAEMPFER 0
#define BREMS_SHYAM_MOSEL      1

class PBremsstrahlung : public PChannelModel
{
 public:

    PBremsstrahlung(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

//    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);
    Bool_t Init(void);
    using PChannelModel::SampleMass;
    Bool_t SampleMass(void);
    
    using PChannelModel::GetWeight;
    Double_t GetWeight(void);
    
    Double_t    Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper to use GetRandom of ROOT
    
    void SetSqrtS(Double_t s) {sqrt_s=s;};
    void SetMode(Char_t mode);
    void SetAuthor(Int_t a) {author=a;};       //0=KK, 1=SM
    void SetFunc(TGraph * gr) {graph=gr;};     //overwrites author by handmade array
    void SetFunc(TGraph2D * gr) {graph2d=gr;}; //overwrites author by handmade array
    
    //This are variables needed for the standalone Draw()
    void SetNeutron(Int_t n) {neutron_position=n;};
    void SetP2E(Double_t s) {p2_energy=s;};


 private:

    enum Mode {gNN=0,gDN,cSum,gN1520,FSI,VMD,cFsiVmd};
    Char_t model;

    Double_t mn, mp;
    Double_t sqrt_s,threshold;
    PParticle *dilepton,*parent;
    PParticle *p1,*p2,*p3,*p4;
    
    Int_t dilepton_position,neutron_position,author;
    Double_t p2_energy;
    PAdaptiveMesh * bin;

    Double_t    EvalSM(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;

    Double_t Pol1(Double_t z, Double_t p1, Double_t p2) const;
    Double_t Pol2(Double_t z, Double_t p1, Double_t p2, Double_t p3) const;
    //Double_t Pol3(Double_t z, Double_t p1, Double_t p2, Double_t p3, Double_t p4);*/
    Double_t AGauss(  Double_t x,   Double_t *par) const;
    Double_t DLines(Double_t x,   Double_t par0, Double_t par1,Double_t par2,Double_t par3) const;

    TGraph   *graph;
    TGraph2D *graph2d;

    ClassDef(PBremsstrahlung,0)  // pn/pp Bremsstrahlung

};

#endif
