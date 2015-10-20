// Author: I. Froehlich
// Written: 7.05.2007
// Revised: 
// PMesh
// Linear Mesh

#ifndef _PMESH_H_
#define _PMESH_H_

#include "TObject.h"
#include "TF1.h"

class PMesh : public TF1 {

 private:

    Double_t *td;
    Double_t max, min;
    Int_t size;

 public:

    //constructor
    PMesh(Int_t size, const Char_t * name);
    ~PMesh();
 
    Int_t    GetSize(void){return size;}; //meshsize
    void     SetMax(Double_t pmax){max=pmax;fXmax=pmax;};
    void     SetMin(Double_t pmin){min=pmin;fXmin=pmin;};
    Double_t GetMax(void){return max;};
    Double_t GetMin(void){return min;};
    
    void     SetNode(Int_t node,Double_t v);
    Double_t GetNode(Int_t node);

    Double_t GetLinearIP(Double_t m) const ;
    void     Print(const Option_t* = NULL) const;
  
    Double_t Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    Double_t EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper
    
    ClassDef(PMesh,0)  //The linear mesh array
};

#endif
