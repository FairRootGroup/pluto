// Author: I. Froehlich
// Written: 1.08.2009
// Revised: 
// PArray
// 1 or 2dim array from file(s)

#ifndef _PARRAY_H_
#define _PARRAY_H_

#include "TObject.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TArrayD.h"




#define PARRAY_GRANULARITY 10
#define PARRAY_MAX_COLUMNS 10

class PArray : public TObject {

 private:

    Int_t dim;        //! Dimension of the function
    FILE* fp;         //! input file pointer

    TArrayD *vals_1d[PARRAY_MAX_COLUMNS];
    TArrayD *vals_2d;
    Int_t real_size_1d;
    Double_t scaling;
   

 public:

    //constructor
    PArray(Int_t dimension);
    ~PArray();
 
    void Scaling(Double_t sc) {scaling=sc;};
    Bool_t OpenFile(const char * filename, Int_t syntax, Int_t num_columns, Double_t y=0);
    TGraph*   GetTGraph(Int_t xcol, Int_t ycol);
    TGraph2D* GetTGraph2D(Int_t xcol, Int_t ycol);


    ClassDef(PArray,0)  //Tool class to convert txt-files to TGraph(s)
};

#endif
