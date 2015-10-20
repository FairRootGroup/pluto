// Author: I. Froehlich
// Written: 9.7.2007
// Revised: 
// PValues Class Header

#ifndef _PVALUES_H_
#define _PVALUES_H_

#include "TObject.h"

#define MAX_VALUES 10

#define T_MATRIX 1
#define U_MATRIX 2
#define TU_MATRIX 3
#define CHANNEL_POS 4

#define IS_BREAKUP 10
#define P_SCATTER 11

class PValues: public TObject {

 public:
    PValues();
    
    PValues(const PValues  & p);


    void Print(const Option_t* delme = NULL) const;
    bool SetValue(int id , double val);
    bool GetValue(int id , double * val);
    int  StringToValueID(char * st);


 protected:

    int array_id[MAX_VALUES];
    double array_val[MAX_VALUES];

    int pointer;

    ClassDef(PValues,0)  // User value container

};

#endif
