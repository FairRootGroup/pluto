// Author: I. Froehlich
// Written: 18.9.2008
// Revised: 

#ifndef _PDELTADALITZKRIVORUCHENKO_H_
#define _PDELTADALITZKRIVORUCHENKO_H_

#include "PDalitzDecay.h"

class PDeltaDalitzKrivoruchenko : public PDalitzDecay  {
  
 public:

    Bool_t FreezeOut(void);

    PDeltaDalitzKrivoruchenko(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    double dGdM(const int& id, const double& m, const double& ecm);

    
    ClassDef(PDeltaDalitzKrivoruchenko,0)  //Another Delta Dalitz description
};

#endif
