// Author: M.A. Kagarlis
// Written: 23.02.99
// Revised: 29.06.00
// PFilter Class Header

#ifndef _PFILTER_H_
#define _PFILTER_H_

#include "TFormula.h"
class PReaction;

class PFilter : public TObject {

 public:

  PFilter(PReaction *, char *);
  // PFilter constructor by pointer to the PReaction where the filter should apply,

ClassDef(PFilter,0)//Pluto Filter Class
};

#endif // _PFILTER_H_

