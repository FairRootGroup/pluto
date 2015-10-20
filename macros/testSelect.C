#ifndef __CINT__
#include "PParticle.h"
#endif

Int_t testSelect(PParticle* p) {  // test user selection function

  if (p==-1) {
    printf("= ok <<<\n");
    return -1;
  }

  if (p->Charge() > 0.) {  // +charged
        printf("testing %i has %f %f %f\n",p->ID(),p->Px(),p->Py(),p->Pz());
    return 1;
  }
  else return 0;

}




