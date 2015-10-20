// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 
// 
// (copied in part from PData) IF

#ifndef _PSTDMODELS_H_
#define _PSTDMODELS_H_

#include "TROOT.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include <iostream>
#include "TObjArray.h"

#include "PDistributionManagerUtil.h"

#include "PChannelModel.h"
#include "PDalitzDecay.h"
#include "PAngularDistribution.h"
#include "PDeltaAngularDistribution.h"
#include "PDalitzDistribution.h"
#include "PScatterDistribution.h"
#include "PPiOmegaAngularDistribution.h"
#include "PChannel.h"
#include "PEEDirectDecay.h"
#include "PComplexBreitWigner.h"

using namespace std;

class PStdModels : public TObject {
    
 public:

    //constructor
    PStdModels();
    //destructor
    ~PStdModels();

    
    void Add(PDistributionManagerUtil * pdmutil);
    

 private:

    void GenericPhysics(PDistributionManagerUtil * pdmutil);

    TString *id,*de;
    PChannelModel *model;
    TObjArray * arr;
    TObjArray * GetModels(void);

    Bool_t generic_physics_done;

    ClassDef(PStdModels,0) //Pluto Std Models Class
};

#endif // _PSTDMODELS_H_







