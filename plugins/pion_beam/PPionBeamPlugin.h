// Author:  Schuldes/Froehlich
// Written: Jul 2009


#ifndef _PPionBeamAmplitude_PLUGIN_H_
#define _PPionBeamAmplitude_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"   

#include "PPionBeamAmplitude.h"
#include "PPropagator.h"
#include "PInclusiveModel.h"


using namespace std;

class PPionBeamPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PPionBeamPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PPionBeamPlugin();


    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    PPionBeamAmplitude *Pi_minusBeamAmplitude, *Pi_plusBeamAmplitude; 
    PPropagator        *Rho0Propagator,*OmegaPropagator; 
    PInclusiveModel    *Pi_minusBeamAmplitude_gen, *Pi_plusBeamAmplitude_gen;
   
    ClassDef(PPionBeamPlugin,0) 
};

#endif //_PPionBeam_PLUGIN_H_
