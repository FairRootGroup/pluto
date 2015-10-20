// Author: I. Froehlich
// Written: 16.9.2009
// Revised: 
// 


#ifndef _PETADECAYSPLUGIN_H_
#define _PETADECAYSPLUGIN_H_

#define ETA_DOUBLE_DALITZ_BR 0.00005
#define ETA_EE_PIPI_BR       0.00001

#include "TROOT.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include <iostream>
#include "TObjArray.h"

#include "PChannelModel.h"
#include "PDalitzDecay.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

#include "PEtaDoubleDalitz.h"
#include "PEtaDoubleDalitzEnv.h"
#include "PEtaDoubleDalitzFF.h"
#include "PEtaPiPiGamma.h"


using namespace std;

class PEtaDecaysPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PEtaDecaysPlugin();
    PEtaDecaysPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PEtaDecaysPlugin();

    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    PEtaDoubleDalitz    *eta_dd_simple;
    PEtaDoubleDalitzFF  *eta_dd_ff;
    PEtaDoubleDalitzEnv *eta_dd_complex;

    PEtaPiPiGamma       *eta_pipi_gamma;



    ClassDef(PEtaDecaysPlugin,0) // A plugin for (rare) eta decays
};

#endif // _PETADECAYSPLUGIN_H_







