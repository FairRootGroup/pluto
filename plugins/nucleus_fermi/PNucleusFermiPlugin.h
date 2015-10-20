// Author:
// Written: 20.5.2009
// Revised: 
// 


#ifndef _PNUCLEUS_FERMI_PLUGIN_H_
#define _PNUCLEUS_FERMI_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

#include "PFermiMomentumGA.h"

using namespace std;

class PNucleusFermiPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PNucleusFermiPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PNucleusFermiPlugin();


    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:
    Bool_t gamma_active,proton_active;

    ClassDef(PNucleusFermiPlugin,0) //Plugin to activate fermi motion for various nuclei
};

#endif //_PNUCLEUS_FERMI_PLUGIN_H_







