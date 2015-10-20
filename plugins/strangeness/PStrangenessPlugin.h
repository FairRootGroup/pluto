// Author: I. Froehlich
// Written: 12.10.2008
// Revised: 
// 


#ifndef _PSTRANGENESS_PLUGIN_H_
#define _PSTRANGENESS_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

using namespace std;

class PStrangenessPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PStrangenessPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PStrangenessPlugin();


    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    Int_t is_initialized;

    ClassDef(PStrangenessPlugin,0) //Plugin to add new particles with strangeness and their decays
};

#endif //_PSTRANGENESS_PLUGIN_H_







