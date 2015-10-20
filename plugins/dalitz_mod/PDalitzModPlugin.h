// Author: I. Froehlich
// Written: 17.9.2008
// Revised: 
// 


#ifndef _PDALITZMOD_PLUGIN_H_
#define _PDALITZMOD_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"
#include "PDeltaDalitzKrivoruchenko.h"

#include "PInclusiveModel.h"
//#include "PNNFSI.h"

using namespace std;

class PDalitzModPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PDalitzModPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PDalitzModPlugin();


    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    // PInclusiveModel *generators;

    Double_t static_br_thresh;   //Threshold for disabling static br in GeV
       
    PDeltaDalitzKrivoruchenko * kriv1, * kriv2;

    ClassDef(PDalitzModPlugin,0) //Plugin to modify Dalitz decays
};

#endif //_PDALITZMOD_PLUGIN_H_







