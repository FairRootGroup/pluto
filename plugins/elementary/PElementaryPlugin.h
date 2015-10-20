// Author: I. Froehlich
// Written: 21.8.2008
// Revised: 
// 


#ifndef _PTEISPLUGIN_H_
#define _PTEISPLUGIN_H_

#include "TROOT.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include <iostream>
#include "TObjArray.h"

#include "PChannelModel.h"
#include "PDalitzDecay.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

using namespace std;

class PElementaryPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PElementaryPlugin();
    PElementaryPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PElementaryPlugin();

    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    


    ClassDef(PElementaryPlugin,0) // A plugin for elemenaty collisions, cross sections etc...
};

#endif // _PELEMENTARYPLUGIN_H_







