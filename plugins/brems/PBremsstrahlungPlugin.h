// Author: I. Froehlich
// Written: 27.7.2008
// Revised: 
// 


#ifndef _PBREMSSTRAHLUNG_PLUGIN_H_
#define _PBREMSSTRAHLUNG_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

#include "PBremsstrahlung.h"
#include "PInclusiveModel.h"
#include "PNNFSI.h"
#include "PArray.h"

using namespace std;

class PBremsstrahlungPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PBremsstrahlungPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PBremsstrahlungPlugin();


    Bool_t ExecCommand(const char * command, Double_t value); 

    Bool_t Activate(void);

 private:

    PBremsstrahlung *pn,*pp,*pn_kk,*pp_kk;
    PBremsstrahlung *pn_sm_total,*pp_sm_total;     //Shyam sum
    PBremsstrahlung *pn_sm_elastic,*pp_sm_elastic; //Shyam elastic
    PBremsstrahlung *pn_sm_delta,*pp_sm_delta;     //Shyam delta
    PBremsstrahlung *pn_sm_n1520,*pp_sm_n1520;     //Shyam n1520

    PInclusiveModel *pn_gen, *pp_gen;
    PNNFSI          *pn_fsi, *pp_fsi, *np_fsi;

    PArray          *array_pp_sm_125, *array_pn_sm_125;

    Double_t kin_max;
    Int_t    current_author;

    ClassDef(PBremsstrahlungPlugin,0) //Plugin to activate KK or SM Bremsstrahlung
};

#endif //_PBREMSSTRAHLUNG_PLUGIN_H_







