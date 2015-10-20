// Author: I. Froehlich
// Written: 3.7.2006
// Revised: 
// Base class for common distributions

#ifndef _PDISTRIBUTIONMANAGER_H_
#define _PDISTRIBUTIONMANAGER_H_

#include "TObjArray.h"
#include "TF1.h"
#include "TF2.h"
#include "PParticle.h"
#include "PDistribution.h"
#include "PDistributionCollection.h"
#include "PStdModels.h"
#include "PDistributionManagerUtil.h"
#include "PBatch.h"
#include "PProjector.h"
#include "PCommandList.h"

#define PDISTRIBUTIONMANAGER_MAX_COLLECT 100

class PDistributionManager;
//R__EXTERN PDistributionManager *gDM;
PDistributionManager& fDistributionManager();
PDistributionManager* makeDistributionManager();

class PDistributionManager : public TObject {
  
 public:

  PDistributionManager();

  int  Add(PDistribution * dist) {
      //Adds user-defined distribution
      return pdmutil->Add(dist);
  };

  int  Add(PDistribution * dist, const Char_t *gr) {
      //Adds user-defined distribution
      return pdmutil->Add(dist,gr);
  };

  void AlternativeTo(const char * a, const char *b) {
      pdmutil->AlternativeTo(a,b);
  };

  void Print(const Option_t* delme= "root") const {
      pdmutil->Print(delme);
  };

  int  Attach(PChannel * ch);
// {
      //Attach the PChannel to the known distributions
//      return pdmutil->Attach(ch);
//  };

  void Disable(const Char_t * id) {
      //Disable id, if possible
      pdmutil->Disable(id);
  };

  void Enable(const Char_t * id) {
      //Enable id, if possible
      pdmutil->Enable(id);
  };

  PDistribution * GetDistribution(const Char_t * id) {
      return pdmutil->GetDistribution(id);
  };

  void SetVerbose(Int_t v) {
      pdmutil->SetVerbose(v);
  };

  void DisableAddWarning(void) {
      makeDistributionManagerUtil()->no_warning=kTRUE;
  }

  void LinkDB(void) {
      // Link the channel models to data base
      pdmutil->LinkDB();
  };

  void ExpandGroup(const char * gr, Int_t ex=1) {
      pdmutil->ExpandGroup(gr,ex);
  }

  int from_pdecaymanager;

  Int_t   AddGroup(const Char_t *id, const Char_t *de) {
      return pdmutil->AddGroup(id,de);
  };

  void    SetGroup(const Char_t *id) {
      //set standard group for the next Add(PDistribution * dist)
      pdmutil->SetGroup(id);
  };

  Int_t   AddSubGroup(const Char_t *id, const Char_t *de, const Char_t *sub) {
      return pdmutil->AddSubGroup(id,de,sub);
  };

  Bool_t AddPlugin(PDistributionCollection *plugin);

  Bool_t Exec(const char * command); 
  Bool_t ExecAll(const char * command = "init");
  Bool_t Activate(const char * name);

  Bool_t Startup(void); //called before event loop
  Bool_t Startup(char * command); //add commands for startup

  Bool_t Unpack(const char * filename);  //unpack and execute command list

  PProjector *GetLoopFilter(void) {
      return loop_start;
  }

 private:
  

  PDistributionManagerUtil * pdmutil;
  PStdModels * std_models;

  PDistributionCollection *collect[PDISTRIBUTIONMANAGER_MAX_COLLECT];
  Int_t collect_pointer;

  void PluginInfo(const char * info);

  PBatch *batch;
  PProjector *loop_start;

   ClassDef(PDistributionManager,0)  //The main manager for model and distribution handling
};

#endif


