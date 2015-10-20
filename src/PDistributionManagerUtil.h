// Author: I. Froehlich
// Written: 26.7.2008
// Revised: 
// Util class for common distributions

#ifndef _PDISTRIBUTIONMANAGERUTIL_H_
#define _PDISTRIBUTIONMANAGERUTIL_H_



#include "TObjArray.h"
#include "TF1.h"
#include "TF2.h"
#include "PParticle.h"
#include "PDistribution.h"
#include "PDynamicData.h"
#include "PChannel.h"

#define MAX_DISTRIBUTION_LIST 1000
#define MAX_GROUP_LIST 100
#define MAX_CORR_LIST 10000

class PDistributionManagerUtil;

PDistributionManagerUtil& fDistributionManagerUtil();
PDistributionManagerUtil* makeDistributionManagerUtil();

class PDistributionManagerUtil : public TObject {
  
 public:

  PDistributionManagerUtil();
  int  Add(PDistribution * dist); //Adds user-defined distribution
  int  Add(PDistribution * dist, const Char_t *gr); //Adds user-defined distribution
  int  Add(TObjArray * arr);      //Adds a bulk of user-defined distribution
  void AlternativeTo(const char * a, const char *b);

  void Print(const Option_t* delme= "root") const;
  Int_t PrintGroup(Int_t group_id, Int_t width,Int_t indent,
		   const Char_t* name,
		   Int_t *num_enabled_mods,Int_t *num_total_mods,
		   Int_t *num_subs, Int_t *will_print)const;
  


  int  Attach(PChannel * ch);     //Attach the PChannel to the known distributions

  Bool_t Disable(const Char_t * id);     //Disable id, if possible
  Bool_t Enable(const Char_t * id);     //Enable id, if possible

  PDistribution * GetDistribution(const Char_t * id);
  void SetVerbose(Int_t v) {verbosity=v;};
  void LinkDB(void);             // Link the channel models to data base

  void ExpandGroup(const char *, Int_t ex=1);


  Int_t   AddGroup(const Char_t *id, const Char_t *de);
  void    SetGroup(const Char_t *id);  //set standard group for the next Add(PDistribution * dist)
  Int_t   AddSubGroup(const Char_t *id, const Char_t *de, const Char_t *sub);

  Bool_t no_warning;

 private:
  
   Int_t position;
   PDistribution * distribution[MAX_DISTRIBUTION_LIST];
   int alt_distribution[MAX_DISTRIBUTION_LIST];

  
   Int_t   group_position;
   const Char_t *group_identifier[MAX_GROUP_LIST], *group_description[MAX_GROUP_LIST];
   Int_t   group_expanded[MAX_GROUP_LIST];
   Int_t   group_corr[MAX_GROUP_LIST];
  
   Int_t   GetGroup(const Char_t *id, Int_t warning=1);  //get the group number from string
   //return a "none" (=0) when id not found
   Int_t current_group;
   Int_t AddCorrelation(Int_t gr, Int_t dis);  //correlates a certain group and distribution
   Int_t corr_position;
   Int_t corr_gr[MAX_CORR_LIST], corr_dis[MAX_CORR_LIST];
   void DisableAlts(int id);
   Int_t verbosity;
   Int_t linkdb_done; //flag for linkdb freeze out BUGBUG: Needs to be checked in Add()

   ClassDef(PDistributionManagerUtil,0)  //Util class for the distribution manager
};

#endif


