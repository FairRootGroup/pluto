// Author: I. Froehlich
// Written: 13.8.2008
// Revised: 

#ifndef _PDISTRIBUTIONCOLLECTION_H_
#define _PDISTRIBUTIONCOLLECTION_H_

#define PDISTRIBUTIONCOLLECTION_MAX_PLUGIN 100

#include "PDistribution.h"


class PDistributionCollection : public PDistribution {
  
 public:

    PDistributionCollection();
    PDistributionCollection(const Char_t *id, const Char_t *de);
    PDistribution* Clone(const char*delme=NULL) const;

    virtual Bool_t Activate(void);
    const char * GetDependency(Int_t * pointer);
    

 protected:
  
    Bool_t RequiresPlugin(const char * name);

    Int_t  plugin_pointer;
    const char *plugin_name[PDISTRIBUTIONCOLLECTION_MAX_PLUGIN];


    ClassDef(PDistributionCollection,0)  // Base class for distribution collections (plugins)
};

#endif


