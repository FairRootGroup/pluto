/////////////////////////////////////////////////////////////////////
//
// This is the base class for distribution collections
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PDistributionCollection.h"


ClassImp(PDistributionCollection)

PDistributionCollection::PDistributionCollection()  {
    Fatal("PDistributionCollection()","Wrong ctor");
} ;

PDistributionCollection::PDistributionCollection(const Char_t *id, const  Char_t *de) :
    PDistribution(id, de) {
    plugin_pointer=0;
} ;

PDistribution* PDistributionCollection::Clone(const char*delme) const {
    return new PDistributionCollection((const PDistributionCollection &)* this);
};

Bool_t PDistributionCollection::Activate(void) {
    return kFALSE;
}


Bool_t PDistributionCollection::RequiresPlugin(const char * name) {

    if (plugin_pointer == PDISTRIBUTIONCOLLECTION_MAX_PLUGIN) {
	Warning("RequiresPlugin","PDISTRIBUTIONCOLLECTION_MAX_PLUGIN reached");
	return kFALSE;
    }
    plugin_name[plugin_pointer]=name;
    plugin_pointer++;
    return kTRUE;

}

const char * PDistributionCollection::GetDependency(Int_t * pointer) {
    if (plugin_pointer == *pointer) return NULL;
    (*pointer)++;
    return plugin_name[(*pointer)-1];
}
