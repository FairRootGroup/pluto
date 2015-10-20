////////////////////////////////////////////////////////
//  Pluto bulk interface base class implementation file
//
//  This class serves as a base class for bulk modifications
//  It can be used for bulk decays (e.g. Pythia, Pluto)
//  or to add array of particles to the particle stream
// 
//                    Author:  Ingo Froehlich
//                    Written: 14/05/2007
//                    Revised:
//
////////////////////////////////////////////////////////

#include "PBulkInterface.h"



PBulkInterface::PBulkInterface() {
    
    bulk_id = (++gBulkCounter);
    fPriority = -1;
    
}

bool PBulkInterface::Modify(PParticle ** array, int *decay_done, int * num, int maxnum) {
    //modify the particle array "array"
    //input: particle array with "*num" members
    //new particles may be added (they are already instantiated) and *num increased
    //maxnum is the size of the particle array
    //setting the decay_done prevent further modifiers to decay this particle

    return kFALSE;
}

Int_t PBulkInterface::gBulkCounter=0;

ClassImp(PBulkInterface) 
