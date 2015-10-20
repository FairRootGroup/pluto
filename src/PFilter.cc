////////////////////////////////////////////////////////////////////
//  PFilter Class 
//
//  Removed -> Features are covered now by the PProjector
//
//                                  Author:  M.A. Kagarlis
//                                  Written: 23.02.99
//                                  Revised: 29.06.00
////////////////////////////////////////////////////////////////////


#include "PReaction.h"
#include "PFilter.h"

ClassImp(PFilter)

PFilter:: PFilter(PReaction *r, char * expression) {

    Error("PFilter","The PFilter class has been removed. Use PProjector instead");

}

