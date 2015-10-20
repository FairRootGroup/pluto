////////////////////////////////////////////////////////
//  File input interface implementation file
//
//  This class serves as template for file output classes
//
//                    Author:  Ingo Froehlich
//                    Written: 14/05/2007
//                    Revised:
//
////////////////////////////////////////////////////////

#include "PDebugFileOutput.h"



PDebugFileOutput::PDebugFileOutput() {
    fp = NULL;
}

PDebugFileOutput::~PDebugFileOutput() {
    if (fp) fclose(fp);
};

bool PDebugFileOutput::OpenFile(char * filename) {
    sprintf(filename_app,"%s.dbg",filename);
    fp = fopen(filename_app,"w");
    if (fp) return kTRUE;
    return kFALSE;
}


bool PDebugFileOutput::CloseFile(void) {
    if (fp) fclose(fp);
    fp=NULL;
    return kTRUE; 
}


bool PDebugFileOutput::WriteEvent(void) {
    //next event
    if (fp) {
	fprintf(fp,"next event\n");
	return kTRUE; 
    }
    return kFALSE; 
}


bool PDebugFileOutput::WriteParticle(PParticle *par) {
    //write one particle
    if (fp) {
	fprintf(fp,"name: %s\n",par->Name());
	fprintf(fp,"debug: %s\n",par->GetDebugString()->Data());
	
	return kTRUE; 
    }
    return kFALSE; 
}


ClassImp(PDebugFileOutput) 
