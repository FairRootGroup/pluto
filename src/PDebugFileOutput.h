// Author: Ingo Froehlich
// Written: 14/05/2007
// Modified: 
// PDebugFileOutput Class Header

#ifndef _PFILEDEBUGOUTPUT_H_
#define _PFILEDEBUGOUTPUT_H_

#include "PFileOutput.h"

class PDebugFileOutput: public PFileOutput {

private:

 protected:
    FILE *fp;
    
 public:
  PDebugFileOutput();
  ~PDebugFileOutput();

  bool OpenFile(char * filename);      //filename
  bool CloseFile(void);                //
  bool WriteEvent(void);               //next event
  bool WriteParticle(PParticle *par);  //write one particle


ClassDef(PDebugFileOutput,0) // Pluto debug file output

};

#endif 

















