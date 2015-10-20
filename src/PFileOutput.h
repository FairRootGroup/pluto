// Author: Ingo Froehlich
// Written: 14/05/2007
// Modified: 
// PFileOutput Class Header

#ifndef _PFILEOUTPUT_H_
#define _PFILEOUTPUT_H_

#include "PParticle.h"
#include "PBulkInterface.h"
#include "PChannel.h"

class PFileOutput: public PBulkInterface {

 private:

 protected:

     Int_t cnt;              //! number of particles in Event
     Int_t getVERTEX;        //! transport getVERTEX switch from PReaction
     Int_t writeINDEX;       //! transport writeINDEX switch from PReaction
     Int_t allPARTICLES;     //! transport allPARTICLES switch from PReaction
     Int_t asciiOUTPUT;      //! transport asciiOUTPUT switch from PReaction
     PChannel** channel;     //! transport channel array from PReaction
     Int_t      nChannel;    //! transport number of channels in array

     char filename_app[1024];
  public:
    PFileOutput();

    void SetHeader(Int_t my_cnt, 
		   Int_t my_allPARTICLES,
		   Int_t my_getVERTEX,
		   Int_t my_asciiOUTPUT,
		   Int_t my_writeINDEX,
		   PChannel** my_channel,
                   Int_t my_nCahnnel
		  ) {
	//number of expected WriteParticles:
	cnt=my_cnt;
	//PReaction flags:
	allPARTICLES = my_allPARTICLES;
	getVERTEX    = my_getVERTEX;
	asciiOUTPUT  = my_asciiOUTPUT;
	writeINDEX   = my_writeINDEX;
	channel      = my_channel;
        nChannel     = my_nCahnnel;
    };

    virtual bool OpenFile(char * filename);      //filename
    virtual bool CloseFile(void);                //
    virtual bool WriteEvent(void);               //next event
    virtual bool WriteParticle(PParticle *par);  //write one particle


    ClassDef(PFileOutput,0) // Pluto file output template
};
#endif 

















