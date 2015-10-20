////////////////////////////////////////////////////////
//  File output interface implementation file
//
//  This class serves as ASCII file output class.
//  This IO class supports all standard ascii io formats
//  of PReaction plus additional formats exclusivly used
//  by HGeant of HADES. The output flags getVERTEX and writeINDEX
//  are transported from the settings of PReaction.
//  Below the coding of the header flag for the Geant IO are listed.
//
//  getVERTEX        -> write out the vertex coordinates of the particles
//  writeINDEX       -> write out parent index
//  writeSEQNUMBER   -> write out event sequence number for each particle (used for embedding)
//
//    //----------------------------------------------
//    // construct the header flag : normal flag + 10
//    Int_t shift = 0;
//    if(writeSEQNUMBER) shift = 10;
//    if(!getVERTEX && writeINDEX == 0) flag =   (2 + shift);
//    if(!getVERTEX && writeINDEX != 0) flag = - (2 + shift);
//    if( getVERTEX && writeINDEX == 0) flag =   (4 + shift);
//    if( getVERTEX && writeINDEX != 0) flag = - (4 + shift);
//    //----------------------------------------------
//
//  The above header flags are analyzed by HGeant during ascii file io
//  to set the input mode.
//
//
//----------------------------------------------
//  USAGE:
//
//  {
//      PReaction my_reaction(..settings...);      // create reaction
//                                                
//      ...                                        // add all channels , bulk etc
//
//      PHGeantOutput* output = new PHGeantOutput();  // created output
//      output->OpenFile("filename");                 // open file
//      my_reaction.AddFileOutput(output);            // connect output to raection
//
//      ...                                        // do all other actions
//
//      my_reaction.Loop(nEvents);                 // run event loop
//      output->CloseFile();                       // close output
//  }
//
//  Author:  Jochen Markert
//  Written: 05/12/2007
//  Revised:
//
////////////////////////////////////////////////////////

#include "PHGeantOutput.h"
#include "PReaction.h"
#include "PData.h"
#include "TString.h"



PHGeantOutput::PHGeantOutput() {
    // init variables
    asciiFile        = NULL;
    ctEvt            = 0;
    writeSEQNUMBER   = kTRUE;

    //globals from data base
    event_impact_param = makeStaticData()->GetBatchValue("_event_impact_param");
    seqnr              = makeStaticData()->GetBatchValue("_event_seqnr");

}

bool PHGeantOutput::OpenFile(char * filename)
{
    // open the ascii output file. returns kFALSE
    // if not successful.

    if(filename){
	TString name = filename;
	if(name.CompareTo("") != 0){
	    asciiFile = fopen(name.Data(),"w");
	    if(!asciiFile){
		Error("PHGeantOutput::OpenFile()","NULL pointer retrieved for output!");
                return kFALSE;
	    } else {
                Info("PHGeantOutput::OpenFile()","Opened ascii file IO = %s!",name.Data());
		return kTRUE;
	    }
	} else {
               Error("PHGeantOutput::OpenFile()","No File name set for output!");
               return kFALSE;
	}
    }
    return kFALSE;
}


bool PHGeantOutput::CloseFile(void) {
    // close ascii output file and put the
    // pointer to Null.
    if(asciiFile){
	fclose(asciiFile);
        asciiFile = NULL;
        return kTRUE;
    }

    return kFALSE;
}

bool PHGeantOutput::WriteEvent(void)
{
    // Write the event header. The io flag is coded
    // by the diffent switches. The blast parameter
    // is taken from the first channel (if it exists)
    // otherwise blast will be 0.

    Int_t flag    = 0;
    Float_t blast = 0;

    if(channel) {
        blast = channel[0]->GetBT();
    }

    //----------------------------------------------
    // construct the header flag : normal flag + 10
    Int_t shift = 0;
    if(writeSEQNUMBER) shift = 10;
    if(!getVERTEX && writeINDEX == 0) flag =   (2 + shift);
    if(!getVERTEX && writeINDEX != 0) flag = - (2 + shift);
    if( getVERTEX && writeINDEX == 0) flag =   (4 + shift);
    if( getVERTEX && writeINDEX != 0) flag = - (4 + shift);
    //----------------------------------------------


    if(asciiFile){
	fprintf(asciiFile," %i %i %f %f %i\n",ctEvt,cnt,blast,*event_impact_param,flag);

        ctEvt ++;
        return kTRUE;
    }


    return kFALSE;
}


bool PHGeantOutput::WriteParticle(PParticle *par)
{
    // write one particle to the ascii output file.
    // the format of the line depends on the settings of
    // difefrent switches.

    if(asciiFile)
    {
	if (!getVERTEX) {
	    if(writeINDEX == 0) {
		if(!writeSEQNUMBER){
		    fprintf(asciiFile," %e %e %e %e %i %i %i %e\n",                            //  8 vars
			    par->E() , par->Px(), par->Py(), par->Pz(),
			    par->ID(), par->GetSourceId()  , par->GetParentId(),
			    par->W());
		    return kTRUE;
		} else {
		    fprintf(asciiFile," %e %e %e %e %i %i %i %e %i\n",                         //  9 vars
			    par->E() , par->Px(), par->Py(), par->Pz(),
			    par->ID(), par->GetSourceId()  , par->GetParentId(),
			    par->W(), (Int_t) *seqnr );
                   return kTRUE;
		}
	    } else {
		if(!writeSEQNUMBER){
		    fprintf(asciiFile," %e %e %e %e %i %i %i %i %e\n",                         //  9 vars
			    par->E() , par->Px()         , par->Py()         , par->Pz(),
			    par->ID(), par->GetSourceId(), par->GetParentId(), par->GetParentIndex(),
			    par->W());
		    return kTRUE;
		} else {
		    fprintf(asciiFile," %e %e %e %e %i %i %i %i %e %i\n",                       // 10 vars
			    par->E() , par->Px()         , par->Py()         , par->Pz(),
			    par->ID(), par->GetSourceId(), par->GetParentId(), par->GetParentIndex(),
			    par->W(), (Int_t)*seqnr);
                    return kTRUE;
		}
	    }
	} else {

	    if (writeINDEX == 0) {
		if(!writeSEQNUMBER){
		    fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %e\n",                 // 12 vars
			    par->E()     , par->Px()         , par->Py()         , par->Pz(),
			    par->T()/300., par->X()          , par->Y()          , par->Z(),
			    par->ID()    , par->GetSourceId(), par->GetParentId(), par->W());
		    return kTRUE;
		} else {
		    fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %e %i\n",              // 13 vars
			    par->E()     , par->Px()         , par->Py()         , par->Pz(),
			    par->T()/300., par->X()          , par->Y()          , par->Z(),
			    par->ID()    , par->GetSourceId(), par->GetParentId(), par->W(),
                            (Int_t)*seqnr);
		    return kTRUE;
	    }
	    } else {
		if(!writeSEQNUMBER){
		    fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %i %e\n",              // 13 vars
			    par->E()     , par->Px()         , par->Py()         , par->Pz(),
			    par->T()/300., par->X()          , par->Y()          , par->Z(),
			    par->ID()    , par->GetSourceId(), par->GetParentId(),
			    par->GetParentIndex(), par->W());
		    return kTRUE;
		} else {
		    fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %i %e %i\n",           // 14 vars
			    par->E()     , par->Px()         , par->Py()         , par->Pz(),
			    par->T()/300., par->X()          , par->Y()          , par->Z(),
			    par->ID()    , par->GetSourceId(), par->GetParentId(),
			    par->GetParentIndex(), par->W(), (Int_t)*seqnr);
		    return kTRUE;
		}
	    }
	}
    }
    return kFALSE;
}


ClassImp(PHGeantOutput) 
