////////////////////////////////////////////////////////
//  Author:  Jochen Markert
//  Written: 11/28/2007
//  Revised:
//
//  vertex information in a file. 
//  format : ntuple name "vertex" leafs "vX:vY:vZ".
//      or : ntuple name "vertex" leafs "vX:vY:vZ:seqNr".
//  Can be used e.g. for embedded particles. The
//  eventloop will be stopped latest at the last entry of the ntuple
//  no matter if a larger number was given to the loop.
////////////////////////////////////////////////////////

#include "PVertexFile.h"
#include "PChannel.h"
#include "PStaticData.h"

PVertexFile::PVertexFile()
{
    nEvent       = 0;
    input        = NULL;
    vertex_tuple = NULL; 

    vertex_x = makeStaticData()->GetBatchValue("_event_vertex_x");
    vertex_y = makeStaticData()->GetBatchValue("_event_vertex_y");
    vertex_z = makeStaticData()->GetBatchValue("_event_vertex_z");
    seqnr    = makeStaticData()->GetBatchValue("_event_seqnr");

}

PVertexFile::~PVertexFile()
{
    if(input) {
	input->Close();
	input        = NULL;
        vertex_tuple = NULL;
    }
}


Bool_t PVertexFile::OpenFile(TString inputfile) {
    // Inputfile : root ntuple file for the vertex point.
    // If no file is provided vertex will be 0,0,0.


    if(inputfile.CompareTo("") != 0 ){
	input = new TFile(inputfile.Data(),"READ");
	if(!input){
	    Error("OpenFile","Could not open root file for vertex input! Vertex will be 0,0,0!");
	}
        vertex_tuple = (TNtuple*)input->Get("vertex");
	if(!vertex_tuple){
	    Error("OpenFile","Could not read vertex ntuple from root file! Vertex will be 0,0,0!");
	    input -> Close();
            input = NULL;
	}

    }
    return kTRUE;

}

Bool_t PVertexFile::Modify(PParticle ** mstack, int *decay_done, int * num, int stacksize)
{
    // Read the particles from the defined stack and copy this to the official
    // particle stream

    Float_t xVertex,yVertex,zVertex;
    xVertex = yVertex = zVertex = 0.0;
    Int_t   seqNr=-1;

    //-----------------------------------------------------
    // reading vertex coordinates from ntuple for embedding
    if(vertex_tuple)
    {
	if(!vertex_tuple->GetEntry(nEvent)) {
	    return kFALSE;
	} else {
	    Float_t* vertex = vertex_tuple->GetArgs();

	    if(vertex_tuple->GetNvar() == 4) {
		xVertex   = vertex[0];
		yVertex   = vertex[1];
		zVertex   = vertex[2];
		seqNr = (Int_t) vertex[3];
		nEvent ++;
	    } else if(vertex_tuple->GetNvar() == 3) {
		xVertex   = vertex[0];
		yVertex   = vertex[1];
		zVertex   = vertex[2];
		seqNr = -1;
		nEvent ++;
	    } else {
		Error("Modify()","Ntuple has wrong number of variables!");
                return kFALSE;
	    }
	}
    }

    *vertex_x = xVertex;
    *vertex_y = yVertex;
    *vertex_z = zVertex;
    *seqnr    = (Double_t) seqNr;

    return kTRUE;
};

ClassImp(PVertexFile) 
