// Author: Jochen Markert
// Written: 11/28/2007
// Modified:
// PVertexFile Class Header

#ifndef _PVERTEXFILE_H_
#define _PVERTEXFILE_H_

#include "PBulkInterface.h"
#include "PStaticData.h"
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"

class PVertexFile: public PBulkInterface {

private:

    Long64_t nEvent;        // entry counter from read from ntuple
    TFile* input;           // root input for vertex ntuple
    TNtuple* vertex_tuple;  // ntuple which contains the vertex point

    Double_t *vertex_x,*vertex_y,*vertex_z, *seqnr;

protected:


public:

    PVertexFile();
    ~PVertexFile();

    Bool_t OpenFile(TString inputfile = "");
    Bool_t Modify(PParticle ** mstack, int *decay_done, int * num, int stacksize);

    ClassDef(PVertexFile,0) // Read vertex list from file
};
#endif


















