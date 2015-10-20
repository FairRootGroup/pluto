// Author: Ingo Froehlich
// Written: 14/02/2007
// Modified: 

#ifndef _PPROJECTOR_H_
#define _PPROJECTOR_H_

#include "PBulkInterface.h"
#include "PBatch.h"
#include "TH2.h"
#include "TH3.h"

#define PROJECTOR_MAX_BATCH 100
#include "TNtuple.h"

#define PROJECTOR_MAX_BRANCHES 50
#define PROJECTOR_MAX_STACK   100


class PProjector: public PBulkInterface {

 private:

    PBatch *batch[PROJECTOR_MAX_BATCH];
    TH3 * hist3[PROJECTOR_MAX_BATCH];
    TH2 * hist2[PROJECTOR_MAX_BATCH];
    TH1 * hist1[PROJECTOR_MAX_BATCH];
    Int_t fFillFlag[PROJECTOR_MAX_BATCH];
    int batch_pos, force_weight;

    int pid_param, link_param,batch_particle_param,batch_value_param,
	stream_default_pos_param,stream_max_pos_param;

    int key,key_x,key_y,key_z;
    Double_t *x,*y,*z,*w;

    Double_t *proj_nr;

    //for NTuple part:
    TNtuple *fp_out[PROJECTOR_MAX_BATCH],*fp_in[PROJECTOR_MAX_BATCH];
    Int_t key_list_out[PROJECTOR_MAX_BATCH][PROJECTOR_MAX_BRANCHES];
    Int_t key_pos_out[PROJECTOR_MAX_BATCH];
    Int_t key_list_in[PROJECTOR_MAX_BATCH][PROJECTOR_MAX_BRANCHES];
    Int_t key_pos_in[PROJECTOR_MAX_BATCH];
    Float_t values[PROJECTOR_MAX_BRANCHES];
    Float_t values_in[PROJECTOR_MAX_BATCH][PROJECTOR_MAX_BRANCHES];
    Int_t num_events_in[PROJECTOR_MAX_BATCH];
    Int_t num_events_in_c[PROJECTOR_MAX_BATCH];

    //for embedding particles
    Int_t stackpointer;
    PParticle stack[PROJECTOR_MAX_STACK];
    Int_t  SetParticles(PParticle ** mstack, int *decay_done, int * num, int stacksize, Int_t first_time);  //set the particles from the stream

 protected:
    
    
 public:
    
    PProjector();
    ~PProjector();


    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface
    
    Bool_t AddCommand(const char * command); //adds a command line to batch
    Bool_t Do(const char * command) {return AddCommand(command);};
    Bool_t AddHistogram(TH3 * histo, const char * command = "", 
			Int_t fillflag = 1); //Add a command and a histogram as a tool object
    Bool_t AddHistogram(TH2 * histo, const char * command = "", 
			Int_t fillflag = 1); //Add a command and a histogram as a tool object
    Bool_t AddHistogram(TH1 * histo, const char * command = "", 
			Int_t fillflag = 1); //Add a command and a histogram as a tool object

    Bool_t AddInputTNtuple(TNtuple * n, const char * command = "");
    Bool_t Input(TNtuple * n) {return AddInputTNtuple(n);};

    Bool_t AddOutputTNtuple(TNtuple * n, const char * command = "");

    void Print(const Option_t* delme=NULL) const ;

    ClassDef(PProjector,0) // Project particle array to histograms
};
#endif 

















