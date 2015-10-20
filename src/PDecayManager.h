// PDecayManager Class Header
//
//  Author:   Volker Hejny
//  Written:  27.08.99
//  Revised:  21.06.00 MK
//  Revised:  21.03.05 R. Holzmann
//
// This header file contains 
//   template class PNextList
//   template class PStack
//   class PReactionList
//   class PDecayManager
//
// While all these classes are global, the first three originally
// planned as to be private to PDecayChannel. Due to some (maybe
// trivial) problems with nested classes using rootcint, they were
// put in a global scope.
// --------------------------------------------------------------------------

#ifndef _PDECAYLIST_H_
#define _PDECAYLIST_H_

#include "TObject.h"
class TTree;
class TClonesArray;
class PReaction;
class PChannel;
class PParticle;
class TPythia6;

#include "PDistributionManager.h"
#include "PReaction.h"
#include "PBulkInterface.h"
#include "PProjector.h"

// the following templates and the class PReactionList are
// private members of PDecayManager, but I due to problems using
// ClassImp() and ClassDef() with nested classes, they are
// not derived from TObject. The members are all public for better 
// use inside PDecayManager.

// template PNextList:
// implements a linked list of pointers(!) to object T. A delete does 
// not imply a delete of the object T itself. If this is wanted -> use
// member function Delete().
template<class T>
class PNextList {
 public:
    T*          Curr;	    			// Pointer to object T
    PNextList*  Next;				// Pointer to next list node
  
    PNextList()  : Curr(NULL), Next(NULL) {};
    PNextList(T *t) : Curr(t), Next(NULL) {};
    ~PNextList() { if (Next) delete Next; };

    void Add(T *t) {
	if (!Curr) Curr = t;			// if node still empty
	else if (!Next) Next = new PNextList<T>(t); // if node has no successor
	else Next->Add(t);				
    }
  
    void Delete() {				// delete all T objects
	delete Curr;	
	if (Next) Next->Delete();
    }
}; 
// end of template PNextList

// template PStack:
// implements a last-in-first-out standard stack via a single linked 
// list. Additionally the number of objects is counted.
template<class T>
class PStack {
 public:
    PNextList<T>*   top;				// Pointer to top node
    Int_t           Count;
  
    PStack() : top(NULL), Count(0) {};
    ~PStack() { delete top; };			// delete the nodes, not
                                                // the object itself
    void Push(T* t) {
	PNextList<T>* temp = new PNextList<T>(t);   // new node with object T
	temp->Next = top;                           // position before top node
	top = temp;                                 // make new node the top node
	Count++;
    }
  
    T* Pop() {
	if (!top) return NULL;
	T* t = top->Curr;				// this object is to be returned
	PNextList<T>* t2 = top;			// store top to delete
	top = top->Next;                            
	t2->Next = NULL;                            // unlink
	delete t2;                                   
	Count--;
	return t;
    }

    // Clone doubles the whole stack. Either an empty stack is
    // used (giving UseThis as an argument) or a new stack is
    // produced.  
    PStack<T>* Clone(PStack<T>* UseThis = NULL) {
	PStack<T>* t1;
	if (UseThis) t1 = UseThis;
	else  	 t1 = new PStack();
	if (!top) return t1;			// return an empty stack
	t1->Push(top->Curr);			// push top node object
	PNextList<T>* t2 = top;
	while (t2->Next) {				// add all node in t2 to t1
	    t2 = t2->Next;
	    t1->top->Add(t2->Curr);
	}
	t1->Count = Count;				// synchronize Count
	return t1;
    }
};
// end of template PStack

// class PReactionList:
// implements two connected stacks of PChannel to organize the
// already processed decay channels and the work to do. 
// The variable ReactionWeight contains the overall weight of
// this specific channel.
class PReactionList {
 public:
    PStack<PChannel>*	Finished;
    PStack<PChannel>*	ToDo;
    Double_t		ReactionWeight;
    Int_t			ID;
    static Int_t		maxID;
  
    PReactionList() {
	Finished       = new PStack<PChannel>;
	ToDo           = new PStack<PChannel>;
	ReactionWeight = 1.;
	ID             = ++maxID;
    };
    ~PReactionList() {
	delete Finished;
	delete ToDo;
    }
    PReactionList* Clone() {
	PReactionList* newRL = new PReactionList;
	Finished->Clone(newRL->Finished);
	ToDo->Clone(newRL->ToDo);
	newRL->ReactionWeight = ReactionWeight;
	return newRL;
    }
};
// end of class PReactionList

// class PDecayManager:
// PDecayManager manages the decay modes of all available particles
// in Pluto and builds all possible decay chains from one starting
// particle.
// The list of decay modes is empty by default. The user has to take
// care about the filling of these list.

class PDecayManager : public TObject {

 private:

    Int_t       verbose;		// verbose flag (0/1)
    Double_t    CurrentWeight;		// weight of currently active reaction
    Int_t       CurrentReactionNumber;	// serial number of currently active reaction
  
    PNextList<PParticle>   *UsedParticles;	// collection of used particles to delete them later on
    PNextList<PParticle*>  *UsedParticleArrays; // same for particle arrays
    PNextList<PChannel>	   *UsedChannels;       // same for used channels
						 
    PNextList<PReactionList> *ReactionList;       // collection of all possible decay branches
    PNextList<PReactionList> *CurrentReactionListPointer; // pointer to the currently active decay branch
    PChannel**	      ListForReaction;     // used for PReaction
    PReaction*	      CurrentReaction;     // pointer to currently used PReaction
    Int_t             NumberOfReactions;   // nb of reactions in ReactionList

    Int_t             decaychannel_param;  // param for all known particle decay modes
    Bool_t            fHGeant;             // set if PLUTO runs in HGeant
    void*             userSelection;       // selection function 
    Int_t             nTrigCond;           // trigger multiplicity

    Bool_t            fWriteIndex;         // write parent indices out, if set

    TPythia6          *fPythia;            // pointer to Pythia object
 

    // ContructPChannel is used internally to construct a new PChannel from
    // parent particle *p and the decay channel *c1 and store it in the
    // ToDo stack of PReactionList *RL.
    // If CopyFlag is set, the particle *p is copied before using. This
    // is necessary for the top channel.
    void ConstructPChannel(PParticle* p, PDecayChannel *c1, 
			   PReactionList* RL, Int_t CopyFlag=0);
  
    // utility function for PrintReactionList
    void PrintReactionListEntry(PReactionList*, ostream &os ) const;
    void PrintChain(PParticle* p, PChannel** l, Int_t c, ostream &os) const;

    // wrapper for PList.getName() to avoid error for composed 
    // particles
    using TObject::GetName;
    const char* GetName(Int_t id) const;
    Int_t maxFileSize;
    Float_t tauMax;

    PDistributionManager * pdist;

    int fileoutput_pos;
    PFileOutput *files[MAX_FILEOUTPUT];


    int bulkdecay_pos, pro_bulkdecay_pos;
    
    PBulkInterface *bulk[MAX_BULKDECAY];
    PBulkInterface *pro_bulk[MAX_BULKDECAY];
    PProjector *current_projector;



 public:
  
    PDecayManager();
    ~PDecayManager();

    void SetVerbose(Int_t v = 1);			// set verbose on/off 
  
    // add/assign/get a specific decay channel to one particle either by
    // particle id, by an instance of PParticle or just by name 
    void AddChannel(Int_t id, PDecayChannel* n);
    void AddChannel(PParticle* p, PDecayChannel* n);       
    void AddChannel(char* p, PDecayChannel* n);       
    PDecayChannel* GetChannel(Int_t id) const;
    PDecayChannel* GetChannel(PParticle* p)const;    
    PDecayChannel* GetChannel(char* n)const;    
  
    // set the default decay channel for this particle (according to PDG)
    void SetDefault(Int_t id,Int_t recursive=0);
    void SetDefault(PParticle* p,Int_t recursive=0);       
    void SetDefault(char* p,Int_t recursive=0);       

    // clear the channel
    void Clear(Int_t id);
    void Clear(PParticle* p);       
    void MyClear(char *);

    void Clear(const Option_t* delme=NULL) {
	if (delme) MyClear((char*) delme);
    };
    
       

    // by specifying one start particle InitReaction does look for all
    // possible decay chains and lists them in ReactionList. The initial
    // decay channel can be given either by a special decay channel or,
    // if it is a standard particle, the standard decay channel can be used.
    void InitReaction(PParticle *start, PDecayChannel *c1 = NULL);
  
    // GetNextReaction calls the constructor of PReaction for the next
    // reaction in ReactionList. The given parameters are forwarded
    // to the constructor. If wf is set, the weight stored in the reaction 
    // list is transfered to the initial particle; each reaction is simulated
    // with the same nb. of events, product particles have a proper weight set.
    // The second form is obsolete (and implies wf=0)
    PReaction* GetNextReaction(int wf, char* name, int f0=0, int f1=0, 
			       int f2=0, int f3=0, TTree *tt=NULL);
    PReaction* GetNextReaction(char* name, int f0=0, int f1=0, int f2=0, 
			       int f3=0, TTree *tt=NULL);
  
    // loop does a loop over all possible reaction channels and write
    // the output in a file name.root. The arguments are the same as
    // in GetNextReaction.
    Int_t Loop(int num, int wf, char* name, int f0=0, 
	       int f1=0, int f2=0, int f3=0, int rf=0);
    Int_t loop(int num, int wf, char* name, int f0=0, 
	       int f1=0, int f2=0, int f3=0, int rf=0) {
	return Loop(num, wf, name, f0,f1,f2,f3,rf);
    }

    // GetCurrentWeight provides the weight factor of the actual
    // reaction
    Double_t    GetCurrentWeight();
  
    // prints out the decay channel list
    void Print(const Option_t* delme=NULL) const {
	if (delme && strlen(delme)) {
	    MyPrint((char*) delme);
	}
	else MyPrint();
    };

    void MyPrint() const;

    // print out information for one specific particle
    void Print(Int_t id) const;
    void Print(PParticle *p) const;
    void MyPrint(char* name) const;

    void SetHGeant(Int_t fH) {fHGeant = fH;}
    void SetUserSelection(void *f) {userSelection = f;}
    void SetUserSelection(Int_t (*f)(PParticle*)) {userSelection = (void*)f;}
    void SetTrigCond(Int_t n) {nTrigCond = n;}



    // visualize the decay chains
    void PrintReactionList() const;
  
    void SetPythia(TPythia6 *p) {fPythia = p;}

    void SetMaxFileSize(Int_t bytes) {maxFileSize=bytes;}

    void SetDecayAll(Float_t tau=1.) {tauMax = tau;}

    void DisableHelicityAngle(void) {pdist->Disable("helicity_angles");};
    PDistributionManager * GetDistributionManager(void){return makeDistributionManager();};

    void SetWriteIndex(Bool_t flag) { fWriteIndex=flag; }

    Bool_t AddFileOutput(PFileOutput * file) {
	if (fileoutput_pos == MAX_FILEOUTPUT ) {
	    cout << "MAX_FILEOUTPUT reached" << endl;
	    return kFALSE;
	}
	files[fileoutput_pos++]=file;
	return kTRUE;
    }

    //These interfaces are analogue to PReaction:
    Bool_t AddBulk(PBulkInterface * mybulk);
    Bool_t AddPrologueBulk(PBulkInterface * mybulk); //Bulk IO before any decay

    Bool_t Do(char * command) {
	return GetCurrentProjector()->AddCommand(command);
    }
    Bool_t Do(TH1F * f, char * command) {
	return GetCurrentProjector()->AddHistogram(f,command);
    }
    Bool_t Do(TH2F * f, char * command) {
	return GetCurrentProjector()->AddHistogram(f,command);
    }
    Bool_t Output(TNtuple * f, char * command = (char *)"") {
	return GetCurrentProjector()->AddOutputTNtuple(f,command);
    }
    Bool_t Input(TNtuple * f) {
	return GetCurrentProjector()->AddInputTNtuple(f);
    }
    PProjector *GetCurrentProjector(void) {
	if (!bulkdecay_pos) {
	    current_projector = new PProjector();
	    AddBulk(current_projector);
	} else {
	    if (strcmp("PProjector",bulk[bulkdecay_pos-1] -> GetName())==0) {
		current_projector = (PProjector *)bulk[bulkdecay_pos-1];
	    } else {
		current_projector = new PProjector();
		AddBulk(current_projector);
	    }
	}


	return current_projector;
    }

    ClassDef(PDecayManager,0)//Pluto Decay Manager Class
	};
// end of class PDecayManager

#endif // _PDECAYLIST_H_









