////////////////////////////////////////////////////////////////////
//  PDecayManager
//
//  Author:   Volker Hejny
//  Written:  Sep. 99
//  Revised:  21.03.05 RH (writing parent indices to ascii file implemented)
//  Revised:  10.12.04 RH (beam smearing/skewing implemented)
//  Revised:  28.08.03 RH (random decays implemented)
//  Revised:  04.07.03 RH (adjusted for multiple Root files: 1.9Gb limit!)
//  Revised:  15.12.00 RH (adjustments for Pythia decay routines)
//  Revised:  26.09.00 R. Holzmann (adjusted to gcc 2.95.2 20000220)
//  Revised:  08.08.00 R. Holzmann (weight treatment fixed)
//  Revised:  21.06.00 MK (adjustments for decay modes from PData)
//  Revised:  14.04.00 R. Holzmann (minor adjustments for 
//                                  PFireball and HGeant operation)
//  Revised:  27.07.07 IF: Changes for the new framework,
//                         Removed beam smearing as covered by the Distr.Manager
// 
//  The purpose of PDecayManager is to produce a set of possible
//  decay branches of a initial particle using a list of particle
//  decay modes. The list of decay modes is empty by default. The 
//  user has to take care about the filling of this list. This can
//  be done using the defaults supplied by this file (at end) or
//  by hand.
//
//  ----
//
//  The following templates and the class PReactionList are
//  private members of PDecayManager, but I due to problems using
//  ClassImp() and ClassDef() with nested classes, they are not
//  derived from TObject and not known inside ROOT. The members are 
//  all public for better use inside PDecayManager.
//  
//  template<class T> class PNextList
//    It implements a linked list of pointers to object T. All elements
//    of this templates are public. The destructor does not imply a
//    delete of the object T itself. If this is wanted -> use member
//    function Delete().
//
//    T 	 *Curr;		// pointer to object T
//    PNextList  *Next;		// pointer to next list node
//
//    PNextList();		// standard constructor
//    PNextList(T*);		// initalized with object T
//    ~PNextList();		// standard destructor
//				
//    void Add(T*);		// add object T to list
//    void Delete();            // delete the objects T in the list
//
//
//  template<class T> class PStack
//    It implements a last-in-first-out standard stack via a single linked
//    list. Additionally the number of objects is counted. All elements
//    are public.
//
//    PNextList<T>	*top;	// pointer to top node
//    Int_t		Count;	// number of elements in the stack
//
//    PStack();			// standarad constructor
//    ~PStack();		// the destructor deletes the nodes
//				// starting at top node, but not the
//				// objects T itself.
//
//    void Push(T*);
//    T *Pop();			// standard stack operators
//
//    PStack<T> *Clone(PStack<T> *UseThis);
//				// returns a pointer to a clone of the stack.
// 				// either an existing, empty stack is used 
//				// (giving UseThis as an optional argument)
//				// or a new stack is constructed internally.
//
//
//  class PReactionList
//    It implements two connected stacks of PChannel objects to organize 
//    the already processed decay channels and the work to do.
//    'Processed' means that all possible decay channels of the daughter
//    particles have been put on the ToDo stack.  
//    The variable ReactionWeight contains the overall weight of
//    this specific channel. All elements of this class are public.
//
//    PStack<PChannel> 	*Finished;	// channels already looked at
//    PStack<PChannel> 	*ToDo;		// decays to be done
//    Double_t		ReactionWeight;	// overall branching ratio
//    Int_t		ID;		// serial number of this object
//    static Int_t	maxID;		// maximum serial number
//  
//    PReactionList(); 
//    ~PReactionList(); 		// again: deletes only the stacks,
//					// not the objects in it.
//    PReactionList *Clone(); 		// doubles the whole PReactionList
//					// returning a pointer to the 
//					// new one.
//
// 
////////////////////////////////////////////////////////////////////
using namespace std;
//#include <strstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "TClonesArray.h"
#include "PData.h"
#include "PUtils.h"
#include "PReaction.h"
#include "PDecayChannel.h"
#include "PDecayManager.h"

ClassImp(PDecayManager)
    static char ret[200];
Int_t PReactionList::maxID = 0;



// --------------------------------------------------------------------------
PDecayManager::PDecayManager() {

    makeStaticData();

    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins

    Int_t key=makeDataBase()->GetEntry("std_set");
    if (key<0) cout << "PDecayManager: std_set not found" << endl;
    Int_t listkey=-1;
    while (makeDataBase()->MakeListIterator(key, "snpart","slink" , &listkey)) {
	PDecayChannel *ch=new PDecayChannel();
	//loop over all particles
	makeDataBase()->SetParamTObj (listkey, "decaychannel", ch);   
    }
    decaychannel_param=makeDataBase()->GetParamTObj("decaychannel");

    UsedParticles 	= new PNextList<PParticle>;
    UsedParticleArrays 	= new PNextList<PParticle*>;
    UsedChannels 		= new PNextList<PChannel>;
    ReactionList		= new PNextList<PReactionList>;
    CurrentReactionListPointer = NULL;
    verbose = 0;
    userSelection = NULL;
    nTrigCond = 0;
    fHGeant = kFALSE;
    fWriteIndex = kFALSE;
    fPythia = NULL;
    maxFileSize=0;
    tauMax=-1.;

    pdist =makeDistributionManager(); //create static object, if not yet done

    Info("PDecayManager()","(%s)", PRINT_CTOR);

    fileoutput_pos = 0;
    bulkdecay_pos = pro_bulkdecay_pos = 0;

}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayManager::~PDecayManager() {
    printf("decay manager destructor ok\n");

    UsedParticles->Delete(); 	delete UsedParticles;
    UsedParticleArrays->Delete(); delete UsedParticleArrays;
    UsedChannels->Delete(); 	delete UsedChannels;
    ReactionList->Delete();	delete ReactionList;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::SetVerbose(Int_t v) {
    // Setting verbose flag to optional parameter v. Omitting 'v' switches
    // on verbose option.
    verbose = v;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::AddChannel(Int_t id, PDecayChannel* n) {
    // Adds a specific PDecayChannel 'n' to the particle with ID 'id'. 
    if (!makeStaticData()->IsParticleValid(id)) { 
	Warning("AddChannel","id %i s unvalid",id);
	return;
    }
    else {
	Int_t key=makeDataBase()->GetEntryInt("pid",id);
	PDecayChannel*ch;
	TObject * o;
	if (key<0) return;
	if (!makeDataBase()->GetParamTObj (key, "decaychannel", &o)) return;
	ch=(PDecayChannel*)o;
	if (ch==NULL) {
	    makeDataBase()->SetParamTObj (key, "decaychannel", ch);   
	} else
	    ch->AddChannel(n);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::AddChannel(PParticle* p, PDecayChannel* n) {
    // Adds a specific PDecayChannel 'n' to the particle represented by a
    // pointer to a PParticle object. 
    if (p) AddChannel(p->ID(),n);  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::AddChannel(char* p, PDecayChannel* n) {
    // Adds a specific PDecayChannel 'n' to the particle represented just
    // by its name. 
    AddChannel(makeStaticData()->GetParticleID(p),n);  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayManager::GetChannel(Int_t id)const {
    // Returns a pointer to PDecayChannel for the particle with the ID 'id.'

    Int_t key=makeStaticData()->GetParticleKey(id); 
    if (key<0) {
	Warning("GetChannel","id %i not found in data base",id);
	return NULL;
    }
    TObject *ch;
    if (!makeDataBase()->GetParamTObj (key, decaychannel_param, &ch)) {
	Warning("GetChannel","no decaychannel found for ID: %i",id);
	return NULL;
    }
    return (PDecayChannel*)ch;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayManager::GetChannel(PParticle* p)const {
    // Returns a pointer to PDecayChannel for the particle represented by a
    // pointer to a PParticle object.

    if (p) return GetChannel(p->ID()); //TODO: This would be much faster if PParticles would know their key
    return NULL;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayManager::GetChannel(char* p) const{
    // Returns a pointer to PDecayChannel for the particle represented
    // by its name.

    if (p) return GetChannel(makeStaticData()->GetParticleID(p));
    return NULL;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::SetDefault(Int_t id,Int_t recursive) {
    // Sets the default decay branches for the particle with ID 'id'.
    if (!makeStaticData()->IsParticleValid(id)) return;
    
    Clear(id);
    PDecayChannel* c = GetChannel(id);
    if (!c) {
	Warning("SetDefault","Channel not present");
    }

    const char* name = GetName(id);
    if (verbose) cout << endl << "Setting defaults for particle #" << id << ", " 
		      << name << endl;
    
    if (verbose) {
	cout << "Information found:" << endl;
	makeStaticData()->PrintParticle(id); 
    }
  
    //now loop over decay modes
    Int_t key=makeDataBase()->GetEntryInt("pid", id);
    Int_t listkey=-1;
    Int_t tid[11];
    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	tid[0]=10; 
	makeStaticData()->GetDecayModeByKey(listkey,tid); // retrieve current mode info
	
	int *pos;
	makeDataBase()->GetParamInt (listkey, "didx" , &pos); //small workaround
	Double_t weight = makeStaticData()->GetDecayBR(*pos);        // static branching ratio
	c->AddChannel(weight, tid[0], &tid[1]);    // constructor by product pid-array pointer
	if (recursive) {
	    for (int i=0;i<tid[0];i++)
		SetDefault(tid[i+1]);
	}
    }

}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::SetDefault(PParticle* p,Int_t recursive) {
    // Sets the default decay branches for the particle represented by a
    // a pointer to a PParticle object.
    if (p) SetDefault(p->ID(),recursive);  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::SetDefault(char* p,Int_t recursive) {
    // Sets the default decay branches for the particle represented by its
    // name.
    SetDefault(makeStaticData()->GetParticleID(p),recursive);  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::Clear(Int_t id) {
    // Deletes the current list of decay branches for the particle with the
    // ID 'id'. This is necessary when one wants to introduce a new set
    // of branches to a particle. The SetDefault() functions are calling
    // this function automatically.
    
    Int_t key=makeDataBase()->GetEntryInt("pid",id); //TODO: make this faster
    if (key<0);
    TObject *ch;
    if (!makeDataBase()->GetParamTObj (key, "decaychannel", &ch)) return;
    if (ch==NULL) return;
    PDecayChannel *ch2=new PDecayChannel();
    if (!makeDataBase()->SetParamTObj (key, "decaychannel", ch2)) return;
    delete (PDecayChannel*)ch;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::Clear(PParticle* p) {
    // Deletes the current list of decay branches for the particle represented
    // by a pointer to a PParticle object.
    Clear(p->ID());  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::MyClear(char* p) {
    // Deletes the current list of decay branches for the particle represented
    // by its name.
    Clear(makeStaticData()->GetParticleID(p));  
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::InitReaction(PParticle *start, PDecayChannel *CurrentChannel) {
    // This function calculates the whole tree of decay branches of particle
    // 'start'. If there is a 'CurrentChannel' given, the initial branches
    // are fetched from this PDecayChannel. Otherwise the information is
    // used, which has been added by AddChannel() or SetDefault() calls.
    // Calling this function is mandatory as an initialization step
    // before producing any data.

    if (UsedParticles) {  UsedParticles->Delete();
    delete UsedParticles;
    UsedParticles = new PNextList<PParticle>;
    }
    if (UsedParticleArrays) { UsedParticleArrays->Delete();  
    delete UsedParticleArrays;
    UsedParticleArrays = new PNextList<PParticle*>;
    }
    if (UsedChannels) {   UsedChannels->Delete();
    delete UsedChannels;
    UsedChannels = new PNextList<PChannel>;
    }
    if (ReactionList) {   ReactionList->Delete();
    delete ReactionList;
    ReactionList = new PNextList<PReactionList>;
    }
    CurrentReactionNumber = 0;
    CurrentReactionListPointer = NULL;
  

    // build the initial channels to start with
  
    if (verbose) cout << endl << "processing decay:" << endl; 
  
    // if no channel given use the standard one.
    if (!CurrentChannel) CurrentChannel = GetChannel(start);
    
    Int_t copyFlag = ( start->IsFireball() || start->IsFileInput() )
	? 0 : 1; // don't want to make clones
    Int_t nReac1 = 0, nReac2 = 0;
    while (CurrentChannel) {
	//CurrentChannel->Print();
	//start->Print();
	PReactionList* tempRL = new PReactionList; 
	ConstructPChannel(start, CurrentChannel, tempRL, copyFlag);
	ReactionList->Add(tempRL);
	nReac1++;
	CurrentChannel = CurrentChannel->GetNext();
    }

    // ok, now we are processing ReactionList as long as there are some
    // channels to process
   
    Int_t MoreToDo = 1;
    while (MoreToDo) {
	MoreToDo = 0;
	PNextList<PReactionList>* x = ReactionList;
	// process the whole reaction list
	while (x) {	 
	    PReactionList* tempRL = x->Curr;
	    PChannel* tempPC = tempRL->ToDo->Pop();
	    // process the whole ToDo stack
	    while (tempPC) {
		MoreToDo = 1;
		tempRL->Finished->Push(tempPC);
 
		Int_t           i,loop;
		Int_t           NOP = tempPC->GetNumPar();	// number of products
		PParticle**     LOP = tempPC->GetParticles();   // list of products
		PDecayChannel** pdc = new PDecayChannel*[NOP];  // array of decay channels
		// one for each product
      

		for(i=0;i<NOP;i++) pdc[i] = GetChannel(LOP[i+1]);
      
		// like a multi digit counter, the decay channels for one particle
		// are processed via the Next entry one by one. If the last has been
		// processed, it is restored and loop is set to 1. Then the channels
		// of the next particle will be switched. The decay channel of the
		// first particle will always change. The loop has finished when the
		// decay channel of the last particle reaches its last entry.
		loop = 0;
		while (pdc[NOP-1]) {
		    // first lets clone the current reaction list
		    PReactionList* newRL = tempRL->Clone();
		    // collect the current decay channels
		    for (i=0;i<NOP;i++) {
			ConstructPChannel(LOP[i+1], pdc[i], tempRL);
			// if it is the first one or loop has been set -> Get the next entry
			if ( (!i) || (loop) ) {
			    loop = 0;
			    pdc[i] = pdc[i]->GetNext();
			    // if it is not the last one and there is no next channel -> reset
			    // and set loop to increment the next channel
			    if ( (i<(NOP-1)) && !pdc[i]) {
				pdc[i] =  GetChannel(LOP[i+1]);
				loop = 1;
			    }
			}
		    } // end of collecting channels
		    // if the last channel has not reached the end, store the newly
		    // created reaction list and switch to it, otherwise delete it.
		    if (pdc[NOP-1]) { ReactionList->Add(newRL);
		    nReac2++;
		    tempRL = newRL; }
		    else delete newRL;
		} // end of loop over all channel combinations
     
		// Get the next parent channel from the ToDo stack
		tempPC = tempRL->ToDo->Pop();

	    } // end of ToDo stack.
      
	    // Get next entry in reaction list.
	    x = x->Next;
	} // end of ReactionList    
    } // there is nothing left to do
    NumberOfReactions = nReac1 + nReac2;  // number of reactions in ReactionList 
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::ConstructPChannel(PParticle* p, PDecayChannel *c1,
				      PReactionList* RL, Int_t CopyFlag) {
    // This function is used internally inside InitReaction(). For one
    // given parent particle 'p' and PDecayChannel 'c1' it constructs
    // a PChannel and puts it in the ToDo stack of the PReactionList *RL.    

    //    cout << "***enter for  " << p->ID() << endl;

    // multiplying all the channel weights
    Double_t CurrentWeight = c1->GetWeight();
    if (CurrentWeight<1.0e-20) return;  
    RL->ReactionWeight *= CurrentWeight;
  
    // Get daughter particles, create a new array (for PChannel) and store
    // the original particle in it.
    Int_t  NumOfProducts;
    Int_t* ListOfProducts = c1->GetDaughters(NumOfProducts);
    PParticle** array = new PParticle*[NumOfProducts+1];
    UsedParticleArrays->Add(array);
    if (CopyFlag)     array[0] = p->Clone();
    else 		    array[0] = p;
//   if (verbose) {
//     cout << "(" << setw(2) << RL->ID << ", " << setw(8) << setprecision(2);
//     cout.setf(ios::scientific);
//     cout << RL->ReactionWeight << setprecision(0) << ") ";
// #warning: cout.setf(0, ios::scientific) commented for the moment!
//     cout.setf(ios::left, ios::scientific);
//     cout << GetName(array[0]->ID()) << " --> ";
//   }
  
    // Process the daughter particles.
    for (Int_t i=0; i<NumOfProducts; i++) {
	Int_t id = ListOfProducts[i];
	// Store this particle as a new one in the array
	array[i+1] = new PParticle(id);
	UsedParticles->Add(array[i+1]);
//    if (verbose) cout << GetName(id) << " ";

    } // out of daughter particles

    // again, create new channel with array and put in on the stack
    PChannel* tempPC = new PChannel(array, NumOfProducts);
    UsedChannels->Add(tempPC);
    RL->ToDo->Push(tempPC);
//  if (verbose) cout << endl;

    if (verbose) {
	cout << "Next channel for " << p->Name() << ", weight: " << setw(8) << setprecision(2);
	cout.setf(ios::scientific);
	cout << RL->ReactionWeight << setprecision(0) << endl;
	//cout << " is: " << tempPC->GetName() 
	//   << endl;
	tempPC->Print();
    }

    pdist->from_pdecaymanager = 1;
    pdist->Attach(tempPC); //Attach this channel to the known physics
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PReaction* PDecayManager::GetNextReaction(int wf, char* name,  int f0,
					  int f1, int f2, int f3,
					  TTree *tt) {
    // GetNextReaction() looks for the next decay tree in 'ReactionList'
    // and constructs a new PReaction. If the weight flag 'wf' is set, 
    // the channel weight is copied to the parent particle of the first
    // PChannel. This results in equal statistics of all channels. 
    // The options are passed to the appropriate PReaction constructor.
    //
    // Normally this routine is called internally from the loop() function.
  
    char* filename;
  
    if (verbose) cout << "Selecting next reaction:" << endl;
    if (ListForReaction) delete [] ListForReaction;
  
    if (CurrentReaction) delete CurrentReaction;
  
    if (!ReactionList) {
	if (verbose) cout << " No initialization done." << endl;
	return NULL;
    }

    if (!CurrentReactionListPointer) CurrentReactionListPointer = ReactionList;
    else CurrentReactionListPointer = CurrentReactionListPointer->Next;
  
    if (!CurrentReactionListPointer) {
	if (verbose) cout << " No reaction left." << endl;
	return NULL;
    }
    
    CurrentReactionNumber++;
    if (verbose) cout << endl << "This is reaction #" << CurrentReactionNumber << endl;

    filename = new char[1024];
    if (f3==2) sprintf(filename,"%s",name);  // all in 1 file 
    else sprintf(filename,"%s_%03d",name,CurrentReactionNumber);
  
    PReactionList* tempRL = CurrentReactionListPointer->Curr->Clone();
    Int_t NumOfChannels = tempRL->Finished->Count;
    ListForReaction = new PChannel*[NumOfChannels];
    for (Int_t i = NumOfChannels-1; i>=0 ; i--) {
	ListForReaction[i] = tempRL->Finished->Pop();
    }
    CurrentWeight = tempRL->ReactionWeight;
    if (verbose) cout << " Weight of reaction is " << CurrentWeight << endl;
    if (wf)
	(ListForReaction[0]->GetParticles())[0]->SetW(CurrentWeight);
    CurrentReaction = new PReaction(ListForReaction,filename,NumOfChannels,
				    f0,f1,f2,f3,tt);
    CurrentReaction->SetHGeant(fHGeant);
    CurrentReaction->SetWriteIndex(fWriteIndex);
    if (maxFileSize>0) CurrentReaction->SetMaxFileSize(maxFileSize);
    if (fPythia) CurrentReaction->SetPythia(fPythia);
    if (userSelection) {
	CurrentReaction->SetUserSelection(userSelection);
	CurrentReaction->SetTrigCond(nTrigCond);
    }
    if (tauMax>=0.) CurrentReaction->SetDecayAll(tauMax);

    //Copy bulk pointers to the PReaction
    for (int i=0;i<bulkdecay_pos;i++)
	CurrentReaction->AddBulk(bulk[i]);
    for (int i=0;i<pro_bulkdecay_pos;i++)
	CurrentReaction->AddPrologueBulk(pro_bulk[i]);



    delete tempRL;

    return CurrentReaction;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PReaction* PDecayManager::GetNextReaction(char* name,  int f0,
					  int f1, int f2, int f3,
					  TTree *tt) {
    // This call sets the weight flag to zero and is obsolete.
    return GetNextReaction(0,name,f0,f1,f2,f3,tt);
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
Int_t PDecayManager::Loop(int num, int wf, char* name, int f0,
			  int f1, int f2, int f3, int rf) {
    // This is the standard call to process all possible reaction channels
    // of the reaction set up before.
    // The total number of events given as the first argument are used
    // as a 'equivalent number'. It reflects the sum of weights to
    // get a proper normalization between the single channels. If the
    // weight flag 'wf' is set, the channel probability is also put
    // into this weight. This results in equal statistics for each channel.
    // The other parameters are the same as for the PReaction::loop() call.
    // The real number of events is returned.
  
    if (rf==1 && f3==1) f3=2;   // for random decays, we want common ascii file

    Int_t rStart, rStop, rEntries;
    char* rDescription = new char[1024];
    Double_t rWeight=0.;

    Int_t SumOfEvents = 0;
    Double_t WeightSum = 0;
    char* filename = new char[1024];
    sprintf(filename,"%s.root",name);
    TFile *f=NULL;
    TTree *td=NULL, *ti=NULL;  
    CurrentReaction = NULL;
    ListForReaction = NULL;
  


    if(!fHGeant) { // do not open Root output if run from HGeant
	f  = new TFile(filename,"RECREATE");
	f->SetCompressionLevel(1);
    }

    td = new TTree("data","Event data");

    sprintf(filename,"%s.evt",name);
    if (f3==2) PReaction::asciif=fopen(filename,"w");   // common ascii file
    PReaction::globalEventCounter = 0;

    if (rf==1) { // random decay mix of all reactions
	PReaction *pReactionArray[NumberOfReactions];
	double rWeights[NumberOfReactions];  

	Double_t wsave[NumberOfReactions];
	for (int i=0;i<NumberOfReactions;i++) {
	    CurrentReaction = NULL;
	    ListForReaction = NULL;
	    pReactionArray[i] = GetNextReaction(wf,name,f0,f1,f2,f3,td);
	    wsave[i] = GetCurrentWeight();
	    WeightSum += wsave[i];
	    rWeights[i] = WeightSum;    // cumulative weigth
	}


	for (int i=0;i<NumberOfReactions;i++) {
	    cout << endl << "Reaction " << i+1 << " -------------- with weight: "
		 << wsave[i] << endl;
	    pReactionArray[i]->Print();
	}

	for (int i=0;i<NumberOfReactions;i++) rWeights[i]/=WeightSum; // normalize

	cout << endl << "Do " << num << " events now..." << endl;

	for (int i=0;i<num;i++) { // do now num events randomly from reaction list
//	    cout << i << endl
	    double rand=PUtils::sampleFlat();
	    int ind = PUtils::FindIndex(NumberOfReactions,rWeights,rand);
	    pReactionArray[ind]->InitChannels();
//      pReactionArray[ind]->Print();
	    pReactionArray[ind]->DisableWeightReset();
	    pReactionArray[ind]->loop(1,0,0); //BUGBUG: global weight has to be set
	}
	SumOfEvents = num;
    }
    else {  // sequential treatment of reactions

	ti = new TTree("info","Run information");
	ti->Branch("rStart",&rStart,"rStart/I");
	ti->Branch("rStop",&rStop,"rStop/I");
	ti->Branch("rEntries",&rEntries,"rEntries/I");
	ti->Branch("rWeight",&rWeight,"rWeight/D");
	ti->Branch("rDescription",(void*)rDescription,"rDescription/C");

	PReaction *r=GetNextReaction(wf,name,f0,f1,f2,f3,td);
	while (r) {
//	std::ostringstream os(rDescription); //BUGBUG: Is this length-save?
//	std::ostringstream os(rDescription,1024);
//      ostrstream os(rDescription,1024);
//      PReactionList* RL = CurrentReactionListPointer->Curr->Clone();
//      PrintReactionListEntry(RL,os);
//      os << ends; //do not use it at the moment!!!
	    rStart = SumOfEvents + 1;  
	    if (verbose) r->Print();
	    rWeight = GetCurrentWeight();
	    WeightSum += rWeight;
	    if (verbose) cout << " Weight of this channel is " << rWeight << endl;
	    if (f3==1) PReaction::globalEventCounter=0;
	    if (wf) rEntries = r->loop(num,1);
	    else rEntries = r->loop((Int_t)(num*rWeight + 0.5),0);
	    SumOfEvents += rEntries;
	    rStop = SumOfEvents;

	    r=GetNextReaction(wf,name,f0,f1,f2,f3,td);
	    if (ti) ti->Fill();
	}
    }
    if (f) {
	f=td->GetCurrentFile();
	f->Write();
	f->Close();
	delete f;
    }
    if (f3==2) fclose(PReaction::asciif);
    delete [] filename;
    delete [] rDescription;
    if (verbose) cout << "Sum of weights is " << WeightSum << endl;
    return SumOfEvents;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
Double_t PDecayManager::GetCurrentWeight() {
    // Returns the weight of the current reaction channel.
    return CurrentWeight;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::MyPrint() const {
    // Prints the particle decay list.
    cout << "Particle Decay List:" << endl;
    cout << "--------------------" << endl;
//  for (Int_t i=0;i<PData::maxnumpar;Print(i++)); 
  
    Int_t key=makeDataBase()->GetEntry("std_set");
    if (key<0) cout << "PDecayManager: std_set not found" << endl;
    Int_t listkey=-1;
    while (makeDataBase()->MakeListIterator(key, "snpart","slink" , &listkey)) {
	//loop over all particles
	int id=makeStaticData()->GetParticleIDByKey(listkey);
	if (id<1000) Print(id); //Skip Quasi-Particles
    }
    PrintReactionList();
} 
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::Print(Int_t id) const  {
    // Prints the decay channels of the given particle (specified by its id).
    if (!makeStaticData()->IsParticleValid(id)) return;
    PDecayChannel* DC=GetChannel(id);
    if (!DC) return;

    
    if (DC->GetWeight()) {

	cout.setf(ios::left);
	cout.width(20);
	cout << GetName(id);
	DC->Print();
	cout << endl;
    }
    return;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::Print(PParticle* p) const  {
    // Prints the decay channels of the given particle (specified by a pointer
    // to a PParticle).
    if (p) Print(p->ID());
    return;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::MyPrint(char* p) const  {
    // Prints the decay channels of the given particle (specified by its name).
    if (p) Print(makeStaticData()->GetParticleID(p));
    return;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::PrintReactionList() const {
    // Prints the current reaction list
    PNextList<PReactionList>* x = ReactionList;
    Int_t Entry = 0;

    while(x) {
	if ( x->Curr ) {
	    PReactionList* RL = x->Curr->Clone();
	    Entry++;
	    cout << "-----" << endl;
	    cout << "Decay chain #" << Entry << ", probability " 
		 << RL->ReactionWeight << endl;
	    PrintReactionListEntry(RL,cout);
	    delete RL;
	}
	x = x->Next;
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::PrintReactionListEntry(PReactionList* RL, 
					   ostream &os) const{
    // Prints one entry of the reaction list (used internally).
    Int_t count = RL->Finished->Count;
    PChannel** temp = new PChannel*[count];
    for (Int_t i = count-1; i>=0 ; i--) {
	temp[i] = RL->Finished->Pop();
    }
    PrintChain((temp[0]->GetParticles())[0],temp,count,os);  
    delete [] temp; 
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayManager::PrintChain(PParticle *p, PChannel** l, Int_t c, 
			       ostream &os) const {
    // Prints one reaction chain (used internally).
    static Int_t indent = 0;
    indent += 2;
    for (Int_t i = 0; i<c; i++) {
	Int_t pcount      = l[i]->GetNumPar();
	PParticle** plist = l[i]->GetParticles();
	if (p==plist[0]) {
	    if (i==0) indent = 2;
	    os.width(indent) ; os << " ";
	    os.setf(ios::left);
	    os << setw(12) << GetName(p->ID()) 
	       << setw(0) << " --> ";
	    for (Int_t k=1; k<=pcount; k++) {
		os << GetName((plist[k])->ID()) << " ";
	    }
	    os << endl;
	    for (Int_t k=1; k<=pcount; k++) {
		PrintChain(plist[k],l,c,os);
	    }
	}
    }
    indent -= 2;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
const char*	PDecayManager::GetName(Int_t id) const {
    // Returns the name of a particle (specified by its id).
    if (id<=0) 	return "???";
    if (id<1000)	return (char*)makeStaticData()->GetParticleName(id);
    char temp[200];
    Int_t Upper = id/1000;
    Int_t Lower = id - 1000*Upper;
    sprintf(ret,"%s ",makeStaticData()->GetParticleName(Lower));
    id = Upper;
    while (id) {
	Upper = id/1000;
	Lower = id - 1000*Upper;
	sprintf(temp,"%s + %s",makeStaticData()->GetParticleName(Lower),ret);
	strcpy(ret,temp);
	id = Upper;
    }
    return ret;
}
// --------------------------------------------------------------------------


Bool_t PDecayManager::AddBulk(PBulkInterface * mybulk) {
    //Add a bulk interface to the list
    //Each bulk object will be executed during the event loop
    //after the normal decay
  
    if (bulkdecay_pos == MAX_BULKDECAY ) {
	Error("AddBulk","MAX_BULKDECAY reached");
	return kFALSE;
    }
    if (bulkdecay_pos && (bulk[bulkdecay_pos-1]->GetPriority() > mybulk->GetPriority())) {
	bulk[bulkdecay_pos]=bulk[bulkdecay_pos-1];
	bulk[bulkdecay_pos-1]=mybulk;
	bulkdecay_pos++;
    }
    else
	bulk[bulkdecay_pos++]=mybulk;

    return kTRUE;
}

Bool_t PDecayManager::AddPrologueBulk(PBulkInterface * mybulk) {
    //Add a bulk interface to the list
    //Each bulk object will be executed during the event loop
    //before the normal decay

    if (pro_bulkdecay_pos == MAX_BULKDECAY ) {
	Error("AddPrologueBulk","MAX_BULKDECAY reached");
	return kFALSE;
    }
    if (pro_bulkdecay_pos && (pro_bulk[pro_bulkdecay_pos-1]->GetPriority() > mybulk->GetPriority())) {
	pro_bulk[pro_bulkdecay_pos]=pro_bulk[pro_bulkdecay_pos-1];
	pro_bulk[pro_bulkdecay_pos-1]=mybulk;
	pro_bulkdecay_pos++;
    }
    else
	pro_bulk[pro_bulkdecay_pos++]=mybulk;

    return kTRUE;
}






