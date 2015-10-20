//////////////////////////////////////////////////////////////////////////
//  PReaction Class implementation file
//
//  A PReaction object represents a complete reaction process,
//  made up of a sequence of PChannels, namely a succession of
//  particle decays. All the participating PParticles and PChannels 
//  must be instantiated before the PReaction is set up.
//
//                                  Author:  M.A. Kagarlis
//                                  Written: 03.02.99
//                                  Revised: 02/11/2005 by R. Holzmann
//                                  Revised: 2007-2009 I. Froehlich
////////////////////////////////////////////////////////////////////////
//
//  A few important variables explained:
//
//  nchan    = number of channels in reaction
//  ndpar    = number of produced particles
//  ntpar    = number of stable (tracked) particles in final state
//  nclones  = number of particles in reaction, incl. beam and target
//
//  cindex[nclones] = array holding for a given particle m the index j+1
//                    of its channel: cindex[m] = j+1 (note the +1!)
//
//  dindex[nchan-1] = array holding for a given decay channel k
//                    the index m of the parent, i.e. decaying particle:
//                    dindex[k] = m
//
////////////////////////////////////////////////////////////////////////

const char *EType[4] = {
    "PReaction: ok",
    "PReaction: invalid option (must be 0 or 1)",
    "PReaction: invalid channel address",
    "PReaction: insufficient energy"
};

const char *RMessage[2] = {
    "PReaction: control external - user functions disabled",
    "PReaction: calculating widths in PData..."
};

const char *OType[2] = {
    "tracked particles on file",
    "all particles on file"
};

const char * percent="\%";

#include "Api.h"
#include "PData.h"
#include "PReaction.h"
#include "PUtils.h"
#include "TStopwatch.h"
#include "PFilter.h"

#ifdef USE_PYTHIA6
#include "TPythia6.h"
#endif



Int_t passEvent(float eb, float bp, float ph, int n,  // passes event to HGEANT
		int* id, int* src, int* par, int* parindex, float* px, float* py, float* pz, 
		float* vx, float* vy, float* vz, float* vt, float* w);

Int_t CalledSelection(PParticle* p);
Int_t CalledAnalysis(PParticle** p, Int_t n);

TMethodCall* gMethodCall1;  // global pointer needed for user selection
TMethodCall* gMethodCall2;  // global pointer needed for user analysis

FILE *PReaction::asciif=NULL;
Int_t PReaction::globalEventCounter=0;

TClonesArray* PReaction::evt=NULL;
Int_t PReaction::activeCnt=0;

ClassImp(PReaction)

// define option flags
    const unsigned int PReaction::ff0=1;
const unsigned int PReaction::ff1=1<<1;
const unsigned int PReaction::ff2=1<<2;
const unsigned int PReaction::ff3=1<<3;
const unsigned int PReaction::ff4=1<<4;

PReaction:: PReaction(PChannel **pchannel, char *file_name, int n,
		      int f0, int f1, int f2, int f3, TTree *ttree) {
    // Reaction constructor by pointer to array of pointers to Channels,
    // output file name (see below), options (see below), pointer to external event tree
    // (if called from the decay manager).
    // Options: 
    //   f0: 0=Tracked, 1=All particles in ROOT output file_name.root
    //   f1: unused
    //   f2: 0=Do not calculate, 1=Calculate production vertices for product particles
    //   f3: 0=Do not produce, 1=Produce GEANT ascii output file_name.evt (tracked particles always)
    //   _________________________________________________________________
    //   Format of the ASCII output in file_name.evt:
    //   EvSeqNb NumPart Ebeam bPar phiEv flag <== number of particles in event
    //   E1 px1  py1  pz1  ID1  weight1      <-- 3-momentum components in GeV/c,
    //   E2 px2  py2  pz2  ID2  weight2        <-- GEANT particle ID,
    //   E3 px3  py3  pz3  ID3  weight3          <-- weight
    //   .    .    .    .   .    .   .   .   .   .  .   . 
    //   _________________________________________________________________
    // See also PProjector/PBatch for additional conventions if filters are used.
    //   _________________________________________________________________

    makeStdData()->fillDataBase();//init physics, if not yet done
    makeDistributionManager()->DisableAddWarning();
    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins

    if (!pchannel) n=0;
    sub_reaction = NULL;

    status=0; tree=ttree;
    int check_options = (f0!=0&&f0!=1) + (f1!=0&&f1!=1)
	+ (f2!=0&&f2!=1) + (f3!=0&&f3!=1&&f3!=2);
    if (check_options) {
	status=1;
	printf("%s\n",EType[status]);
	return;
    } else {
	allPARTICLES = f0;
	resetCHANNELS = 1-f1;
	getVERTEX = f2;
	asciiOUTPUT = f3;
	extTREE = tree!=NULL;
	ropt = 0 | ff0*allPARTICLES | ff1*resetCHANNELS | ff2*getVERTEX 
	    | ff3*(asciiOUTPUT>0)  | ff4*extTREE;
    }
    nchan=n;
    SetUp(pchannel);
    SetName(file_name);
}


PReaction:: PReaction() {
    makeStdData()->fillDataBase();//init physics, if not yet done
    makeDistributionManager()->DisableAddWarning();
    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins
    sub_reaction = NULL;

    status=0;
    tree=NULL;
    allPARTICLES = 1;
    resetCHANNELS = 0;
    getVERTEX = 1;
    asciiOUTPUT = 0;
    extTREE = 0;
    ropt = 0;
    nchan=0;
    SetUp(NULL);
    SetName(NULL);
}


PReaction:: PReaction(char *filename) {
    makeStdData()->fillDataBase();//init physics, if not yet done
    makeDistributionManager()->DisableAddWarning();
    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins
    sub_reaction = NULL;

    status=0;
    tree=NULL;
    allPARTICLES = 1;
    resetCHANNELS = 0;
    getVERTEX = 1;
    asciiOUTPUT = 0;
    extTREE = 0;
    ropt = 0;
    nchan=0;
    SetUp(NULL);
    SetName(filename);
};


PReaction:: PReaction(PChannel **pchannel, int n, unsigned int ff, 
		      TTree *ttree, char *file_name) {
    // same as above, but passing one unsigned int as flag
    // position of arguments shuffled to break ambiguity with previous constructor

    makeStdData()->fillDataBase();//init physics, if not yet done
    sub_reaction = NULL;

    status=0;
    tree=ttree;
    allPARTICLES = ((ff&ff0)==ff0);
    resetCHANNELS = 1-((ff&ff1)==ff1);
    getVERTEX = ((ff&ff2)==ff2);
    asciiOUTPUT = ((ff&ff3)==ff3);
    extTREE = (tree!=NULL);
    ropt = ff&(ff4*extTREE);
    nchan=n;
    SetUp(pchannel);
    SetName(file_name);
}




PReaction:: PReaction(Double_t momentum,
                      char* beam, char* target,
                      char* reaction, char* file_name,
                      Int_t f0, Int_t f1, Int_t f2, Int_t f3, TTree *ttree) {
    // build reaction completely from descriptor string
    
    char dummy[100];
    sprintf(dummy,"_P1=%lf",momentum);
    parse_script(dummy,
		 beam, target,
		 reaction, file_name,
		 f0, f1,f2, f3, ttree);
} 

PReaction:: PReaction(char* e,
                      char* beam, char* target,
                      char* reaction, char* file_name,
                      Int_t f0, Int_t f1, Int_t f2, Int_t f3, TTree *ttree) {
    // build reaction completely from descriptor string

    parse_script(e,
		 beam, target,
		 reaction, file_name,
		 f0, f1,f2, f3, ttree);
} 


Bool_t PReaction:: parse_script (char * command,
                      char* beam, char* target,
                      char* reaction, char* file_name,
                      Int_t f0, Int_t f1, Int_t f2, Int_t f3, TTree *ttree) {
    // Build "reaction" completely from descriptor string  (author: V. Hejny)
    // For "command, the batch script syntax is supported
    // to set the beam parameters. The following variables
    // can be used:
    // _T1, _T2, _P1, _P2 as kinetic beam energy/momentum 
    // for the beam/target (in GeV), and
    // _theta1, _theta2, _phi as beam inclination (in rad). 
    // It should be mentioned, that for the target, the
    // notation is opposite, i.e. theta2 has the opposite
    // rotation as theta1, and T2 goes into the opposite 
    // direction as T1
    //
    // If no variable is used, only T1 as the beam energy is
    // assumed
    //
    // Examples:
    // "2.2"               -> Beam energy is 2.2 GeV
    // "_T1=1.5; _T2=1.2"  -> Collider experiment
    // "_T1=1.5; _theta=2*TMath::DegToRad();"
    //                     -> Beam inclination of 2 deg

  
    r_beam=beam,r_target=target;
    reaction_string=reaction;
    makeDistributionManager()->DisableAddWarning();
    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins
    sub_reaction = NULL;

    PUtils::remove_spaces(&command);
    if (strlen(command)>0 && isdigit(command[0])) {
	char * newe = new char[strlen(command)+5];
	sprintf(newe,"_T1=%s",(char *)command);
	command=newe;
    }

    makeGlobalBatch()->SetVarList((char *)
				  "_T1;_T2;_P1;_P2;_theta1;_theta2;_phi;");
    makeGlobalBatch()->Execute(command);
    makeGlobalBatch()->SetVarList(NULL);
    
    Double_t beam_energy1=0.;
    if (makeStaticData()->GetBatchValue("_T1",0)) {
	beam_energy1=*(makeStaticData()->GetBatchValue("_T1",0));
    }
    Double_t beam_energy2=-0.;
    if (makeStaticData()->GetBatchValue("_T2",0)) {
	beam_energy2=*(makeStaticData()->GetBatchValue("_T2",0));
    }
    Double_t beam_momentum1=-1.;
    if (makeStaticData()->GetBatchValue("_P1",0)) {
	beam_momentum1=*(makeStaticData()->GetBatchValue("_P1",0));
    }
    Double_t beam_momentum2=-1.;
    if (makeStaticData()->GetBatchValue("_P2",0)) {
	beam_momentum2=*(makeStaticData()->GetBatchValue("_P2",0));
    }
    Double_t beam_theta1=0.;
    if (makeStaticData()->GetBatchValue("_theta1",0)) {
	beam_theta1=*(makeStaticData()->GetBatchValue("_theta1",0));
    }
    Double_t beam_theta2=0.;
    if (makeStaticData()->GetBatchValue("_theta2",0)) {
	beam_theta2=*(makeStaticData()->GetBatchValue("_theta2",0));
    }
    Double_t beam_phi1=0.;
    if (makeStaticData()->GetBatchValue("_phi",0)) {
	beam_phi1=*(makeStaticData()->GetBatchValue("_phi",0));
    }
    // Double_t beam_phi2=0.;
    // if (makeStaticData()->GetBatchValue("_phi2",0)) {
    //   beam_phi2=*(makeStaticData()->GetBatchValue("_phi2",0));
    // }

    char *array[200];
    Int_t array_s=200; 
    PUtils::Tokenize(reaction, ";", array, &array_s);

    Int_t total_channels = 0;
    TList plutoList;
   
    for (int i=0; i<array_s; i++) {
       
	PParticle pb,pt;

	Int_t n = 0;
	if (beam_momentum1>0) {
	    pb = PParticle (beam, 0,0,beam_momentum1);  
	} else {
	    pb = PParticle (beam, beam_energy1);  
	}
	if (beam_momentum2>0) {
	    pt = PParticle (target, 0,0,-beam_momentum2);  
	} else {
	    pt = PParticle (target, 0,0,-sqrt(beam_energy2*beam_energy2+2*beam_energy2*
					      makeStaticData()->GetParticleMass(target)));  
	}
	pb.RotateY(beam_theta1);
	pb.RotateZ(beam_phi1);
	pt.RotateY(-beam_theta2);
	pt.RotateZ(beam_phi1);

	cout << "<Beam>" << endl;
	pb.Print();	
	cout << "<Target>" << endl;
	pt.Print();

	PParticle *q = new PParticle(pb+pt);

	plutoList.AddLast(q);
	ParseChannel(q,array[i] , plutoList, n);

	total_channels+=n;
    }

    PChannel **pchannel = new PChannel*[total_channels];
    TIter next(&plutoList);
    
    Int_t pos = 0;
    while (TObject *t = next()) {
	if (t->IsA() == PChannel::Class()) { 
	    pchannel[pos++] = (PChannel*) t;
	}
    }


    status=0;
    tree=ttree;
    int check_options = (f0!=0&&f0!=1) + (f1!=0&&f1!=1)
	+ (f2!=0&&f2!=1) + (f3!=0&&f3!=1&&f3!=2);
    if (check_options) {
	status=1;
	printf("%s\n",EType[status]);
	return kFALSE;
    } else {
	allPARTICLES = f0;
	resetCHANNELS = 1-f1;
	getVERTEX = f2;
	asciiOUTPUT = f3;
	extTREE = tree!=NULL;
	ropt = 0 | ff0*allPARTICLES | ff1*resetCHANNELS | ff2*getVERTEX
	    | ff3*(asciiOUTPUT>0)  | ff4*extTREE;
    }
    nchan=total_channels;
    SetUp(pchannel);
    SetName(file_name); 
    return kTRUE;
} 


PReaction:: PReaction(PParticle * q,
                      char* reaction, char* file_name,
                      Int_t f0, Int_t f1, Int_t f2, Int_t f3, TTree *ttree) {
    // build reaction completely from descriptor string  (author: V. Hejny)

    reaction_string=reaction;
    makeDistributionManager()->DisableAddWarning();
    makeDistributionManager()->ExecAll("init"); //init physics, if not yet done, and allow for didx-plugins
    sub_reaction = NULL;

    char *array[200];
    Int_t array_s=200; 
    PUtils::Tokenize(reaction, ";", array, &array_s);

    Int_t total_channels = 0;
    TList plutoList;
   
    for (int i=0; i<array_s; i++) {
       
      Int_t n = 0;

      plutoList.AddLast(q);
      ParseChannel(q,array[i] , plutoList, n);

      total_channels+=n;
    }

    PChannel **pchannel = new PChannel*[total_channels];
    TIter next(&plutoList);
    
    Int_t pos = 0;
    while (TObject *t = next()) {
      if (t->IsA() == PChannel::Class()) { 
	pchannel[pos++] = (PChannel*) t;
      }
    }


    status=0;
    tree=ttree;
    int check_options = (f0!=0&&f0!=1) + (f1!=0&&f1!=1)
	+ (f2!=0&&f2!=1) + (f3!=0&&f3!=1&&f3!=2);
    if (check_options) {
	status=1;
	printf("%s\n",EType[status]);
	return;
    } else {
	allPARTICLES = f0;
	resetCHANNELS = 1-f1;
	getVERTEX = f2;
	asciiOUTPUT = f3;
	extTREE = tree!=NULL;
	ropt = 0 | ff0*allPARTICLES | ff1*resetCHANNELS | ff2*getVERTEX
	    | ff3*(asciiOUTPUT>0)  | ff4*extTREE;
    }
    nchan=total_channels;
    SetUp(pchannel);
    SetName(file_name); 
} 


void PReaction::AddReaction(char* reaction) {
    if (sub_reaction) sub_reaction->AddReaction(reaction);
    else {

	if (original_filename=="")
	    ConvertFilename();

	if (strlen(filename) > 0)
	  sub_reaction = new PReaction((char*)"",r_beam,r_target,
					 reaction,(char*)original_filename.Data(),
					 allPARTICLES,1-resetCHANNELS,
					 getVERTEX,asciiOUTPUT);	
	else
	    sub_reaction = new PReaction((char*)"",r_beam,r_target,
					 reaction,NULL,
					 allPARTICLES,1-resetCHANNELS,
					 getVERTEX,asciiOUTPUT);
	sub_reaction->ConvertFilename();
      
    }
}

void PReaction::ConvertFilename(void) {
    if (strlen(filename) == 0) return;
    TString original_filename2=filename;

    filename.Append("_to_");
    filename.Append(reaction_string);
    filename.ReplaceAll(" [ ",".");
    filename.ReplaceAll(" ] ","");
    filename.ReplaceAll(" [",".");
    filename.ReplaceAll(" ]","");
    filename.ReplaceAll("[ ",".");
    filename.ReplaceAll("] ","");

    filename.ReplaceAll("[",".");
    filename.ReplaceAll("]","");
    filename.ReplaceAll(" ","_");

    //spectator
    filename.ReplaceAll("(","_");
    filename.ReplaceAll(")","");

    SetName((char*)filename.Data());
    original_filename=original_filename2;
    
    cout << "New name: " << filename << endl;
}

PReaction:: ~PReaction() {
    // Reaction destructor

    cindex.~TArrayI(); 
    if (particle_stack)
	delete [] particle_stack;

    if (evt) {
	delete evt;
	evt = NULL;
    }
    if (dindex.GetArray()) dindex.~TArrayI();
    if (ftrack.GetArray()) ftrack.~TArrayI();
    if (!extTREE&&tree) delete tree;

    if (rootfile && !extTREE && (strlen(filename) > 0)) {
	rootfile=tree->GetCurrentFile(); // needed if a new file was opened by Root
	rootfile->Write();
	rootfile->Close();
	delete rootfile;
    }

    if (gMethodCall1) delete gMethodCall1;
    if (gMethodCall2) delete gMethodCall2;
}


PParticle *PReaction::MakeParticle(char * name) {
    //    cout << name << endl;
    int is_spec =-1;
    int len = strlen(name)-1;
    if ((name[0]=='(') && (name[len]==')'))  { 
	is_spec = 1;
	PUtils::remove_brackets(&name,'(',')');
    }
    if (name[0]=='(') { 
	is_spec = 2;
	name++;
    }
    if (name[len]==')') { 
	is_spec = 3;
	name[len]='\0';
    }

    PParticle *dummy=new PParticle(name);
    dummy->SetSpectator(is_spec);
    return dummy;
};

Int_t PReaction:: ParseChannel(PParticle *parent, const char* channel,
			       TList &plutoList, Int_t &numChannels) {
    // parse reaction descriptor string  (author: V. Hejny)
    TList plist;
    Int_t num_parts = 0;
    PParticle *current = NULL;
    Int_t i=0;
    Int_t start = i;
    while(1) {
	switch (channel[i]) {
	    case '[':
		if (!current) {
		    while(1) {
			if (channel[i]==']') { i++; break; }
			else if (channel[i]==0) break;
			else i++;
		    }
		    break;
		}
		i++;
		i += ParseChannel(current, &channel[i], plutoList, numChannels);
		if (channel[i]) i++;
		start = i;
		break; 
	    case 0:
	    case ']': {
		if (start!=i) {
		    TString name(&channel[start], i-start);
//		    current = new PParticle( (char*) name.Data());
		    current = MakeParticle( (char*) name.Data());
		    plist.AddLast(current);
		    plutoList.AddLast(current);
		    num_parts++;
		}
		if (!num_parts) return i;
		PParticle **ppl = new PParticle*[num_parts+1];
		ppl[0] = parent;
		Int_t pos = 1;
		TIter next(&plist);
		while (PParticle *p = (PParticle*) next()) {
		    ppl[pos++] = p;
		}
		PChannel *ch = new PChannel(ppl,num_parts);
		plutoList.AddFirst(ch);
		numChannels++;
		return i;
		break;
	    } 
	    case ' ': {
		if (start==i) { start = ++i; break; }
		TString name(&channel[start], i-start);
//		current = new PParticle( (char*) name.Data());
		current = MakeParticle( (char*) name.Data());
		plist.AddLast(current);
		plutoList.AddLast(current);
		num_parts++;
		i++;
		start = i;
		break;
	    }
	    default:
		i++;
	}
    }
}

//************Interfaces ***************


Bool_t PReaction::AddBulk(PBulkInterface * mybulk) {
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
    if (sub_reaction) sub_reaction->AddBulk(mybulk);
    return kTRUE;
}

Bool_t PReaction::AddPrologueBulk(PBulkInterface * mybulk) {
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
    if (sub_reaction) sub_reaction->AddPrologueBulk(mybulk);
    return kTRUE;
}

void PReaction::SetReactionId() {
    // compute the reaction id from the list of product ids of 1st channel
    // (to be used in SetUp() as sourceId)

    PChannel *pch0 = channel[0];  // 1st channel in reaction
    Int_t n = pch0->GetNumPar();  // nb of decay products
    Int_t *ids = pch0->GetPids(); // array of product ids
    if (ids[0]/100 != 5) {
	Int_t fac = 1;
	reactionId = 0;
	if(n>5) n = 5;              // cannot handle more in 32 bits
	for (Int_t i=n;i>0;i--) {
	    reactionId += ids[i]*fac; // code reaction id
	    fac *= 100;
	}
	(pch0->GetParticles())[0]->SetSourceId(reactionId); // source id of 1st part. in 1st ch.
    } else reactionId = ids[0];   // this is a thermal source, dilepton generator, etc. 
}

void PReaction:: SetUp(PChannel **pchannel) {
    // get the channels and particles, identify the physics, set up the defaults.

//    cout << "PReaction:: setUp called" << endl;
    
    int j, i, cnew, l=-1, m=-1, k, pcount;

    num_filters=0;
    reset_count=0;
    HGeant = 0;
    evt=NULL;
    rootfile=NULL;
    userSelection=NULL;
    gMethodCall1=NULL;
    userAnalysis=NULL;
    gMethodCall2=NULL;
    decayALL=0;
    fPythia=NULL;
    nMaxBytes=0;
    nTrigCond=0;
    // thetaBeam=0.;
    // phiBeam=0.;
    // sigmaBeam=0.;
    writeINDEX=1;
    fileoutput_pos=0;
    bulkdecay_pos=0;
    pro_bulkdecay_pos=0;
    makeDistributionManager()->LinkDB();
    weight_reset=1;
    pre_heating=0;
    current_projector=NULL;
    is_inline=0;
    if (makeDistributionManager()->GetLoopFilter()) {
	AddBulk(makeDistributionManager()->GetLoopFilter());
	//Add projector from command file
	//TODO: add dummy bulk for separation, otherwise
	//GetCurrentPorjector gets the global one
    }
    
    event_impact_param = (makeStaticData()->GetBatchValue("_event_impact_param"));
    event_plane        = (makeStaticData()->GetBatchValue("_event_plane"));
    vertex_x = makeStaticData()->GetBatchValue("_event_vertex_x");
    vertex_y = makeStaticData()->GetBatchValue("_event_vertex_y");
    vertex_z = makeStaticData()->GetBatchValue("_event_vertex_z");
    weight_version = makeStaticData()->GetBatchValue("_system_weight_version");

    PParticle **pnew=NULL, **poth=NULL;

    particle_stack = new PParticle*[MIN_PARTICLE_STACKSIZE];
    stacksize = MIN_PARTICLE_STACKSIZE; //Default. BUGBUG: Can be extended later

    ndpar=0;                                // initialize total number of product particles
    int quasi=0;
    channel = NULL;
    if (pchannel) {
	if (pchannel[0]->GetQuasi()) {
	    //Workaround for d+p quasi-construction 
	    //Simply move everything by one and fill in the scattering
	    channel = new PChannel*[nchan+1];
	    for (int j=0;j<nchan;++j) {
		channel[j+1]=pchannel[j];
	    }
	    channel[0]=pchannel[0]->GetQuasi();
	    pchannel[0]->ClearQuasi();
	    nchan++;
	    quasi=1;
	}
	else {
	    channel = pchannel;                     // address of the PChannel array
	}
    }
    dindex.Set(nchan-1);                    // indices of parents to channels after the 1st

    for (j=0;j<nchan;++j) {                 // get number of product particles
	if (!channel[j]) {
	    status=2;                           // error
	    printf("%s\n",EType[status]);
	    return;
	}

	ndpar += channel[j]->GetNumPar();     // update for this channel
    }
    // this assumes all sequential channels

    Int_t is_fireball=0;

    if (nchan && (channel[0]->GetParticles()[0]->IsFireball())) {
	//cancel first PFireball, that it does not appear as an "unstable"
	//particle in the loop
	nclones = ndpar;
	is_fireball=1;
    } else
	nclones = ndpar + 1;                    // initialize number of clones for root tree
  
    for (j=0;j<nchan;++j) {
	if (makeDistributionManager()->from_pdecaymanager == 0) {
//	  cout << "ataching " << j << endl;
//	    channel[j]->Print();
	    makeDistributionManager()->Attach(channel[j]); 

	} //if no DecayManager, do this by hand
    }


    for (j=0;j<nchan;++j) {                 // loop over channels
	cnew = channel[j]->GetNumPar();       // number of particles produced in current channel

	pnew = channel[j]->GetParticles();    // address of the array of these particles
	pcount = 0;                           // reset decayed photon counter

	if (j==0) {                           // for the entry channel:

	    particle=new PParticle*[nclones];
	    cindex.Set(nclones);                //    channel index of each reaction particle
	    ftrack.Set(nclones);                //    identifies whether particle is tracked or not
	    if (!is_fireball) {
		m=0;
		particle[0] = pnew[0];       // the parent particle is fixed: first reaction particle
	    }

	    cindex[0] = 1;                      // channel index of 1st particle (or beam)
	    ftrack[0] = 0;
	}

	for (i=1;i<=cnew;++i) {               // now loop over particles of current channel
	    ++m;                                // count of particles in reaction so far
	    particle[m] = pnew[i];              // the current particle
	    particle[m]->SetParent(pnew[0]);    // set pointer to parent
	    cindex[m] = j+1;                    // its channel index
	    int dcount=0;                       // reset decayed-particle count
	    if (j!=nchan-1) {                   // if this is not the last chanel:
		for (k=j+1;k<nchan;++k) {         // loop over subsequent channels
		    poth=channel[k]->GetParticles();// address of that channels particle array
		    if (pnew[i]==poth[0]) {         // current particle decays at a later channel
			if (pnew[i]->Is("dilepton")) ++pcount;
			++dcount;	    
			dindex[k-1] = m;              // keep track which particle decays at which channel
		    }
		}
	    }
	    if (!dcount) {
//		tparticle[++l] = pnew[i];         // the present particle survives decay
		
		ftrack[m]=1;                      // identify as tracked particle

	    } else {
		ftrack[m]=0;                 // not a tracked particle

	    }
	}

    }
    ntpar = l+1;

    PParticle* pParent=NULL;

    if (nchan>0) SetReactionId();
    else nclones=0;

    for (j=0; j<nclones; j++) {
//	particle[j]->Print();

	pParent = particle[j]->GetParent();
	//if (pParent==NULL) cout << "ids= " << particle[j]->ID() << endl;
	//else cout << "ids= " << particle[j]->ID() <<" " << pParent->ID()<< endl;
	particle[j]->SetIndex(j);             // set index in array
	particle[j]->SetParentIndex(-1);      // parent not known yet
	particle[j]->SetDaughterIndex(-1);    // daughter(s) not known yet
	particle[j]->SetSiblingIndex(-1);     // sibling(s) not known yet
	particle[j]->SetSibling(NULL);
    }
  

    //initialize PChannels
    for (j=0;j<nchan;++j) {
	channel[j]->SetDaughters();
    }
    for (j=0;j<nchan;++j) {
	channel[j]->SetPrintTentative(0);
	channel[j]->Init();
	channel[j]->SetPrintTentative(1);
    }

    

}

void PReaction::InitChannels() {
    for (int j=0;j<nchan;++j) {
	channel[j]->SetDaughters();
    }
    for (int j=0;j<nchan;++j) {
	channel[j]->SetPrintTentative(0);
	channel[j]->Init();
	channel[j]->SetPrintTentative(1);
    }

    
}

void PReaction:: SetName(char * name) {
    // Sets up output file names

    if (extTREE&&reset_count) {
	printf("%s\n",RMessage[0]);
	return;
    }
    if (name) {
	filename=name;
	file2=name;
	if (extTREE) file2=tree->GetCurrentFile()->GetName();
	else file2+=".root";  // root file name
    } else {
	//Do NOT write the event file to disc
	//This can be useful when projecting directly to
	//the histogram
	filename="";
	file2="";
    }
    if (rootfile) {
	rootfile->Close();
	delete rootfile;
    }
    if (asciiOUTPUT && (strlen(filename) > 0) ) {      
	// ascii output required
	file1=filename+".evt";
    }
//if (!extTREE && !HGeant) rootfile=new TFile(file2,"RECREATE",file2);
    loop_count=0;              // reset loop counter
    original_filename="";  //not yet needed
}


void PReaction:: InitLoop() {
    // resets the dynamical objects of the reaction, creates tree and opens output files

    int size=ntpar;
    static int count=0;

    if (!evt) evt=new TClonesArray("PParticle",size);   // create on 1st time
    else return;

    reset_count=++count;

    if (!extTREE && !HGeant && (strlen(filename) > 0)) 
	rootfile=new TFile(file2,"RECREATE",file2);
    if (!extTREE && (strlen(filename) > 0))
	tree = new TTree("data","event branch");

    if (allPARTICLES) size=nclones;

    if (tree) {
	if (extTREE && tree->GetBranch("Particles")) {     // reconnect branches
	    tree->SetBranchAddress("Npart",&activeCnt);
	    tree->SetBranchAddress("Impact",event_impact_param);
	    tree->SetBranchAddress("Phi",event_plane);
	    tree->SetBranchAddress("Particles",&evt);
	}
	else {                                             // create branches 
	    tree->Branch("Npart",&activeCnt,"Npart/I");
	    tree->Branch("Impact",event_impact_param,"Impact/D");
	    tree->Branch("Phi",event_plane,"Phi/D");
	    tree->Branch("Particles",&evt,32000,99);
	}
    }

    if (asciiOUTPUT==1) asciif=fopen(file1,"w");

    if (nMaxBytes>0 && tree) tree->SetMaxTreeSize(nMaxBytes);


    Int_t listkey=-1;

    num_filters=0;
    
    Int_t primary_key= makeDataBase()->GetEntry("batch_objects");
    
    if (primary_key>-1)
	while (makeDataBase()->MakeListIterator(primary_key, NBATCH_NAME, LBATCH_NAME,
						&listkey)) {
	    //cout << makeDataBase()->GetName(listkey);
	    if ((*makeDataBase()->GetName(listkey) == '#') && 
		PUtils::ValidVariableName(makeDataBase()->GetName(listkey))) {
		//Filter found
		cout << "Found filter variable " << makeDataBase()->GetName(listkey) << endl;
		if (num_filters < MAX_REACTION_FILTERS) {
		    filter_keys[num_filters]=listkey;
		    filter_counter[num_filters]=0;
		    if (!makeDataBase()->GetParamDouble (listkey,"batch_value", &filter_values[num_filters]))
			Warning("InitLoop", "Filter content variable not found");
		    num_filters++;
		} else {
		    Warning("InitLoop", "Too many filters (MAX_REACTION_FILTERS reached)");
		}
	    }
	    
	} //end iterator
    
}




int PReaction::Loop(int nevents, int wf, int verbose) {
    // Simulate "nevents" events

    // double initialWeight = 1.;
    // if (channel) (channel[0]->GetParticles())[0]->W();

    if (weight_reset) PChannel::SetGlobalWeight(1./nevents);

    double t0=-1;
    if (extTREE&&reset_count&&verbose) printf("%s\n",RMessage[0]);

    if (!is_inline && !extTREE && reset_count>0) {
	Error("Loop","Called a second time -> abort");
	return 0;  // prevent executing loop() more than once
    }
    // from one and the same PReaction object

    int i, j, k, l1, error_count=0, empty_count=0, good_event=0, size,
	size_tracked, twentypercent=nevents/5, cpc=1, ipc=cpc*twentypercent;
    int error_count_array[100];
    for (int i=0;i<100;i++) error_count_array[i]=0;


    // PParticle **file_particle, **ascii_particle, **tree_particle; 
    // particles to be kept on file

    //Create additional particles for the bulk decay
    PParticle *p_array[stacksize];
    PParticle *stable_particle[stacksize];


//    int decay_done[stacksize];

    if (bulkdecay_pos) {
	for (k=0;k<stacksize;k++) { 
	    p_array[k] = new PParticle("dummy");
	    particle_stack[k] = p_array[k];
	}
    } else {
	//Fix stack for the rest of reaction lifetime to save time
	stacksize = nclones;
    }

    ++loop_count;                            // initialized in setNAME
    InitLoop();                              // rebuild dynamic objects
    TClonesArray &pclone=*evt;
    TStopwatch timer;                        // time loop
    timer.Start();
    Int_t print_welcome=1;
    Int_t last_nonempty = -1;

/////////////////////////////    Event loop    ////////////////////////////////

    if (verbose) printf("%s\n",RMessage[1]); //calculating widths in PData...

    makeDistributionManager()->Startup();

    for (i=0;(i<nevents) || (nevents<0);++i) {                // number of events to generate
	l1=0;
	//CRASH;
	//In the prologue all particles are undecayed
	for (k=0; k<stacksize;k++) {
	    decay_done[k]=0;
	}

	//Fill our stack with the particles from
	//the PChannel list
	for (k=0; k<nclones;k++) {
	    particle_stack[k] = particle[k];
	}	

	//First of all, the size is the filled size
	//(can also be 0)

	size = nclones;

	if (nevents>0) {
	    // Timing:
	    if (print_welcome) {
		if (good_event==1) {
		    t0=timer.RealTime();
		    if (verbose) printf("...widths calculated in %f sec\n",t0);
		    if (verbose) printf("New loop: %i events\n",nevents);
		    timer.Continue();
		}
		
		if (verbose && i==1000) {   // estimate total exec time from 1000 events
		    Double_t elaps = timer.RealTime()-t0;
		    printf("Expected execution time = %.3f sec\n",nevents*0.001*elaps);
		    timer.Continue();
		}
		print_welcome=0;
	    }
	    
	    if (i==ipc) {
		if (verbose) printf(" %i%s done in %f sec\n",cpc*20,percent,timer.RealTime()-t0);
		++cpc;
		ipc=cpc*twentypercent;
		timer.Continue();
	    }
	}

	if (channel) {
	    (channel[0]->GetParticles())[0]->SetVertex(*vertex_x,
						       *vertex_y,
						       *vertex_z,0.);
	    (channel[0]->GetParticles())[0]->clearDebugString(); 
	}

	*event_impact_param = 0.0;             // reset reaction impact parameter
	*event_plane = 0.0;              // and event-plane angle

////////////////////////    Decay all particles (Prologue)  ////////////////////////////

	//BUGBUG: If Prologue modifier decay daughter of a PChannel, this should be skipped later

	Bool_t statusOfModify = kTRUE;
	for (int bu =0; bu < pro_bulkdecay_pos ; bu++) {
	    int old_size=size;
	    if(!(statusOfModify = pro_bulk[bu]->Modify(particle_stack, 
						       decay_done, &size, stacksize))) { 
		break; 
	    }

	    for (int ii=old_size;ii<size;ii++) {
		//set index for "new" particles. 
		particle_stack[ii]->SetIndex(ii);		
	    }	
	    for (int ii=0;ii<size;ii++) {
		//check for decayed particles and set them inactive
		if (!allPARTICLES && decay_done[ii]) 
		    particle_stack[ii]->SetInActive();  // decayed particle
	    }
	}
	
        // if Modify failed we will stop the eventloop
	if(statusOfModify == kFALSE){
	    if (nevents>0)
		Warning("Loop()",
			"Bulk return kFALSE. Not full number of events calculated! nEvents = %i !",i);
	    break;
	}

////////////////////////    End of Prologue  ////////////////////////////



	//Loop over all Particles and set the initial weight
	//This includes the particles from the prologue
	for (int ii=0;ii<size;ii++) {
	    if (*weight_version) 
		particle_stack[ii]->SetW(PChannel::GetGlobalWeight() 
					 *particle_stack[ii]->GetMultiplicity());
	    particle_stack[ii]->SetStatus(STATUS_NOT_DECAYED);
	    particle_stack[ii]->SetActive(); //active by default, otherwise it is 
	    //un-initialized in PChannel::decay()
	}
	for (j=0;j<nchan;++j) {
	    PParticle * parent = (channel[j]->GetParticles())[0];
	    if (*weight_version) 
		parent->SetW(PChannel::GetGlobalWeight() * parent->GetMultiplicity());
	    parent->SetGenW(1.);
	    parent->SetInvGenW(1.);
	    //parent->Print();
	}
	


    repeat:


	if (error_count==10000&&!good_event) {
	    status=3;
	    //BUGBUG: Here we should print a more complete message
	    printf("%s\n",EType[status]); //insufficient energy
	    break;
	}

	Int_t ret=0;
	Int_t start=0;

    repeat2:

	for (j=start;j<nchan;++j) {
	    if ( (ret=channel[j]->decay()) ) {    // execute steps in the reaction
		++error_count;
		error_count_array[ret]++;

		if (ret!=8) {
		    last_nonempty = i;
		    goto repeat;            // FAILED, repeat event
		}

	    } 
	    
	    //Clear status flag for the parent particle
	    (channel[j]->GetParticles())[0]->SetStatus(ret);

	} // end of channel loop

	for (j=0;j<nchan;++j) {
	  if (channel[j]->GetNumNotFinalized()) {
	      // cout << channel[j]->GetNumNotFinalized() << endl;
	    for (int k=0; k < channel[j]->GetNumNotFinalized(); k++) {
	      if (!channel[j]->GetDistributionNotFinalized(k)->EndOfChain()) {
		  //		start=j;
		  start=0;
		++error_count;
		error_count_array[STATUS_REDO_CHAIN]++;
		goto repeat2;
	      }
	    }
	  }
	}

	if (pre_heating) {
	    pre_heating--;
	    goto repeat;            // dummy, repeat event
	}
	
	if (ret==8) break;                // EOF in PChannel::readFileInput()

	Double_t my_global_weight=1.;	

	if (*weight_version) {
	    
	    //BUGBUG: The following lines should go into a separate class:
#if 1
	    //Now starting to combine the individual weights
	    //to (at least) a chain weight
	    for (k=0;k<size;k++) { 
		particle_stack[k]->weight_divisor=1.;
		//cout << particle_stack[k]->W() << endl;
	    }
	    
	    //Step1: Collect the divisors 
	    for (k=size-1;k>=0;k--) {
		if (particle_stack[k]->GetParent() && (particle_stack[k]->W()>0.)) {
		    particle_stack[k]->GetParent()->weight_divisor =
			particle_stack[k]->GetParent()->W() /
			particle_stack[k]->W();
		}
	    }
	    
	    //Step2: Combine the divisors 
	    for (k=size-1;k>=0;k--) {
		if (particle_stack[k]->GetParent())
		    particle_stack[k]->GetParent()->weight_divisor *=
			particle_stack[k]->weight_divisor;
	    }
	
	    //Step3: InitWeight
	    for (k=size-1;k>=0;k--) {
		particle_stack[k]->SetW(
		    particle_stack[k]->W() / (particle_stack[k]->weight_divisor));
	    }
	    
	    //Step4: Now we assume that the "Adam" has the correct weight
	    //       for all daughters. This will be copied, thus a
	    //       chain has the same weight
	    //
	    //       At the same time we collect already the "global" weight
	    
	    
	    for (k=size-1;k>=0;k--) {
		if (particle_stack[k]->GetParent()) {
		    particle_stack[k]->SetW (particle_stack[k]->GetParent()->W());
		} else //adam
		    my_global_weight *= particle_stack[k]->W();
	    }
	    
	    //cout << "W: "<< my_global_weight << endl;
	    
	    //Step5: (optional)
	    //Here, we give all particles the same weight (=event weight)
	    //BUGBUG: Problem when we have different chains and the default weight
	    //is the number of events???
	    for (k=size-1;k>=0;k--) 
		particle_stack[k]->SetW (my_global_weight);
#endif

	}

	////////////////////////    Decay all particles (Epilogue)  ////////////////////////////

	//Set the decay flag of the original particles	
	for (k=0; k<size;k++) {
	    if ((particle_stack[k]->GetStatus()) == STATUS_OK)
		decay_done[k]=1;
	    else
		decay_done[k]=0;
	}
	
	//set the remaining particles to 0
	for (k=size; k<stacksize;k++) {
	    decay_done[k]=0;
	}

	statusOfModify = kTRUE;
	for (int ii=0;ii<size;ii++) {
	    //check for decayed particles and set them inactive ... to be consistent with below
	    if (!allPARTICLES && decay_done[ii]) {		   
		particle_stack[ii]->SetInActive();  // decayed particle
	    }
	}
	for (int bu =0; bu < bulkdecay_pos ; bu++) {
	    int old_size=size;

	    //BUGBUG: Preliminary (see above)
	    bulk[bu]->SetWeight(my_global_weight);
// 	    for (k=0;k<size;++k) {
// 		cout << "act: " << particle_stack[k]->IsActive() << endl;}
	    if(!(statusOfModify = bulk[bu]->Modify(particle_stack, 
						   decay_done, &size, stacksize))) { break; }
	    for (int ii=old_size;ii<size;ii++) {
		//set index for "new" particles. 
		particle_stack[ii]->SetIndex(ii);
	    }
	    for (int ii=0;ii<size;ii++) {
		//check for decayed particles and set them inactive
		if (!allPARTICLES && decay_done[ii]) {		   
		    particle_stack[ii]->SetInActive();  // decayed particle
		}
	    }

	}
// 	cout << "after EPI" << endl;
// 	for (k=0;k<size;++k) if (particle_stack[k]->IsActive()) particle_stack[k]->Print();

        // if Modify failed we will stop the eventloop
	if(statusOfModify == kFALSE) {
	    if (nevents>0)
		Warning("Loop()",
			"Bulk return kFALSE. Not full number of events calculated! nEvents = %i !",i);
	    break;
	}
    
	Int_t cnt0 = 0;         // check for empty event
	for (k=0;k<size;++k) {
//     	    cout << "act: " << particle_stack[k]->IsActive() << "done: " <<decay_done[k]<< 
//     		":" << particle_stack[k]->ID() << endl;
	    if (particle_stack[k]->IsActive()==kTRUE) cnt0++;
	}
	if (cnt0==0 && nchan != 0) {
	    empty_count++;
	    error_count++;
	    if (last_nonempty < i) {
		//firts time here....
		error_count_array[STATUS_EMPTY_EVENT]++;
		
	    } else {
		//passed already
		error_count_array[STATUS_PSEUDO_EMPTY_EVENT]++;
	    }
	    goto repeat; 
	} else {
	    last_nonempty = i;
	}



//-------------------------------------------------------------------------------------

// 	Double_t xSave = tree_particle[0]->Px(); // save momenta of 1st particle 
// 	Double_t ySave = tree_particle[0]->Py();
// 	Double_t zSave = tree_particle[0]->Pz();

// This should go into BeamSmearing: (BUGBUG)
// 	if (thetaBeam>0. || sigmaBeam>0.) {    // we have a skewed and/or smeared beam axis
// 	    Double_t thB=0.;
// 	    Double_t phB=0.;
// 	    if (sigmaBeam>0.) {   // gaussian smearing of beam axis
// 		//  thB = TMath::Abs(PUtils::sampleGaus(0.,sigmaBeam));
// 		//  phB = 6.283185308*PUtils::sampleFlat();
// 		Double_t thx = PUtils::sampleGaus(0.,sigmaBeam);   // this gives better results
// 		Double_t thy = PUtils::sampleGaus(0.,sigmaBeam);
// 		thB = sqrt(thx*thx+thy*thy);
// 		phB = TMath::ACos(thx/thB);
// 		if (thy<0.) phB = 6.283185308-phB;
// 	    }
// 	    TVector3 beamAxis(TMath::Sin(thB)*TMath::Cos(phB),
// 			      TMath::Sin(thB)*TMath::Sin(phB),
// 			      TMath::Cos(thB));
// 	    if (thetaBeam>0.) {   // skewing of beam axis (i.e. average beam is off z axis)
// 		TVector3 skew(TMath::Sin(thetaBeam)*TMath::Cos(phiBeam),
// 			      TMath::Sin(thetaBeam)*TMath::Sin(phiBeam),
// 			      TMath::Cos(thetaBeam));
// 		beamAxis.RotateUz(skew);
// 	    }

// 	    for (k=0;k<size;k++) {  // rotate now all active particles
// 		if (treepar[k]->IsActive()) treepar[k]->RotateUz(beamAxis);
// 	    }
// 	}

/////////////////////////    User selection    ////////////////////////////////
//
// A user selection is implemented via a C++ function
// Int_t userSelect(PParticle*) returning 0 on failure and >0 on success.
//
// From a macro, the user selection (and analysis) can be used in 2 ways:
//
// 1) CINT interpreted function
//
// .
// .
// .
// gROOT->ProcessLine(".L userSelect.C");
// r->SetUserSelection(userSelect);
// r->loop(1000);
//
//
// 2) or ACLIC compiled function
//
// .
// .
// .
// gROOT->ProcessLine(".L userSelect.C+");  // or .L userSelect.C++
// r->SetUserSelection(userSelect);
// r->loop(1000); 
//
//
// 3) true compiled function (i.e. embedded into Pluto code)
//
// .
// .
// .
// r->SetUserSelection(userSelect);
// r->loop(1000); 
// 
// Remark:
// A select function explicitely compiled into the code (e.g. in PReaction.cc)
// has to be Set directly with userSelection = select in SetUserSelection(f),
// and not with via TMethodCall. In addition, to be set from the interpreter, select
// has to made known to Cint via an entry in the LinkDef.h file.
// 
///////////////////////////////////////////////////////////////////////////////

	if (userSelection) {  // test user selection function
	    Int_t nTrig = 0;
	    Int_t ret = 0;
	    for (k=0;k<size;k++) {  // loop over all particles
		if (particle_stack[k]->IsActive()) {
		    if ( (ret=userSelection(particle_stack[k])) > 0) nTrig+=ret;
		    else particle_stack[k]->SetInActive();
		}
	    }
	    if (nTrig < nTrigCond) continue;   // skip event
	}

///////////////////////////    Analyze data    ////////////////////////////////
	Int_t retVal=1;
	if (userAnalysis) retVal = userAnalysis(particle_stack,size);  // call analysis
	if (!retVal) continue;

///////////////////////////    Analyze #filters ////////////////////////////////

	Int_t filters_sum=0;
	for (k=0; k<num_filters;k++) {
	    //	    cout << makeDataBase()->GetName(filter_keys[k]) << ":" << 
	    //		*(filter_values[k]) << endl;
	    
	    if (*(filter_values[k])<PUtils::sampleFlat()) {
		//Random variable failed
		filter_counter[k]++;
		filters_sum++;
	    }

	}

	if (filters_sum) continue;


	good_event++;
	PReaction::globalEventCounter++;  // count events across several PReaction::loop() calls



///////////////////////////    Output data    /////////////////////////////////

	size_tracked=0;

	//Write size_tracked
	for (k=0;k<size;k++) {  // loop over all particles
	    if (!decay_done[k]) {
		stable_particle[size_tracked]=particle_stack[k];
		size_tracked++;
//		cout << "found stable: " << particle_stack[k]->ID() << endl;
	    }
	}


	// Output:
	Int_t cnt = 0;
	if (HGeant) {                           // pass event to HGeant and execute
	    int id_tmp[1000], src_tmp[1000], par_tmp[1000], ind_tmp[1000];
	    float px_tmp[1000], py_tmp[1000], pz_tmp[1000], w_tmp[1000];
	    float vx_tmp[1000], vy_tmp[1000], vz_tmp[1000], vt_tmp[1000];
	    cnt = 0;
	    for(k=0;k<size_tracked;++k) {
		PParticle *pt=stable_particle[k];
		if (pt->IsActive()==kFALSE) continue; // skip inactive particles
		id_tmp[cnt] = pt->ID();
		px_tmp[cnt] = pt->Px();
		py_tmp[cnt] = pt->Py();
		pz_tmp[cnt] = pt->Pz();

		if (getVERTEX) {
		    vx_tmp[cnt] = pt->X();        // vertex in mm
		    vy_tmp[cnt] = pt->Y();
		    vz_tmp[cnt] = pt->Z();
		    vt_tmp[cnt] = pt->T()/300.;  // decay time in nanoseconds
		}
		else vx_tmp[cnt] = vy_tmp[cnt] = vz_tmp[cnt] = vt_tmp[cnt] = 0.;

		w_tmp[cnt] = pt->W();
		src_tmp[cnt] = pt->GetSourceId();
		par_tmp[cnt] = pt->GetParentId();
		ind_tmp[cnt] = pt->GetParentIndex();
		cnt++;
		if(cnt==1000) break;
	    }
	    float Ebeam = channel[0]->GetBT();
	    float bpar = *event_impact_param;
	    float phiEv = 57.29578 * (*event_plane);   // go to degrees
	    passEvent(Ebeam,bpar,phiEv,cnt,id_tmp,src_tmp,par_tmp,ind_tmp,
		      px_tmp,py_tmp,pz_tmp,vx_tmp,vy_tmp,vz_tmp,vt_tmp,w_tmp);
	    gROOT->ProcessLine("doGeant(\"trigger 1\");"); 
	} 

	if (asciiOUTPUT) {                      // write event to ASCII file
      
	    if (!getVERTEX) {
		cnt = 0;
		for(k=0;k<size_tracked;++k)  
		    if (stable_particle[k]->IsActive()==kTRUE) cnt++;
		if (writeINDEX==0) {
		    fprintf(asciif," %i %i %f %f 2\n",
			    PReaction::globalEventCounter,cnt,
			    channel[0]->GetBT(),*event_impact_param);
		} else { 
		    if (channel) fprintf(asciif," %i %i %f %f -2\n",
					 PReaction::globalEventCounter,cnt,
					 channel[0]->GetBT(),*event_impact_param);
		    else fprintf(asciif," %i %i %f %f -2\n",
				 PReaction::globalEventCounter,cnt,
				 0.,*event_impact_param);
		}
		for (k=0;k<size_tracked;++k) {               // tracked particles only
		    PParticle *pt=stable_particle[k];
		    if (pt->IsActive()==kFALSE) continue; // skip inactive particles
		    if(writeINDEX==0) {
			fprintf(asciif," %e %e %e %e %i %i %i %e\n",
				pt->E(), pt->Px(), pt->Py(), pt->Pz(),
				pt->ID(), pt->GetSourceId(), pt->GetParentId(), pt->W());
		    } else {
			fprintf(asciif," %e %e %e %e %i %i %i %i %e\n",
				pt->E(), pt->Px(), pt->Py(), pt->Pz(),
				pt->ID(), pt->GetSourceId(), pt->GetParentId(), pt->GetParentIndex(), pt->W());
		    }
		}
	    } else {
		cnt = 0;
		for(k=0;k<size_tracked;++k)  
		    if (stable_particle[k]->IsActive()==kTRUE) cnt++;
		if (writeINDEX==0) {
		    fprintf(asciif," %i %i %f %f 4\n",
			    PReaction::globalEventCounter,cnt,
			    channel[0]->GetBT(),*event_impact_param);
		} else {
		    if(channel) {
		    fprintf(asciif," %i %i %f %f -4\n",
			    PReaction::globalEventCounter,cnt,
			    channel[0]->GetBT(),*event_impact_param);
		    } else {
			fprintf(asciif," %i %i %f %f -4\n",
				PReaction::globalEventCounter,cnt,0.,*event_impact_param);
		    }
		}
		for (k=0;k<size_tracked;++k) {               // tracked particles only
		    PParticle *pt=stable_particle[k];
		    if (pt->IsActive()==kFALSE) continue; // skip inactive particles
		    if (writeINDEX==0) {
			fprintf(asciif," %e %e %e %e %e %e %e %e %i %i %i %e\n",
				pt->E(), pt->Px(), pt->Py(), pt->Pz(),
				pt->T()/300., pt->X(), pt->Y(), pt->Z(),
				pt->ID(), pt->GetSourceId(), pt->GetParentId(), pt->W());
		    } else {
			fprintf(asciif," %e %e %e %e %e %e %e %e %i %i %i %i %e\n",
				pt->E(), pt->Px(), pt->Py(), pt->Pz(),
				pt->T()/300., pt->X(), pt->Y(), pt->Z(),
				pt->ID(), pt->GetSourceId(), pt->GetParentId(),
				pt->GetParentIndex(), pt->W());
		    }
		}
	    }
	}


	if (!HGeant) {
	    cnt = 0;
	    pclone.Delete();

	    //First count the number of particles
	    for (k=0;k<size;++k) {
		if (particle_stack[k]->IsActive()==kFALSE) continue; //skip inactive prtcls.
		cnt++;
	    }

	    Int_t final_num = 0;

	    for (int z=0;z<fileoutput_pos;z++) {
		final_num = size; //reset number
		files[z]->SetHeader(cnt, allPARTICLES, getVERTEX, asciiOUTPUT, writeINDEX , channel, nchan);
		files[z]->Modify(particle_stack, 
				 decay_done, &final_num, 
				 stacksize);
		files[z]->WriteEvent();
	    }
	    cnt = 0; //reset cnt

	    for (k=0;k<size;++k) {               // update TClonesArrays

		if (particle_stack[k]->IsActive()==kFALSE) 
		    continue; //skip inactive prtcls.
		if (!allPARTICLES && decay_done[k]) 
		    continue; //skip decayed prtcls.

		particle_stack[k]->SetT(particle_stack[k]->T()/300.);// go from mm/c to ns
		int save_sclone=particle_stack[k]->GetScatterClone();
		particle_stack[k]->SetScatterClone(0);
		pclone[cnt] = new(pclone[cnt]) PParticle(*particle_stack[k]);
		particle_stack[k]->SetScatterClone(save_sclone);

		for (int z=0;z<fileoutput_pos;z++)
		    files[z]->WriteParticle(particle_stack[k]);
		cnt++;
	    }
	    activeCnt = cnt;

	    if ((tree) &&  ((strlen(filename) > 0) || extTREE)) //otherwise memory leak
		tree->Fill();
	}
    
// 	if (thetaBeam>0. || sigmaBeam>0.) {    // Set 1st particle back to be ready for next event.
// 	    // This is needed because it is not resampled by
// 	    // PChannel::decay() and would cumulate successive rotations.
// 	    treepar[0]->SetPx(xSave);
// 	    treepar[0]->SetPy(ySave);
// 	    treepar[0]->SetPz(zSave);
// 	}


    } // <====================== End of event loop

    if (bulkdecay_pos) {
	for (k=0;k<stacksize;k++) delete p_array[k];
    }

    timer.Stop();
    if (verbose) printf(" 100%s done in %f sec\n CPU time %f sec\n",
			percent,timer.RealTime(),timer.CpuTime());
    if (wf && verbose) printf("\n Total of %i events processed (incl. repeated empty evts.)\n\n", PReaction::globalEventCounter+empty_count);
    //____________________________________________________________________________

    if (error_count && verbose) {
	printf(" %i aborted events were repeated, error codes:\n ",error_count);
	for (int i=0;i<100;i++)
	    if (error_count_array[i])
		printf("%i=%i ",i,error_count_array[i]);
	printf("\n");
    }

    for (i=0; i<num_filters;i++) 
	if (filter_counter[i]) printf(" %i events failed filter %s\n",filter_counter[i],
				      makeDataBase()->GetName(filter_keys[i]));
	
    if (asciiOUTPUT==1) fclose(asciif);
    if (tree) {
	if (!extTREE && (strlen(filename) > 0)) {
	    rootfile = tree->GetCurrentFile();
	    rootfile->Write();
	    rootfile->Close();
	    rootfile=NULL; //otherwise destructor crashes
	    //  delete tree;
	    tree=NULL;
	}
    }


    if (sub_reaction) {
	//Delete evt, otherwise next loop cannot start
	if (evt) evt=NULL;
	good_event+=sub_reaction->Loop(nevents, wf,  verbose);
    }
    return good_event;
}

void PReaction::SetFilter(int ch_id, PFilter * filter) {
    //outdated

    Error("SetFilter","The PFilter class has been removed. Use PProjector instead");

}

void PReaction:: Close() {
    // close root output file

    if (!extTREE && (strlen(filename) > 0)) 
	rootfile->Close();
    else printf("%s\n",RMessage[0]);
}

void PReaction:: PrintReport() const {
    printf("   Reaction Channels:\n");
    for (int j=0;j<nchan;++j) {
	printf("     %i. ",j+1);
	channel[j]->PrintReport();
    }
}

void PReaction:: Print(const Option_t* delme) const {

    int ii=0, esize=nclones-ndpar, j=esize, iid;
    printf("\n Reaction of %i Particles interacting via %i Channels\n",nclones,nchan);
    printf("   Reaction Particles:\n");
    if (esize==2) {                                    // composite particle at entry channel
	printf("     0. %s (beam)\n",particle[0]->Name());    // beam
	printf("     1. %s (target)\n",particle[1]->Name());  // target
    } else if (esize==1) {                             // decay of elementary particle at entry channel
	printf("     0.");
	iid=particle[0]->ID();
	if (iid<1000) {
	    if (particle[0]->IsFireball() )
		printf(" Fireball (%s) \n",makeStaticData()->GetParticleName(particle[0]->ID()-500));
	    else if (particle[0]->IsFileInput() )
		printf(" File Input (%s) \n",particle[0]->Name());
	    else if (particle[0]->IsDilepton() )
		printf(" Source of %ss \n",makeStaticData()->GetParticleName(particle[0]->ID()-500));
	    else
		printf(" %s (decay)\n",particle[0]->Name());
	} else {
	    printf(" quasi-particle (%s beam and %s target)\n",
		   makeStaticData()->GetParticleName(iid%1000),makeStaticData()->GetParticleName(iid/1000));
	}
    }	   
    for (;j<nclones;++j) {
	printf("     %i. %s",j,makeStaticData()->GetParticleName(particle[j]->ID()));
	if ((!allPARTICLES||asciiOUTPUT)*ftrack[j]) {
	    printf(" (tracked particle %i)\n",ii);
	    ++ii;
	} else { printf("\n"); }
    }
    printf("   Reaction Channels:\n");
    for (j=0;j<nchan;++j) {
	printf("     %i. ",j+1);
	channel[j]->Print(delme);
    }
    if (bulkdecay_pos || pro_bulkdecay_pos) {
	printf("   Bulk Classes:\n");
	if (pro_bulkdecay_pos) {
	    printf("     Prologue: ");
	    for (j=0;j<pro_bulkdecay_pos;j++) 
		cout << "<" << pro_bulk[j]->ClassName() << "> ";
	    cout << endl;
	}
	if (bulkdecay_pos) {
	    printf("     Epilogue: ");
	    for (j=0;j<bulkdecay_pos;j++) 
		cout << "<" << bulk[j]->ClassName() << "> ";
	    cout << endl;
	}
    }
 

    if (strlen(filename)>0) {
	printf("   Output Files:\n");
	printf("     Root : %s, %s",(const char*)file2,OType[allPARTICLES]);
	if (getVERTEX) printf(" including vertices.\n");
	else printf(".\n");
	if (asciiOUTPUT) printf("     Ascii: %s, %s\n",(const char*)file1,OType[0]);
    }
    else
	printf("   *NO* output File\n");

    if (sub_reaction) {
	printf("++++++++++++++++++++");
	sub_reaction->Print();

    }
}

void PReaction::SetUserSelection(void *f) {  // set user selection function
    userSelection=NULL;
    if (!f) return ;
    Int_t type =  G__isinterpretedp2f(f);     // determine type of funntion
    if (type==0 || type==3) {  // this is a precompiled global function
	userSelection = (int (*)(PParticle*))f;
	printf("\n>>> User selection ");
    } else { // this is an interpreted function or a method call
	char *funcname = G__p2f2funcname(f);      // modeled on TMinuit constructor 
	TMethodCall* fMethodCall=NULL;
	if (funcname) {
	    fMethodCall = new TMethodCall();
	    fMethodCall->InitWithPrototype(funcname,"PParticle*");
	} else {
	    printf("Function: %s cannot be compiled\n",funcname);
	    return;
	}
	gMethodCall1 = fMethodCall;
	userSelection = CalledSelection;
	printf("\n>>> User selection %s() ",funcname);
    }
    (*userSelection)((PParticle*)-1);       // call once for eventual init work
}

Int_t CalledSelection(PParticle* p) { // interface for interactive user
                                      // selection function
    Long_t args[1];
    args[0] = (Long_t)p;
    gMethodCall1->SetParamPtrs(args);  // pass parameter...
    Long_t result;
    gMethodCall1->Execute(result);     // ... and execute
    return (Int_t)result;
} 

void PReaction::SetUserAnalysis(void *f) { // set user analysis function
    userAnalysis=NULL;
    if (!f) return ;
    Int_t type =  G__isinterpretedp2f(f);     // determine type of function
    if (type==0 || type==3) {  // this is a precompiled global function
	userAnalysis = (int (*)(PParticle**,int))f;
	printf("\n>>> User analysis ");
    } else { // this is an interpreted function or a method call
	char *funcname = G__p2f2funcname(f);      // modeled on TMinuit constructor 
	TMethodCall* fMethodCall=NULL;
	if (funcname) {
	    fMethodCall = new TMethodCall();
	    fMethodCall->InitWithPrototype(funcname,"PParticle**,Int_t");
	} else {
	    printf("Function: %s cannot be compiled\n",funcname);
	    return;
	}
	gMethodCall2 = fMethodCall;
	userAnalysis = CalledAnalysis;
	printf("\n>>> User analysis %s() ",funcname);
    }
    (*userAnalysis)((PParticle**)-1,0);      // call once for eventual init work
}

Int_t CalledAnalysis(PParticle** p, Int_t n) { // interface for user
                                               // analysis function
    Long_t args[2];
    args[0] = (Long_t)p;
    args[1] = (Long_t)n;
    gMethodCall2->SetParamPtrs(args);  // pass parameter...
    Long_t result;
    gMethodCall2->Execute(result);     // ... and execute
    return (Int_t)result;
} 


//Int_t select(PParticle* p) {  // example of a user selection function
//
//  printf("Entering select %x\n",(UInt_t)p);
//  if ((int)p==-1) {
//    printf("= ok <<<\n");
//    return -1;
//  }
//
//  if (p->Charge() != 0.) {  // charged
//    printf("testing %i has %f %f %f\n",p->ID(),p->Px(),p->Py(),p->Pz());
//    return 1;
//  }
//  else return 0;
//
//}







