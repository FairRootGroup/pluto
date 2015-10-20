////////////////////////////////////////////////////////
//  Pluto projector interface
//
//  The main idea of this class is to provide a very
//  simple, fast and easy-to-use analysis tool, which
//  make the writing of an analysis macro unnecessary.
//  It uses the Bulk interface and can be plugged into
//  the reaction before or after the decay
//
//  The syntax of the commands are based on the PBatch syntax,
//  for details look into the documentation of this class.
//
//  Commands for the batch can be added in 2 ways: Either
//  with AddCommand(char * command) or
//  in one step together with a histogram which will 
//  be filled after the command has been executed.
//
//  Input: The input for the commands must be particles
//  from the data stream in the fillowing format:
//  [pid,num] where pid is the pid-string (e.g. "pi0")
//  and num is the number of the particle (can be omitted if only
//  one particle of this type present)
//
//  Output: _x,_y as doubles. The latter one only for 2dim histograms.
//
//                    Author:  Ingo Froehlich
//                    Written: 14/02/2008
//                    Revised: 
//
////////////////////////////////////////////////////////



#include "PProjector.h"
#include "PChannel.h"


PProjector::PProjector() {
    batch_pos=0;
    pid_param= makeDataBase()->GetParamInt("pid");
    link_param= makeDataBase()->GetParamInt("link");
    batch_particle_param = makeDataBase()->GetParamTObj("batch_particle");
    w = makeStaticData()->GetBatchValue("_w");
    
    if (batch_particle_param<0) 
	batch_particle_param = makeDataBase()->MakeParamTObj("batch_particle", "PParticle storage for batch");

    batch_value_param = makeDataBase()->GetParamDouble("batch_value");

    stream_default_pos_param = makeDataBase()->GetParamInt(STREAM_DEFAULT_POS);
    if (stream_default_pos_param<0)  
	stream_default_pos_param = makeDataBase()->MakeParamInt(STREAM_DEFAULT_POS,"Default position");
    stream_max_pos_param = makeDataBase()->GetParamInt(STREAM_MAX_POS);
    if (stream_max_pos_param<0)  
	stream_max_pos_param = makeDataBase()->MakeParamInt(STREAM_MAX_POS,"Max position in stream");


    for (int i=0;i<PROJECTOR_MAX_BATCH;i++) { 
	hist2[i]=NULL;
	hist1[i]=NULL;
	fp_out[i]=NULL;
	fp_in[i]=NULL;
	key_pos_out[i]=0;
	key_pos_in[i]=0;
    }

    batch_pos=0;
    force_weight=1;

    fPriority=PPROJECTOR_PRIORITY;
    proj_nr = makeStaticData()->GetBatchValue("_system_embedded_particle_projector");
}

PProjector::~PProjector() {
    for (int i=0;i<batch_pos;i++) delete batch[i];
}

Bool_t PProjector::AddCommand(const char * command) {

    //adds a command line to batch
    if (batch_pos==PROJECTOR_MAX_BATCH) {
	Error("AddCommand","PROJECTOR_MAX_BATCH reached");
	return kFALSE;
    }
    
    batch[batch_pos]=new PBatch();
    batch[batch_pos]->SetPosition(batch_pos,bulk_id);  //Set absolute adress
    if (!batch[batch_pos]->AddCommand((char*) command))  {
	delete batch[batch_pos];
	return kFALSE;
    }

    key=makeDataBase()->GetEntry("batch_objects");
    if (batch[batch_pos]->Status()) {
	//[+]
	if (*(proj_nr)) {
	    Error("AddCommand","Embedded particle ([+]) used 2 times");
	}
	*(proj_nr) = bulk_id;
    }

    batch_pos++;
    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH3 * histo, const char * command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist3[batch_pos-1]=histo;

    //get the result
    key=makeDataBase()->GetEntry("batch_objects");

    key_x=makeDataBase()->GetEntry("_x");
    key_y=makeDataBase()->GetEntry("_y");
    key_z=makeDataBase()->GetEntry("_z");

    if (key_x < 0) {
	Error ("AddHistogram","result _x not found");
    }
    if (key_y < 0) {
	Error ("AddHistogram","result _y not found");
    }
    if (key_z < 0) {
	Error ("AddHistogram","result _z not found");
    }
    if (!makeDataBase()->GetParamDouble (key_x,batch_value_param,&x)) {
	Error ("AddHistogram","Double _x not found");
    }
    if (!makeDataBase()->GetParamDouble (key_y,batch_value_param,&y)) {
	Error ("AddHistogram","Double _y not found");
    }
    if (!makeDataBase()->GetParamDouble (key_z,batch_value_param,&z)) {
	Error ("AddHistogram","Double _z not found");
    }

    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1]=fillflag;

    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH2 * histo, const char * command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist2[batch_pos-1]=histo;

    //  batch[batch_pos-1]->Print();

    //get the result
    key=makeDataBase()->GetEntry("batch_objects");

    key_x=makeDataBase()->GetEntry("_x");
    key_y=makeDataBase()->GetEntry("_y");

    if (key_x < 0) {
	Error ("AddHistogram","result _x not found");
    }
    if (key_y < 0) {
	Error ("AddHistogram","result _y not found");
    }
    if (!makeDataBase()->GetParamDouble (key_x,batch_value_param,&x)) {
	Error ("AddHistogram","Double _x not found");
    }
    if (!makeDataBase()->GetParamDouble (key_y,batch_value_param,&y)) {
	Error ("AddHistogram","Double _y not found");
    }

    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1]=fillflag;

    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH1 * histo, const char * command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist1[batch_pos-1]=histo;

    //  batch[batch_pos-1]->Print();

    //get the result
    key=makeDataBase()->GetEntry("batch_objects");

    key_x=makeDataBase()->GetEntry("_x");

    if (key_x < 0) {
	Error ("AddHistogram","result _x not found");
    }

    if (!makeDataBase()->GetParamDouble (key_x,batch_value_param,&x)) {
	Error ("AddHistogram","Double _x not found");
    }
    //  cout << x << endl;

    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1]=fillflag;

    return kTRUE;
}

Bool_t PProjector::AddOutputTNtuple(TNtuple * n, const char * command) {

    if (!AddCommand(command)) {
	return kFALSE;
    }

    fp_out[batch_pos-1] = n;

    //get the result
    key=makeDataBase()->GetEntry("batch_objects");
    //Examine NTuple file and get the branches

    TIter iter(fp_out[batch_pos-1]->GetListOfBranches());
    while(TBranch *br = (TBranch *)iter.Next()) {
	const char * name = br->GetName();

	//Each branch should be correlated to the batch key
	
	if (key_pos_out[batch_pos-1]==PROJECTOR_MAX_BRANCHES) {
	    Error("AddNTuple","Too many branches in NTuple");
	    return kFALSE;
	}

	key_list_out[batch_pos-1][key_pos_out[batch_pos-1]]=makeDataBase()->GetEntry((char *) name);

	if (key_list_out[batch_pos-1][key_pos_out[batch_pos-1]] < 0) {
	    Warning("AddOutputNTuple","Branch %s found in NTuple but not defined as a batch value",name);
	}

	key_pos_out[batch_pos-1]++;
    }
    return kTRUE;    
}

Bool_t PProjector::AddInputTNtuple(TNtuple * n,const  char * command) {

    if (!AddCommand(command)) {
	return kFALSE;
    }

    fp_in[batch_pos-1] = n;
    num_events_in[batch_pos-1] = n->GetEntries();
    num_events_in_c[batch_pos-1] = 0; //counted events

    //get the result
    key=makeDataBase()->GetEntry("batch_objects");
    //Examine NTuple file and get the branches

    TIter iter(fp_in[batch_pos-1]->GetListOfBranches());
    while(TBranch *br = (TBranch *)iter.Next()) {
	const char * name = br->GetName();

	//Each branch should be correlated to the batch key
	
	if (key_pos_in[batch_pos-1]==PROJECTOR_MAX_BRANCHES) {
	    Error("AddNTuple","Too many branches in NTuple");
	    return kFALSE;
	}

	key_list_in[batch_pos-1][key_pos_in[batch_pos-1]]=makeDataBase()->GetEntry((char *) name);

	if (key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] < 0) {
	    key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);
	    Info("AddInputTNtuple","Created variable %s for the TNtuple branch",name);
	    //cout << "Created " << name << endl;
	}

	//Check if double is existing
	Double_t *val;
	if (!makeDataBase()->GetParamDouble(key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] ,batch_value_param,&val)) {
	    Double_t * delme =  new Double_t(0.);
	    makeDataBase()->SetParamDouble(key_list_in[batch_pos-1][key_pos_in[batch_pos-1]],"batch_value", delme);
	    //cout << "Added double to " << name << endl;
	}

	//Set the branch adress to the floats
	n->SetBranchAddress(name,&(values_in[batch_pos-1][key_pos_in[batch_pos-1]]));

	key_pos_in[batch_pos-1]++;
    }

    fPriority=FILEINPUT_PRIORITY;

    return kTRUE;    
}

Int_t  PProjector::SetParticles(PParticle ** mstack, int *decay_done, int * num, int stacksize, Int_t first_time) {
    //loop over batch object and see what I can do

    Int_t listkey=-1,*i_result,new_particles=0,counter_all=0,particle_key;

    if (first_time) {
 	for (int i=0;i<*num;i++) {
	    //cout << "FT:" << i << ":" << mstack[i]->ID() << endl;
 	    //Delete old default pos list
	    // Reset value on 1st call -> BUGBUG: can cause trouble when jumping over 2 projectors
 	    particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
 	    if (makeDataBase()->GetParamInt (particle_key, stream_default_pos_param ,&i_result)) {
 		(*i_result)=0;
 	    } 
 	}
	particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
	if (makeDataBase()->GetParamInt (particle_key, stream_default_pos_param ,&i_result)) {
	    (*i_result)=0;
	} 
	for (int i=0;i<*num;i++) {
 	    //Delete old max pos list
 	    particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
 	    if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param ,&i_result)) {
 		(*i_result)=0;
 	    } 
 	}
	particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
	if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param ,&i_result)) {
	    (*i_result)=0;
	} 
	for (int i=0;i<*num;i++) {
	    //Create max list
	    if (mstack[i]->IsActive()) { //count only active particles
		particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
		if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param ,&i_result)) {
		    (*i_result)++;
		} else {
		    Int_t * dummy = new Int_t(1);
		    makeDataBase()->SetParamInt(particle_key ,STREAM_MAX_POS, dummy);
		}
		particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
		if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param ,&i_result)) {
		    (*i_result)++;
		} else {
		    Int_t * dummy = new Int_t(1);
		    makeDataBase()->SetParamInt(particle_key ,STREAM_MAX_POS, dummy);
		}
		counter_all++;
	    }	    
	} 
	//cout << "counted max " << counter_all << endl;
    }//END first_time



    while (makeDataBase()->MakeListIterator(key, NBATCH_NAME, LBATCH_NAME , &listkey)) {
        //loop over all particles
	
	if (makeDataBase()->GetParamInt (listkey, pid_param,&i_result)) {
	    Int_t pid=*i_result;

	    //fill object
	    Int_t pos = -1;
	    if (makeDataBase()->GetParamInt (listkey, link_param,&i_result)) {
		pos=*i_result-1;
	    }
	    
	    //First clear the entry to avoid the use of old objects
	    if (pos>=0) {
		makeDataBase()->SetParamTObj (listkey ,batch_particle_param, NULL);
	    } else if ((pos == -112)  &&  first_time && (*(proj_nr) == bulk_id)) { 
		//stumbled over "+"
		//cout << "CALLED + in "<< bulk_id << endl;
		new_particles++;
		if (new_particles == PROJECTOR_MAX_STACK) {
		    Warning("Modify","PROJECTOR_MAX_STACK reached");
		    return kFALSE;
		}
		if (new_particles>stackpointer) {
		    //create new pparticle
		    stack[stackpointer] = new PParticle(0,0,0,0);
		    stackpointer++;
		    Info("SetParticles","New particle created");

		}
		mstack[*num]=&(stack[new_particles-1]);
		stack[new_particles-1].SetID(pid);
		stack[new_particles-1].SetW(1.0);
		makeDataBase()->SetParamTObj (listkey ,batch_particle_param, &(stack[new_particles-1]));
		(*num)++;
	    } else if (pos == -1) {
		//No pos, default 
		makeDataBase()->SetParamTObj (listkey ,batch_particle_param, NULL);
		Int_t particle_key = makeStaticData()->GetParticleKey(pid);
		if (makeDataBase()->GetParamInt (particle_key, stream_default_pos_param ,&i_result)) {
		    pos = (*i_result);
		} else
		    pos = 0;
	    } else if (pos < -1000) { //found link to variable
		//cout << "pos  is now:" << pos << edl;
		Double_t *res;
		if (makeDataBase()->GetParamDouble ((-(pos+1))-1000, batch_value_param ,&res)) {
		    //cout << "key: " << ((-(pos+1))-1000) << " res: " << *res << endl;
		    //makeDataBase()->ListEntries(-1,1,"*name,batch_value,*num_batch,*pid,*link");
		    pos = ((Int_t) *res)-1;
		} else {
		    pos = -1000;
		}
	    }

	    if (pos == -1000) {
		Error("SetParticles","Unkown particle position for %s",makeDataBase()->GetName(listkey));
	    }

	    //cout << "pos  is now:" << pos << " for " << makeDataBase()->GetName(listkey) <<  endl;

	    for (int i=0;i<*num;i++) {
		//cout << mstack[i]->IsActive() << endl;
		if (mstack[i]->IsActive()) {
		    
		    //		cout << "stack_pid:"<< mstack[i]->ID() << endl;
		    //mstack[i]->Print();
		    if ((mstack[i]->ID() == pid) || (pid==0)) { //0=DUMMY
			//cout << "match " << pid << " at " << pos << endl;
			if (pos==0) {
			    //cout << "match2" << endl;
			    TObject * delme =  (TObject *) mstack[i];
			    makeDataBase()->SetParamTObj (listkey ,batch_particle_param, delme);
			    //makeDataBase()->ListEntries(listkey,1,"name,*pid,*batch_particle");
			} 
			pos--;
		    }

		}

	    }

	}
	


    }
    return 0;

}

Bool_t PProjector::Modify(PParticle ** mstack, int *decay_done, int * num, int stacksize) {
    //cout << "Modify " << endl;

    *w = current_weight;

    SetParticles(mstack, decay_done, num, stacksize, 1);
    //cout << "num:  " << *num << endl;
    //excuting batch
    Int_t startcommand=0;
    for (int i=0;i<batch_pos;i++) {
	
	Int_t retval = batch[i]->Execute(startcommand);
	startcommand=0;

	if ((retval == kTRUE) || (retval == kFOREACH)) {

	    current_weight = *w;
	    
	    //if (hist1[i]) cout << fFillFlag[i] << endl;
	    if (hist2[i] && fFillFlag[i]) {
		if (hist2[i]->GetSumw2()->GetSize() || force_weight)
		    hist2[i]->Fill((*x),(*y),current_weight);
		else
		    hist2[i]->Fill((*x),(*y));
	    }
	    if (hist1[i] && fFillFlag[i]) {
		if (hist1[i]->GetSumw2()->GetSize()|| force_weight)
		    hist1[i]->Fill((*x),current_weight);
		else {
		    hist1[i]->Fill((*x));
		}
		//hist1[i]->Fill((*x));
		//cout << *x << ":"<< current_weight << endl;
		//cout << PChannel::GetGlobalWeight() << endl;
	    }
	    
	    if (fp_out[i]) {
		//fill the ntuple
		Double_t *val;
		
		for (int j=0;j<key_pos_out[i];j++) {
		    if (key_list_out[i][j]>-1) {
			if (makeDataBase()->GetParamDouble (key_list_out[i][j],batch_value_param,&val))
			    values[j]=(Float_t) (*val);
		    }
		}
		
		fp_out[i]->Fill(values);
	    }

	    if (fp_in[i]) {
		//read the ntuple
		if (num_events_in_c[i] == num_events_in[i]) {
		    Info("Modify","NTuple <%s>: number of events reached",fp_in[i]->GetTitle());
		    return kFALSE;
		}

		fp_in[i]->GetEntry(num_events_in_c[i]);
		Double_t *val;

		for (int j=0;j<key_pos_in[i];j++) {
		    if (key_list_in[i][j]>-1) {
			if (makeDataBase()->GetParamDouble (key_list_in[i][j],batch_value_param,&val))
			    *val = values_in[i][j];
		    }
		}
		num_events_in_c[i]++;
	    }
	    

	} else if (retval == kGOTO) {
	    if (batch[i]->GetNewBulk() == bulk_id) {
		//stay in the same PProjector
		SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles for "formore"
		Int_t new_batch = batch[i]->GetNewBatch();
		startcommand = batch[i]->GetNewCommand();
		i =  new_batch -1;
	    } else {
		Error("Modify","Jumping with a GOTO over different Projector not yet implemented");
	    }
	}

	if (retval == kFOREACH) {
	    startcommand = batch[i]->GetNewCommand();
	    i -= 1; //redo current loop
	    SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	}

	if (retval == kFOREACHEND) {
	    SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	    startcommand = 0;
	}

	if (retval == kUPDATE) {
	    SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	    retval = kTRUE;
	    startcommand = batch[i]->GetNewCommand();
	    i -=1; //stay in same batch
	}

    }

    Int_t *i_result;
    //Before leaving clean the max_list again, because next time we could have a different configuration of PIDS
    for (int i=0;i<*num;i++) {
	//Delete old max pos list
	Int_t particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
	if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param ,&i_result)) {
	    (*i_result)=0;
	} 
    }


    return kTRUE;
}


void PProjector::Print(const Option_t* ) const {
    for (int i=0;i<batch_pos;i++) 
	batch[i]->Print();

}


ClassImp(PProjector)
