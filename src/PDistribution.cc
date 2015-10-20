/////////////////////////////////////////////////////////////////////
//  PDistribution Class implementation file
//
//  PDistribution keeps information about generic distributions
//  This is only the base class
//  Algorithms implemented in the inherited classes
//  Distributions are inherited from TF1, this means that
//  the distributions may use Eval to make drawing very easy
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>

#include "PDistribution.h"


ClassImp(PDistribution)

PDistribution::PDistribution() {
    identifier=new char[1];
    description=new char[1];
    version_flag=VERSION_SAMPLING; 
    no_daughters=0;
    enable=1;
    primary_key=-1;
    group=NULL;
    position=0;
    linkdbdone=0;
    debug_flag=0;
    for (int i=0; i<MAX_PARTICLE_LIST; i++) {
	particle[i]=NULL;
	names[i]=NULL;
	particle_flag[i]=0;
	pid[i]=0;
	private_flag[i]=NULL;
	private_flag_int[i]=NULL;
    }
}

PDistribution::PDistribution(const Char_t *id, const Char_t *de):
    TF1(id,"0", 0, 1) {

    identifier=new char[strlen( id ) + 1];
    description=new char[strlen( de ) + 1];
    strcpy((char *)identifier,id);
    strcpy((char *)description,de);

    no_daughters=0;

    w_max=-1.;
    w_num=0;
    w_sum=0;
    exp_w_mean=1.;
    dynamic_range=1.;
    opt_string=NULL;
    
    version_flag=VERSION_SAMPLING; 
    relative_warning=1;

    enable=1;
    primary_key=model_def_key=sec_key=-1;

    group=NULL;
    position=0;
    linkdbdone=0;
    debug_flag=0;
    preliminary_particles=0;
    for (int i=0; i<MAX_PARTICLE_LIST; i++) {
	particle[i]=NULL;
	names[i]=NULL;
	particle_flag[i]=0;
	pid[i]=0;
	private_flag[i]=NULL;
	private_flag_int[i]=NULL;
    }

    // Copied from TF1 constructor... 
    SetName(id);
    fNpar = 0;

    fXmin      = 0.;
    fXmax      = 0.;
    
    fNpx       = 1;
    fType      = 0;
    //    fFunction  = 0;
    fNdim = 1;
    //... to be overwritten by PChannelModel

}

PDistribution::~PDistribution() {
    delete[] identifier;
    delete[] description;
}

PDistribution* PDistribution::Clone(const char*delme) const {
    Fatal("Clone","Clone calling virt. function for %s",description);
    return NULL;
};

Bool_t PDistribution::FreezeOut(void) {
    //called once after the link to the data base has been done
    //LinkDB()

    
    return kTRUE;    
};


Bool_t PDistribution::Init(void) {
    //reading the private flags e.g. 
    //called directly after attaching the distribution to the channel
    //particle pointers complete at this point

    
    return kTRUE;    
};



Bool_t PDistribution::IsValid(void) {
    //check after genbod if the generated daughters are matching to the distributions
    //Each inherited class must decide on its own, what should be done here
    //Can be used to cut-out distributions from phase space


    return kTRUE;   
};

Bool_t PDistribution::Prepare(void) {
    //make the preparation is front of genbod here
    //e.g. rotations
    return kTRUE;   
};

Bool_t PDistribution::SampleMass(void) {
    //Can be user-defined mass sampling

    return kFALSE; //Switched off by default
};

Bool_t PDistribution::SampleMomentum(void) {
    //Can be user-defined momentum sampling (e.g. genbod)

    return kFALSE; //Switched off by default
};

Bool_t PDistribution::SampleAngle(void) {
    //Can be user-defined angular sampling

    return kFALSE; //Switched off by default
};

Bool_t PDistribution::Finalize(void) {
    //rotate back, if needed
    //return kFALSE if EndOfChain should be called
    return kTRUE;   
};


Bool_t PDistribution::EndOfChain(void) {
  return kTRUE; 
};

Bool_t PDistribution::CheckAbort(void) {
    //Moreover, the class must decide if the complete chain must be aborted
    //e.g. jump over a resonance, and look to the angular distributions
    //of the decay products
    //If true, resample the complete reaction chain
     return kFALSE;
};

Double_t PDistribution::GetWeight(void) {
    //Support for weights in the PChannel::Decay()
    return 1.;
}

Bool_t PDistribution::WriteDebugInfo(PParticle *par) {
    par->addDebugString((char *)GetDescription());

    return kTRUE;
};

void PDistribution::Reset(void) {
    for (int i=0; i<MAX_PARTICLE_LIST; i++) {
	particle[i]=NULL;
    }
    ResetStatus();
};

Int_t PDistribution::GetStatus(void) {
    for (int i=0; i<position; i++) {
	if (particle[i]==NULL) {
	    //cout << "NOT SET:" << pid[i]<< endl;
	    return -1;
	}
    }
    return 0;
}

Int_t PDistribution::Add(const Char_t * opt) {
    parse_s=4;
    PUtils::Tokenize(opt, ",", parse, &parse_s);
    if (parse_s==1) {
	Warning("Add","Not enough options in %s",opt);
	return -1;
    }
    if (parse_s==2) return Add(parse[0],parse[1]);
    if (parse_s==3) return Add(parse[0],parse[1],parse[2]);
    if (parse_s==4) return Add(parse[0],parse[1],parse[2],parse[3]);
    return -1;
}

Int_t PDistribution::Set(const Char_t * opt) {
    parse_s=2;
    PUtils::Tokenize(opt, ",", parse, &parse_s);
    if (parse_s==1) {
	Warning("Set","Not enough options in %s",opt);
	return -1;
    }
    if (parse_s==2) return Set(parse[0],parse[1]);
  
    return -1;
}

Int_t PDistribution::Add(const Char_t * name, int flag, const Char_t * pflag) {
    if (position == MAX_PARTICLE_LIST) {
	Error("Add","MAX_PARTICLE_LIST reached");
	return -1;
    }

    if (GetVersionFlag(VERSION_IS_PRIMARY) 
	&& !(flag==PARTICLE_LIST_DAUGHTER) 
	&& !(flag==PARTICLE_LIST_PARENT) 
	&& relative_warning ) {
	Warning("Add","The primary model (%s) must not have other relatives then parent or daughter",identifier);
    }

    if (preliminary_particles) {
 	position=0; //Overwrite everything
 	preliminary_particles=0;
    }

    particle[position]=NULL;
    names[position]=name;
    particle_flag[position]=flag;
    private_flag[position]=pflag;
    if (strcmp(name, "q")==0) { //particle is quasi
	pid[position]=DISTRIBUTION_QUASI_PID;
	relative_warning = 0; //switch off for grandparents
    }
    else if (strcmp(name, "?")==0) //can be anything
	pid[position]=DISTRIBUTION_SOMETHING_PID;
    else if (strcmp(name, "N")==0) //any nucleon
	pid[position]=DISTRIBUTION_NUCLEON_PID;
    else
	pid[position]=makeStaticData()->GetParticleID(name);
    position++;
    return 0;
};

Int_t PDistribution::Set(const Char_t * name, const Char_t * pflag) {
    for (int i=0;i<position;i++) {
//	cout << names[i] << private_flag[i] << endl;
	if ((strcmp(name, names[i])==0)) {
	    private_flag[i]=pflag;
	    return 0;
	}
    }
    return -1;
}

Int_t PDistribution::Add(const Char_t * name, const Char_t *  flag1, 
			 const Char_t *  flag2, const Char_t *  flag3) {

    Int_t tmp = GetFlag(flag1);
    if (tmp == 0) {
	cout << "PDistribution::Add: Could not get flag. First flag must be *not* private" << endl;
	return -1;
    }

    tmp |= (GetFlag(flag2) | GetFlag(flag3));
    
    const Char_t * pflag = NULL;
    if ((GetFlag(flag2) == 0) && (flag2  != NULL)) pflag=flag2;
    if ((GetFlag(flag3) == 0) && (flag3  != NULL)) pflag=flag3;
    if (!flag2 || GetFlag(flag2)) {  //catch the "parent,sibling" case
	TString * delme = new TString(flag1);
	delme->ToLower(); //memory leak but avoids GetDaughter bug
	pflag=(Char_t*)delme->Data(); //copy the official flag to private one to make GetParticle("parent");
    }

    return Add(name, tmp, pflag);
}

Int_t PDistribution::GetFlag(const Char_t *  flag) {
    if (flag == NULL) return 0;
    if (strcmp(flag, "daughter") == 0) return PARTICLE_LIST_DAUGHTER;
    if (strcmp(flag, "DAUGHTER") == 0) return PARTICLE_LIST_DAUGHTER;
    if (strcmp(flag, "parent") == 0) return PARTICLE_LIST_PARENT;
    if (strcmp(flag, "GRANDDAUGHTER") == 0) return PARTICLE_LIST_GRANDDAUGHTER;
    if (strcmp(flag, "granddaughter") == 0) return PARTICLE_LIST_GRANDDAUGHTER;    
    if (strcmp(flag, "PARENT") == 0) return PARTICLE_LIST_PARENT;
    if (strcmp(flag, "grandparent") == 0) return PARTICLE_LIST_GRANDPARENT;
    if (strcmp(flag, "GRANDPARENT") == 0) return PARTICLE_LIST_GRANDPARENT;
    if (strcmp(flag, "grandgrandparent") == 0) return PARTICLE_LIST_GRANDGRANDPARENT ;
    if (strcmp(flag, "GRANDGRANDPARENT") == 0) return PARTICLE_LIST_GRANDGRANDPARENT ;
    if (strcmp(flag, "sibling") == 0) return PARTICLE_LIST_SIBLING;   
    if (strcmp(flag, "SIBLING") == 0) return PARTICLE_LIST_SIBLING;   
    return 0;

}

void  PDistribution::ResetRelatives(Int_t flag) {
    //Cleans all particles which are not a direct parent or daughter
    //This is used to reset the PChannel and make distributions again tentative

    for (int i=0; i<position; i++) {
	if ((particle_flag[i] != PARTICLE_LIST_DAUGHTER) &&
	    (particle_flag[i] != PARTICLE_LIST_PARENT)) {
	    if (flag && (particle_flag[i]==flag)) {
		particle[i] = NULL;
	    } else if (!flag) {
		particle[i] = NULL;
	    }
	}
    }
    ResetStatus();
}

Int_t PDistribution::SetParticle(PParticle *part, int mypid, int flag) {

    if (no_daughters && (flag == PARTICLE_LIST_DAUGHTER)) 
	return 0; //Do not check daughters if not required

     for (int i=0; i<position; i++) {
	if (((pid[i] == mypid ) && (particle_flag[i] == flag) && (particle[i] == NULL))	    
	    || ((pid[i] == DISTRIBUTION_SOMETHING_PID)  
		&& (particle[i] == NULL) && (particle_flag[i] == flag))
	    || ((pid[i] == DISTRIBUTION_NUCLEON_PID) &&
		(makeStaticData()->GetParticleID("n") == mypid 
		 || makeStaticData()->GetParticleID("p") == mypid ) && 
		(particle[i] == NULL) && (particle_flag[i] == flag))
	    ||
	    ((flag == PARTICLE_LIST_DAUGHTER) && (mypid>1000) && (pid[i]==DISTRIBUTION_QUASI_PID) 
	     && (particle[i] == NULL))) {
	    particle[i] = part;

	    /* BUGBUG: If "anything" already set, and the template is matching
	       then anything has the be shifted, or ? */


	    //Add granddaughters
	    if ((flag == PARTICLE_LIST_DAUGHTER) && (mypid>1000) && (pid[i]==DISTRIBUTION_QUASI_PID)) {
		int pid1=mypid/1000,pid2=mypid%1000;
		PParticle *beam= part->GetScattering(0);
		PParticle *target= part->GetScattering(1);
		if (!beam || !target) return -1; //have not found anything
		int xpid1=beam->ID(),xpid2=target->ID();

		for (int j=0; j<position; j++) {
		    if ((particle_flag[j] == PARTICLE_LIST_GRANDDAUGHTER) && (particle[j] == NULL)) {
			if (pid[j] == pid1 && xpid1 == pid[j]) {
			    particle[j] = beam;
			    xpid1=-1;
			}
			else if (pid[j] == pid2 && xpid1 == pid[j]) {
			    particle[j] = beam;
			    xpid1=-1;
			}
			else if (pid[j] == pid1 && xpid2 == pid[j]) {
			    particle[j] = target;
			    xpid2=-1;
			}
			else if (pid[j] == pid2 && xpid2 == pid[j]) {
			    particle[j] = target;
			    xpid2=-1;
			}
		    }
		}
	    }
	    return 0;
	}

	    

    }
    return -1; //have not found anything
};

Int_t PDistribution::CheckDaughters(void) {
    if (no_daughters)  {
	return 0;
    }

    //check if all daughters are complete, return -1 if not
    for (int i=0; i<position; i++) {
	if ((particle_flag[i] == PARTICLE_LIST_DAUGHTER) && (particle[i] == NULL)) {
	    return -1;
	}
    }
    return 0; //all daughters are set
};

void PDistribution::ResetStatus(void) {
    for (int i=0; i<position; i++) {
	if (private_flag_int[i]) {
	    private_flag[i] = private_flag_int[i];
	    private_flag_int[i] = NULL; //remove saved flag
	}
    }
};

PParticle *PDistribution::GetParticle(const Char_t * pflag) {
    //Read the particle in the Init() method of the inherited distribution
    //"pflag" can be a private identifier string (as used in the Add() method)
    //or one of the official strings like daughter, parent etc...
    //Each particle object is been read out only one time


    for (int i=0; i<position; i++) {
	if (private_flag[i]) {
	    if (pflag) {
		if ((strcmp(pflag, private_flag[i]) == 0) ||
		(strcmp(pflag, makeStaticData()->GetParticleName((pid[i] < 1000)? pid[i] : 0)) == 0)) {
		    private_flag_int[i] = private_flag[i]; //save flag for later use
		    current_flag = particle_flag[i];
		    private_flag[i] = NULL; //remove used flag

		    return particle[i];
		}
	    } else {
		opt_string = private_flag[i];
		return particle[i];
	    }	    
	}   
    }
    return NULL;
};

void PDistribution::BasePrint(void) const  {
    cout << "[" << identifier; 
    if (GetVersionFlag() & VERSION_INVERT_WEIGHTING)
	cout << ",IW";
    if (GetVersionFlag() & VERSION_SAMPLING)
	cout << ",S";
    cout << "] " << description << " <" << ClassName() << ">" << endl;
    for (int i=0; i<position; i++) {
	cout << "    " << names[i] << ", PID: " << pid[i];
	if (particle_flag[i] == PARTICLE_LIST_DAUGHTER)
	    cout << " is a Daughter";
	else if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    cout << " is a Parent";
	if (private_flag[i]) 
	    cout << " [" << private_flag[i] << "]" << endl;
	else cout << endl;
	if (particle[i]) { //particle object are already defined
	    particle[i]->Print();
	}
    }
}

void PDistribution::Print(const Option_t* delme) const {  
    //Debug info
    BasePrint();
}

void PDistribution::SubPrint(Int_t opt) const {
//Print sub-models
    return;
}

const Char_t *PDistribution::Path() const {
    if (GetVersionFlag(VERSION_IS_PRIMARY)) return "\0";
    if (model_def_key < 0) return NULL;
    return makeDataBase()->GetName(model_def_key);
}

int PDistribution::GetDepth(int i) {

    return -1; //nothing
}


Bool_t PDistribution::ExecCommand(const char * command, Double_t value) {
    //Pure virtual
    return kFALSE;
}

Bool_t PDistribution::Exec(const char * command) {

    if (strlen(command) == 0) {
	return ExecCommand("",0);
    };

    char *array[200];
    Int_t array_s=200; 
    PUtils::Tokenize(command, ";", array, &array_s);

    
    for (int i=0; i<array_s; i++) {
	char *array2[200];
	Int_t array2_s=200; 
	PUtils::Tokenize(array[i], "=", array2, &array2_s);

	if (array2_s > 2) {
	    Warning("Exec","Syntax error: Too many ='s");
	    return kFALSE;
	} else if (array2_s == 2){
	    Double_t arg;
	    sscanf(array2[1],"%lf",&arg);
	    if (!ExecCommand(array2[0],arg)) {
		Warning("Exec","Syntax error: %s",array2[1]);
	    }
	} else {
	    if (!ExecCommand(array[i],0.)) {
		Warning("Exec","Syntax error: %s",array[i]);
	    }
	}
    }

    return kTRUE;
}

