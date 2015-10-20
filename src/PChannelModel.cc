/////////////////////////////////////////////////////////////////////
//
// This is the base class for coupled channel sampling distributions
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PChannelModel.h"


ClassImp(PChannelModel)

PChannelModel::PChannelModel()  {
} ;

PChannelModel::PChannelModel(const Char_t *id, const Char_t *de, Int_t key) :
    PDistribution(id, de) {
    // Constructor for coupled-channel calculations
    // The PChannelModel must be correlated to a key in the PDataBase
    // If "key" < 0, we try to figure out the key by parsing the descriptor

    model_def_key = is_pid = is_channel = -1;
    sec_key=-1;

    //cout << id << endl;

    makeStaticData();
    SetVersionFlag(VERSION_IS_PRIMARY);
    if (key<0) {
	//Needs some parser to find out what we need
	char parser_delim[2]="_";
	char path_delim='/';

	if (key==-2) { // new syntax
	    path_delim='&';
	    strcpy(parser_delim,"#");
	}

	char *arr1[10+2];
	Int_t arr1_s=12; //max 10 decay products + parent + id

	Int_t alt_position=0,path_position=0;
	for (Int_t i=strlen(id);i>=0;i--) {
	    if (id[i]=='@') alt_position=i;
	}

	if (alt_position) { 
	    alt_position++; //skip "@"
	}

	char * db_id = new char[strlen(id)+1-alt_position];
	strcpy(db_id,&(id[alt_position]));
	//	cout << "db_id " << db_id << endl;
	for (Int_t i=strlen(id);i>=alt_position;i--) {
	    if (id[i]==path_delim) path_position=i;
	}

	if (path_position) { 
	    path_position++; //skip "/"
	    model_def_key = makeStaticData()->
		MakeDirectoryEntry("modeldef",NMODEL_NAME,LMODEL_NAME,&(id[path_position]));
	    //cout << &(id[path_position-1]) << endl;
	    char * nid = new char[strlen(id)+1];
	    strcpy(nid,id);
	    nid[path_position-1]='\0';
	    id=nid;
	}

	PUtils::Tokenize(&(id[alt_position]), parser_delim, arr1, &arr1_s);

	if (path_position) { 
	    char * nid = new char[strlen(id)+1];
	    strcpy(nid,id);
	    nid[path_position-1]=path_delim;
	    id=nid;
	}
	//cout << arr1_s << ":" << arr1[0] << endl;
	if ((arr1_s <= 2)) { //Particle
	    //key=makeDataBase()->GetEntry(arr1[0]);
	    key=makeStaticData()->GetParticleKey(arr1[0]);
 	    if (path_position) {
		sec_key = makeStaticData()->GetSecondaryKey(key,model_def_key);
		is_pid=makeStaticData()->GetParticleIDByKey(key);
		if ( sec_key>0 ) {
		    key = sec_key;
		}
		else {
		    key = makeStaticData()
			->AddAlias(makeStaticData()->
				   GetParticleName(is_pid), db_id);
		    sec_key = key;
		    makeDataBase()->SetParamInt (key, "defkey", &model_def_key);
		}
		ClearVersionFlag(VERSION_IS_PRIMARY);
		//cout << key << endl;
 	    }
	} else if (arr1_s > 2) { //Decay

	    Int_t arr1_id[10+2];
	    arr1_id[0]=makeStaticData()->GetParticleID(arr1[0]);
	    
	   
	    for (int i=1;i<=(arr1_s-2);i++) 
		arr1_id[i] = makeStaticData()->GetParticleID(arr1[i+1]);

	    key = makeStaticData()->GetDecayKey(arr1_id,arr1_s-2);
	    if (key<0) key = makeStaticData()->AddDecay(arr1_id,arr1_s-2);
	    is_channel = makeStaticData()->GetDecayIdxByKey(key);

	    //cout << "parse decay got " << key << endl;

 	    if (path_position) {
		//before we continue here, we have to check if alias
		//for this decay already exists
		sec_key = makeStaticData()->GetSecondaryKey(key,model_def_key);
		if ( sec_key>0 ) {
		    key = sec_key;
		}
		else {		    
		    key = makeStaticData()->AddAlias(
			makeStaticData()->GetDecayNameByKey(key), db_id);
		    //defkey has to be stored for future use!
		    makeDataBase()->SetParamInt (key, "defkey", &model_def_key);
		    sec_key = key;
		}
		ClearVersionFlag(VERSION_IS_PRIMARY);
 	    }

	    //cout << "after sec " << key << endl;

	    //To make life more easy, add "default" track template
 	    Add(arr1[0],"parent");
 	    for (int i=1;i<=(arr1_s-2);i++) 
 		Add(arr1[i+1],"daughter");
 	    preliminary_particles=1;


	} // end "decay"

	if (alt_position) { 
	    //reset ID:
	    char * newid = new char[alt_position];
	    strncpy(newid,id,alt_position);
	    newid[alt_position-1]='\0';
	    //cout << id << ":" << newid << endl;
	    identifier=newid;
	    
	}
//	cout << id <<" key: " << key << endl;

    }

    

    if (key<0) {
	Warning("PChannelModel","Primary key not identified");
    }

    primary_key=key;
    width_init=0;
    mesh = NULL;
    loop_flag=0;

    if (sec_key<0) { //No sec. model
	is_channel=makeStaticData()->GetDecayIdxByKey(key);
	is_pid=makeStaticData()->GetParticleIDByKey(key);
    } 

    if (is_pid<0 && !is_channel) {
	Warning("PChannelModel","Undetermined particle for key %i", key);
    }

    mmin=0;
    mmax=0;
    Int_t *maxmesh_db;
    if (makeDataBase()->GetParamInt (key, "maxmesh",&maxmesh_db))
	maxmesh=*maxmesh_db;
    else
	maxmesh=302;

    draw_option=0;
    didx_option=-1;
    //some TF1 defaults...
    fNpx       = (maxmesh-1)*3; //3 interpolation points per meshpoint should be enough

    fNpar = 4; //Parameter1: Draw mass/width, Parameter2: Draw _specific_ didx for partial decay
    if (fNpar) {
	fNames      = new TString[fNpar];
	fParams     = new Double_t[fNpar];
	fParErrors  = new Double_t[fNpar];
	fParMin     = new Double_t[fNpar];
	fParMax     = new Double_t[fNpar];
	for (int i = 0; i < fNpar; i++) {
	    fParams[i]     = 0;
	    fParErrors[i]  = 0;
	    fParMin[i]     = 0;
	    fParMax[i]     = 0;
	}
	fNames[0] = "Drawing option: 0: Weigth, 1:Width, 2: BR";
	fNames[1] = "0: Total weight, 1:Partial width(didx) for Weight";
	fNames[2] = "Amplitude(abs) for 1st propagator term for didx";
	fNames[3] = "Amplitude(phase) for 1st propagator term for didx";
    }
    fParams[0]=0;
    fParams[1]=-1;
    fParams[2]=0;
    fParams[3]=0;
    mc_max = 1000; //for the integration methods in Width-mesh calculation

    if (is_pid>=0) SetDynamicRange(PData::LMass(is_pid),PData::UMass(is_pid));

    if ((is_channel>=0) && (sec_key<0)) {
	SetDynamicRange(makeStaticData()->GetDecayEmin(is_channel),
			PData::UMass(makeStaticData()->GetDecayParentByKey(key)
				     ));
	exp_w_mean=makeStaticData()->GetDecayBR(is_channel);

    }
    didx_param     = makeDataBase()->GetParamInt("didx");
    scfactor_param = makeDataBase()->GetParamDouble("scfactor");
    unstable_width = makeStaticData()->GetBatchValue("_system_unstable_width");
} ;

PDistribution* PChannelModel::Clone(const char*delme) const {
    return new PChannelModel((const PChannelModel &)* this);
};

PChannelModel * PChannelModel::GetSecondaryModel(const char * name) {
    sec_key = makeStaticData()->MakeDirectoryEntry("modeldef",NMODEL_NAME,LMODEL_NAME,name);
    if (sec_key<0) return NULL;

    Int_t listkey = makeStaticData()->GetSecondaryKey(primary_key,sec_key);
    if (listkey<0) return NULL;

    Int_t model_param = makeDataBase()->GetParamTObj("model");
    TObject *t_result = NULL;
    if (makeDataBase()->GetParamTObj (listkey , model_param, &t_result)) {
	if ((((PChannelModel *) t_result) -> GetDef())  == sec_key)
	    return (PChannelModel *) t_result;
	else {
	    Warning("GetSecondaryModel","Consistency check failed");
	    return NULL;
	}
    }
    return NULL;
}

Bool_t PChannelModel::SampleMass(Double_t *mass,  Int_t *didx) {
    //Samples the masses of the known decay products (or single particle)
    //Order of particles is the same as in data base
    
    return kFALSE;
}

Double_t PChannelModel::EvalPar(const Double_t *x, const Double_t *params) {
    if (params) {
	draw_option=(int)params[0];
	didx_option=(int)params[1];
    }

    return Eval(x[0]);
}
 
Double_t PChannelModel::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t my_x[11];
    Int_t my_didx_option[11];
    my_x[0]=x;
    my_didx_option[0]=didx_option;

    if (draw_option==0) {
//	cout << my_didx_option[0] << endl;
	return ((PChannelModel*)this)->GetWeight(my_x,my_didx_option);
	//return res;
    }
    if (draw_option==1) {
	((PChannelModel*)this)->GetWidth(x,&res,didx_option);
	return res;
    }
    if (draw_option==2) {
	((PChannelModel*)this)->GetBR(x,&res);
	return res;
    }
    if (draw_option==3) {
	return ((PChannelModel*)this)->GetAmplitude(my_x,my_didx_option).Re();
    }
    if (draw_option==4) {
	return ((PChannelModel*)this)->GetAmplitude(my_x,my_didx_option).Im();
    }
    return 0;
}


Double_t PChannelModel::GetWeight(Double_t *mass, Int_t *didx) {
    //Get the weight for the masses of the known decay products 
    //(or single particle)
    //Decay: Order of particles is the same as in data base,
    //mass[0] is the parent mass
    //Single particle: Use only mass[0]

    if (loop_flag>0) {
	loop_flag=0;
	return 1.;
    }
    loop_flag++;
    Double_t r=GetAmplitude(mass, didx).Rho2();
    if (loop_flag>0)
	loop_flag--;
    return r;
}


TComplex PChannelModel::GetAmplitude(Double_t *mass, Int_t *didx) {
    //Get the amplitude for the masses of the known decay products 
    //(or single particle)
    //Decay: Order of particles is the same as in data base,
    //mass[0] is the parent mass
    //Single particle: Use only mass[0]

    if (loop_flag>0) {
	loop_flag=0;
	return TComplex (1.,0,kTRUE);
    }
    loop_flag++;
    TComplex r=TComplex (GetWeight(mass,didx),0,kTRUE);
    if (loop_flag>0)
	loop_flag--;
    return r;
};

Bool_t PChannelModel::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    //Calculates the mass-dependent width (either partial for the decay or
    //total for the particles)
    //If didx is set, use the partial width for a particular target channel

    // Use by default only static values
    if (didx>=0) {
	//Do not use the local width
	//but use the partial width only
	*width = makeStaticData()->GetDecayPartialWidth(didx);
	return kTRUE;
    }
    if (is_channel) {
	*width = makeStaticData()->GetDecayPartialWidth(is_channel);
	return kTRUE;
    }
    if (is_pid) {
	*width = makeStaticData()->GetParticleTotalWidth(is_pid);
	return kTRUE;
    }
    return kFALSE;
}


Bool_t PChannelModel::GetBR(Double_t mass, Double_t *br, Double_t totalwidth) {
    //Calculates the mass-dependent br
    //This method is meaningless for particles
    //the totalwidth may be set by the user, otherwise we use the static one

    if (is_pid>=0)
	return kFALSE;
    
    Double_t lwidth;

    //test if the GetWidth was overloaded
    //if not, simply use the static values
    if (!GetWidth(mass,&lwidth)) {
 	*br = makeStaticData()->GetDecayBR(is_channel);
 	return kTRUE;
    }
    
    //Make everything coherent to Dyn.Data
    Double_t sc=1.;
    Double_t *scfactor=&sc;
//    cout << didx_param <<","<< is_channel <<","<< scfactor_param<<endl;
    makeDataBase()->GetParamDouble (didx_param, is_channel , scfactor_param, &scfactor);

    lwidth*=*scfactor;
    //if total width is below unstable width-> static BR
    if (totalwidth<0.)
	totalwidth=makeStaticData()->GetParticleTotalWidth(
	    makeStaticData()->GetDecayParent(is_channel));

    if (totalwidth <= (*unstable_width)) 
	*br = makeStaticData()->GetDecayBR(is_channel);
    else  {
	//if we have a local width, the BR is the partial width divided by total
//	cout << lwidth << ":" << totalwidth << endl;
	*br = lwidth/totalwidth;
    }
    return kTRUE; 
	
}


int PChannelModel::GetDepth(int i) {
    //By default don't care
    return -1;
}
