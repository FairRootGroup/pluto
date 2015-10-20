////////////////////////////////////////////////////////
// Batch commands for particle and value operations
//
// This class allows for particle operations for
// analysing angles, masses, or making boost and
// rotations of TLorentzVectors/PParticles
//
// It is used for the PProjector to make histograms
// _inside_ the event loop without the need of an
// analysis macro
//
// The syntax is inspired by C++, but only a few
// rules of C++ have been implemented. Furthermore, the syntax
// has been simplified to make it more readable, and an interface
// to the particle input (and maybe to the database in a later step)
// has been foreseen
//  
// As usual we have assignments, operators, functions and methods
// Each operation can be seperated by a semicolon (;)
//
// A simple case is outlined in the following command
// _x = [p,1]->M();
// This means that the double _x (which is automatically constructed)
// is assigned to the mass of the particle [p,1]
// Particle in brackets are the only allowed use without instantiantion
// as they are (hopefully) filled by the input. The additional number 
// can be used for ordering information
//
// + and - operation
// This adds/subtracts PParticles or doubles:
// pp = [p,1] + [p,2];
// ...nothing else then the pp composite
// _x = pp->M();
// gives the invariant mass of the pp pair
// CAVEAT: Do _not_ use object name which are
// already assigned to particles in the data base, e.g.
// pi0=[pi0] is not allowed, use _pi0=[pi0] instead
// 
// Brackets: Objects can be nested:
// _x = ([p,1] + [p,2])->M();
// ...makes the same as the 2 lines above
//
// Build-in methods:
// M()              : Invariant mass of the object
// M2()             : Invariant mass^2 of the object
// Boost(obj)       : Boost into rest frame of "obj"
//   Example:
//   pp = p1->Boost([p,1] + [p,2]);
//   ...boost p1 in the rest frame of [p,1]+[p,2]
// Rot(obj)         : Rotate object such that obj
//                    would point to z-Direction
// Angle(obj)       : Opening angle between obj
// Theta()          : Theta of momentum
// Print()          : Dump 4momentum or double
//
// In addition, all browsable methods of the class PParticle
// can be used (functions without arguments)
//
// Example:
// obj->Px();
//
// Build-in functions:
// cos()
// fabs()
//
// In addition, the TFormula syntax can be used
//
// Example:
// val_new = (val + 1.2);
// val = (val * TMath::Pi());
//
// Conditions:
// if(arg)          : Interrupts the current chain if arg == 0
//
// A very complex example: Look for the helicity distribution
// of the e+e- pair in the eta Dalitz decay
// _eta=[eta]; _ep=[e+]; _em=[e-]; 
// _eta->Boost([p + p]); _em->Boost([p + p]); _ep->Boost([p + p]); 
// _ep->Rot(_eta); _em->Rot(_eta); _eta->Rot(_eta) ; 
// _ep->Boost(_eta); _em->Boost(_eta); dil=_ep+_em; 
// _ep->Rot(dil); _em->Rot(dil); dil->Rot(dil); _ep->Boost(dil); _em->Boost(dil) ; 
// s1= _ep->Theta(); _x = cos(s1)
//
//
//                    Author: I. Froehlich
//                    Written: 14.02.2008
//                    Revised: 
//
////////////////////////////////////////////////////////

//makeDataBase()->ListEntries(-1,1,"name,*batch_particle,*pid,*batch_value");

#include "PBatch.h"
#include "TString.h"
#include "PUtils.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "PParticle.h"
#include "PStaticData.h"

#include "TClass.h"
#include "TMethod.h"

#include <cmath>

Int_t PBatch::stack_num_pos=0;
Int_t PBatch::stack_num_batch[MAX_STACK_GOSUB], 
    PBatch::stack_num_bulk[MAX_STACK_GOSUB], 
    PBatch::stack_num_command[MAX_STACK_GOSUB];

PBatch& fBatch()
 {
   static PBatch* ans = new PBatch();
   return *ans;
 }

PBatch * makeGlobalBatch()
 {
   return &fBatch();
 }

PBatch::PBatch() {

    makeStaticData();

    command_pointer=last_command_pointer=0;
    method_pointer=0;

    fHisto1=NULL;
    fHisto2=NULL;
    fHisto3=NULL;

    batch_particle_param = makeDataBase()->GetParamTObj("batch_particle");
    batch_value_param= makeDataBase()->GetParamDouble("batch_value");
    pid_param= makeDataBase()->GetParamInt("pid");

    if (batch_particle_param<0) 
	batch_particle_param = 
	    makeDataBase()->MakeParamTObj("batch_particle", "PParticle storage for batch");

    if (batch_value_param<0) 
	batch_value_param = 
	    makeDataBase()->MakeParamDouble("batch_value", "Value storage for batch");

    //This is used for the "goto"
    num_command_param= makeDataBase()->GetParamInt("num_command");    
    if (num_command_param<0) 
	num_command_param= 
	    makeDataBase()->MakeParamInt("num_command", "Number of command for a label");
    num_batch_param= makeDataBase()->GetParamInt("num_batch");    
    if (num_batch_param<0) 
	num_batch_param= 
	    makeDataBase()->MakeParamInt("num_batch", "Number of batch object for a label");
    num_bulk_param= makeDataBase()->GetParamInt("num_bulk");    
    if (num_bulk_param<0) 
	num_bulk_param= 
	    makeDataBase()->MakeParamInt("num_bulk", "Number of bulk for a label");
    num_batch=num_bulk=-1;

    //Used for "formore"
    stream_default_pos_param = makeDataBase()->GetParamInt(STREAM_DEFAULT_POS);
    if (stream_default_pos_param<0)  
	stream_default_pos_param = makeDataBase()->MakeParamInt(STREAM_DEFAULT_POS,"Default position");
    stream_max_pos_param = makeDataBase()->GetParamInt(STREAM_MAX_POS);
    if (stream_max_pos_param<0)  
	stream_max_pos_param = makeDataBase()->MakeParamInt(STREAM_MAX_POS,"Max position in stream");

    batch_update_param = makeDataBase()->GetParamInt("batch_update");
    if (batch_update_param<0)
	makeDataBase()->MakeParamInt("batch_update", 
				     "If set this is a variable which must trigger an update"); 

    for (int i= 0;i< MAX_COMMAND_POINTER; i++) error_flag[i]=0;

    x = makeStaticData()->GetBatchValue("_x");
    y = makeStaticData()->GetBatchValue("_y");
    z = makeStaticData()->GetBatchValue("_z");
    eval_err_dumped=status=0;
    varlist = NULL;
    //pid_param=base->GetParamInt("pid");
}


Int_t PBatch::Execute(Int_t command_pos) {

    Int_t retval = kFALSE; //standard, can be overwritten for a foreach
    //cout << command_pos << ":" << command_pointer << endl;

    for (int i=command_pos;i<command_pointer;i++) {
	TObject * res=NULL;
	TObject * res2=NULL;
	TObject * res3=NULL;
	Double_t *val=NULL;
	Double_t *val2=NULL;
	Double_t *val3=NULL;
	Double_t *val4=NULL;
	//cout << "command:" << lst_command[i] << " num "<< i <<  endl;

	if (lst_command[i] == COMMAND_PLUS) {
	    //
	    //+ of pparticles or values
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    makeDataBase()->GetParamTObj (lst_key[2][i],batch_particle_param, &res2);
	    makeDataBase()->GetParamTObj (lst_key_a[i],batch_particle_param, &res3);
	
	    if ( res && res2 && res3) {
		PParticle * a1  = (PParticle * ) res;
		PParticle * a2  = (PParticle * ) res2;
		PParticle * r  = (PParticle * ) res3;
		*r = *a1;
		r->AddTmp(*a2);
		found = kTRUE;
	    } 
	
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    makeDataBase()->GetParamDouble (lst_key[2][i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val3);
	
	    if ( val && val2 && val3) {
		*val3=*val2 + *val;
		found = kTRUE;
	    }
	    if (!found) return retval;
	} else if (lst_command[i] == COMMAND_MINUS) {
	    //
	    //+ of pparticles or values
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    makeDataBase()->GetParamTObj (lst_key[2][i],batch_particle_param, &res2);
	    makeDataBase()->GetParamTObj (lst_key_a[i],batch_particle_param, &res3);
	
	    if ( res && res2 && res3) {
		PParticle * a1  = (PParticle * ) res;
		PParticle * a2  = (PParticle * ) res2;
		PParticle * r  = (PParticle * ) res3;
		*r = *a1 - *a2;
		found = kTRUE;
	    } 
	
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    makeDataBase()->GetParamDouble (lst_key[2][i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val3);
	
	    if ( val && val2 && val3) {
		*val3=*val - *val2;
		found = kTRUE;
	    }
	    if (!found) return found;
	} else if (lst_command[i] == COMMAND_EQUAL) {
	    //
	    // compare values (+/- 0.5)
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    makeDataBase()->GetParamDouble (lst_key[2][i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val3);
	
	    if ( val && val2 && val3) {
		if (fabs(*val - *val2)<0.5) *val3=1.;
		else *val3=0.;
		found = kTRUE;
	    }
	    if (!found) return retval;
	} else if (lst_command[i] == COMMAND_BOOST) { 
	    //
	    //boost pparticles
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a1  = (PParticle * ) res;
	    makeDataBase()->GetParamTObj (lst_key[2][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a2  = (PParticle * ) res;
	
	    a1->Boost(-a2->BoostVector()) ;
	} else if (lst_command[i] == COMMAND_ANGLE ) { 
	    //
	    //angle between particle tracks
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a1  = (PParticle * ) res;
	    makeDataBase()->GetParamTObj (lst_key[2][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a2  = (PParticle * ) res;
	    
	    Double_t angle= a1->Vect().Angle(a2->Vect());
	    
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) {
		Warning("Execute","->Angle(): Result value not found");
		return kFALSE;
	    }
	    *val= angle;
	} else if (lst_command[i] == COMMAND_ROT) { 
	    //
	    //rotate
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a1  = (PParticle * ) res;
	    makeDataBase()->GetParamTObj (lst_key[2][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a2  = (PParticle * ) res;
	    double tmp_phi=a2->Phi();
	    double tmp_theta=a2->Theta();

	    a1->RotateZ(-tmp_phi);
	    a1->RotateY(-tmp_theta);	 
	} else if (lst_command[i] == COMMAND_IS) {
	    //
	    // '='
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    makeDataBase()->GetParamTObj (lst_key_a[i],batch_particle_param, &res2);
	    
	    if ( res && res2) {
		PParticle * r  = (PParticle * ) res2;
		PParticle * a1  = (PParticle * ) res;
		//store PID, if old PParticle had already an pid, keep it
		Int_t pid = r->ID();
		Double_t w = r->W();
		*r=*a1;
		if (pid) {
		    r->SetID(pid);
		    r->SetW(w);
		}
		found = kTRUE;
	    } 

	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);

	    if ( val && val2) {
		*val2=*val;
		found = kTRUE;
		//		cout << *val2 << endl;
	    }
	    
	    if (!found) return retval;
	    Int_t *update;
	    //cout << "checking :" << lst_key_a[i] << endl;
	    if (makeDataBase()->GetParamInt (lst_key_a[i],batch_update_param, &update)) {
		//cout << "found :" << lst_key_a[i] << endl;
		if (*update == 1) {
		    //cout << "Update forced" << endl;
		    locnum_command = i + 1;
		    return kUPDATE;
		}
	    }

	} else if (lst_command[i] == COMMAND_MASS2) {
	    //
	    // M2()
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a1  = (PParticle * ) res;
	    
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) return retval;
	    *val= a1 -> M2();
	} else if (lst_command[i] == COMMAND_MASS) {
	    //
	    // M()
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) {
		return retval;
	    }
	    PParticle * a1  = (PParticle * ) res;
	    
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) {
		Warning("Execute","->M(): Result value %s not found",
			makeDataBase()->GetName(lst_key_a[i]));
		return retval;
	    }
	    *val= a1 -> M();
	} else if (lst_command[i] == COMMAND_THETA) {
	    //
	    // Theta()
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle * a1  = (PParticle * ) res;
	    
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) return retval;
	    *val= (a1 -> Theta());
	} else if (lst_command[i] == COMMAND_INTERNAL) {
	    //
	    // Pointer to PUtils
	    //
	    TObject *a1;
	    if (flag_command_int[i] == 0) {
		makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
		if (!res) return retval;
		a1  = (PParticle * ) res;
	    }
	    else {
		a1 = makePUtilsREngine();
	    }

	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) return retval;

 	    Double_t *argval;

	    //First, we have to set the variables
	    methods[lst_command_int[i]]->ResetParam();
	    for (int j=0;j<4;j++) {
		if (methods_arg_flags[j][lst_command_int[i]] == METHOD_RETURN_DOUBLE) {
		    makeDataBase()->GetParamDouble (lst_key[j+2][i],batch_value_param, &argval); //j=1 is object
		    if (!argval) {
			Warning("Execute","Double argument for key %i not found",lst_key[j+2][i]);
			return retval;
		    }
		    methods[lst_command_int[i]]->SetParam(*argval);
		    //cout << j << ":" << *argval << endl;
		}  else if (methods_arg_flags[j][lst_command_int[i]] == METHOD_RETURN_INT) {
		    makeDataBase()->GetParamDouble (lst_key[j+2][i],batch_value_param, &argval); //j=1 is object
		    if (!argval) {
			Warning("Execute","Int argument for key %i not found",lst_key[j+2][i]);
			return retval;
		    }
		    methods[lst_command_int[i]]->SetParam( (Long_t)*argval);
		} 
	    }

	    //Execute internal
	    if (methods_flags[lst_command_int[i]] == METHOD_RETURN_DOUBLE) {
		methods[lst_command_int[i]]->Execute(a1,*val);
		//cout << "Double_t meth " << method_name[lst_command_int[i]] << " called, result " << *val << endl;
	    }
	    else if (methods_flags[lst_command_int[i]] == METHOD_RETURN_INT ) { //INT
		Long_t ret;
		methods[lst_command_int[i]]->Execute(a1,ret);
		//cout << "Int_t meth " << method_name[i] << " called, result " << ret << endl;
		(*val)=(Double_t) ret;
	    } else { //VOID
		methods[lst_command_int[i]]->Execute(a1);
		//cout << "void meth " << method_name[lst_command_int[i]] << " called" << endl;
		(*val)=0;
	    }

	} else if (lst_command[i] == COMMAND_PVALUE) {
	    //
	    // Reads the PValue ("PParticle.val")
	    // or the database entry
	    //
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) {
		if (!error_flag[i]) {
		    error_flag[i]=1;
		    Error("Execute","Result object %s not found",
			      makeDataBase()->GetName(lst_key_a[i]));
		}
		return retval;
	    }	    
	    if (lst_key[3][i] == 0) { //PValue
		makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
		if (!res) return retval;
		PParticle * a1  = (PParticle * ) res;
		

		if (!a1->GetValue(lst_key[2][i],val)) {
		    if (!error_flag[i]) {
			error_flag[i]=1;
			Error("Execute","PValue %i not set",lst_key[2][i]);
		    }
		    return retval;
		} 
	    } else { //database entry
		if (lst_key[3][i]<0) {
		    makeDataBase()->GetParamDouble (lst_key[1][i],(- lst_key[3][i]) - 1, &val2);
		    if (!val2) {
			Error("Execute","Connection to data base for double %i failed",(- lst_key[3][i]) - 1);
			return retval;
		    }
		    *val = *val2;		    
		} else {
		    Int_t *intval;
		    makeDataBase()->GetParamInt (lst_key[1][i],lst_key[3][i] - 1, &intval);
		    *val = (Double_t)*intval;
		} 

	    }
	} else if (lst_command[i] == COMMAND_PFORMULA) {
	    //
	    // Pointer to PFormula
	    //
	    Double_t pars[MAX_COMMAND_OPTIONS];
	    for (int j=0;j<lst_options_counter[i];j++) {
		makeDataBase()->GetParamDouble (lst_key[j+1][i],batch_value_param, &val);
		if (!val) return retval;
		pars[j]=*val;
	    }
	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);
	    if (!val) return retval;
	    lst_form[i]->SetParameters(pars);
	    *val = lst_form[i]->Eval(1.);
	} else if (lst_command[i] == COMMAND_PRINT) {
	    //
	    // Dumps the value(s)
	    //
	    makeDataBase()->GetParamTObj (lst_key[1][i],batch_particle_param, &res);
	    if (res) {
		cout << "*******: " << makeDataBase()->GetName(lst_key[1][i]) << endl;
		PParticle * a1  = (PParticle * ) res;
		//a1->Dump();
		a1->Print();
	    }
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    if (val) {
		cout << "*******: " << makeDataBase()->GetName(lst_key[1][i]) << endl;
		cout << *val << endl;
	    }

	} else if (lst_command[i] == COMMAND_COS) {
	    //
	    // Internal cos()
	    //
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    if (!val) {
		return retval;
	    }
 	    Double_t myres=cos(*val);
	    
 	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);

 	    if (!val) return retval;
 	    *val= myres;
	} else if (lst_command[i] == COMMAND_FABS) {
	    //
	    // Internal fabs()
	    //
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    if (!val) {
		return retval;
	    }
 	    Double_t myres=fabs(*val);
	    
 	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);

 	    if (!val) return retval;
 	    *val= myres;
	} else if (lst_command[i] == COMMAND_EVAL) {
	    //
	    // Eval the attached histogram and returns the value
	    //
	    if (!x) {
		return retval;
	    }
	    Double_t myres=0.;
	    if (fHisto1) {
		int bin = fHisto1->FindBin(*x);
		myres = fHisto1->GetBinContent(bin);
                //dw= hist3D->GetBinError(bin);
	    } else if (fHisto2) {
		int bin = fHisto2->FindBin(*x,*y);
		myres = fHisto2->GetBinContent(bin);
	    } else if (fHisto3) {
		int bin = fHisto3->FindBin(*x,*y,*z);
		myres = fHisto3->GetBinContent(bin);
	    } else {
		if (!eval_err_dumped) {
		    eval_err_dumped=1;
		    Error("Execute","Eval() called, but no object present");
		}
		return retval;
	    }

 	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);

 	    if (!val) return retval;
 	    *val= myres;
	} else if (lst_command[i] == COMMAND_IF) {
	    //
	    // if (....)
	    //
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    if (!val) {
		Warning("Execute","if: argument not found");
		return retval;
	    }
 	    Double_t myres=*val;
	    
 	    makeDataBase()->GetParamDouble (lst_key_a[i],batch_value_param, &val);

 	    if (!val) return retval;

	    //	    cout << *val << endl;
 	    *val= myres;
	    if (fabs(*val) == 0) return retval;

	} else if (lst_command[i] == COMMAND_P3M) {
	    //
	    // Constructor for PParticles (with mass)
	    //
	    makeDataBase()->GetParamTObj (lst_key_a[i],batch_particle_param, &res);
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    makeDataBase()->GetParamDouble (lst_key[2][i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key[3][i],batch_value_param, &val3);
	    makeDataBase()->GetParamDouble (lst_key[4][i],batch_value_param, &val4);

	    if (!res || !val || !val2 || !val3 || !val4) {
		return kTRUE;
	    }
	    
	    ((PParticle * ) res)->
		SetPxPyPzE((*val),(*val2),(*val3),
			   sqrt((*val)*(*val)+(*val2)*(*val2)+(*val3)*(*val3)+(*val4)*(*val4)));
	    
	} else if (lst_command[i] == COMMAND_P3E) {
	    //
	    // Constructor for PParticles (with energy)
	    //
	    makeDataBase()->GetParamTObj (lst_key_a[i],batch_particle_param, &res);
	    makeDataBase()->GetParamDouble (lst_key[1][i],batch_value_param, &val);
	    makeDataBase()->GetParamDouble (lst_key[2][i],batch_value_param, &val2);
	    makeDataBase()->GetParamDouble (lst_key[3][i],batch_value_param, &val3);
	    makeDataBase()->GetParamDouble (lst_key[4][i],batch_value_param, &val4);

	    if (!res || !val || !val2 || !val3 || !val4) {
		return kTRUE;
	    }
	    
	    ((PParticle * ) res)->
		SetPxPyPzE((*val),(*val2),(*val3),(*val4));
	    
	} else if (lst_command[i] == COMMAND_GOTO) {
	    //
	    // Goto
	    //
	    Bool_t labelfound = kTRUE;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i],num_command_param, &locnum_command)) 
		labelfound = kFALSE;
	    makeDataBase()->GetParamInt(lst_key[1][i],num_batch_param, &locnum_batch);
	    makeDataBase()->GetParamInt(lst_key[1][i],num_bulk_param, &locnum_bulk);
	    if ((locnum_command>=0) && (labelfound) && (locnum_batch < 0)) { //local batch
		i = locnum_command-1;
	    }
	    else if((locnum_command>=0) && (labelfound)) { //go back to PProjector
		return kGOTO;
	    }
	    if (!labelfound)
		Warning("Execute","Label '%s' not found in goto command",makeDataBase()->GetName(lst_key[1][i]));
	} else if (lst_command[i] == COMMAND_GOSUB) {
	    //
	    // Gosub
	    //
	    //space left on stack?
	    if (stack_num_pos == MAX_STACK_GOSUB) {
		Warning("Execute","Cannot call '%s': stack full (loop?)",makeDataBase()->GetName(lst_key[1][i]));
		return kFALSE;
	    }

	    Bool_t labelfound = kTRUE;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i],num_command_param, &locnum_command)) 
		labelfound = kFALSE;
	    makeDataBase()->GetParamInt(lst_key[1][i],num_batch_param, &locnum_batch);
	    makeDataBase()->GetParamInt(lst_key[1][i],num_bulk_param, &locnum_bulk);
	    if((locnum_command>=0) && (labelfound)) {
		stack_num_batch[stack_num_pos]=num_batch;
		stack_num_bulk[stack_num_pos]=num_bulk;
		stack_num_command[stack_num_pos]=i;
		stack_num_pos++;
		if (locnum_batch < 0) { //local batch
		    i = locnum_command-1;
		}
		else  { //go back to PProjector
		    return kGOTO;
		}
	    }
	    if (!labelfound)
		Warning("Execute","Label '%s' not found in goto command",makeDataBase()->GetName(lst_key[1][i]));
	} else if (lst_command[i] == COMMAND_RETURN) {
	    //
	    // Returns from gosub
	    //
	    //something on stack?
	    if (stack_num_pos == 0) {
		Warning("Execute","Cannot return: stack empty");
		return kFALSE;
	    }
	    stack_num_pos--;
	    locnum_command=stack_num_command[stack_num_pos]+1;
	    locnum_batch=stack_num_batch[stack_num_pos];
	    locnum_bulk=stack_num_bulk[stack_num_pos];
	    if (locnum_batch < 0) { //local batch
		i = locnum_command-1;
	    }
	    else  { //go back to PProjector
		return kGOTO;
	    }
	} else if (lst_command[i] == COMMAND_EXIT) {
	    //
	    // Exit
	    //
	    locnum_command=999999999;
	    locnum_batch=999999999;
	    locnum_bulk=num_bulk;
	    return kGOTO;
	} else if ((lst_command[i] == COMMAND_FORMORE) || (lst_command[i] == COMMAND_FOREACH)) {
	    //
	    // Build-in loops for PProjector
	    //
	    Int_t  stream_max;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i],stream_max_pos_param, &stream_max)) {
		//cout << "particle not in stream" << endl;
		return kFALSE; //particle not in stream, abort
	    }

	    Int_t  *stream_default;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i],stream_default_pos_param, &stream_default)) {
		Int_t * dummy = new Int_t(1);
		stream_default = dummy;
		makeDataBase()->SetParamInt(lst_key[1][i] ,STREAM_DEFAULT_POS, dummy);
	    } else {
		(*stream_default)++;
	    }
	    //cout << lst_key[1][i] << ":" << stream_max << ":" << (*stream_default) << endl;
	    if (lst_command[i] == COMMAND_FORMORE) {
		if ( (*stream_default) >= stream_max) { //Maximum reached, abort
		    //reset value:
		    (*stream_default)=0;
		    //cout << "max reached" << endl;
		    return kFALSE; 
		} 
	    } else { //FOREACH
		if ( (*stream_default) > stream_max) { //Maximum reached, abort
 		    //reset value:
 		    (*stream_default)=0;
 		    //cout << "max reached" << endl;
 		    return kFOREACHEND; 
 		} else {
 		    retval = kFOREACH;
 		}
		locnum_command = i; //Jump to the "foreach" command
	    }

	} else if (lst_command[i] == COMMAND_ECHO) {
	    //
	    // echo blabla
	    //
	    char puffer[1000]; //I hope there will be never a string longer...
	    unsigned int puffer_pointer=0;
	    cout << "<PBatch> "; 
	    Int_t seek_mode=0;
	    for (unsigned int j=0;j<=strlen(echo_string[i]); j++) {
		char current = (*(echo_string[i]+j));
		if (current == '$')
		    seek_mode =1;
		else {
		    //if (seek_mode && (((*(echo_string[i]+j)) == ' ')  ||  (*(echo_string[i]+j)) == '\0') ) {
		    if (seek_mode && ((!(current=='#') and !(current=='_') and !isalnum(current)))) {
			//end
			seek_mode =0;
			puffer[puffer_pointer]='\0';
			//cout << ":"<< puffer << ":" << endl;
			Double_t *x=makeStaticData()->GetBatchValue(puffer,0);
			if (x) cout << *x;
			else cout << "[Error: " << puffer << " not found]";
			    
			puffer_pointer=0;
		    }
		    if (seek_mode) {
			puffer[puffer_pointer]=current;
			//cout << (int) puffer[puffer_pointer] << endl;
			if (puffer_pointer < 1000) puffer_pointer++;
			else {
			    cout << "[Error: variable too large]";
			    seek_mode = 0;
			}
		    }
		    
		    if (!seek_mode) cout << current;
		}
	    }
	    cout << endl;
	}
    }

    if (retval == kFOREACH) return kFOREACH;
    return kTRUE;

}

Bool_t PBatch::AddCommand(char * command) {

    PUtils::remove_spaces(&command);
    if (strlen(command)==0) return kTRUE;
    Bool_t has_something=kFALSE;

    for (unsigned int i=0;i<strlen(command);i++) {
	if ((command[i] != ' ') && (command[i] != ';'))
	    has_something=kTRUE;
    }
    if (!has_something) return kTRUE;

    //cout << "AddCommand:" << command << endl;
    //First check if we have a composite command

    int is_composite=0;

    for (UInt_t i=0; i<strlen(command); i++) {
	if (command[i]==';') is_composite=1;
    }

    if (is_composite) {
	char *array[200];
	Int_t array_s=200; //max products
	PUtils::Tokenize(command, ";", array, &array_s);

	for (int i=0; i<array_s; i++) {
	    AddCommand(array[i]);
	}
	return kTRUE;
    }

    //single command


    //echo has the highest priority:
	if (!strncmp(command,"echo",4)) {
// 	    key_a = makeStaticData()->
// 		MakeDirectoryEntry("batch_objects",command3);
	    if (strlen(command) > 4) {
		char * dummy = new char[strlen(command)-2];
		strncpy(dummy,command+4,strlen(command)-4);
		dummy[strlen(command)-4]='\0';
		//cout << "echo:" << dummy << endl;
		PUtils::remove_spaces(&dummy);
		echo_string[command_pointer]=dummy;
	    } else echo_string[command_pointer]=new char('\0');
	    AddCommand(COMMAND_ECHO,-1,-1,-1);
	    return kTRUE;
	}    

    //first, check for a label
    for (UInt_t i=1; i<strlen(command); i++) {
	
	Bool_t ret=kFALSE;
	if (command[i]==':') {
	    if (((i == (strlen(command))) ||  (command[i+1] != ':'))
		&&  ((i == 0) ||  (command[i-1] != ':')) ){
		//add new label
		char *label_name = new char[i+1];
		strncpy(label_name,command,i);
		label_name[i] = '\0';
		//cout << "found label: " << label_name  << endl;
		int key_a = makeStaticData()->
		    MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command);
		if (key_a>=0) {
		    Int_t * delme =  new Int_t(command_pointer);
		    AddCommand(COMMAND_LABEL,key_a,-1,-1);
		    ret = AddCommand(command+i+1);
		    int key_l = GetKey(label_name,1,1);
		    if (key_l  > 0) {
			//cout << "Label key is : " << key_l << endl;
			makeDataBase()->SetParamInt(key_l ,"num_command", delme);
			//cout <<  "num_command:  " << *delme << endl;
			delme =  new Int_t(num_batch);
			makeDataBase()->SetParamInt(key_l ,"num_batch", delme);
			//cout <<  "num_batch:  " << *delme << endl;
			delme =  new Int_t(num_bulk);
			makeDataBase()->SetParamInt(key_l ,"num_bulk", delme);
			//cout <<  "num_bulk:  " << *delme << endl;
		    }
		}
		return ret;
	    }
	}
    }


    //second thing: I check for a "=" operator

    int is_operator=0;
    Bool_t found=kFALSE;
    char *prod[2];
    Int_t prod_s=2; //max 2 products
    int key_a,key1,key2;

    for (UInt_t i=0; i<strlen(command); i++) {
	if (command[i]=='=') is_operator++;
    }

    //let us copy the command first, because it might be modified
    //after the pointer is given to the data base
    char * command2=new char[strlen(command)+1];
    strcpy(command2,command);
    command=command2;
    char * command3=new char[strlen(command)+1];
    strcpy(command3,command);
    Int_t copy_ctor = 0;
    if (!is_operator) {
	//first the functions
	if (!strncmp(command,"goto",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    char * dummy = new char[strlen(command)-3];
	    strncpy(dummy,command+4,strlen(command)-4);
	    dummy[strlen(command)-4]='\0';
	    key1 = GetKey(dummy,1,1);
	    //cout << "found goto to key " << key1 << endl;
	    if (AddCommand(COMMAND_GOTO,key_a,key1,-1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"gosub",5)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    char * dummy = new char[strlen(command)-3];
	    strncpy(dummy,command+5,strlen(command)-5);
	    dummy[strlen(command)-5]='\0';
	    key1 = GetKey(dummy,1,1);
	    //cout << "found goto to key " << key1 << endl;
	    if (AddCommand(COMMAND_GOSUB,key_a,key1,-1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"return",6)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    if (AddCommand(COMMAND_RETURN,key_a,key1,-1)) {
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"exit",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    if (AddCommand(COMMAND_EXIT,key_a,key1,-1)) {
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"formore",7)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    char * dummy = new char[strlen(command)-3];
	    strncpy(dummy,command+7,strlen(command)-4);
	    dummy[strlen(command)-7]='\0';
	    if (strcmp(dummy,"(*)"))
		key1 = GetKey(dummy,1,1);
	    else
		key1 = GetKey((char *)"dummy",1,1);
	    //cout << "found formore to key " << key1 << endl;
	    if (AddCommand(COMMAND_FORMORE,key_a,key1,-1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"foreach",7)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    char * dummy = new char[strlen(command)-3];
	    strncpy(dummy,command+7,strlen(command)-4);
	    dummy[strlen(command)-7]='\0';
	    if (strcmp(dummy,"(*)"))
		key1 = GetKey(dummy,1,1);
	    else
		key1 = GetKey((char *)"dummy",1,1);
	    if (AddCommand(COMMAND_FOREACH,key_a,key1,-1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"Eval()",6)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    if (AddCommand(COMMAND_EVAL,key_a,key1,-1)) {
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}

	if (!strncmp(command,"cos(",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    key1 = GetKey(command+3,1,-1);
	    if (AddCommand(COMMAND_COS,key_a,key1,-1)) {
		
		//final thing is to create the result
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"fabs(",5)) {

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);

	    key1 = GetKey(command+4,1,-1);
	    if (AddCommand(COMMAND_FABS,key_a,key1,-1)) {
		//final thing is to create the result
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}

	if (!strncmp(command,"if(",3) || !strncmp(command,"if ",3)) {

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	   
	    key1 = GetKey(command+2,1,-1);
	    if (AddCommand(COMMAND_IF,key_a,key1,-1)) {
		//final thing is to create the result
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}

	if (!strncmp(command,"P3M(",4)) {
	    char *prodx[4];
	    Int_t prodx_s=4; //max 4 products

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);

	    PUtils::Tokenize(command+4, ",", prodx, &prodx_s);
	    if (prodx_s<4) {
		Warning("AddCommand","P3M needs 4 arguments");
	    }

	    *(prodx[3]+strlen(prodx[3])-1) = '\0'; //remove trailing )

	    key1 = GetKey(prodx[0],1,-1);
	    key2 = GetKey(prodx[1],1,-1);
	    int key3 = GetKey(prodx[2],1,-1);
	    int key4 = GetKey(prodx[3],1,-1);

	    if (AddCommand(COMMAND_P3M,key_a,key1,key2,key3,key4)) {
		//adding objects, if needed
		TObject * delme = NULL;
		makeDataBase()->GetParamTObj (key_a,batch_particle_param, &delme);
		if (!delme) {
		    delme =  (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj (key_a ,"batch_particle", delme);
		}
		Double_t *val=NULL;
		makeDataBase()->GetParamDouble (key1,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key1 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key2,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key2 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key3,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key3 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key4,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key4 ,"batch_value", val);
		}
		
		return kTRUE;
	    }
	    return kFALSE;
	}


	if (!strncmp(command,"P3E(",4)) {
	    char *prodx[4];
	    Int_t prodx_s=4; //max 4 products

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);

	    PUtils::Tokenize(command+4, ",", prodx, &prodx_s);
	    if (prodx_s<4) {
		Warning("AddCommand","P3E needs 4 arguments");
	    }

	    *(prodx[3]+strlen(prodx[3])-1) = '\0'; //remove trailing )

	    key1 = GetKey(prodx[0],1,-1);
	    key2 = GetKey(prodx[1],1,-1);
	    int key3 = GetKey(prodx[2],1,-1);
	    int key4 = GetKey(prodx[3],1,-1);

	    if (AddCommand(COMMAND_P3E,key_a,key1,key2,key3,key4)) {
		//adding objects, if needed
		TObject * delme = NULL;
		makeDataBase()->GetParamTObj (key_a,batch_particle_param, &delme);
		if (!delme) {
		    delme =  (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj (key_a ,"batch_particle", delme);
		}
		Double_t *val=NULL;
		makeDataBase()->GetParamDouble (key1,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key1 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key2,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key2 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key3,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key3 ,"batch_value", val);
		}
		val=NULL;
		makeDataBase()->GetParamDouble (key4,batch_value_param, &val);
		if (!val) {
		    val= new Double_t(0.);
		    makeDataBase()->SetParamDouble (key4 ,"batch_value", val);
		}

		return kTRUE;
	    }
	    return kFALSE;
	}





	//Methods below here...
	//cout << "adding a method" << endl;

	//we start from right to left. Only methods which are
	//not embedded in brackets are taken into account
	Int_t dot_version=0,prev_dot_version=0;
	Int_t method_position=-1;
	Int_t brack_counter=0,total_brack_counter=0;
	for (unsigned int i=0;i<(strlen(command)-1);i++) {
	    if (command[i]=='(') {
		brack_counter++;
		if (prev_dot_version)
		    total_brack_counter++;
	    }
	    if (command[i]==')') brack_counter--;
	    if (brack_counter==0) {
		if (command[i]=='.' && isalpha(command[i+1])) {
		    method_position=i;
		    dot_version=1;
		}
		if (command[i]=='-' && command[i+1]=='>') {
		    method_position=i;
		    dot_version=2;
		}
		if (total_brack_counter && prev_dot_version && ((command[i]!=' ') || (command[i]!=';')))
		    dot_version=0; //must be something different
		prev_dot_version = dot_version;
	    }
	}

	if ((dot_version>0) && (!found)) {
	    //catch the cases where comparators are used
	    for (unsigned int i=0;i< (unsigned int)method_position;i++) {
		if ((command[i]=='~') || (command[i]=='=') || (command[i]=='<')  || (command[i]=='>'))
		    dot_version=0;
	    }
	}

	if ((dot_version>0) && (!found)) {
	    //cout << "adding a method" << endl;
	    char *prodx[2];
	    prodx[0] = new char[method_position+1];
	    strncpy(prodx[0], command, method_position);
	    prodx[0][method_position]='\0';
	    
 	    key_a = makeStaticData()->
 		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command);

	    if (dot_version==1) {
		prodx[1] = new char[strlen(command)-method_position];
		strncpy(prodx[1], command+method_position+1, strlen(command)-method_position-1);
		prodx[1][strlen(command)-method_position-1]='\0';
		
	    } else {
		prodx[1] = new char[strlen(command)-method_position-1];
		strncpy(prodx[1], command+method_position+2, strlen(command)-method_position-2);
		prodx[1][strlen(command)-method_position-2]='\0';
		
	    } 

	    if (!strcmp(prodx[1],"M2()")) {
		key1 = GetKey(prodx[0],1,-1);
		if (key1 < 0) return kFALSE;
		if (AddCommand(COMMAND_MASS2,key_a,key1,-1)) 
		    found=kTRUE;
	    }
	    if (!strcmp(prodx[1],"M()")) {
		key1 = GetKey(prodx[0],1,-1);
		if (AddCommand(COMMAND_MASS,key_a,key1,-1)) 
		    found=kTRUE;
	    }
	    if (!strcmp(prodx[1],"Theta()")) {
		key1 = GetKey(prodx[0],1,-1);
		if (AddCommand(COMMAND_THETA,key_a,key1,-1)) 
		    found=kTRUE;
	    }
	    if (!strcmp(prodx[1],"Print()")) {
		key1 = GetKey(prodx[0],1,-1);
		if (AddCommand(COMMAND_PRINT,key_a,key1,-1)) 
		    found=kTRUE;
	    }
// 	    if (!strcmp(prodx[1],"Parent()")) {
// 		key1 = GetKey(prodx[0],1,0);
// 		if (AddCommand(COMMAND_PARENT,key_a,key1,-1)) 
// 		    found=kTRUE;
// 	    }
// 	    if (!strcmp(prodx[1],"ID()")) {
// 		key1 = GetKey(prodx[0],1,0);
// 		if (AddCommand(COMMAND_ID,key_a,key1,-1)) 
// 		    found=kTRUE;
// 	    }

	    //with arguments:
	    if (!strncmp(prodx[1],"Boost(",6)) {
		key1 = GetKey(prodx[0],1,-1);
		key2 = GetKey(prodx[1]+5,1,-1);
		if (AddCommand(COMMAND_BOOST,key_a,key1,key2)) 
		    found=kTRUE;
	    }
	    if (!strncmp(prodx[1],"Angle(",6)) {
		key1 = GetKey(prodx[0],1,-1);
		key2 = GetKey(prodx[1]+5,1,-1);
		if (AddCommand(COMMAND_ANGLE,key_a,key1,key2)) 
		    found=kTRUE;
	    }
	    if (!strncmp(prodx[1],"Rot(",4)) {
		key1 = GetKey(prodx[0],1,-1);
		key2 = GetKey(prodx[1]+3,1,-1);
		if (AddCommand(COMMAND_ROT,key_a,key1,key2)) 
		    found=kTRUE;
	    }
	  
	    //If there are no (), it could be an internal PValue
	    Int_t found_brackets=0;
	    for (unsigned int j=0;j<strlen(prodx[1]);j++) {
		if (((prodx[1])[j] == ')') || ((prodx[1])[j] == '('))
		    found_brackets++;
	    }

	    if (!found_brackets) {
		int val_id = pdummy.StringToValueID(prodx[1]);
		int database_id = 0;
		if (val_id<0) {
		    //check for database entry
		    //positive: param_int (+1)
		    //negative: param_double (-1)
		    database_id = makeDataBase()->GetParamInt(prodx[1]) -
			makeDataBase()->GetParamDouble(prodx[1]);
		    if (database_id == 0)
			Error("AddCommand","[%s] The value %s is unknown",command,prodx[1]);
		} 
		if ((val_id>=0) || database_id){		    
		    if (strcmp(prodx[0],"*")==0) { //dummy *
			key1 = GetKey((char *)"dummy",1,0);
		    } else {
			key1 = GetKey(prodx[0],1,0);
		    }
		    if (key1 >=0) {
			AddCommand(COMMAND_PVALUE,key_a,key1,val_id,database_id);
			found = kTRUE;
		    } 
		}
	    }

	    //Try to get build-in method
	    //only without arguments

	    if (!found) {
		Int_t handle = GetMethodHandle(prodx[1]);
		if (handle>=0) {		
		    key1 = GetKey(prodx[0],1,-1);
		    if (AddCommand(COMMAND_INTERNAL,key_a,key1,arg1,arg2,arg3,arg4))  {
			found=kTRUE;
			lst_command_int[command_pointer-1]=handle;
			flag_command_int[command_pointer-1]=0;
		    }
		}
	    }
	    //Make another try using the PFormula

	    if (!found) {
		found = EvalPFormula(command);
		if (found) return kTRUE;
	    }
	    
	    
	    if (!found) {
		Error("AddCommand","[%s] The method ->%s is unknown",command,prodx[1]);
		return kFALSE;
	    } else {
		
		
// 		Int_t type=CheckObjectType(key1);
		
// 		if (type<0) {
// 		    Error("AddCommand","[%s] Object not found: %s",command,
// 			  makeDataBase()->GetName(key1));
// 		    return kFALSE;
// 		}
		//can be done later (e.g. file input)
		
		//final thing is to create the result
		Double_t * delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
		
		return kTRUE;
		
	    }
	}
	else {
	    //Warning("AddCommand","[%s] Unknown single command",command);
	    key_a=makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command);
	    prod[1]=command; //continue and try the best...
	    //return kFALSE; //not yet done
	}
    }
    else if (is_operator>1) {
	//equal
	if (strstr(command, "==")) {
	    for (UInt_t i=0; i<(strlen(command)-1); i++) {
		if (command[i]=='=' && command[i+1]=='=' ) {
		    Info("AddCommand","[%s] the '==' is not exact for doubles - lets use '~'",
			 command);
		    command[i]= '~';
		    command[i+1]= ' ';

		}
	    }
	} else {
	    
	    Error("AddCommand","[%s] Too many ='s",command);
	    return kFALSE;
	}
    } else {
	PUtils::Tokenize(command, "=", prod, &prod_s);

	key_a=GetKey(prod[0],2,2);
	copy_ctor = 1;
    }
    
    found=kFALSE;
    
    if (!is_operator) {
	//look for the arguments
	if (CheckAndSplit(prod[1],'+',&key1,&key2)) {
	    if (AddCommand('+',key_a,key1,key2))
		found=kTRUE;
	}
	else if (CheckAndSplit(prod[1],'-',&key1,&key2)) {
	    if (AddCommand('-',key_a,key1,key2))
		found=kTRUE;
	}
	else if (CheckAndSplit(prod[1],'~',&key1,&key2)) {
	    if (AddCommand(COMMAND_EQUAL,key_a,key1,key2))
		found=kTRUE;
	}
    }

    if (!found && copy_ctor) {
	//copy ctor
	//cout << "copy ctor:" << prod[1] << endl;

	key1 = GetKey(prod[1],1,-1);
	

	key2 = -1;
	if (is_operator == 1) {
	    if (AddCommand(COMMAND_IS,key_a,key1,key2))
		found=kTRUE;
	} else
	    if (AddCommand(COMMAND_EQUAL,key_a,key1,key2))
		found=kTRUE;
//	GetKey(prod[0],0,1);


	Bool_t makenew=kFALSE;

	//is lvalue one of the "switching" functions?
	for (UInt_t i=0;i<strlen(prod[0]);i++) 
	    if (*(prod[0]+i) == ')' || *(prod[0]+i) == '(') {
		makenew=kTRUE;
	    }

	if (makenew) AddCommand(prod[0]);
    } else if (!found) {
	//if nothing helps, try to use the PFormula....
	found = EvalPFormula(command);
	if (found) return kTRUE;
	//....or the wrapper to PUtils
	Int_t handle = GetMethodHandle(command, 1);
	if (handle>=0) {		
	    //cout << "found PUtils with " << command3<< endl;
	    key1 = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command3);
	    Double_t * delme =  new Double_t(0.);
	    makeDataBase()->SetParamDouble(key1 ,"batch_value", delme);
	    //key1 = GetKey(command3,1,1);
	    if (AddCommand(COMMAND_INTERNAL,key_a,key1,arg1,arg2,arg3,arg4))  {
		found=kTRUE;
		lst_command_int[command_pointer-1]=handle;
		flag_command_int[command_pointer-1]=1;
	    }
	    return kTRUE;
	}
	
    }

    if (found) {
	//further look to the arguments, are they file-input???

	Int_t result_type=0;

	Int_t num_args=1;
	if (key2>-1) num_args=2;

	for (int pos=0;pos<num_args;pos++) {
	    //loop over 2 args
	    Int_t key=key1;
	    if (pos==1) key=key2;
	    
	    Int_t *ii;

	    if (key>-1) {

		if (makeDataBase()->GetParamInt (key, pid_param, &ii)) {
		    
		    result_type++;
		    
		} else {
		    //have to check if objects are there and what kind they are
		    //4momentum?
		    
		    Int_t type=CheckObjectType(key);
		    
		    if (type<0) {
			Error("AddCommand","[%s] Object not identified: %s",command,
			      makeDataBase()->GetName(key));
			return kFALSE;
			//		    AddCommand(makeDataBase()->GetName(key));
		    }
		    
		    if (type==IS_OBJECT) result_type++;
		    else result_type--;
		    
		}
	    }
 	    
	}

	if (result_type >0) {
	    TObject * delme;
	    if (!makeDataBase()->GetParamTObj (key_a,"batch_particle",&delme)) {
		TObject * delme =  (TObject *) (new PParticle(0,0,0,0));
		makeDataBase()->SetParamTObj (key_a ,"batch_particle", delme);
	    }
	    
	} else {
	    Double_t * delme;
	    if (!makeDataBase()->GetParamDouble (key_a,"batch_value",&delme)) {
		delme = new Double_t (0.);
		makeDataBase()->SetParamDouble (key_a ,"batch_value", delme);
		//Error("AddCommand","[%s] Result double not yet impl.",command);
	    }
	}
	

	return kTRUE;

    } //end found operator


    Error("AddCommand","[%s] Syntax error",command);

    return kFALSE;
}

void PBatch::AddSpacePlaceholder(char * command) {
    int numbra=0;
    for (unsigned int i=0;i<strlen(command);i++) {
	if (command[i] == '[') numbra++;
	if (command[i] == ']') numbra--;
	if (numbra && (command[i] == ' ') ) command[i] ='#';
    }
}

void PBatch::RemoveSpacePlaceholder(char * command) {
    int numbra=0;
    for (unsigned int i=0;i<strlen(command);i++) {
	if (command[i] == '[') numbra++;
	if (command[i] == ']') numbra--;
	if (numbra && (command[i] == '#') ) command[i] =' ';
    }
}

Int_t PBatch::EvalPFormula(char * command) {
    //This helper function evaluates a (possible)
    //PFormula command, and put the arguments on the stack
    //if this did not worked out, the stack is resetted
    AddSpacePlaceholder(command);
    //    cout << "PFormula:" << command << endl;
    PFormula * tmp = new PFormula(command,command);    
    Int_t worked_out=1, num_params=0;

    const char * mod_command = command;
    Int_t lst_key_tmp[MAX_COMMAND_OPTIONS]; //can be overwritten by "daughter" commands
    Int_t lst_options_counter_tmp=0;

    for (int i=0;i<MAX_COMMAND_OPTIONS;i++) lst_key_tmp[i] = -1;

    while (tmp->error_code && num_params<(MAX_COMMAND_OPTIONS-1)) {
	
	//try to get the ugly guy
	//copy the error string first
	char * internal_command = new char[strlen(tmp->error_string.Data())+1];
	strcpy(internal_command,tmp->error_string.Data());

	//if the error string is *exactly* the input command, something is wrong....
	if (strcmp(internal_command,command)==0) {
	    delete(tmp);
	    RemoveSpacePlaceholder(command);
	    return 0;
	}

	//	mod_command = ...aus PFormula, after replacing	
	//	TString *op = new TString(tmp->GetTitle());
	TString *op = new TString(tmp->chaine);
	char opt[5];
	sprintf(opt,"[%i]",num_params);
	//op->ReplaceAll(tmp->error_string,opt);
	ReplaceAll(op,tmp->error_string,opt);
	mod_command = op->Data();
	delete(tmp);
	RemoveSpacePlaceholder(internal_command);
	Int_t key = GetKey(internal_command,1,-1);
	if (key>=0) lst_key_tmp[num_params+1]=key;
	lst_options_counter_tmp++;
	tmp = new PFormula(mod_command,mod_command);   
	num_params++;

    }

    if (tmp->error_code) worked_out=0;

    if (worked_out) {
	lst_form[command_pointer]=tmp;
	lst_options_counter[command_pointer]=lst_options_counter_tmp;
	for (int i=0;i<MAX_COMMAND_OPTIONS;i++)
	    lst_key[i][command_pointer]=lst_key_tmp[i];
	Int_t key_a = makeStaticData()->
	    MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,command);
	Double_t * delme =  new Double_t(0.);
	makeDataBase()->SetParamDouble(key_a ,"batch_value", delme);
	lst_key_a[command_pointer]=key_a;
	//	lst_key[1][command_pointer]=key1;
	lst_command[command_pointer]=COMMAND_PFORMULA;
	command_pointer++;
    } 

    RemoveSpacePlaceholder(command);
    return worked_out;
}

void  PBatch::ReplaceAll(TString *op,const char * oldstring,const char * newstring) {
    //Similar to TString::ReplaceAll, but checks also the following
    //char

    bool doloop = kTRUE;


    while (doloop) {
	const char * mystring = op->Data();
	doloop = kFALSE;

	for (int i=0;i<((int)strlen(mystring)-(int)strlen(oldstring)+1);i++) {
	    //loop over possible matching positions

	    if (strncmp(oldstring,mystring+i,strlen(oldstring))==0) {
		//match
		char nextchar = mystring[i+strlen(oldstring)];

		if ((nextchar == '\0') || (!isalpha(nextchar) && 
					   (nextchar != '.') &&
					   (nextchar != '-') &&
					   (nextchar != '_')
					   )) {
		    op->Replace(i,strlen(oldstring),newstring,strlen(newstring));
		    i=strlen(mystring)+1;
		    doloop = kTRUE;  //yet another try
		}
	    }
	}
    }
}

Int_t PBatch::GetMethodHandle(char * name, Int_t flag) {
    //seeks for the internal method
    //returns a -1 if not found or invalid

    TList *list;
    if (flag==0)
	list=  gROOT->GetClass("PParticle")->GetListOfAllPublicMethods();
    else if (flag==1)
	list=  gROOT->GetClass("PUtilsREngine")->GetListOfAllPublicMethods();
    else Fatal("GetMethodHandle","Unsupported flag");

    TIterator *iter = list->MakeIterator();
    TMethod *meth=NULL;
    
    while ((meth=(TMethod *) iter->Next())) {

	UInt_t realname_len = 0;
	for (;realname_len<strlen(name);realname_len++) if (name[realname_len]=='(') break;

	if ((strncmp(meth->GetName(),name,realname_len) == 0) && 
	    (strlen(meth->GetName()) == realname_len)) {
	    //found name

	    if ((strcmp(meth->GetReturnTypeName(),"Double_t") ==0 ) || 
		(strcmp(meth->GetReturnTypeName(),"Int_t") ==0 ) || 
		(strcmp(meth->GetReturnTypeName(),"void") ==0 )) {

		//Set handle and prepare the pointers		 
		//is MethodCall existing?

		for (int i=0;i<method_pointer;i++) {

		    if (strcmp(method_name[i],meth->GetName()) == 0) {
			//found old entry			 
			CrackMethodArgs(name+realname_len+1);
			return i;
		    }
		}

		if (method_pointer == MAX_COMMAND_TMETHODS) {
		    Error("GetMethodHandle","MAX_COMMAND_TMETHODS reached");
		    return -1;
		}

		arg1=arg2=arg3=arg4=-1;

		TString argstring("");
		if (meth->GetNargs() != 0) {
		    TList *arg_list = meth->GetListOfMethodArgs();
		    TIterator *arg_iter = arg_list->MakeIterator();
		    TMethodArg *arg_meth;

		    methods_arg_flags[0][method_pointer]
			= methods_arg_flags[1][method_pointer]
			= methods_arg_flags[2][method_pointer]
			= methods_arg_flags[3][method_pointer] = -1;
		    int j=0;
		    while ((arg_meth=(TMethodArg *) arg_iter->Next())) {
			//set up the argument list
			//cout << arg_meth->GetTypeName() << endl;
			if (j==4) {
			    Error("GetMethodHandle","More then 4 argument, not yet supported");
			    return -1;
			}
			if (strlen(argstring.Data()) > 0) argstring += TString(",");
			if (strcmp(arg_meth->GetTypeName(),"Double_t") ==0 ) {
			    methods_arg_flags[j][method_pointer] = METHOD_RETURN_DOUBLE;
			    argstring += TString("Double_t");
			}
			else if (strcmp(arg_meth->GetTypeName(),"Int_t") ==0 ) {
			    methods_arg_flags[j][method_pointer] = METHOD_RETURN_INT;
			    argstring += TString("Int_t");
			}
			else if (strcmp(arg_meth->GetTypeName(),"void") ==0 )
			    methods_arg_flags[j][method_pointer] = METHOD_RETURN_VOID;
			else {
			    Error("GetMethodHandle","Method %s has a %s as an argument, not yet supported",
				  name,arg_meth->GetTypeName());
			    return -1;
			}
			j++;
		    }
		     
		    CrackMethodArgs(name+realname_len+1);
		     

		}

		methods[method_pointer] = new TMethodCall();


		 
		if (strcmp(meth->GetReturnTypeName(),"Double_t") ==0 )
		    methods_flags[method_pointer] = METHOD_RETURN_DOUBLE;
		else if (strcmp(meth->GetReturnTypeName(),"Int_t") ==0 )
		    methods_flags[method_pointer] = METHOD_RETURN_INT;
		else 
		    methods_flags[method_pointer] = METHOD_RETURN_VOID;

		char * new_name = new char[strlen(argstring.Data())+1];
		strcpy(new_name,argstring.Data());

		if (flag==0)
		    methods[method_pointer]->InitWithPrototype(gROOT->GetClass("PParticle"),
							       meth->GetName(),new_name);
		else if (flag==1)
		    methods[method_pointer]->InitWithPrototype(gROOT->GetClass("PUtilsREngine"),
							       meth->GetName(),new_name);
		new_name = new char[strlen(meth->GetName())+1];
		strcpy(new_name,meth->GetName());
		method_name[method_pointer]= new_name;
		method_pointer++;
		return method_pointer-1;
	    }  else {
		//unsupported return type
		Error("GetMethodHandle","Method %s has a return type %s, not yet supported",
		      name,meth->GetReturnTypeName());
		return -1;
	    }
	}

    }
    return -1;
}

void  PBatch::CrackMethodArgs(char * name) {
    //Input: the arg string including the trailing ")"
    arg1=arg2=arg3=arg4=-1;
    //cout << "CrackMethodArgs:" << name << endl;
    //nested objects should not be cracked!
    //workaround: replace , by "

    Int_t numbrack=0;
    for (UInt_t i=0;i<strlen(name);i++) {
	if ((name[i]=='(') || (name[i]=='[')) numbrack++;
	if ((name[i]==')') || (name[i]==']')) numbrack--;
	if ((name[i]==',') && (numbrack>0)) name[i]='"';
    }

    char *prodx[4];
    Int_t prodx_s=4; //max 4 args
    
    PUtils::Tokenize(name, ",", prodx, &prodx_s);
    if (!prodx_s) return;    
    *(prodx[prodx_s-1]+strlen(prodx[prodx_s-1])-1) = '\0'; //remove trailing )

    numbrack=0;
    for (UInt_t i=0;i<strlen(name);i++) {
	if ((name[i]=='(') || (name[i]=='[')) numbrack++;
	if ((name[i]==')') || (name[i]==']')) numbrack--;
	if ((name[i]=='"') && (numbrack>0)) name[i]=',';
    }
    
    arg1=GetKey(prodx[0],1,-1);
    if (prodx_s>1) arg2=GetKey(prodx[1],1,-1);
    if (prodx_s>2) arg3=GetKey(prodx[2],1,-1);
    if (prodx_s>3) arg4=GetKey(prodx[3],1,-1);

}

Int_t PBatch::CheckObjectType(Int_t key) {
    
    TObject *obj=NULL;
    
    makeDataBase()->GetParamTObj (key,"batch_particle", &obj);
    
    //anything else?
    Double_t *res=NULL;
    makeDataBase()->GetParamDouble (key,"batch_value", &res);
    
    
    if (!obj && !res) {
	return -1;
    }
    if (obj) return IS_OBJECT;
    if (res) return IS_DOUBLE;

    return 0;

}

Bool_t PBatch::AddCommand(char command,int key_a,int key1,int key2, int key3, int key4, int key5) {

    if (command == '=' && key_a==key1)
	return kTRUE; //filter nonsense

    if (command_pointer == MAX_COMMAND_POINTER) {
	Error ("AddCommand","MAX_COMMAND_POINTER reached");
	return kFALSE;
    }

    lst_command[command_pointer]=command;
    lst_key_a[command_pointer]=key_a;
    lst_key[1][command_pointer]=key1;
    lst_key[2][command_pointer]=key2;
    lst_key[3][command_pointer]=key3;
    lst_key[4][command_pointer]=key4;
    lst_key[5][command_pointer]=key5;

    command_pointer++;

    return kTRUE;

}

Bool_t PBatch::GetArguments(const char *a, const char *b, 
			    char * name,char ** function,char ** arg1,char ** arg2) {
    //looks for syntax like f(a,b);
    
    char *prod[2];
    Int_t prod_s=2; //max 2 products

    PUtils::Tokenize(name,a, prod, &prod_s);

    if (prod[1]==NULL) {
	*function=NULL;
	prod[1]=prod[0];	
    }

    *(prod[1]+strlen(prod[1])-1) = '\0';
    //kill ")"

    int komma = 0;

    for (UInt_t i=0; i<strlen(prod[1]); i++) {
	if (*(prod[1]+i)==',') komma++;
    }
    
    if (komma>1) {
	Error("AddCommand","[%s] Too many kommas",prod[1]);
	return kFALSE;
    }

    if (komma) {
	char *prodx[2];
	Int_t prodx_s=2; //max 2 products
	
	PUtils::Tokenize(prod[1],",", prodx, &prodx_s);
	*arg1=prodx[0];
	*arg2=prodx[1];
    } else {
	*arg1=prod[1];
	*arg2=NULL;
    }


    return kTRUE;

}

Int_t PBatch::GetDelimPosition(char *arg,char delim, Int_t *yes) {
    
    Int_t brackets=0, split_pos=-1;
    Int_t found_something=0;
    for (UInt_t i=0; i<strlen(arg); i++) {
	if ((arg[i]=='(') || (arg[i]=='[')) brackets++;
	if ((arg[i]==')') || (arg[i]==']')) brackets--;
	if ((arg[i]==delim) && (brackets==0)) {
	    split_pos=i;
	    if (arg[i+1]=='>' || !found_something) split_pos=-1;
	    //...skip -> and - in the beginning
	    else
		if (yes) (*yes)++;
	}
	if (arg[i]!=' ')  found_something++;
    }

    return split_pos;
}

Bool_t PBatch::CheckAndSplit(char * arg,char delim,int * key1, int * key2) {

    int yes = 0, split_pos=-1;
    
    char *arg2 = new char[strlen(arg)+1];
    strcpy(arg2,arg);
    arg=arg2;

    if (strstr(arg,"->")) return kFALSE;


    split_pos=GetDelimPosition(arg,delim,&yes);

    if (!yes) {
	return kFALSE; 
    }

    char *prod[2];
    
    prod[0]=new char [split_pos+1];
    prod[1]=new char [strlen(arg)-split_pos+1];
    strncpy(prod[0],arg,split_pos);
    *(prod[0]+split_pos) = '\0';
    strncpy(prod[1],arg+split_pos+1,strlen(arg)-split_pos-1);
    *(prod[1]+strlen(arg)-split_pos-1) = '\0';

    *key1 = GetKey(prod[0],1,-1);
    *key2 = GetKey(prod[1],1,-1);

    return kTRUE;
}




Int_t PBatch::GetKey(char * name, int fl, int makeflag) {

    //makeflag=0  : Pure GetKey, just take the key from DB
    //              no AddCommand()
    //makeflag=1  : put our object into the DB in any case
    //              this is the case for lvalues
    //              e.g. obj = ....
    //              (but see also SetVarList)
    //makeflag=2  : Same as above, but it must be a valid variable name
    //makeflag=-1 : put into DB only if we have a composite object
    //              This is the case if we find () or + or - 
    //              e.g. obj->...
    //              make AddCommand
    //fl=2        : Take SetVarList into account

    //cout << "getKey" << name << ":" << makeflag <<  endl;
  
    //do we contain ()?
    //int br = 0;
    int found_br = 0;

    //remove brackets if fl==1
    PUtils::remove_spaces(&name);
    if (fl) {
	found_br =  PUtils::remove_brackets(&name,'(',')' );
    }
    PUtils::remove_spaces(&name);

    Int_t key = -1; 

    if (strlen(name)==0) return -1;

    //first check for file-input
    if (name[0] == '[' && name[strlen(name)-1]==']' && strlen(name)>2) {
	//found [] at the very begin & end
		
	//further look to the content
	//if we find any additional brackets, it is a composite object!
	Int_t num_brackets =0;
	for (UInt_t i=1;i<(strlen(name)-1);i++) {
	    if (name[i] == '[' || name[i]==']') num_brackets++;
	}

	if (!num_brackets) {

	    key = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);
	    
	    char * function,*arg1,*arg2;
	    GetArguments("[","]",name,&function,&arg1,&arg2);
	   	    
	    if (!arg1) { 
		Error("AddCommand","[%s] Argument not found",name);
		return kFALSE;
	    }
	   
	    Int_t pid = 0; //DUMMY

	    if (strcmp(arg1,"*")) {
		//no Joker
		pid = makeStaticData()->GetParticleID(arg1);
		if (pid==0)  { 
		    Error("AddCommand","[%s] Unknown particle %s",name,arg1);
		    return kFALSE;
		}
	    }
	    
	    Int_t *ii=new int(pid);  //never destructed, but called only once!
	   	
	    if (!makeDataBase()->SetParamInt (key, "pid", ii)) {
		delete ii;
		return kFALSE;
	    }
	    
	    Int_t number=-999;
	    
	    if (arg2) { //try to get number
		//cout << "arg: " << arg2 << endl;
		Int_t * delme = new Int_t (number);
		Int_t numkey = -999;
		//First, let's check for a variable
		if (arg2[0] == '$') {
		    if (PUtils::ValidVariableName(arg2+1)) {
			numkey = makeDataBase()->GetEntry(arg2+1);
			if (numkey>0) {
			    Double_t *val=NULL;
			    if (!makeDataBase()->GetParamDouble (numkey,batch_value_param, &val)) {
				numkey = -1;
			    } else {
				Int_t * result = new Int_t (1);
				makeDataBase()->SetParamInt (numkey, "batch_update", result);   
			    }
			}
			
		    } 
		    if (numkey<0) {
			Error("AddCommand","[%s] Unknown variable '%s'",name,arg2);
		    } 
		    *delme = -1000 - numkey; //-999 means not found!
		} else if (PUtils::IsInt(arg2)) {
		    sscanf(arg2,"%i",&number);
		    //cout << "found number: " << number << endl;
		    *delme = number;
		    if (strcmp(arg2,"+")==0) { 
			*delme = -111;
			status=1;
		    }
		} else Error("AddCommand","[%s] Unknown value '%s'",name,arg2);
		//mis-use link:
		//cout << "delme: " << *delme << endl;
		makeDataBase()->SetParamInt (key ,"link", delme);
	    }
	    makeflag = 0;
	}
    } //END file input

    if (makeflag==-1 && PUtils::ValidVariableName(name)) {
	//cancel AddCommand
	makeflag=0;
	//cout << name << " cancelled" << endl;
    }

    if (makeflag==2 && !PUtils::ValidVariableName(name)) {
	//cancel AddCommand
	makeflag=0;
	Error("GetKey","[%s] is not a valid variable name",name);
    }
    if (makeflag==2) makeflag=1;

    if (makeflag==-1) AddCommand(name);

    
    if (makeflag) {
	if ((varlist == NULL) || (fl!=2)) {
	    key = makeStaticData()->
		MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);
	} else {
	    for (unsigned int j=0; j< ((strlen(varlist) - strlen(name)) ); j++) {
		if ((strncmp(varlist+j,name,strlen(name))==0) && 
		    varlist[j+strlen(name)]==';') {
		    key = makeStaticData()->
			MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);
		}
	    }
	    if (key<0) Error("AddCommand","[%s] is not an allowed object name, the list in this context is [%s] ",name,varlist);
	}
    } else if (key<0) {
	key = makeDataBase()->GetEntry(name);
    }



    if (key<0) {
	Error ("AddCommand","[%s] Unknown object",name);
    }
    
    return key;
}




void PBatch::Print(const Option_t* delme) const {

    cout << "Command list:" << endl;

    for (int i=0;i<command_pointer;i++) {

	if (lst_command[i] == COMMAND_PFORMULA) {

	    cout << makeDataBase()->GetName(lst_key_a[i]) 
		 << " <" << lst_command[i] << "> ";
	    for (int j=0;j<lst_options_counter[i];j++) {
		cout << makeDataBase()->GetName(lst_key[j+1][i])  << " ";
	    }
	    cout << endl;
	} else {
	    if (lst_key[2][i]>-1) {
		cout << makeDataBase()->GetName(lst_key_a[i]) 
		     << " <" << lst_command[i] << "> " <<  makeDataBase()->GetName(lst_key[1][i])  
		     << "," << makeDataBase()->GetName(lst_key[2][i]) << endl;
	    } else if (lst_key[1][i]>-1) {
		cout << makeDataBase()->GetName(lst_key_a[i]) 
		     << " <" << lst_command[i] << "> " <<  makeDataBase()->GetName(lst_key[1][i])  
		     << endl;
	    } else if (lst_key_a[i]>-1) {
		cout << makeDataBase()->GetName(lst_key_a[i]) 
		     << " <" << lst_command[i] << "> " 
		     << endl;
	    } else { cout 
		     << " <" << lst_command[i] << "> " 
		     << endl;
	    }
	}
    }

}


ClassImp(PBatch)

