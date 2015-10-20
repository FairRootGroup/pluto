/////////////////////////////////////////////////////////////////////
//
// Generic angular distribution
// This algorithm evaluates the angular distribution of
// the particle "primary" which must be a "daughter"
//
// The angle to be evaluated is the relative angle between
// the "primary" particle momentum in the "reference" frame
// and the "reference" particle momentum. This momentum, however
// is determined in the "base_reference" frame.
// If "ang_reference" is used, the angle is determined
// relative to this particle
//
//
// If "base_reference" is not set, the default is the lab frame
// If "reference" is not set, the default is the parent particle
// If "ang_reference" is not set, but "beam", ang_reference
// is set to beam
//
//
// If "align" is set, the reference will be the composite particle
// of the primary + align 
// This can be used for the alignment of two particles
//
// The 2nd parameter for TF2 is the mass of the reference
// If "mass_reference" is set, this particle (e.g. q) will be used
// Such a feature is useful for energy-dependent parameterizations
//
//
//                                  Author:  I. Froehlich
//                                  Written: 01.10.2006
//                                  Revised: 24.10.2010
//                                  
//                                  
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PAngularDistribution.h"


ClassImp(PAngularDistribution)

PAngularDistribution::PAngularDistribution() {
    Fatal("PAngularDistribution()","Wrong constructor called");
};

PAngularDistribution::PAngularDistribution(const Char_t *id,const  Char_t *de) :
    PDistribution(id, de) {

    angles1 = NULL;
    angles2 = NULL;
    anglesg = NULL;
    reference = NULL;
    base_reference = NULL;
    ang_reference = NULL;
    parent = NULL;
    align = NULL;
    check_abort = never_abort = always_reject = kFALSE;
    rotate = kTRUE;
    mass_reference = NULL;
    beam = target = NULL;
    align_is_daughter = kFALSE;
    direct_sampling_possible=kFALSE;
    reflection_symmetry=kFALSE;
    spline=kFALSE;
    g_spline=NULL;
    n_daughters=0;

    fNpar = 1 ; //Parameter1: mass_reference
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
	fNames[0] = "Mass reference for 2-dim angular distributions";
    }
    fParams[0]=0;

    fNpx       = 1000;   
    fXmin      = -1;
    fXmax      = 1;

    for (int i=0;i<MAX_ANG_NUM;i++) {
	daughter[i]=NULL;
    }


    weight_max = 1.;
};

PDistribution* PAngularDistribution::Clone(const char*delme) const {
    return new PAngularDistribution((const PAngularDistribution &)* this);
};

Bool_t PAngularDistribution::Init(void) {

    if ((angles1 == NULL ) && (angles2 == NULL ) && (anglesg == NULL )) {
	//At least one of the parameterized functions must be present
	Error("Init", "Angular distribution not found");
	return kFALSE;
    }

    //looking for primary. This is also mandatory
    primary = GetParticle("primary");
    if (!primary) {
	Error("Init", "Primary particle not found");
	return kFALSE;
    }


    //primary->Print();
    reference = GetParticle("reference");
    base_reference = GetParticle("base_reference");
    mass_reference = GetParticle("mass_reference");
    ang_reference  = GetParticle("ang_reference");
    align = GetParticle("align");
    
    if (align) {//check if align is a daughter or parent
	if ((current_flag == PARTICLE_LIST_PARENT) ||
	    (current_flag == PARTICLE_LIST_DAUGHTER))
	    align_is_daughter = kTRUE;
    }

    //In addition we need the beam and target for some models
    beam = GetParticle("beam");
    target = GetParticle("target");

    if (!ang_reference && beam) ang_reference=beam;

    //Now get the parent
    for (int i=0; i<position; i++) {
	if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    parent=particle[i];
    }

    if (!parent) {
	Error("Init", "Parent not found");
	return kFALSE;
    }

    if (!reference) reference = parent;
    
    for (int i=0; i<position; i++) {
	if ((particle_flag[i] & PARTICLE_LIST_GRANDPARENT) 
	    | (particle_flag[i] & PARTICLE_LIST_GRANDGRANDPARENT) 
	    | (particle_flag[i] & PARTICLE_LIST_SIBLING)) 
	    check_abort = kTRUE;
	
    }

    //get ADDITIONAL daughters
    bool myloop=1;
    for (n_daughters=0;n_daughters<MAX_ANG_NUM && myloop;n_daughters++) {
	daughter[n_daughters]=GetParticle("daughter");
	if (!daughter[n_daughters]) {
	    myloop=0;
	} 
    }
    n_daughters--;
    

    direct_sampling_possible=kFALSE; 

    if (!align &&             // no composite
	!mass_reference &&    // too complicated
	reference==parent &&  // direct sampling in parent frame
	n_daughters &&        // No mistake in template
	!always_reject &&
	angles1)              // Only 1dim-sampling
	direct_sampling_possible=kTRUE;
    //N.B. this is overwritten by inherited models (e.g.PDeltaAngularDistribution)

    return kTRUE;    
};

Bool_t PAngularDistribution::Prepare(void) {

    direct_sampling_done=kFALSE;
    return kTRUE;
};

Bool_t PAngularDistribution::Finalize(void) {
    return kTRUE;
};

Bool_t PAngularDistribution::CheckAbort(void) {
    //    cout << check_abort << endl;
    if (never_abort) return kFALSE;
    return check_abort;
};


Bool_t PAngularDistribution::SampleAngle(void) {
    //Try first to sample the angle, rather to use
    //the rejection method. This works
    //in the simple cases, and only if there
    //are no correlated angular distributions

    if (!direct_sampling_possible) return kFALSE;

     // primary->Print();
     // for (int i=0;i<n_daughters;i++) daughter[i]->Print();

    if (!Rotate(1)) return kFALSE;    

      // cout << "---" << endl;
      // primary_tmp.Print();
      // for (int i=0;i<n_daughters;i++) daughter[i]->Print();
      // cout << "*********************" << endl;

    //now we rotate all daughters such that primary_tmp is pointing to the z-axis
    double prim_phi   = primary_tmp.Phi();
    double prim_theta = primary_tmp.Theta();

#if 1

    //primary->Print();
    // primary_tmp.Print();
    // for (int i=0;i<n_daughters;i++) daughter[i]->Print();

    primary_tmp.RotateZ(-prim_phi);
    primary_tmp.RotateY(-prim_theta);
    for (int i=0;i<n_daughters;i++) {
	daughter[i]->RotateZ(-prim_phi);
	daughter[i]->RotateY(-prim_theta);
    }



    //sample the polar angle
    double polar_angle=acos(SamplePolarAngle(0.5*(cos(prim_theta)+1)));   
    //2r=cos_theta-1 was the definition in the baryon_cos alg of the original PChannel
    //polar_angle=prim_theta;

    if ((polar_angle < -2.*TMath::Pi()) || (polar_angle > 2.*TMath::Pi())) {
	Warning("SampleAngle","Wrong polar angle");
	return kFALSE;
    }

    //rotate the particles such that primary has the correct polar angle
    primary_tmp.RotateY(polar_angle);
    for (int i=0;i<n_daughters;i++)
	daughter[i]->RotateY(polar_angle);

    //restore the old azimuthal angle
    primary_tmp.RotateZ(prim_phi);
    for (int i=0;i<n_daughters;i++)
	daughter[i]->RotateZ(prim_phi);
#endif
    //rotate back to the old frame
    if (!RotateBack(1)) return kFALSE;        

    //if all worked out we can copy the primary :
    *primary = primary_tmp;
    

    
    direct_sampling_done=kTRUE;
    
    return kTRUE;
}

Bool_t PAngularDistribution::Rotate(Int_t rotate_daughters) {
    //Tool method to rotate all tmp particles into
    //the defined frames
    //if rotate_daughters=1 the daughters (the original one!)
    //are rotated as well. This is needed for direct samplings

    Int_t loc_n_daughters = n_daughters * rotate_daughters;
    
    primary_tmp = primary;  //particle under investigation. Make better a copy
    primary_tmp.Boost(parent->BoostVector());  // go back to lab frame
    for (int i=0;i<loc_n_daughters;i++) 
	daughter[i]->Boost(parent->BoostVector());
    PParticle compound(primary_tmp);
    PParticle atmp;
    if (ang_reference) ang_tmp = ang_reference;
    //if (ang_reference) ang_tmp . Print();

    if (align) {
	atmp=align;
	if (align_is_daughter)
	    atmp.Boost(parent->BoostVector());  // go back to lab frame
	compound.AddTmp(atmp); //Do not use "+" in evtloop
	reference = &compound;
    }

    if (base_reference==NULL) {  // no 2nd reference frame
	if (rotate) {
	    primary_tmp.RotateZ(-reference->Phi());
	    primary_tmp.RotateY(-reference->Theta());
	    primary_tmp.Boost(0,0,-reference->Beta());
	    ang_tmp.RotateZ(-reference->Phi());
	    ang_tmp.RotateY(-reference->Theta());
	    ang_tmp.Boost(0,0,-reference->Beta());
	    for (int i=0;i<loc_n_daughters;i++) {
		daughter[i]->RotateZ(-reference->Phi());
		daughter[i]->RotateY(-reference->Theta());
		daughter[i]->Boost(0,0,-reference->Beta());
	    }
	} else {
	    primary_tmp.Boost(-reference->BoostVector());
	    for (int i=0;i<loc_n_daughters;i++) {
		daughter[i]->Boost(-reference->BoostVector());	
	    }
	}
    } else { // first go to 2nd reference (e.g. base_reference=eta), then
	// rotate and boost to 1st reference (e.g. reference=dilepton)
	reference_tmp = reference;

	primary_tmp.Boost(-base_reference->BoostVector());  // daughter in base_reference
	ang_tmp.Boost(-base_reference->BoostVector());  
	reference_tmp.Boost(-base_reference->BoostVector()); // reference in base_reference
	for (int i=0;i<loc_n_daughters;i++) {
	    daughter[i]->Boost(-base_reference->BoostVector());
	}
	if (rotate) {
	    primary_tmp.RotateZ(-reference_tmp.Phi());  // rotate daughter reference onto z-axis
	    primary_tmp.RotateY(-reference_tmp.Theta());
	    primary_tmp.Boost(0,0,-reference_tmp.Beta());
	    ang_tmp.RotateZ(-reference_tmp.Phi());
	    ang_tmp.RotateY(-reference_tmp.Theta());	 
	    ang_tmp.Boost(0,0,-reference_tmp.Beta());
	    for (int i=0;i<loc_n_daughters;i++) {
		daughter[i]->RotateZ(-reference_tmp.Phi());
		daughter[i]->RotateY(-reference_tmp.Theta());
		daughter[i]->Boost(0,0,-reference_tmp.Beta());
	    }
	} else {
	    primary_tmp.Boost(-reference_tmp.BoostVector()); // go to reference
	    for (int i=0;i<loc_n_daughters;i++) 
		daughter[i]->Boost(-reference_tmp.BoostVector());
	}
    }
    
    //finally take ang_reference into account
    if (ang_reference && rotate) {
	primary_tmp.RotateZ(-ang_tmp.Phi());  // rotate daughter reference onto z-axis
	primary_tmp.RotateY(-ang_tmp.Theta());
	for (int i=0;i<loc_n_daughters;i++) {
	    daughter[i]->RotateZ(-ang_tmp.Phi());
	    daughter[i]->RotateY(-ang_tmp.Theta());
	}	    
    }

    return kTRUE;
}

Bool_t PAngularDistribution::RotateBack(Int_t rotate_daughters) {
    //restore the all tmp particles and the daughters

    Int_t loc_n_daughters = n_daughters * rotate_daughters;

    if (ang_reference && rotate) {
	primary_tmp.RotateY(ang_tmp.Theta());
	primary_tmp.RotateZ(ang_tmp.Phi());  // rotate daughter reference onto z-axis
	for (int i=0;i<loc_n_daughters;i++) {
	    daughter[i]->RotateY(ang_tmp.Theta());
	    daughter[i]->RotateZ(ang_tmp.Phi());
	}	    
    }

    if (base_reference==NULL) {  // no 2nd reference frame
	if (rotate) {
	    primary_tmp.Boost(0,0, reference->Beta()); // go to 1st reference frame	    
	    primary_tmp.RotateY( reference->Theta());
	    primary_tmp.RotateZ( reference->Phi());
	    for (int i=0;i<loc_n_daughters;i++) {
		daughter[i]->Boost(0,0,reference->Beta());
		daughter[i]->RotateY(reference->Theta());
		daughter[i]->RotateZ(reference->Phi());
	    }
	} else {
	    primary_tmp.Boost(reference->BoostVector());
	    for (int i=0;i<loc_n_daughters;i++) {
		daughter[i]->Boost(reference->BoostVector());
	    }
	}
    } else { // first go to 2nd reference (e.g. base_reference=eta), then
	// rotate and boost to 1st reference (e.g. reference=dilepton)
	primary_tmp.Boost(0,0,reference_tmp.Beta()); // go to reference
	primary_tmp.RotateY(reference_tmp.Theta());
	primary_tmp.RotateZ(reference_tmp.Phi());  // rotate daughter reference onto z-axis

	for (int i=0;i<loc_n_daughters;i++) {
	    daughter[i]->Boost(0,0,reference_tmp.Beta());
	    daughter[i]->RotateY(reference_tmp.Theta());
	    daughter[i]->RotateZ(reference_tmp.Phi());
	}	

	primary_tmp.Boost(base_reference->BoostVector());  // daughter in base_reference
	reference_tmp.Boost(base_reference->BoostVector()); // reference in base_reference
	for (int i=0;i<loc_n_daughters;i++) {
	    daughter[i]->Boost(base_reference->BoostVector()); 
	}
    }
    primary_tmp.Boost(-parent->BoostVector());  // go back to parent frame
    for (int i=0;i<loc_n_daughters;i++) 
	daughter[i]->Boost(-parent->BoostVector());
    return kTRUE;
}



Bool_t PAngularDistribution::IsValid(void) {

    if (direct_sampling_done) return kTRUE;

    Double_t tmp_c0,f;
 
    if (!Rotate(0)) return kFALSE;

    tmp_c0=cos(primary_tmp.Theta());
    
    if ((reflection_symmetry) && (tmp_c0<0.)) {
        tmp_c0=-tmp_c0;
    }

    if (angles1) {
	f=angles1->Eval(tmp_c0);
    }
    else if (angles2) {
	if (mass_reference)
	    f=angles2->Eval(tmp_c0,mass_reference->M());
	else
	    f=angles2->Eval(tmp_c0,reference->M());
    } else {
        f=anglesg->Eval(tmp_c0,g_spline);
    }

    if (f > weight_max) {
	weight_max = f*1.1;
	Warning("IsValid","[%s] Weight > max, new max is %lf",GetName(),weight_max);
    }

    if ((f/weight_max)>PUtils::sampleFlat()) return kTRUE; // sample now distribution
    
    
    return kFALSE;

};

double PAngularDistribution::SamplePolarAngle(double r) {
    if (angles1) {
	return angles1->GetRandom();
    }
    else { //No GetRandom for TGraph, wrapper for TF2

	if (mass_reference) {
	    q_value=mass_reference->M();
	    SetParameter(0,mass_reference->M());
	}
	else {
	    q_value=reference->M();
	    SetParameter(0,reference->M());
	}
	return this->GetRandom();
    }
}

Double_t PAngularDistribution::EvalPar(const Double_t *x, const Double_t *params) {
    if (params) {
	q_value=(int)params[0];
    }
    if (Eval(x[0])<0) {
	Warning("EvalPar", "q_value is negative: %lf",q_value);
    }
    return Eval(x[0]);
}
 
Double_t PAngularDistribution::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    if ((reflection_symmetry) && (x<0.)) x=-x;

    if (angles1) return angles1->Eval(x);
    if (angles2) return angles2->Eval(x,q_value);
    if (anglesg) return anglesg->Eval(x,g_spline,"");
    return 0;
}

void PAngularDistribution::Print(const Option_t* delme)  const{
    BasePrint();

    cout << "    Formula used: ";
    if (angles1) { 
	if (angles1->GetExpFormula() != TString("")) cout << angles1->GetExpFormula();
	else cout << "<compiled>";
    }
    else if (angles2) {
	if (angles2->GetExpFormula() != TString("")) cout << angles2->GetExpFormula();
	else cout << "<compiled>";
    }
    else if (anglesg) {
        cout << "<TGraph>";
    }

    else cout << "NONE";
    cout << endl;

}


