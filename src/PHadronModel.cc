/////////////////////////////////////////////////////////////////////
//
// Class for calculating the total mass-dependent width
// from first principles
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PHadronModel.h"


ClassImp(PHadronModel)

PHadronModel::PHadronModel()  {
} ;

PHadronModel::PHadronModel(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_pid<0)
	Warning("PHadronModel","The model (%s) should be bound to PARTICLES only",de);

    global_weight_scaling=1.;
} ;

PDistribution* PHadronModel::Clone(const char*delme) const {
    return new PHadronModel((const PHadronModel &)* this);
};



Bool_t PHadronModel::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
//Calculates the mass-dependent width (either partial for the decy or
//total for the particles)
//
//If didx is set, return the partial width of the decay channel didx
//This can be used for mass sampling in PChannel
//(else set didx to -1)
//
//Returns the mass-dependent width for a hadronic resonance id of mass m.
//It is based on first principles and is valid for arbitrary resonances
//given the pole mass, vacuum width, and decay modes (units GeV/c**2).
//
 
    if (didx>=0) {
	//Do not calculate the width from first principles
	//but use the partial width only
	*width = makeDynamicData()->GetDecayPartialWidth(mass,didx);
	return kTRUE;
    }

    int twidx=makeStaticData()->GetTWidx(is_pid);

    double g0=
	makeStaticData()->GetParticleTotalWidth(is_pid); // vacuum width

    if (!width_init) // if 1st call, initialize flags
	makeDynamicData()->GetParticleDepth(is_pid);   

    if (twidx==-1 || g0 < (*unstable_width) || mass==0. || 
	makeStaticData()->GetParticleNChannels(is_pid)==0) { 
	*width=g0; // no mass-dependent width
	return kTRUE;
    }

    double dm;            // mass range
    int j, didx2, np;
    double mm, g_tmp;

    Bool_t self_consistency_loop=kTRUE;
    Bool_t next_self_consistency_loop=kFALSE;

    if (!width_init) {   
	// Below: enter once on the first call
	width_init++;
	
	//cout << "Width 1st call for "  << description << ":" << is_pid << ":" << width_init<< endl;
	global_weight_scaling = makeDynamicData()->GetParticleScalingFactor(is_pid);

	mmin=PData::LMass(is_pid);   // mass threshold
	mmax=PData::UMass(is_pid);   // mass ceiling
	
	dm=(mmax-mmin)/(maxmesh-3.); // mass increment for the mesh
 
	Info("GetWidth","Width 1st call for %s, mass range %f GeV to %f GeV",
	     description,mmin,mmax);

	while (self_consistency_loop) {
	    // The mass-dependent width is evaluated on a mesh as a function of mass.
	    

	    

	    
	    //now loop over decay modes
	    Int_t listkey=-1;
	    Int_t tid[11];
	    while (makeDataBase()->MakeListIterator(primary_key, "pnmodes", 
						    "link", &listkey)) {
		tid[0]=10; 
		makeStaticData()->GetDecayModeByKey(listkey,tid); 
		// retrieve current mode info
		didx2=makeStaticData()->GetDecayIdxByKey(listkey);
		np=tid[0];
		if (makeStaticData()->GetPWidx(didx2)!=-1 ) {  
		    // Note: The next call is required 
		    // if the current particle has nested decays via 
		    // any recognized modes.
		    // The innermost widths must be available by the time 
		    // the local Width sampling 
		    // is called and their calculation is forced here.
		    makeDynamicData()->GetDecayPartialWidth(mass,didx2);
		}
	    } //end iterator
	    // _____________________________________________________________________
	    // Store the total width: Unitarity is assumed for normalization, i.e. the
	    // available decay modes are taken to exhaust the total transition strength.
	    // This means:  Sum(all branching ratios)=1

	    if (mesh) delete mesh;

	    mesh = new PMesh(maxmesh-2,"mesh");
	    for (j=0;j<maxmesh-2;++j) {  // loop over the mesh points [0-200]
		g_tmp=0.;                // reset the sum of mass-dependent partial widths
		mm=mmin+j*dm;            // update mass on the mesh
		while (makeDataBase()->MakeListIterator(
			   primary_key, "pnmodes", "link", &listkey)) {
		    didx2=makeStaticData()->GetDecayIdxByKey(listkey);
		    g_tmp+=makeDynamicData()->GetDecayPartialWidth(mm,didx2);
		}
		mesh->SetNode(j,g_tmp);
	    }
	    mesh->SetMin(mmin); // store mass threshold
	    mesh->SetMax(mmax); // store mass ceiling

	    // End of first call: the mass-dependent width is stored for future use


	    //******** Scaling
	    //By the way... Normalize Int(weight) to =1
	    Double_t global_weight=0;
	    for (j=0;j<maxmesh-2;++j) {  // loop over the mesh points [0-200]
		mm=mmin+j*dm;            // update mass on the mesh
		global_weight+=dm*makeDynamicData()->GetParticleTotalWeight(mm,is_pid);		
	    }
	    
	    if (global_weight)
		makeDynamicData()->SetParticleScalingFactor(is_pid,1/(global_weight));
	    global_weight_scaling = makeDynamicData()->GetParticleScalingFactor(is_pid);

	    //*********Self-consistency
	    //It is important that width_init is already defined, otherwise we go in endless loop here
	    listkey=-1;
	    self_consistency_loop=kFALSE;
	    while (makeDataBase()->MakeListIterator(
		       primary_key, "pnmodes", "link", &listkey)) {
		//Loop again over decay modes
		didx2=makeStaticData()->GetDecayIdxByKey(listkey);
		Double_t partial_weight=0.; // reset the sum of mass-dependent branching ratio
		global_weight=0.; 

		int brflag=makeStaticData()->GetDecayBRFlag(didx2);
		if (brflag) { //Do integration 
		    for (j=0;j<maxmesh-2;++j) {  // loop over the mesh points [0-200]
			mm=mmin+j*dm;            // update mass on the mesh
			
			if (makeDynamicData()->GetParticleTotalWidth(mm,is_pid)>0)

			    partial_weight+=dm*(makeDynamicData()->GetDecayPartialWidth(mm,didx2)*
						makeDynamicData()->GetParticleTotalWeight(mm,is_pid)*
						makeStaticData()->GetParticleTotalWidth(is_pid)/
						makeDynamicData()->GetParticleTotalWidth(mm,is_pid));
			//Obtain the partial fold
			//PWidth/TWidth is the mass-dependent BR
			//*TWeight will be the Partial mass shape
			//Int(Partial mass shape)/Int(Total mass shape) must be the static BR
			
			global_weight+=dm*makeDynamicData()->GetParticleTotalWeight(mm,is_pid);
			
			
		    }
//		    cout << didx2 << ":" << partial_weight << ":" << global_weight << endl;
		} else { //pole width
		    mm=makeStaticData()->GetParticleMass(is_pid);
		    partial_weight=(makeDynamicData()->GetDecayPartialWidth(mm,didx2)*
				     makeDynamicData()->GetParticleTotalWeight(mm,is_pid)*
				     makeStaticData()->GetParticleTotalWidth(is_pid)/
				     makeDynamicData()->GetParticleTotalWidth(mm,is_pid));
		    global_weight=makeDynamicData()->GetParticleTotalWeight(mm,is_pid);
		}
		
		//Global weight should be 1 after Scaling

		Double_t calculated_br=0;
		Double_t tot_sc=1.;
		//take into account scfactor at mass pole
		tot_sc=makeDynamicData()->GetParticleTotalWidth(makeStaticData()->GetParticleMass(is_pid),is_pid)/
		    makeStaticData()->GetParticleTotalWidth(is_pid);

		if (global_weight) //Compare the calculated partial width to the static one
		    calculated_br=tot_sc*partial_weight/global_weight;

								      
		if ((calculated_br/makeStaticData()->GetDecayPartialWidth(didx2)>1.001)
		    || (calculated_br/makeStaticData()->GetDecayPartialWidth(didx2)<0.999)) { 
		    //set SC factor and go on...
		    if (!makeDynamicData()->CheckSCAbort(didx2)) { //check for endless loop
			Warning("GetWidth","Self consistency calculation failed for decay: %s (%i)",
				makeStaticData()->GetDecayName(didx2),didx2);
			Warning("GetWidth","Calculated Width: %f, Static Width: %f", 
				calculated_br, 
				makeStaticData()->GetDecayPartialWidth(didx2));
			next_self_consistency_loop=kTRUE; //break the complete loop			
		    } else {
			if (calculated_br>0) {
			    self_consistency_loop=kTRUE;
			    makeDynamicData()->SetDecaySCFactor(didx2, makeStaticData()->GetDecayPartialWidth(didx2)/calculated_br);
			};
		    }
		}
	    } //END modes iterator
	    if (next_self_consistency_loop) self_consistency_loop=kFALSE; 
	}  //END	while (self_consistency_loop) {
    } //width init

    *width= mesh->GetLinearIP(mass);

    return kTRUE;
}

