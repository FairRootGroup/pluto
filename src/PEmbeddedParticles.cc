////////////////////////////////////////////////////////
//  Pluto embedded particles interface
//
//  With this class one may add additional
//  Particles in the PReaction chain
//
//
//  The embedded particles will be added
//  to the PReaction. If PVertexFile was used used before, the
//  eventloop will be stopped latest at the last entry of the ntuple
//  no matter if a larger number was given to the loop (see PVertexFile).
//  The event sequence number of the real data vertex ntuple (which identifies
//  which event has exacly has been used) is stored in ASCII output.
//  This variables will be later transported to HGeantKine and can be
//  retrieved via HGeantKine->getUserVal().
//  The synchronization in event embeding makes use out
//  of this data.
//                    Author:  Ingo Froehlich
//                    Written: 11/07/2007
//                    Revised: 11/28/2007 (Jochen Markert)
//
////////////////////////////////////////////////////////

#include "PEmbeddedParticles.h"
#include "PChannel.h"


PEmbeddedParticles::PEmbeddedParticles() {
    local_pos = 0;
    for (Int_t i = 0; i < MAX_EMBEDDEDPARTICLES; i ++) {
	local_pMin        [i] = 0;
	local_pMax        [i] = 0;
	local_openingAngle[i] = 0;
	local_theta       [i] = 0;
	local_phi         [i] = 0;
	local_mMax        [i] = 0;
	local_mMin        [i] = 0;

	local_thetaMin    [i] = 0;
	local_thetaMax    [i] = 0;
	local_phiMin      [i] = 0;
	local_phiMax      [i] = 0;

	local_version     [i] = 0;
	downscaling       [i] = 0;
	last_downscaling  [i] = 0;
	local             [i] = NULL;
    }

    vertex_x = makeStaticData()->GetBatchValue("_event_vertex_x");
    vertex_y = makeStaticData()->GetBatchValue("_event_vertex_y");
    vertex_z = makeStaticData()->GetBatchValue("_event_vertex_z");

}

Bool_t  PEmbeddedParticles::SetSampling(Double_t pMin, Double_t pMax,
					Double_t openingAngle, Double_t theta,
					Double_t phi,
					Double_t mMin, Double_t mMax) {
    // The current particle will be re-sampled during the event loop
    // Momentum sampling is done between "pMin" and "pMax" [rad]
    // A region with "openingAngle" [rad] pointing to "theta,phi" will be covered
    // with white mass between "mMin" and "mMax" [GeV/c]. if "mMax" < 0
    // the already set mass of the particle will be used during re-sampling.

    if (local_pos == 0) {
	Warning("SetSampling","no particle defined");
	return kFALSE;
    }
    local_pMin        [local_pos - 1] = pMin;
    local_pMax        [local_pos - 1] = pMax;
    local_mMin        [local_pos - 1] = mMin;
    local_mMax        [local_pos - 1] = mMax;
    local_openingAngle[local_pos - 1] = -cos(openingAngle);

    local_theta       [local_pos - 1] = theta;
    local_phi         [local_pos - 1] = phi;

    local_version     [local_pos - 1] = 0;

    return kTRUE;
};


Bool_t  PEmbeddedParticles::SetSamplingSector(Double_t pMin       , Double_t pMax,
					      Double_t thetaMin   , Double_t thetaMax,
					      Double_t phiMin     , Double_t phiMax,
					      Double_t phiStartVal, Int_t numParticle,
					      Double_t deltaPhi   , Int_t numSectors)
{
    // The current particle will be re-sampled "numParticle" times per sector
    // during the event loop. Momentum sampling is done between "pMin"  "pMax",
    // theta (polar) between "thetaMin" - "thetaMax" and phi between "phiMin" - "phiMax" (0-60 deg).
    // The phi coordinate can be rotated by phiStart. Angles are in degrees, mom in GeV/c.
    // If the mass of the Particle has been set already before adding it to the bulk this
    // mass will be used during re-sampling. Otherwise the particle masses from PStaticData
    // are used.
    if(local_pos == 0){
	Warning("SetSamplingSector","no particle defined");
	return kFALSE;
    }
    startPhi  = phiStartVal;
    nParticle = numParticle;
    Int_t addedParticle = 0;


    if(nParticle > MAX_EMBEDDEDPARTICLES/numSectors){

	nParticle = MAX_EMBEDDEDPARTICLES/numSectors;
	Warning("SetSamplingSector()",
		"Number of particles /sector > MAX_EMBEDDEDPARTICLES/%i, will be trucated to numParticle = %i .",
		numSectors,MAX_EMBEDDEDPARTICLES/numSectors);
    }

    for(Int_t s = 0; s < numSectors; s ++)
    {
	for(Int_t n = 0; n < nParticle; n ++){
	    if(addedParticle > 0){
		// first particle has been already added
		// clone the other one by one
		AddParticle((PParticle*)local[local_pos - 1]->Clone(),
			    downscaling[local_pos - 1]);
		//...local_pos incremented in the line above!
	    }
	    Double_t phi_shift = (s * deltaPhi + startPhi) * TMath::DegToRad();

	    local_pMin    [local_pos - 1] = pMin;
	    local_pMax    [local_pos - 1] = pMax;
	    local_thetaMin[local_pos - 1] = thetaMin * TMath::DegToRad();
	    local_thetaMax[local_pos - 1] = thetaMax * TMath::DegToRad();
	    local_phiMin  [local_pos - 1] = phi_shift + phiMin * TMath::DegToRad();
	    local_phiMax  [local_pos - 1] = phi_shift + phiMax * TMath::DegToRad();
	    local_version [local_pos - 1] = 1;
	    if(local[local_pos - 1]->M() > 0){
		// the mass of the particle has been set
                // outside already
		local_mMax    [local_pos - 1] = -1;
	    }
	    addedParticle ++;
	}
    }

    return kTRUE;
}



Bool_t PEmbeddedParticles::AddParticle(PParticle * particle, int downsc) {
    // Adds the "particle" to the stack to be embedded during the event loop
    // if "downsc" is larger then 1, only each downsc^th event will be filled

    if (local_pos == MAX_EMBEDDEDPARTICLES ) {
	cout << "MAX_EMBEDDEDPARTICLES reached" << endl;
	return kFALSE;
    }
    downscaling     [local_pos]    = downsc;
    last_downscaling[local_pos]    = 1;
    local           [local_pos ++] = particle;
    return kTRUE;
}

Bool_t PEmbeddedParticles::Modify(PParticle ** mstack, int *decay_done, int * num, int stacksize) {
    // Read the particles from the defined stack and copy this to the official
    // particle stream

    int cur = *num;
    Double_t  p,cos_theta,sin_theta,theta,phi,m;


    for (int i = 0; i < local_pos; i ++) {
	if (last_downscaling[i] == downscaling[i])
	{
	    last_downscaling[i] = 1;

	    *(mstack[cur]) = *(local[i]);

	    //set embedded particle to global vertex
	    mstack[cur]->SetVertex(*vertex_x,
				   *vertex_y,
				   *vertex_z,0.);

	    mstack[cur]->SetSourceId(EMBEDDED_SOURCE_ID);
	    mstack[cur]->SetProperTime();

	    if (local_pMax[i] > 0) {
		//do the resample

		if (local_version[i] == 0) {
		    //Cone sampling
		    p         = local_pMin[i] + (PUtils::sampleFlat() * (local_pMax[i] - local_pMin[i]));
		    theta     = acos(1. - ((local_openingAngle[i] + 1) * PUtils::sampleFlat()));
		    phi       = 2. * TMath::Pi() * PUtils::sampleFlat();
		    sin_theta = sin(theta);

		    mstack[cur]->SetPx(p * cos(phi) * sin_theta);
		    mstack[cur]->SetPy(p * sin(phi) * sin_theta);
		    mstack[cur]->SetPz(p * cos(theta));

		    mstack[cur]->RotateY(local_theta[i]);
		    mstack[cur]->RotateZ(local_phi  [i]);

		    if (local_mMax[i] < 0.) {
			mstack[cur]->SetM(local[i]->M());
		    } else {
			m = local_mMin[i] + (PUtils::sampleFlat() * (local_mMax[i] - local_mMin[i])) ;
			mstack[cur]->SetM(m);
		    }
		}
		else if (local_version[i] == 1) {
		    //Sector sampling
		    cos_theta = cos(local_thetaMin[i]) + (PUtils::sampleFlat() * (cos(local_thetaMax[i]) - cos(local_thetaMin[i])));

		    p         = local_pMin  [i] + (PUtils::sampleFlat() * (local_pMax  [i] - local_pMin  [i]));
		    phi       = local_phiMin[i] + (PUtils::sampleFlat() * (local_phiMax[i] - local_phiMin[i]));
		    theta     = acos(cos_theta);
		    sin_theta = sin(theta);

		    mstack[cur]->SetPx(p * cos(phi) * sin_theta);
		    mstack[cur]->SetPy(p * sin(phi) * sin_theta);
		    mstack[cur]->SetPz(p * cos_theta);

		    if (local_mMax[i] < 0.) {
			mstack[cur]->SetM(local[i]->M());
		    } else {
			mstack[cur]->SetM(makeStaticData()->GetParticleMass(mstack[cur]->ID()));
		    }
		    //Vertex set by PParticle ctor

		}
	    }
	    cur ++;

	    if (cur == stacksize) {
		return kFALSE;
	    }

	} else last_downscaling[i] ++; //wait for next time...
    }

    *num = cur;

    return kTRUE;
}

ClassImp(PEmbeddedParticles)
