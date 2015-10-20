/////////////////////////////////////////////////////////////////////
//Plugin for (rare) eta decays
//
//The complex version of the eta Double Dalitz (with the angles)
//is based on a calculation by A. Wirzba and T. Petri (diploma theses)
//
//The mod is automatically enabled, i.e. each PReaction can use it
//without activating it.
//
//However, it can be activated manually via:
//makeDistributionManager()->Exec("eta_decays");
//
//
//eta -> pi+ pi- e+ e-
//--------------------
//This decay uses the "genbod" technique:
//The 4-body decay is done via genbod, and the mass distribution
//  
//  M(q^2) = F(q^2) / (8*q^2)
//
//is folded into the genbod model. For the form factor by default
//the "simple VMD model" based on PSimpleVMDFF is used:
//
//  F(q^2) = m_v^2 / (m_v^2 - q^2) with m_v = 0.77GeV
//
//In total, the following matrix element is used:
//
//  A^2 = M(q^2) * s_pipi * beta_pi^2 sin^2(theta_pi)
//        * lambda(m_eta^2,s_pipi,q^2) (1+beta_e^2*sin^2(theta_e)*sin^2(phi))
//
//The reaction itself can be contructed with (e.g.) the following line:
//PReaction my_reaction("_T1=2.2","p","d","He3 eta [e+ e- pi+ pi-]","eta_pi_pi_dilepton");
//
//Interaction models:
// [eta_ee_pipi_vmd_ff]  VMD form factor (optional)
// [eta_ee_pipi_mass]    Dilepton mass (mandatory)
// [eta_ee_pipi_genbod]  Genbod (mandatory)
// [eta_ee_pipi]         Decay angles  (optional)
//
//
//eta -> e+ e- e+ e-
//------------------
//
//This decay goes via 2 virtual photons:
//
//PReaction my_reaction("_T1=2.2","p","d","He3 eta[dilepton [e+ e-] dilepton [e+ e-]]","etae+e-e+e-");
//
//In the first step, the sampling of the 2 dilepton masses is done. This is
//based on a source code provided by A. Kupsc and implemented in the model [eta_double_dalitz_simple]
//
//In the second step, the angular distributions are sampled via the rejection method
//[eta_double_dalitz_complex]. This distribution is based on the calculation of T. Petri and A. Wirzba
//(Kinematics and matrix elements of the decay eta -> ee mumu)
//
//The angular distributions and the form factors can be disabled with:
//   makeDistributionManager()->Disable("eta_double_dalitz_complex");
//   makeDistributionManager()->Disable("eta_double_dalitz_vmd_ff");
//
//The form factors
//----------------
//
//All form factor models can be obtained via:
//PSimpleVMDFF *ff = (PSimpleVMDFF*)makeDistributionManager()->GetDistribution("...")
//with ... and the distribution name (see above). The equation of the form factor
//can be exchanged via:
//   ff->AddEquation("_ff2 = 1.");
//for details see PSimpleVMDFF
//
//eta -> gamma e+ e-
//------------------
//Not a rare decay, but the distribution of the e+e- decay angle is
//implemented via this class as well. A form factor is not yet
//included by default but van be attached via a compiled class
//(see demo macro)
//
//                             Author:   I. Froehlich
//                             Written:  17.09.2009
//                             Released: 08.12.2010
//                           
//////////////////////////////////////////////////////////////////////

#include "PDataBase.h"
#include "PEtaDecaysPlugin.h"
#include "PSimpleVMDFF.h"
#include "PEtaPiPiDilepton.h"
#include "PEtaPiPiDileptonMass.h"

#include "PGenBod.h"

//#include "PTCrossWeight.h"


PEtaDecaysPlugin::PEtaDecaysPlugin() {
}


PEtaDecaysPlugin::PEtaDecaysPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {

    eta_dd_simple = NULL;
    eta_dd_ff     = NULL;
    
    eta_pipi_gamma = NULL;

}

PEtaDecaysPlugin::~PEtaDecaysPlugin() {
}

Bool_t PEtaDecaysPlugin::Activate(void) {

    return kTRUE;

};


Bool_t PEtaDecaysPlugin::ExecCommand(const char * command, Double_t value) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    //called the 1st time?
    if (!eta_pipi_gamma) {
	//add some groups in the distribution manager
	pdmutil->AddSubGroup("rare_eta_decays", "Rare eta decays", "eta_physics");
	pdmutil->SetGroup("eta_physics");

	//
	// eta -> pipi gamma (not rare)
	//
	
	//Standalone models, up to now nothing special
	eta_pipi_gamma = 
	    new PEtaPiPiGamma("eta_pipi_gamma_matrix_weighting@eta_to_pi+_pi-_g/matrix",
			      "Matrix element eta -> pipi gamma (Weighting)",-1);
	
	eta_pipi_gamma->EnableWeighting();
	makeDistributionManagerUtil()->Add(eta_pipi_gamma);
	eta_pipi_gamma = 
	    new PEtaPiPiGamma("eta_pipi_gamma_matrix@eta_to_pi+_pi-_g/matrix",
			      "Matrix element eta -> pipi gamma",-1);
	makeDistributionManagerUtil()->Add(eta_pipi_gamma);
	
	//
	// eta -> e+e-e+e- (rare)
	//

	Int_t ipid[5];

	//Add dilepton decay, if not yet present:
	ipid[0]=makeStaticData()->GetParticleID("eta");
	ipid[1]=makeStaticData()->GetParticleID("dilepton");
	ipid[2]=makeStaticData()->GetParticleID("dilepton");
	
	if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	    makeStaticData()->AddDecay(-1,"eta -> dilepton + dilepton", 
				       "eta", "dilepton,dilepton",ETA_DOUBLE_DALITZ_BR);
	
	//Create and add simple model:
	eta_dd_simple = new PEtaDoubleDalitz("eta_double_dalitz_simple@eta_to_dilepton_dilepton",
					     "Dilepton generator for eta -> dilepton + dilepton",-1);
	
	makeDistributionManagerUtil()->SetGroup("rare_eta_decays");
	makeDistributionManagerUtil()->Add(eta_dd_simple);

	//the more complex version uses the envelope method
	//of Pluto, the 4 leptons are rejected after the full
	//chain has been sampled

	eta_dd_complex = new PEtaDoubleDalitzEnv("eta_double_dalitz_complex@eta_to_dilepton_dilepton/full",
						 "Full dilepton generator for eta -> e+ e- e+ e-",-1);
	eta_dd_complex->Add("eta,parent");
	eta_dd_complex->Add("dilepton,daughter");
	eta_dd_complex->Add("dilepton,daughter");
	eta_dd_complex->Add("e+,granddaughter"); // <--- the order is very important
	eta_dd_complex->Add("e-,granddaughter");
	eta_dd_complex->Add("e+,granddaughter");
	eta_dd_complex->Add("e-,granddaughter");

	makeDistributionManagerUtil()->Add(eta_dd_complex);

	PSimpleVMDFF *ff = new PSimpleVMDFF("eta_double_dalitz_vmd_ff@eta_to_dilepton_dilepton/formfactor",
					    "Simple VMD form factor for eta Double Dalitz",-1);
	ff->SetVectorMesonMass(0.77);
	ff->SetWeightMax(5.);
	
	makeDistributionManagerUtil()->Add(ff);

	//
	// eta -> e+e-pi+pi- (rare)
	// 
	//
	
	
	//Add  decay, if not yet present:
	ipid[0]=makeStaticData()->GetParticleID("eta");
	ipid[1]=makeStaticData()->GetParticleID("e+");
	ipid[2]=makeStaticData()->GetParticleID("e-");
	ipid[3]=makeStaticData()->GetParticleID("pi+");
	ipid[4]=makeStaticData()->GetParticleID("pi-");
	
	if (makeStaticData()->GetDecayKey(ipid, 4) < 0)
	    makeStaticData()->AddDecay(-1,"eta -> e+ + e- + pi+ + pi-", 
				       "eta", "e+,e-,pi+,pi-",ETA_EE_PIPI_BR);

	//first, the standard ff
	ff = new PSimpleVMDFF("eta_ee_pipi_vmd_ff@eta_to_e+_e-_pi+_pi-/formfactor",
			      "eta -> pi pi dilepton VMD form factor",-1);
	ff->SetVectorMesonMass(0.77);
	makeDistributionManagerUtil()->Add(ff);


	//mass distribution of the di-lepton mass. Here we use the
	//function to be folded into the genbod model
	PEtaPiPiDileptonMass *main = 
	    new PEtaPiPiDileptonMass("eta_ee_pipi_mass@eta_to_e+_e-_pi+_pi-/correlation",
				     "eta -> pi pi dilepton a la Wirzba (dilepton mass)",-1);
	makeDistributionManagerUtil()->Add(main);

	//...and now the adapted genbod:
	PGenBod *genbod = new PGenBod("eta_ee_pipi_genbod@eta_to_e+_e-_pi+_pi-/genbod",
				      "eta -> pi pi dilepton genbod",-1);
	genbod->Add("eta,parent");
	genbod->Add("pi+,daughter");
	genbod->Add("pi-,daughter");
	genbod->Add("e+,daughter,corr1");
	genbod->Add("e-,daughter,corr2");
	makeDistributionManagerUtil()->Add(genbod);

	//main distribution (afterburner to take 4-body correlation into account)
	PEtaPiPiDilepton *final = new PEtaPiPiDilepton("eta_ee_pipi@eta_to_e+_e-_pi+_pi-",
						       "eta -> pi pi dilepton a la Wirzba (decay angles)",-1);
	makeDistributionManagerUtil()->Add(final);

    }

    return kTRUE;
}

ClassImp(PEtaDecaysPlugin)



