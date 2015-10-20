/////////////////////////////////////////////////////////////////////
//  Pluto Standard Models
//
//
//                             Author:  I. Froehlich
//                             Written: 27.6.2007
//                             Revised: 
//////////////////////////////////////////////////////////////////////

#include "PDataBase.h"
#include "PStdModels.h"

#include "PHadronDecayM1.h"
#include "PHadronDecayM2.h"
#include "PHadronDecayM3.h"
#include "PHadronDecayM1N.h"
#include "PEEDirectDecay.h"
#include "PBreitWigner.h"
#include "PFixedDecay.h"
#include "PFixedProduction.h"
#include "PGenBod.h"
#include "PFermiMomentum.h"
#include "PPropagator.h"

PStdModels::PStdModels() {

    generic_physics_done = kFALSE;

    Info("PStdModels()","(%s), Standard model filler", PRINT_CTOR);

}

PStdModels::~PStdModels() {
}


void PStdModels::Add(PDistributionManagerUtil * pdmutil) {

	
    if (!generic_physics_done) {
	generic_physics_done=kTRUE;
	//Include the standard Pluto distributions
	GenericPhysics(pdmutil);
    } else {

    pdmutil->Add(GetModels());

    }



}


double f_eta_decay(double *x, double * par)
{
  // QED helicity angle (Bratkovskaya, Toneev et al.)
    return (1+x[0]*x[0])/2;
}

Double_t f_delta_decay(Double_t *x, Double_t * par)
{
//    return (1+3*x[0]*x[0])/4; 
//    return (1+1.35*x[0]*x[0])/4.; 
    return (1+1.35*x[0]*x[0])/2.35; //faster 
}

Double_t f_delta_decay2(Double_t *x, Double_t * par)
{
    return (5.-3.*(1.35/3.)*x[0]*x[0])/5; //(1.35/3.) is the damping factor
}

void PStdModels::GenericPhysics(PDistributionManagerUtil * pdmutil) {
    //First, add a set of groups, which helps to switch 
    //off certain aspects of physics

    pdmutil->AddGroup("root", "Root group");

    pdmutil->AddSubGroup("eta_physics", "Physics about eta production, and decay", "root");
    //AddGroup("dilepton_physics", "Physics about the decay of dileptons"); //at the moment not more...
    pdmutil->AddSubGroup("helicity_angles", "Helicity angles of dileptons", "root");
    pdmutil->AddSubGroup("resonances_pw", "Partial waves of resonances", "root");
    pdmutil->AddSubGroup("particle_models", "Mass sampling of particles", "root");
    pdmutil->AddSubGroup("decay_models", "Phase space mass sampling & decay partial widths", "root");
    pdmutil->AddSubGroup("polar_angles", "Polar angles in elementary particle production", "root");

    pdmutil->AddSubGroup("genbod_models", "Momentum sampling", "root");
//    AddGroup("mass_sampling", "Sampling of dileptons and resonances");

//    SetGroup("mass_sampling");


//    Add((new PEtaModels())->GetModels());

    pdmutil->SetGroup("particle_models");

    //**** rho/omega propagator used e.g. for pion beam plugin
    PPropagator * Rho0Propagator = new 
	PPropagator("Rho0Propagator@rho0_prop/propagator",
		    "Complex rho0 propagator",-1);
  
    PPropagator * OmegaPropagator = new 
	PPropagator("OmegaPropagator@w_prop/propagator",
		    "Complex omega propagator",-1);

    pdmutil->Add(Rho0Propagator);
    pdmutil->Add(OmegaPropagator);

    //**** Complex BW for rho/omega interference 
    int ipid[11],decaykey;
    ipid[0]=makeStaticData()->GetParticleID("rho0");
    ipid[1]=makeStaticData()->GetParticleID("e+");
    ipid[2]=makeStaticData()->GetParticleID("e-");
    decaykey = makeStaticData()->GetDecayKey(ipid, 2);
    int pkey=makeStaticData()->GetParticleKey("rho0");
    PComplexBreitWigner *pmodel2 = new PComplexBreitWigner("rho0_pionff",
							   "Pion form factor",pkey);
    pmodel2->Add("rho0, parent");
    pmodel2->AddAmplitude(makeStaticData()->GetDecayIdxByKey(decaykey),1,0);

    ipid[0]=makeStaticData()->GetParticleID("w");
    int decaykey_omega = makeStaticData()->GetDecayIdx(ipid, 2);

    pmodel2->AddInterference(makeStaticData()->GetDecayIdxByKey(decaykey),
			     makeStaticData()->GetParticleKey("w"),decaykey_omega,0.0388954,-1.60015);

    pdmutil->Add(pmodel2);


    pkey=makeStaticData()->GetParticleKey("w");
    PComplexBreitWigner *pmodel3 = new PComplexBreitWigner("w_pionff",
							   "Pion form factor",pkey);
    pmodel3->Add("w, parent");
    pdmutil->Add(pmodel3);



    


    //////////////Standard models after THIS point
    pdmutil->Add(GetModels());
    

    pdmutil->SetGroup("decay_models");
    //Alternative models for the rho
    
    ipid[0]=makeStaticData()->GetParticleID("rho0");
    ipid[1]=makeStaticData()->GetParticleID("e+");
    ipid[2]=makeStaticData()->GetParticleID("e-");
//    cout << ipid[0] << ":" << ipid[1] << ":" << ipid[2] << endl;
    decaykey = makeStaticData()->GetDecayKey(ipid, 2);

    PEEDirectDecay *pmodel = new PEEDirectDecay("rho_picutoff_e-_e+",
						"Dilepton direct decay with pion cutoff",decaykey);
    pmodel->SetPiCutoff(1);
    pmodel->Add("rho0, parent");
    pmodel->Add("e+, daughter");
    pmodel->Add("e-, daughter");
    pdmutil->Add(pmodel);

    ipid[1]=makeStaticData()->GetParticleID("mu+");
    ipid[2]=makeStaticData()->GetParticleID("mu-");

    decaykey = makeStaticData()->GetDecayKey(ipid, 2);
    
    PEEDirectDecay *pmodel_mumu = new PEEDirectDecay("rho_picutoff_mu-_mu+",
						     "Dimuon direct decay with pion cutoff",decaykey);
    pmodel_mumu->SetPiCutoff(1);
    pmodel_mumu->Add("rho0, parent");
    pmodel_mumu->Add("mu+, daughter");
    pmodel_mumu->Add("mu-, daughter");
    pdmutil->Add(pmodel_mumu);



//    linkDB();
    
    //
    // Eta Physics
    //

    pdmutil->SetGroup("eta_physics");

    //Eta polar angle
    //Ref.17
    TF2 *eta_angles=new TF2("eta_angles","(1+(3.74041e+01-2.76688e+01*y+5.07488e+00*y*y)*0.5*(3*x*x-1))/2 ",-1,1,0,10);
    gROOT->GetListOfFunctions()->Remove(eta_angles);
    PAngularDistribution * pp_eta_prod_angle = 
	new PAngularDistribution("pp_eta_prod_angle",
				 "Eta polar angles in pp reactions for direct production");
    pp_eta_prod_angle->Add("eta,  daughter, primary");
    pp_eta_prod_angle->Add("p,    daughter");
    pp_eta_prod_angle->Add("p,    daughter");
    pp_eta_prod_angle->Add("q,    parent,   reference");
    pp_eta_prod_angle->SetRotate(kFALSE);
    pp_eta_prod_angle->SetAngleFunction(eta_angles);
    pdmutil->Add(pp_eta_prod_angle);
    

    PAngularDistribution * pp_ns_eta_prod_angle = 
	new PAngularDistribution("pp_ns_eta_prod_angle", 
				 "Eta polar angles in pp reactions via N*");
    pp_ns_eta_prod_angle->Add("eta,  daughter,    primary");
    pp_ns_eta_prod_angle->Add("p,    daughter");
    pp_ns_eta_prod_angle->Add("NS11+,parent");
    pp_ns_eta_prod_angle->Add("p,    parent,      sibling"); //the missing p is the sibling of the N*
    pp_ns_eta_prod_angle->Add("q,    grandparent, reference");
    pp_ns_eta_prod_angle->SetAngleFunction(eta_angles);
    pdmutil->Add(pp_ns_eta_prod_angle);

    //pp alignment for eta production
    TF2 *pp_angles2=new TF2("pp_angles2","(1+(5.037074-4.537543*y+1.010214*y*y)*0.5*(3*x*x-1))/2 ",-1,1,0,10);
    gROOT->GetListOfFunctions()->Remove(pp_angles2);
    PAngularDistribution * pp_eta_pp_align = 
	new PAngularDistribution("pp_eta_pp_align",
				 "pp alignment in pp reactions for direct eta production");
    pp_eta_pp_align->Add("eta,  daughter");
    pp_eta_pp_align->Add("p,    daughter,    primary");
    pp_eta_pp_align->Add("p,    daughter,    align");
    pp_eta_pp_align->Add("q,    parent,      mass_reference");
    pp_eta_pp_align->SetRotate(kFALSE);//DISTO measured just with boost, no rotation
    pp_eta_pp_align->SetAngleFunction(pp_angles2);
    pdmutil->Add(pp_eta_pp_align);
    
    PAngularDistribution * pp_ns_eta_pp_align = 
	new PAngularDistribution("pp_ns_eta_pp_align",
				 "pp alignment in pp reactions for eta production via N*");
    pp_ns_eta_pp_align->Add("eta,    daughter");
    pp_ns_eta_pp_align->Add("p,      daughter,    primary");
    pp_ns_eta_pp_align->Add("NS11+,  parent");
    pp_ns_eta_pp_align->Add("p,      parent,      sibling, align"); //the missing p is the sibling of the N*
    pp_ns_eta_pp_align->Add("q,      grandparent, mass_reference");
    pp_ns_eta_pp_align->SetRotate(kFALSE);//DISTO measured just with boost, no rotation
    pp_ns_eta_pp_align->SetAngleFunction(pp_angles2);
    pdmutil->Add(pp_ns_eta_pp_align);

    //matrix element for eta -> pi+pi-pi0
    PDalitzDistribution * eta_hadronic_decay =
	new PDalitzDistribution("eta_hadronic_decay",
				"Eta matrix element for decay into charged pions");
    eta_hadronic_decay->Add("eta,    parent");
    eta_hadronic_decay->Add("pi0,    daughter,    primary");
    eta_hadronic_decay->Add("pi+,    daughter");
    eta_hadronic_decay->Add("pi-,    daughter");
    eta_hadronic_decay->SetSlopes(-0.94,0.11);
    eta_hadronic_decay->SetMax(2.05);
    pdmutil->Add(eta_hadronic_decay);


    //
    // Anisotropic dilepton decay
    //

    pdmutil->SetGroup("helicity_angles");
    //Ref.7
    TF1 * pseudoscalar_decay = new TF1("helicity",f_eta_decay,-1,1,1);
    gROOT->GetListOfFunctions()->Remove(pseudoscalar_decay);
    PAngularDistribution * eta_dilepton_helicity = 
	new PAngularDistribution("eta_dilepton_helicity",
				 "Helicity angle of the dilepton decay of eta");
    eta_dilepton_helicity->Add("dilepton","PARENT");
    eta_dilepton_helicity->Add("e+",  "DAUGHTER");
    eta_dilepton_helicity->Add("e-",  "DAUGHTER", "primary");
    eta_dilepton_helicity->Add("g", "PARENT", "SIBLING"); //added to distinguish double-dalitz
    eta_dilepton_helicity->Add("eta", "GRANDPARENT", "base_reference");
    eta_dilepton_helicity->SetAngleFunction(pseudoscalar_decay);
    eta_dilepton_helicity->NeverAbort(kTRUE);
    pdmutil->Add(eta_dilepton_helicity);
    
    PAngularDistribution * etaprime_dilepton_helicity = 
	new PAngularDistribution("etaprime_dilepton_helicity",
				 "Helicity angle of the dilepton decay of etaprime");
    etaprime_dilepton_helicity->Add("dilepton","PARENT");
    etaprime_dilepton_helicity->Add("e+",  "DAUGHTER");
    etaprime_dilepton_helicity->Add("e-",  "DAUGHTER", "primary");
    etaprime_dilepton_helicity->Add("eta'", "GRANDPARENT", "base_reference");
    etaprime_dilepton_helicity->SetAngleFunction(pseudoscalar_decay);
    etaprime_dilepton_helicity->NeverAbort(kTRUE);
    pdmutil->Add(etaprime_dilepton_helicity);
    
    PAngularDistribution * pi0_dilepton_helicity =
	new PAngularDistribution("pi0_dilepton_helicity",
				 "Helicity angle of the dilepton decay of pi0");
    pi0_dilepton_helicity->Add("dilepton","PARENT");
    pi0_dilepton_helicity->Add("e+",  "DAUGHTER");
    pi0_dilepton_helicity->Add("e-",  "DAUGHTER", "primary");
    pi0_dilepton_helicity->Add("pi0", "GRANDPARENT", "base_reference");
    pi0_dilepton_helicity->SetAngleFunction(pseudoscalar_decay);
    pi0_dilepton_helicity->NeverAbort(kTRUE);
    pdmutil->Add(pi0_dilepton_helicity);

    //Resonances

    pdmutil->SetGroup("resonances_pw");


    TF1 * pw_decay = new TF1("pw",f_delta_decay,-1,1,1);
    gROOT->GetListOfFunctions()->Remove(pw_decay);
    PScatterDistribution * delta_waves1 = 
	new PScatterDistribution("pp_delta_waves1","Delta+ PW pi0 scattering");
    delta_waves1->Add("D+","PARENT");
    delta_waves1->Add("pi0",  "DAUGHTER", "primary");
    delta_waves1->Add("p",    "DAUGHTER");
    delta_waves1->Add("q"   , "GRANDPARENT");
    delta_waves1->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves1->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves1->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves1);

    PScatterDistribution * delta_waves2 = 
	new PScatterDistribution("pp_delta_waves2","Delta+ PW pi+ scattering");
    delta_waves2->Add("D+","PARENT");
    delta_waves2->Add("pi+",  "DAUGHTER", "primary");
    delta_waves2->Add("n",    "DAUGHTER");
    delta_waves2->Add("q"   , "GRANDPARENT");
    delta_waves2->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves2->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves2->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves2);

    PScatterDistribution * delta_waves3 = 
	new PScatterDistribution("pp_delta_waves3","Delta++ PW pi+ scattering");
    delta_waves3->Add("D++","PARENT");
    delta_waves3->Add("pi+",  "DAUGHTER", "primary");
    delta_waves3->Add("p",    "DAUGHTER");
    delta_waves3->Add("q"   , "GRANDPARENT");
    delta_waves3->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves3->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves3->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves3);

    TF1 * pw_decay2 = new TF1("pw",f_delta_decay2,-1,1,1);
    gROOT->GetListOfFunctions()->Remove(pw_decay2);
    PScatterDistribution * delta_waves1dil = 
	new PScatterDistribution("pp_delta_waves1dil","Delta+ PW DiLepton scattering");
    delta_waves1dil->Add("D+","PARENT");
    delta_waves1dil->Add("dilepton",  "DAUGHTER", "primary");
    delta_waves1dil->Add("p",    "DAUGHTER");
    delta_waves1dil->Add("q"   , "GRANDPARENT");
    delta_waves1dil->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves1dil->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves1dil->SetAngleFunction(pw_decay2);
    pdmutil->Add(delta_waves1dil);


    PScatterDistribution * delta_waves1a = 
	new PScatterDistribution("pp_delta_waves1a","Delta0 PW pi0 scattering");
    delta_waves1a->Add("D0","PARENT");
    delta_waves1a->Add("pi0",  "DAUGHTER", "primary");
    delta_waves1a->Add("n",    "DAUGHTER");
    delta_waves1a->Add("q"   , "GRANDPARENT");
    delta_waves1a->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves1a->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves1a->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves1a);

    PScatterDistribution * delta_waves2a = 
	new PScatterDistribution("pp_delta_waves2a","Delta0 PW pi- scattering");
    delta_waves2a->Add("D0","PARENT");
    delta_waves2a->Add("pi-",  "DAUGHTER", "primary");
    delta_waves2a->Add("p",    "DAUGHTER");
    delta_waves2a->Add("q"   , "GRANDPARENT");
    delta_waves2a->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves2a->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves2a->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves2a);

    PScatterDistribution * delta_waves3a = new PScatterDistribution("pp_delta_waves3a","Delta- PW pi- scattering");
    delta_waves3a->Add("D-","PARENT");
    delta_waves3a->Add("pi-",  "DAUGHTER", "primary");
    delta_waves3a->Add("n",    "DAUGHTER");
    delta_waves3a->Add("q"   , "GRANDPARENT");
    delta_waves3a->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves3a->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves3a->SetAngleFunction(pw_decay);
    pdmutil->Add(delta_waves3a);


    PScatterDistribution * delta_waves2dil = 
	new PScatterDistribution("pp_delta_waves2dil","Delta0 PW DiLepton scattering");
    delta_waves2dil->Add("D0","PARENT");
    delta_waves2dil->Add("dilepton",  "DAUGHTER", "primary");
    delta_waves2dil->Add("n",    "DAUGHTER");
    delta_waves2dil->Add("q"   , "GRANDPARENT");
    delta_waves2dil->Add("N"   , "GRANDGRANDPARENT", "beam");
    delta_waves2dil->Add("N"   , "GRANDGRANDPARENT", "target");
    delta_waves2dil->SetAngleFunction(pw_decay2);
    pdmutil->Add(delta_waves2dil);

    pdmutil->SetGroup("polar_angles");
    TF1 * dummy = new TF1("dummy","1",-1,1);
    gROOT->GetListOfFunctions()->Remove(dummy);
    PDeltaAngularDistribution * delta_production = 
	new PDeltaAngularDistribution("NN_delta+_angle","N+N->N+Delta angular distribution");
    delta_production->Add("N,grandparent,beam");
    delta_production->Add("N,grandparent,target");
    delta_production->Add("q,parent");
    delta_production->Add("N,daughter");
    delta_production->Add("D+,daughter,primary");
    delta_production->SetAngleFunction(dummy);
    pdmutil->Add(delta_production);

    PDeltaAngularDistribution * delta_production2 = 
	new PDeltaAngularDistribution("NN_delta++_angle","N+N->N+Delta angular distribution");
    delta_production2->Add("N"   , "GRANDPARENT", "beam");
    delta_production2->Add("N"   , "GRANDPARENT", "target");
    delta_production2->Add("q"   , "PARENT");
    delta_production2->Add("N,   DAUGHTER");
    delta_production2->Add("D++, daughter, primary");
    delta_production2->SetAngleFunction(dummy);
    pdmutil->Add(delta_production2);

    PDeltaAngularDistribution * delta_production3 = 
	new PDeltaAngularDistribution("NN_delta0_angle","N+N->N+Delta angular distribution");
    delta_production3->Add("N"   , "GRANDPARENT", "beam");
    delta_production3->Add("N"   , "GRANDPARENT", "target");
    delta_production3->Add("q"   , "PARENT");
    delta_production3->Add("N,   DAUGHTER");
    delta_production3->Add("D0, daughter, primary");
    delta_production3->SetAngleFunction(dummy);
    pdmutil->Add(delta_production3);

    PSaid * said = new PSaid ("pp_elastic", "PP elastic scattering with SAID");
    said->Add("q"   , "PARENT");
    said->Add("p"   , "GRANDPARENT", "beam");
    said->Add("p"   , "GRANDPARENT", "target");
    said->Add("p, daughter");
    said->Add("p, daughter, primary");
    said->SetAngleFunction(dummy);
    pdmutil->Add(said);

    PPiOmegaAngularDistribution * omega_production = 
	new PPiOmegaAngularDistribution("pi+n_wp_angle", "pi+ + n --> w + p angular distribution");
    omega_production->Add("pi+"   , "GRANDPARENT", "beam");
    omega_production->Add("n"   , "GRANDPARENT", "target");
    omega_production->Add("q"   , "PARENT");
    omega_production->Add("p",    "DAUGHTER");
    omega_production->Add("w"  ,"DAUGHTER", "primary");
    omega_production->SetAngleFunction(dummy);
    omega_production->SetVersion(PI_OMEGA_piNNw);
    pdmutil->Add(omega_production);

    PPiOmegaAngularDistribution * omega_production2 = 
	new PPiOmegaAngularDistribution("pi-p_wn_angle", "pi- + p --> w + n angular distribution");
    omega_production2->Add("pi-"   , "GRANDPARENT", "beam");
    omega_production2->Add("p"   , "GRANDPARENT", "target");
    omega_production2->Add("q"   , "PARENT");
    omega_production2->Add("n",    "DAUGHTER");
    omega_production2->Add("w"  ,"DAUGHTER", "primary");
    omega_production2->SetAngleFunction(dummy);
    omega_production2->SetVersion(PI_OMEGA_piNNw);
    pdmutil->Add(omega_production2);


    PPiOmegaAngularDistribution * omega_production3 = 
	new PPiOmegaAngularDistribution("pi+p_pi+wp_angle", "pi+ + p -> pi+ + p + w angular distribution");
    omega_production3->Add("pi+"   , "GRANDPARENT", "beam");
    omega_production3->Add("p"   , "GRANDPARENT", "target");
    omega_production3->Add("q"   , "PARENT");
    omega_production3->Add("p",    "DAUGHTER");
    omega_production3->Add("pi+",    "DAUGHTER");
    omega_production3->Add("w"  ,"DAUGHTER", "primary");
    omega_production3->SetAngleFunction(dummy);
    omega_production3->SetVersion(PI_OMEGA_piPPpiw);
    pdmutil->Add(omega_production3);

    
    PPiOmegaAngularDistribution * omega_production4 = 
	new PPiOmegaAngularDistribution("pi+p_D++w_angle", "pi+ + p -> D++ + w angular distribution");
    omega_production4->Add("pi+"   , "GRANDPARENT", "beam");
    omega_production4->Add("p"   , "GRANDPARENT", "target");
    omega_production4->Add("q"   , "PARENT");
    omega_production4->Add("D++",    "DAUGHTER");
    omega_production4->Add("w"  ,"DAUGHTER", "primary");
    omega_production4->SetAngleFunction(dummy);
    omega_production4->SetVersion(PI_OMEGA_piPDw);
    pdmutil->Add(omega_production4);


    //some updates on the vector mesons mass shape:
    pdmutil->Add(new PHadronDecayM3("w_phasespace_pi+_pi-_pi0","Width including 3-body p.s.",-1));
    pdmutil->Add(new PHadronDecayM3("phi_phasespace_pi+_pi-_pi0","Width including 3-body p.s.",-1));

    PEEDirectDecay * ee = (PEEDirectDecay *) 
	pdmutil->GetDistribution("w_ee_e-_e+");
    if (ee) ee->SetHadronicPS(1);
    else Warning("GenericPhysics","[w_ee_e-_e+] not found");
    ee = (PEEDirectDecay *) pdmutil->GetDistribution("phi_ee_e-_e+");
    if (ee) ee->SetHadronicPS(1);
    else Warning("GenericPhysics","[phi_ee_e-_e+] not found");

    pdmutil->ExpandGroup("user");

};





TObjArray * PStdModels::GetModels(void) {

    arr=new TObjArray();

    Info("GetModels","Read std models");

    /***********************************************
     * Loop over the particles                     *
     * for each particle we check one of the       *
     * build-in models                             *
     ***********************************************/
    Int_t header=makeDataBase()->GetEntry("std_set");
    Int_t particlekey=-1;
    Int_t genbod_key=makeStaticData()->MakeDirectoryEntry("modeldef",NMODEL_NAME,LMODEL_NAME,"genbod");

    Double_t unstable_width = 
	*(makeStaticData()->GetBatchValue("_system_unstable_width"));

    while (makeDataBase()->MakeListIterator(header, "snpart","slink" , &particlekey)) {
	Int_t pid = makeStaticData()->GetParticleIDByKey(particlekey);

	//cout << "pid:" << pid << endl;
	Int_t decaykey=-1;
	

	while ((makeStaticData()->GetParticleNChannelsByKey(particlekey)>0) && 
	       (makeDataBase()->MakeListIterator(particlekey, "pnmodes", "link", &decaykey))) {


	    /***********************************************
	     * Loop over the decay modes                   *
	     * for each mode we check one of the           *
	     * build-in models                             *
	     ***********************************************/

	    int decay_has_composite=0;
	    Int_t tid[11];
	    tid[0]=10; 
	    makeStaticData()->GetDecayModeByKey(decaykey,tid); // retrieve current mode info
	    
	    for (int p=1;p<=tid[0];p++) {
		if (tid[p]>1000)
		    decay_has_composite=1;
	    }
	    //cout << decaykey << endl;
	    if (!makeDynamicData()->GetDecayModelByKey(decaykey)) { 
		//Check first if we are alraedy done

		model=NULL;

		if ((tid[0] == 2) && (!decay_has_composite)) {
		    /* 2 decay products */
		    
		    /* Dalitz Decays */
		    if (PData::IsDalitz(pid,tid[1],tid[2])) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_dalitz");
			if (makeStaticData()->IsParticle(tid[1],"dimuon")||
			    makeStaticData()->IsParticle(tid[2],"dimuon"))
			    id->Append("_mumu");
			PDalitzDecay *pmodel = 
			    new PDalitzDecay((char*)id->Data(),"Dalitz decay",decaykey);
			model=(PChannelModel*) pmodel;
			//	makeDataBase()->SetParamInt (decaykey, "maxmesh",new Int_t(2000)); //take into account strong var.
		    }

		    int n1=1,n2=1;
		    if (makeStaticData()->
			GetParticleTotalWidth(tid[1]) < unstable_width) 
			n1 = 0; // too narrow 
		    
		    if (makeStaticData()->
			GetParticleTotalWidth(tid[2]) < unstable_width ) 
			n2 = 0; // too narrow

		    int h=makeStaticData()->IsParticleHadron(tid[1])*
			makeStaticData()->IsParticleHadron(tid[2]);

		    /* Decay of a hadron in 1 stable and 1 unstable product */
		    if ((h && n1 && !n2) || (h && !n1 && n2)) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_m1_");
			id->Append(makeStaticData()->GetParticleName(tid[1]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[2]));
			PHadronDecayM1 *pmodel = 
			    new PHadronDecayM1((char*)id->Data(),
					       "1 unstable hadron (2-body ps)",decaykey);
			model=(PChannelModel*) pmodel;
		    }
		    /* Decay of a hadron in unstable products */
		    if ((h && n1 && n2) || (h && n1 && n2)) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_m2_");
			id->Append(makeStaticData()->GetParticleName(tid[1]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[2]));
			PHadronDecayM2 *pmodel = 
			    new PHadronDecayM2((char*)id->Data(),
					       "2 unstable hadrons (2-body ps)",decaykey);
			model=(PChannelModel*) pmodel;
		    }

		    /* Decay in stable products */
		    if (h && !n1 && !n2) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_fix_");
			id->Append(makeStaticData()->GetParticleName(tid[1]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[2]));
			PHadronDecay *pmodel = 
			    new PHadronDecay((char*)id->Data(),
					     "2-body fixed mass, partial width",decaykey);
			model=(PChannelModel*) pmodel;
		    }
		    /* Decay in ee */
		    if ((PData::IsDirectEE(pid,tid[1],tid[2])||
			PData::IsDirectMuMu(pid,tid[1],tid[2]))
			&& (!makeStaticData()->IsParticle(pid, "eta")
			    && !makeStaticData()->IsParticle(pid, "J/Psi") 
			    && !makeStaticData()->IsParticle(pid, "Psi'") 
			    && !makeStaticData()->IsParticle(pid, "sigma") //Not included  
			)) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_ee_");
			id->Append(makeStaticData()->GetParticleName(tid[1]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[2]));
			PEEDirectDecay *pmodel = 
			    new PEEDirectDecay((char*)id->Data(),
					       "Dilepton direct decay",decaykey);
			model=(PChannelModel*) pmodel;
			makeDataBase()->SetParamInt (particlekey, "maxmesh",new Int_t(2000)); //-> take into accoutn strong variations
			makeDataBase()->SetParamInt (decaykey, "maxmesh",new Int_t(2000)); 
		    }


	    

		} /* END 2 products */

		if ((tid[0] == 3) && (!decay_has_composite)) {
		    int h=makeStaticData()->IsParticleHadron(tid[1])*
			makeStaticData()->IsParticleHadron(tid[2])*
			makeStaticData()->IsParticleHadron(tid[3]);
		    /* 3 decay products */

		    Int_t ust=0;
		    if (makeStaticData()->
			GetParticleTotalWidth(tid[1]) > unstable_width) ust++; 
		    if (makeStaticData()->
			GetParticleTotalWidth(tid[2]) > unstable_width ) ust++; 
		    if (makeStaticData()->
			GetParticleTotalWidth(tid[3]) > unstable_width ) ust++; 

		    if (h && (ust)) {

			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_m3_");
			id->Append(makeStaticData()->GetParticleName(tid[1]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[2]));
			id->Append("_");
			id->Append(makeStaticData()->GetParticleName(tid[3]));
			PHadronDecayM3 *pmodel = 
			    new PHadronDecayM3((char*)id->Data(),"3-body phase space",decaykey);
			model=(PChannelModel*) pmodel;
		    }
		}

		// N>3 models
		if ((tid[0] > 3) && (!decay_has_composite)) {
		    int n1=0;
		    for (int j=0;j<tid[0];j++) {
			if ((makeStaticData()->
			    GetParticleTotalWidth(tid[1+j])> unstable_width) &&
			    makeStaticData()->IsParticleHadron(tid[j+1])) n1++;
		    }
		    if (n1==1) {
			id = new TString(makeStaticData()->GetParticleName(pid));
			id->Append("_m1n_");
			for (int p=1;p<=tid[0];p++) {
			    id->Append(makeStaticData()->GetParticleName(tid[p]));
			    if (p!=tid[0]) id->Append("_");
			}
			PHadronDecayM1N *pmodel = 
			    new PHadronDecayM1N((char*)id->Data(),"N-body phase space with unstable hadron ",decaykey);
			model=(PChannelModel*) pmodel;
		    }

		}
    
		if (!model && !decay_has_composite && tid[0]>1) {
		    //Add dummy model			
		    id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_fixed_");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("_");
		    }

		    PFixedDecay *pmodel = 
			new PFixedDecay((char*)id->Data(),"Fixed product masses",decaykey);
		    model=(PChannelModel*) pmodel;
		    model->SetGroupID("decay_models");
		}

		if (!model && !decay_has_composite && tid[0] == 1) {
		    //beam fusion
		    id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_prod_");
		    
		    id->Append(makeStaticData()->GetParticleName(tid[1]));
		
		    id->Append("@");
		    id->Append(makeStaticData()->GetParticleName(pid));
		    id->Append("#prod#");
		
		    id->Append(makeStaticData()->GetParticleName(tid[1]));
		
		    PFixedProduction *pmodel = 
		      new PFixedProduction((char*)id->Data(),"Particle Production (fusion)",-2);

		    //pmodel->Add(makeStaticData()->GetParticleName(pid),"parent");
		    //pmodel->Add(makeStaticData()->GetParticleName(tid[1]),"daughter");
		    model=(PChannelModel*) pmodel;
		    model->SetGroupID("decay_models");
		}

		if (model) {
		    //If a model has been identified we make it known to the PChannel world
		    if (pid<1000) model->Add(makeStaticData()->GetParticleName(pid),"parent");
		    else model->Add("q","parent");  //"White" decay of composite
		    for (int num=0;num<tid[0];num++)
			model->Add(makeStaticData()->GetParticleName(tid[1+num]),"daughter");
		    model->SetGroupID("decay_models");
		    arr->Add((TObject *)model);
		} 


		if ((tid[0] == 2) && (decay_has_composite) && ( PData::IsN(tid[1]) ||  PData::IsN(tid[2]) )) {
		    //add fermi model
		    id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_fermi_");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("_");
		    }
		    PFermiMomentum *pmodel = new PFermiMomentum(
			(char*)id->Data(),
			"Quasi-free particle production",-1);
		    if (pid<1000) pmodel->Add(makeStaticData()->GetParticleName(pid),"parent");
		    else {
			pmodel->Add("q","parent");  //"White" decay of composite
			pmodel->Add(makeStaticData()->GetParticleName(pid%1000),"grandparent","beam");
			pmodel->Add(makeStaticData()->GetParticleName(pid/1000),"grandparent","target");
		    }
		    for (int num=0;num<tid[0];num++) {
			if (tid[1+num]<1000) 
			    pmodel->Add(makeStaticData()->GetParticleName(tid[1+num]),"daughter","spectator");
			else {
			    //This is the quasi-reaction
			    pmodel->Add("q","daughter","composite"); 
			    pmodel->Add(makeStaticData()->GetParticleName(tid[1+num]/1000),"granddaughter","p1");
			    pmodel->Add(makeStaticData()->GetParticleName(tid[1+num]%1000),"granddaughter","p2");
			}

		    }

		    pmodel->SetGroupID("genbod_models");
		    arr->Add((TObject *)pmodel);
		    
		}

		    
	    } //END check for existing primary model

	    if (!makeDynamicData()->GetDecayModelByKey(decaykey,genbod_key)) { 
		if (!decay_has_composite && tid[0]>1) {
		    //Adding PGENBOD by default
		    id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_genbod_");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("_");
		    }
		    id->Append("@");
		    id->Append(makeStaticData()->GetParticleName(pid));
		    id->Append("#genbod#");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("#");
		    }
		    id->Append("&genbod");
		    PGenBod *pmodel = new PGenBod((char*)id->Data(),"Pluto build-in genbod",-2);
		    //if (pid<1000) 
		    pmodel->Add(makeStaticData()->GetParticleName(pid),"parent");
		    //else pmodel->Add("q","parent");  //"White" decay of composite
		    for (int num=0;num<tid[0];num++)
			pmodel->Add(makeStaticData()->GetParticleName(tid[1+num]),"daughter");
		    pmodel->SetGroupID("genbod_models");
		    arr->Add((TObject *)pmodel);

// 		    //Now check for a possible "reflected" particle
// 		    if (pid > PLUTO_COMPOSITE) {
// 			//parent is composite, add the reflected [clone]
// 			int id1=pid%1000, id2=pid/1000;
// 			Int_t xpid = id1*1000 + id2;
// 			if (makeStaticData()->IsParticleValid(xpid)) {
// 			    id = new TString(makeStaticData()->GetParticleName(xpid));
// 			    id->Append("_genbod_");
// 			    for (int p=1;p<=tid[0];p++) {
// 				id->Append(makeStaticData()->GetParticleName(tid[p]));
// 				if (p!=tid[0]) id->Append("_");
// 			    }
// 			    id->Append("@");
// 			    id->Append(makeStaticData()->GetParticleName(xpid));
// 			    id->Append("#genbod#");
// 			    for (int p=1;p<=tid[0];p++) {
// 				id->Append(makeStaticData()->GetParticleName(tid[p]));
// 				if (p!=tid[0]) id->Append("#");
// 			    }
// 			    id->Append("&genbod");
// 			    PGenBod *pmodel = new PGenBod((char*)id->Data(),"Pluto build-in genbod [clone]",-2);
// 			    pmodel->Add("dummy","parent");
// 			    pmodel->SetGroupID("genbod_models");
// 			    arr->Add((TObject *)pmodel);
// 			}
// 		    }

		} 
	    }

	} //END decay mode iterator

	model=NULL;



	if (!makeDynamicData()->GetParticleModel(pid)) { //check first if we done

	    /* Breit-Wigner model */
	    if (makeStaticData()->IsParticleHadron(pid) &&
		makeStaticData()->GetParticleTotalWidth(pid) > unstable_width) {
		id = new TString(makeStaticData()->GetParticleName(pid));
		id->Append("_bw");
		de = new TString(makeStaticData()->GetParticleName(pid));
		// de->Append(" <PBreitWigner>");
		PBreitWigner *pmodel = new PBreitWigner((char*)id->Data(),(char*)de->Data(),particlekey);
		model=(PChannelModel*) pmodel;
	    }

	    if (model) {
		model->Add(makeStaticData()->GetParticleName(pid),"daughter");
		model->SetGroupID("particle_models");
		arr->Add((TObject *)model);
	    }
	    
	}

    } //END particle iterator
    makeStaticData()->SetFreezeOut();

    return arr;

};

ClassImp(PStdModels)



