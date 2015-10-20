/////////////////////////////////////////////////////////////////////
//Plugin for the Bremsstrahlung from Kaptari/Kaempfer [L1] and
//Shyam/Mosel [L2]
//
//It can be activated with 
//makeDistributionManager()->Exec("brems: kaptari");
//  
//Commands can be formwarded to the plugin via:
//makeDistributionManager()->Exec("brems: command");
//
//The plugin supports the following commands:
//elastic         Elastic term only
//delta           Delta term only
//sum             Coherent sum
//fsi             Add FSI a la Jost function
//
//For the absolute normalization it is suggested to use the weighting option:
//
//kin_max = arg   arg=0...1 maximal kinematic range for di-lepton sampling
//weighting       activate weighting
//
//Example:
//makeDistributionManager()->Exec("brems: kaptari ; kin_max=0.9 ; weighting");
// 
//Author:  I. Froehlich
//Written: 27.7.2008
//
//References:
//[L1] L.P. Kaptari and B. Kaempfer, Nucl. Phys.  A 764 (2006) 338, arXiv:nucl-th/0504072
//[L2] R. Shyam, U. Mosel (Giessen U.) , Phys.Rev.C67:065202,2003, hep-ph/0303035, and priv. comm
//////////////////////////////////////////////////////////////////////

#include "PBremsstrahlungPlugin.h"




PBremsstrahlungPlugin::PBremsstrahlungPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {

    pn_gen = NULL;
    pp_gen = NULL;
    kin_max = 1.;
    pn_fsi = pp_fsi = np_fsi = NULL;
    pn_kk = pp_kk = pn_sm_total = pp_sm_total = NULL;
    pn_sm_elastic = pp_sm_elastic = pn_sm_delta = pp_sm_delta
	= pn_sm_n1520 = pp_sm_n1520 = NULL;
    pn = pn_kk;
    pp = pp_kk;
    array_pp_sm_125 = array_pn_sm_125 = NULL;

    RequiresPlugin("elementary");
}

Bool_t PBremsstrahlungPlugin::Activate(void) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    //First, we check if all composite particles are available

//     if (!makeStaticData()->IsParticleValid("p + p")) {
// 	makeStaticData()->AddParticle(14014, "p + p",2.);
//     }
//     if (!makeStaticData()->IsParticleValid("p + n")) {
// 	makeStaticData()->AddParticle(13014, "p + n",2.);
//     }
//     if (!makeStaticData()->IsParticleValid("n + p")) {
// 	makeStaticData()->AddParticle(14013, "n + p",2.);
//     }
//     if (!makeStaticData()->IsParticleValid("d + p")) {
// 	makeStaticData()->AddParticle(14045, "d + p",2.);
//     }
//     if (!makeStaticData()->IsParticleValid("p + d")) {
// 	makeStaticData()->AddParticle(45014, "p + d",2.);
//     }
    
    Int_t ipid[4];
    ipid[0]=makeStaticData()->GetParticleID("p + n");
    ipid[1]=makeStaticData()->GetParticleID("p");
    ipid[2]=makeStaticData()->GetParticleID("n");
    ipid[3]=makeStaticData()->GetParticleID("dilepton");

    //Add decay
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "p + n -> p + n + dilepton (Bremsstrahlung)", 
				   "p + n","p,n,dilepton", 1.0 );

    ipid[0]=makeStaticData()->GetParticleID("n + p");
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "n + p -> p + n + dilepton (Bremsstrahlung)", 
				   "n + p","p,n,dilepton", 1.0 );

    ipid[0]=makeStaticData()->GetParticleID("p + p");
    ipid[2]=makeStaticData()->GetParticleID("p");
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "p + p -> p + p + dilepton (Bremsstrahlung)", 
				   "p + p","p,p,dilepton", 1.0 );
    //deuteron beam
    ipid[0]=makeStaticData()->GetParticleID("d + p");
    ipid[2]=makeStaticData()->GetParticleID("n + p");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "d + p -> p + (n + p) (Quasi-Free)", 
				   "d + p","p,n + p", 1.0 );
    
    pdmutil->AddSubGroup("bremsstrahlung", "Bremsstrahlung", "root");

    pdmutil->SetGroup("bremsstrahlung");




    


    return kTRUE;
}

PBremsstrahlungPlugin::~PBremsstrahlungPlugin() {
}


Bool_t PBremsstrahlungPlugin::ExecCommand(const char * command, Double_t value) {

    if (strcmp (command,"weighting") == 0) {


	//Monto-Carlo integration
	pn->EnableWeighting();
	pn->SetExpectedWeightMean(-1);
	pn->SetVersionFlag(VERSION_GENERATOR_MC);
	
	pp->EnableWeighting();
	pp->SetExpectedWeightMean(-1);
	pp->SetVersionFlag(VERSION_GENERATOR_MC);


	if (pn_gen) return kTRUE;
	
	//TF1 object representing the di-lepton statistics:
	TF1 *flat = new TF1("flat","1",0,kin_max);
	
	//The "PInclusiveModel" can be used as a generator:
	pn_gen = 
	    new PInclusiveModel("flat@p + n_brems_p_n_dilepton/generator","Dilepton generator",-1);
	
	//The distribution template:
	pn_gen->Add("q,parent");
	pn_gen->Add("p,daughter");
	pn_gen->Add("n,daughter");
	pn_gen->Add("dilepton,daughter,primary");
	pn_gen->SetSampleFunction(flat);
	//Enable distribution as a generator
	pn_gen->EnableGenerator();    
	makeDistributionManagerUtil()->Add(pn_gen);

	//repeat this for the pp case:
	pp_gen = 
	    new PInclusiveModel("flat2@p + p_brems_p_p_dilepton/generator","Dilepton generator",-1);
    
	//The distribution template:
	pp_gen->Add("q,parent");
	pp_gen->Add("p,daughter");
	pp_gen->Add("p,daughter");
	pp_gen->Add("dilepton,daughter,primary");
	pp_gen->SetSampleFunction(flat);
	//Enable distribution as a generator
	pp_gen->EnableGenerator();    
	makeDistributionManagerUtil()->Add(pp_gen);



	return kTRUE;
    } else  if (strcmp (command,"kin_max") == 0) {
	kin_max = value;

	return kTRUE;
    } else  if (strcmp (command,"elastic") == 0) {
	if (current_author) { //sm
	    pp=pp_sm_elastic;
	    pn=pn_sm_elastic;
	    return makeDistributionManagerUtil()->Enable("bremsstrahlung_sm_elastic");
	} else {
	    pn->SetMode(0);
	    pp->SetMode(0);
	}
	return kTRUE;
    } else  if (strcmp (command,"delta") == 0) {
	if (current_author) { //sm
	    pp=pp_sm_delta;
	    pn=pn_sm_delta;
	    return makeDistributionManagerUtil()->Enable("bremsstrahlung_sm_delta");
	} else {
	    pn->SetMode(1);
	    pp->SetMode(1);
	}
	return kTRUE;
    } else  if (strcmp (command,"n1520") == 0) {
	if (current_author) { //sm
	    pp=pp_sm_n1520;
	    pn=pn_sm_n1520;
	    return makeDistributionManagerUtil()->Enable("bremsstrahlung_sm_n1520");
	} else {
	    Warning("ExecCommand","N1520 option not available for KK");
	}
	return kTRUE;
    } else  if (strcmp (command,"sum") == 0) {
	if (current_author) { //sm
	    pp=pp_sm_total;
	    pn=pn_sm_total;
	    return makeDistributionManagerUtil()->Enable("bremsstrahlung_sm_total");
	} else {
	    pn->SetMode(2);
	    pp->SetMode(2);
	}
	return kTRUE;
    } else  if (strcmp (command,"fsi") == 0) {

	if (pn_fsi) return kTRUE;
	makeDistributionManagerUtil()->AddSubGroup("fsi", "Final state effects", "root");
	makeDistributionManagerUtil()->SetGroup("fsi");

	pn_fsi = new PNNFSI("pn_fsi@p + n_brems_p_n_dilepton/nnfsi","PN FSI effects",-1);
	pn_fsi->EnableWeighting();
	pn_fsi->SetExpectedWeightMean(-1);
	makeDistributionManagerUtil()->Add(pn_fsi);

	np_fsi = new PNNFSI("np_fsi@n + p_brems_p_n_dilepton/nnfsi","PN FSI effects",-1);
	np_fsi->EnableWeighting();
	np_fsi->SetExpectedWeightMean(-1);
	makeDistributionManagerUtil()->Add(np_fsi);

	pp_fsi = new PNNFSI("pp_fsi@p + p_brems_p_p_dilepton/nnfsi","PP FSI effects",-1);
	pp_fsi->EnableWeighting();
	pp_fsi->SetExpectedWeightMean(-1);
	pp_fsi->SetA0R0(-7.8098,2.767);

	makeDistributionManagerUtil()->Add(pp_fsi);

	return kTRUE;
    } else  if (strcmp (command,"kaptari") == 0) {
	current_author=0;
	
	if (pp_kk) {
	    pp=pp_kk;
	    pn=pn_kk;
	    makeDistributionManagerUtil()->Enable("bremsstrahlung_kk");
	    return kTRUE;
	}

	PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
	pdmutil->AddSubGroup("bremsstrahlung_kk", "Bremsstrahlung from KK", "bremsstrahlung");
	pdmutil->SetGroup("bremsstrahlung_kk");
	
	pn = new PBremsstrahlung("p + n_brems_p_n_dilepton", "pn Kaptari/Kaempfer bremsstrahlung",-1);
	
	pn->Add("q,parent");
	pn->Add("p,grandparent,p1");
	pn->Add("n,grandparent,p2");
	pn->Add("p,daughter,p3");
	pn->Add("n,daughter,p4");
	pn->Add("dilepton", "daughter");
	
	PBremsstrahlung *newmodel3 = new PBremsstrahlung("n + p_brems_p_n_dilepton", "pn Kaptari/Kaempfer bremsstrahlung [clone]",-1);
	newmodel3->Add("dummy,parent");
	
	
	pp = new PBremsstrahlung("p + p_brems_p_p_dilepton", "pp Kaptari/Kaempfer bremsstrahlung",-1);
	
	pp->Add("q,parent");
	pp->Add("p,grandparent,p1");
	pp->Add("p,grandparent,p2");
	pp->Add("p,daughter,p3");
	pp->Add("p,daughter,p4");
	pp->Add("dilepton", "daughter");
	
	//make it known to the Pluto world:
	pdmutil->Add(pp);
	pdmutil->Add(pn);
	
	pdmutil->Add(newmodel3);
	pp_kk=pp;
	pn_kk=pn;
	return kTRUE;

    } else  if (strcmp (command,"shyam") == 0) {

	char *path = getenv("PLUTOSYS");
	if (!path) {
	    Error ("ExecCommand","$PLUTOSYS not set: dat files cannot be located");
	    return kFALSE;
	}

	current_author=1;

	if (pp_sm_total) {
	    pp=pp_sm_total;
	    pn=pn_sm_total;
	    makeDistributionManagerUtil()->Enable("bremsstrahlung_sm_total");
	    return kTRUE;
	}

	//******* Let us read the array from the files:
	array_pp_sm_125 = new PArray(1); //Dim = 1
	array_pp_sm_125->Scaling(0.001);  //convert to barn
	array_pp_sm_125->OpenFile(path + TString("/plugins/brems/shyam_pp_1.25GeV.dat"), 0, 5);
	

	array_pn_sm_125 = new PArray(2); //Dim = 2, 3 files
	array_pn_sm_125->Scaling(0.001);  //convert to barn
	array_pn_sm_125->OpenFile(path + TString("/plugins/brems/shyam_pn_1GeV.dat"), 0, 5, 1);
	array_pn_sm_125->OpenFile(path + TString("/plugins/brems/shyam_pn_1.25GeV.dat"), 0, 5, 1.25);
	array_pn_sm_125->OpenFile(path + TString("/plugins/brems/shyam_pn_1.5GeV.dat"), 0, 5, 1.5);

	//******* Modify groups
	PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
	pdmutil->AddSubGroup("bremsstrahlung_sm", "Bremsstrahlung from Shyam/Mosel", "bremsstrahlung");
	pdmutil->AddSubGroup("bremsstrahlung_sm_total", 
			     "Bremsstrahlung from Shyam/Mosel (total)", "bremsstrahlung_sm");
	pdmutil->AddSubGroup("bremsstrahlung_sm_elastic", 
			     "Bremsstrahlung from Shyam/Mosel (elastic)", "bremsstrahlung_sm");
	pdmutil->AddSubGroup("bremsstrahlung_sm_delta", 
			     "Bremsstrahlung from Shyam/Mosel (delta)", "bremsstrahlung_sm");
	pdmutil->AddSubGroup("bremsstrahlung_sm_n1520", 
			     "Bremsstrahlung from Shyam/Mosel (N*1520)", "bremsstrahlung_sm");
	
	//******* pp case TOTAL
	pdmutil->SetGroup("bremsstrahlung_sm_total");
	pp_sm_total = new PBremsstrahlung("pp_sm_total@p + p_brems_p_p_dilepton", 
					  "pp Shyam/Mosel bremsstrahlung (total)",-1);
	
	pp_sm_total->Add("q,parent");
	pp_sm_total->Add("p,grandparent,p1");
	pp_sm_total->Add("p,grandparent,p2");
	pp_sm_total->Add("p,daughter,p3");
	pp_sm_total->Add("p,daughter,p4");
	pp_sm_total->Add("dilepton", "daughter");
	pp_sm_total->SetFunc(array_pp_sm_125->GetTGraph(0,1));  //1=total
	pdmutil->Add(pp_sm_total);

	//******* pn case TOTAL
	pn_sm_total = new PBremsstrahlung("pn_sm_total@p + n_brems_p_n_dilepton", 
				    "pn Shyam/Mosel bremsstrahlung (total)",-1);
	
	pn_sm_total->Add("q,parent");
	pn_sm_total->Add("p,grandparent,p1");
	pn_sm_total->Add("n,grandparent,p2");
	pn_sm_total->Add("p,daughter,p3");
	pn_sm_total->Add("n,daughter,p4");
	pn_sm_total->Add("dilepton", "daughter");
	TGraph2D *gr = array_pn_sm_125->GetTGraph2D(0,1);   //1=total
	pn_sm_total->SetFunc(gr);
	pdmutil->Add(pn_sm_total);	
	
	PBremsstrahlung *newmodel3 = new PBremsstrahlung("np_sm_total@n + p_brems_p_n_dilepton", 
							 "pn Shyam/Mosel bremsstrahlung (total) [clone]",-1);
	newmodel3->Add("dummy,parent");
	newmodel3->SetFunc(gr);
	pdmutil->Add(newmodel3);

	//******* pp case ELASTIC
	pdmutil->SetGroup("bremsstrahlung_sm_elastic");
	pp_sm_elastic = new PBremsstrahlung("pp_sm_elastic@p + p_brems_p_p_dilepton", 
					    "pp Shyam/Mosel bremsstrahlung (elastic)",-1);
	
	pp_sm_elastic->Add("p,grandparent,p1");
	pp_sm_elastic->Add("p,grandparent,p2");
	pp_sm_elastic->Add("q,parent");
	pp_sm_elastic->Add("p,daughter,p3");
	pp_sm_elastic->Add("p,daughter,p4");
	pp_sm_elastic->Add("dilepton", "daughter");
	pp_sm_elastic->SetFunc(array_pp_sm_125->GetTGraph(0,2));  //2=elastic
	pdmutil->Add(pp_sm_elastic);

	//******* pn case ELASTIC
	pn_sm_elastic = new PBremsstrahlung("pn_sm_elastic@p + n_brems_p_n_dilepton", 
					    "pn Shyam/Mosel bremsstrahlung (elastic)",-1);
	
	pn_sm_elastic->Add("p,grandparent,p1");
	pn_sm_elastic->Add("n,grandparent,p2");
	pn_sm_elastic->Add("q,parent");
	pn_sm_elastic->Add("p,daughter,p3");
	pn_sm_elastic->Add("n,daughter,p4");
	pn_sm_elastic->Add("dilepton", "daughter");
	gr = array_pn_sm_125->GetTGraph2D(0,2);   //2=elastic
	pn_sm_elastic->SetFunc(gr);
	pdmutil->Add(pn_sm_elastic);	
	
	newmodel3 = new PBremsstrahlung("np_sm_elastic@n + p_brems_p_n_dilepton", 
					"pn Shyam/Mosel bremsstrahlung (elastic) [clone]",-1);
	newmodel3->Add("dummy,parent");
	newmodel3->SetFunc(gr);
	pdmutil->Add(newmodel3);

	//******* pp case DELTA
	pdmutil->SetGroup("bremsstrahlung_sm_delta");
	pp_sm_delta = new PBremsstrahlung("pp_sm_delta@p + p_brems_p_p_dilepton", 
					    "pp Shyam/Mosel bremsstrahlung (delta)",-1);
	
	pp_sm_delta->Add("p,grandparent,p1");
	pp_sm_delta->Add("p,grandparent,p2");
	pp_sm_delta->Add("q,parent");
	pp_sm_delta->Add("p,daughter,p3");
	pp_sm_delta->Add("p,daughter,p4");
	pp_sm_delta->Add("dilepton", "daughter");
	pp_sm_delta->SetFunc(array_pp_sm_125->GetTGraph(0,3));  //3=delta
	pdmutil->Add(pp_sm_delta);

	//******* pn case DELTA
	pn_sm_delta = new PBremsstrahlung("pn_sm_delta@p + n_brems_p_n_dilepton", 
					    "pn Shyam/Mosel bremsstrahlung (delta)",-1);
	
	pn_sm_delta->Add("p,grandparent,p1");
	pn_sm_delta->Add("n,grandparent,p2");
	pn_sm_delta->Add("q,parent");
	pn_sm_delta->Add("p,daughter,p3");
	pn_sm_delta->Add("n,daughter,p4");
	pn_sm_delta->Add("dilepton", "daughter");
	gr = array_pn_sm_125->GetTGraph2D(0,3);   //3=delta
	pn_sm_delta->SetFunc(gr);
	pdmutil->Add(pn_sm_delta);	
	
	newmodel3 = new PBremsstrahlung("np_sm_delta@n + p_brems_p_n_dilepton", 
					"pn Shyam/Mosel bremsstrahlung (delta) [clone]",-1);
	newmodel3->Add("dummy,parent");
	newmodel3->SetFunc(gr);
	pdmutil->Add(newmodel3);

	//******* pp case N1520
	pdmutil->SetGroup("bremsstrahlung_sm_n1520");
	pp_sm_n1520 = new PBremsstrahlung("pp_sm_n1520@p + p_brems_p_p_dilepton", 
					    "pp Shyam/Mosel bremsstrahlung (n1520)",-1);
	
	pp_sm_n1520->Add("p,grandparent,p1");
	pp_sm_n1520->Add("p,grandparent,p2");
	pp_sm_n1520->Add("q,parent");
	pp_sm_n1520->Add("p,daughter,p3");
	pp_sm_n1520->Add("p,daughter,p4");
	pp_sm_n1520->Add("dilepton", "daughter");
	pp_sm_n1520->SetFunc(array_pp_sm_125->GetTGraph(0,4));  //4=n1520
	pdmutil->Add(pp_sm_n1520);

	//******* pn case N1520
	pn_sm_n1520 = new PBremsstrahlung("pn_sm_n1520@p + n_brems_p_n_dilepton", 
					    "pn Shyam/Mosel bremsstrahlung (n1520)",-1);
	
	pn_sm_n1520->Add("p,grandparent,p1");
	pn_sm_n1520->Add("n,grandparent,p2");
	pn_sm_n1520->Add("q,parent");
	pn_sm_n1520->Add("p,daughter,p3");
	pn_sm_n1520->Add("n,daughter,p4");
	pn_sm_n1520->Add("dilepton", "daughter");
	gr = array_pn_sm_125->GetTGraph2D(0,4);   //4=n1520
	pn_sm_n1520->SetFunc(gr);
	pdmutil->Add(pn_sm_n1520);	
	
	newmodel3 = new PBremsstrahlung("np_sm_n1520@n + p_brems_p_n_dilepton", 
					"pn Shyam/Mosel bremsstrahlung (n1520) [clone]",-1);
	newmodel3->Add("dummy,parent");
	newmodel3->SetFunc(gr);
	pdmutil->Add(newmodel3);

	pp=pp_sm_total;
	pn=pn_sm_total;
	

	return kTRUE;

    } 

    return kFALSE;
}



ClassImp(PBremsstrahlungPlugin)



