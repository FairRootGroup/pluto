/////////////////////////////////////////////////////////////////////
//Plugin to modify the Dalitz decays
//
//This plugin is part of the distribution manager
//It can be activated with 
//makeDistributionManager()->Exec("dalitz_mod: command");
//
//where the plugin supports the following commands:
//static_br_thresh=value   Threshold for disabling static br in GeV
//flat_generator           Enables a flat generator for all registered 
//                         Dalitz decays. Instead of the static br we
//                         use the dG/dm directly, if the parent res
//                         width is > static_br_thresh
//krivoruchenko            Use Delta Dalitz from Krivoruchenko [L1]
//
//                         Author:  I. Froehlich
//                         Written: 17.9.2008
//
//References:
//[L1] M.I. Krivoruchenko, A. Faessler (Tubingen U.), Phys.Rev.D65:017502,2002, nucl-th/0104045
//
//////////////////////////////////////////////////////////////////////

#include "PDalitzModPlugin.h"

#include "PDeltaDalitzFF.h"


PDalitzModPlugin::PDalitzModPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {
    
    //    RequiresPlugin("elementary");
}

Bool_t PDalitzModPlugin::Activate(void) {

    //    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
    static_br_thresh = 1111.;  //something very large -> disabled
    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
    kriv1 = kriv2 = NULL;
    pdmutil->SetGroup("decay_models");
    return kTRUE;
}

PDalitzModPlugin::~PDalitzModPlugin() {
}


Bool_t PDalitzModPlugin::ExecCommand(const char * command, Double_t value) {



    if (strcmp (command,"flat_generator") == 0) {

	PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
	pdmutil->LinkDB();

	pdmutil->AddSubGroup("generators", "Generator models", "root");
	pdmutil->SetGroup("generators");

	TF1 *flat = new TF1("flat","1",0,1);

	Int_t header=makeDataBase()->GetEntry("std_set");
	Int_t particlekey=-1;
	//      Int_t generator_key=makeStaticData()->MakeDirectoryEntry("modeldef","generator");

	//loop over particles
	while (makeDataBase()->MakeListIterator(header, "snpart","slink" , &particlekey)) {
	
	    Int_t decaykey=-1;
	    Int_t pid = makeStaticData()->GetParticleIDByKey(particlekey);

	    //loop over decays
	    while ((makeStaticData()->GetParticleNChannelsByKey(particlekey)>0) && 
		   (makeDataBase()->MakeListIterator(particlekey, "pnmodes", "link", &decaykey))) {
		Int_t tid[11];
		tid[0]=10; 
		makeStaticData()->GetDecayModeByKey(decaykey,tid); // retrieve current mode info

		if (PData::IsDalitz(pid,tid[1],tid[2])) {
	    
		    TString * id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_generator_");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("_");
		    }
		    id->Append("@");
		    id->Append(makeStaticData()->GetParticleName(pid));
		    id->Append("#generator#");
		    for (int p=1;p<=tid[0];p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p!=tid[0]) id->Append("#");
		    }
		    id->Append("&generator");

		    PInclusiveModel *dilepton_generator = 
			new PInclusiveModel((char*)id->Data(), "Dilepton generator",-2);    
		    dilepton_generator->Set("dilepton,primary");
		    dilepton_generator->SetSampleFunction(flat);
		    dilepton_generator->EnableGenerator();
	    
		    pdmutil->Add(dilepton_generator);

		    PChannelModel *pmodel=makeDynamicData()->GetDecayModelByKey(decaykey);
		    if (pmodel) {
			pmodel->EnableWeighting();
			pmodel->SetVersionFlag(VERSION_GENERATOR_MC);
			//disable re-scaling if parent width is larger
			//then static_br_thresh
			if (makeStaticData()->GetParticleTotalWidth(pid) > static_br_thresh) {
			    pmodel->SetExpectedWeightMean(-1);
			    Info("ExecCommand","Model <%s> uses dGamma/dM for the branching ratio",
				 pmodel->GetIdentifier());
			} 
		    } else {
			Warning("ExecCommand","Primary model not found");
		    }
		} //isDalitz
	  
	  
	    }
	}
      
	return kTRUE;

    } else if (strcmp (command,"static_br_thresh") == 0) {

	static_br_thresh = value;
	return kTRUE;

    } else if (strcmp (command,"krivoruchenko") == 0) {

	if (kriv1) return kTRUE;
	PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
	pdmutil->LinkDB();
	pdmutil->SetGroup("decay_models");
	kriv1= new PDeltaDalitzKrivoruchenko("D+_krivoruchenko@D+_to_p_dilepton",
					     "dgdm from Krivoruchenko",-1);

	kriv2 = new PDeltaDalitzKrivoruchenko("D0_krivoruchenko@D0_to_n_dilepton",
					      "dgdm from Krivoruchenko",-1);
	pdmutil->Add(kriv1);
	pdmutil->Add(kriv2);

	//Add VMD-FF in the vmd group
	pdmutil->AddSubGroup("vmd", "VMD form factors", "root");
	pdmutil->SetGroup("vmd");

	PDeltaDalitzFF * vmd_newmodel = new PDeltaDalitzFF("D+_iachello_ff@D+_to_p_dilepton/formfactor",
							   "Iachello ff for D+ -> p e+e-",-1);
	pdmutil->Add(vmd_newmodel);
	vmd_newmodel = new PDeltaDalitzFF("D0_iachello_ff@D0_to_n_dilepton/formfactor",
					  "Iachello ff for D0 -> n e+e-",-1);
	pdmutil->Add(vmd_newmodel);
	

	//Add QED-FF in the qed group
	pdmutil->AddSubGroup("qed", "QED form factors", "root");
	pdmutil->SetGroup("qed");

	PDeltaDalitzFF * qed_newmodel = new PDeltaDalitzFF("D+_qed_ff@D+_to_p_dilepton/formfactor",
							   "QED ff for D+ -> p e+e-",-1);
	qed_newmodel->SetQED(1);
	qed_newmodel->SetCC(3.0,0,0);
	pdmutil->Add(qed_newmodel);
	qed_newmodel = new PDeltaDalitzFF("D0_qed_ff@D0_to_n_dilepton/formfactor",
					  "QED ff for D0 -> n e+e-",-1);
	qed_newmodel->SetQED(1);
	qed_newmodel->SetCC(3.0,0,0);
	pdmutil->Add(qed_newmodel);
	return kTRUE;

    } 
	

    return kFALSE;
}



ClassImp(PDalitzModPlugin)



    
