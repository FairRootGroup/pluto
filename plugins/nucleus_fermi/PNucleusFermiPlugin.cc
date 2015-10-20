/////////////////////////////////////////////////////////////////////
//
// This plugin adds various PIDs for nuclear targets, starting with
// ID=601 (Li6), in general (Z-3)*4+2 for the most abundant isotope (*).
//
// The following PIDs have been added so far:
//
// 6Li  Z=3,  ID=601
// 7Li  Z=3,  ID=602 (*)
// 10B  Z=5,  ID=609
// 11B  Z=5,  ID=610 (*)
// 11C  Z=6,  ID=613
// 12C  Z=6,  ID=614 (*)
// 39K  Z=19, ID=666 (*)
// 40K  Z=19, ID=667
// 41K  Z=19, ID=668
// 39Ca Z=20, ID=669
// 40Ca Z=20, ID=670 (*)
// Y would be 746, we will move it
// 90Zr Z=40, ID=745
// 91Zr Z=40, ID=746
// 92Zr Z=40, ID=747
// 93Zr Z=40, ID=748
// 94Zr Z=40, ID=749
// 95Zr Z=40, ID=750 (*,exception)
// 96Zr Z=40, ID=751
// 92Nb Z=41, ID=753
// 93Nb Z=41, ID=754 (*)
// 
//                             Author:  I. Froehlich
//                             Written: 27.7.2008
//                             Revised: 
//////////////////////////////////////////////////////////////////////

#include "PNucleusFermiPlugin.h"

PNucleusFermiPlugin::PNucleusFermiPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {

    gamma_active=proton_active=kFALSE;

    RequiresPlugin("elementary");

}

Bool_t PNucleusFermiPlugin::Activate(void) {

    //To be consistent:
    makeStaticData()->AddAlias("He3","3He");
    makeStaticData()->AddAlias("alpha","4He");
    
    
    //Add our nuclei:
    
    if (!makeStaticData()->IsParticleValid("6He")) { 
	makeStaticData()->AddParticle(600, "6He",6.0188891);
    }

    //7Li  Z=3, ID=602
    if (!makeStaticData()->IsParticleValid("7Li")) { 
	makeStaticData()->AddParticle(602, "7Li",6.533833);
    }

    //6Li  Z=3, ID=601
    if (!makeStaticData()->IsParticleValid("6Li")) { 
	makeStaticData()->AddParticle(601, "6Li",5.601518);
    }

    // 10B  Z=5,  ID=609
    if (!makeStaticData()->IsParticleValid("10B")) { 
	makeStaticData()->AddParticle(609, "10B",10.0129370*0.931494);
    }
    // 11B  Z=5,  ID=610
    if (!makeStaticData()->IsParticleValid("11B")) { 
	makeStaticData()->AddParticle(610, "11B",11.0093054*0.931494);
    } //http://en.wikipedia.org/wiki/Isotopes_of_boron

    //12C  Z=6, ID=614
    if (!makeStaticData()->IsParticleValid("12C")) { 
	makeStaticData()->AddParticle(614, "12C",11.174862);
    }

    //11C  Z=6, ID=613
    if (!makeStaticData()->IsParticleValid("11C")) { 
	makeStaticData()->AddParticle(613, "11C",10.254018);
    }

    // 39K  Z=19, ID=666
    if (!makeStaticData()->IsParticleValid("39K")) { 
	makeStaticData()->AddParticle(666, "39K",38.96370668*0.931494);
    }
    // 40K  Z=19, ID=667
    if (!makeStaticData()->IsParticleValid("40K")) { 
	makeStaticData()->AddParticle(667, "40K",39.96399848*0.931494);
    }
    // 41K  Z=19, ID=668
    if (!makeStaticData()->IsParticleValid("41K")) { 
	makeStaticData()->AddParticle(668, "41K",40.96182576*0.931494);
    } //http://en.wikipedia.org/wiki/Isotopes_of_potassium

    //40Ca Z=20, ID=670
    if (!makeStaticData()->IsParticleValid("40Ca")) { 
	makeStaticData()->AddParticle(670, "40Ca",37.214694);
    }

    //39Ca Z=20, ID=669
    if (!makeStaticData()->IsParticleValid("39Ca")) { 
	makeStaticData()->AddParticle(669, "39Ca",36.290772);
    }

    // 90Zr Z=40, ID=745
    if (!makeStaticData()->IsParticleValid("90Zr")) { 
	makeStaticData()->AddParticle(745, "90Zr",89.9047044*0.931494);
    }
    // 91Zr Z=40, ID=746
    if (!makeStaticData()->IsParticleValid("91Zr")) { 
	makeStaticData()->AddParticle(746, "91Zr",90.9056458*0.931494);
    }
    // 92Zr Z=40, ID=747
    if (!makeStaticData()->IsParticleValid("92Zr")) { 
	makeStaticData()->AddParticle(747, "92Zr",91.9050408*0.931494);
    }
    // 93Zr Z=40, ID=748
    if (!makeStaticData()->IsParticleValid("93Zr")) { 
	makeStaticData()->AddParticle(748, "93Zr",92.9064760*0.931494);
    }
    // 94Zr Z=40, ID=749
    if (!makeStaticData()->IsParticleValid("94Zr")) { 
	makeStaticData()->AddParticle(749, "94Zr",93.9063152*0.931494);
    }
    // 95Zr Z=40, ID=750 (*,exception)
    if (!makeStaticData()->IsParticleValid("95Zr")) { 
	makeStaticData()->AddParticle(750, "95Zr",94.9080426*0.931494);
    }
    // 96Zr Z=40, ID=751
    if (!makeStaticData()->IsParticleValid("96Zr")) { 
	makeStaticData()->AddParticle(751, "96Zr",95.9082734*0.931494);
    } //http://en.wikipedia.org/wiki/Isotopes_of_zirconium

    //92Nb Z=41, ID=753
    if (!makeStaticData()->IsParticleValid("92Nb")) { 
	makeStaticData()->AddParticle(753, "92Nb",85.59005);
    }
    //93Nb Z=41, ID=754
    if (!makeStaticData()->IsParticleValid("93Nb")) { 
	makeStaticData()->AddParticle(754, "93Nb",86.520784);
    }


    //BUGBUG: Have to set Baryon number etc...


    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
    pdmutil->AddSubGroup("fermi", "Fermi momenta", "root");
    pdmutil->SetGroup("fermi");


    PFermiDistributions *model = 
	new PFermiDistributions("3He_fermi@3He/fermi","Fermi momentum of 3He",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);				 
    model = new PFermiDistributions("4He_fermi@4He/fermi","Fermi momentum of 4He",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);				 
    model = new PFermiDistributions("7Li_fermi@7Li/fermi","Fermi momentum of 7Li",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);
    model = new PFermiDistributions("12C_fermi@12C/fermi","Fermi momentum of 12C",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);
    model = new PFermiDistributions("40Ca_fermi@40Ca/fermi","Fermi momentum of 40Ca",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);
    model = new PFermiDistributions("93Nb_fermi@93Nb/fermi","Fermi momentum of 93Nb",-1);
    model->SetRange(0,2.);
    pdmutil->Add(model);

    return kTRUE;

}

PNucleusFermiPlugin::~PNucleusFermiPlugin() {
}


Bool_t PNucleusFermiPlugin::ExecCommand(const char * command, Double_t value) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    if (strcmp (command,"gamma") == 0) {

	if (gamma_active) return kTRUE;

	pdmutil->AddSubGroup("gamma_beam", "Physics with gamma-beam", "root");
	pdmutil->SetGroup("gamma_beam");

	//Li6 (n spectator)
	//*****************
	makeStaticData()->AddParticle(602001,"g + 7Li",6.533833);
	makeStaticData()->AddAlias("g + 7Li","g+7Li");
	makeStaticData()->AddDecay(-1, "g + 7Li -> (g + n) + 6Li (quasi-free)",
				   "g + 7Li","g + n,6Li", 1.0 );
	PFermiMomentumGA *pmodel = 
	    new PFermiMomentumGA("gn_in_7Li@g + 7Li_to_g + n_6Li",
				 "Quasi-free particle production <nucleus_fermi>",-1);	
	// Now add all particles
	// Define spectators and final decay products (the granddaughters)
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("7Li,grandparent,target");
	pmodel->Add("6Li,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");	
	pdmutil->Add(pmodel);

	//(p spectator)
	makeStaticData()->AddDecay(-1, "g + 7Li -> (g + p) + 6He (quasi-free)",
				   "g + 7Li","g + p,6He", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gp_in_7Li@g + 7Li_to_g + p_6He",
				 "Quasi-free particle production <nucleus_fermi>",-1);	
	// Now add all particles
	// Define spectators and final decay products (the granddaughters)
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("7Li,grandparent,target");
	pmodel->Add("6He,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");	
	pdmutil->Add(pmodel);


	//12C (n spectator)
	//*****************
	makeStaticData()->AddParticle(614001,"g + 12C",11.174862); 
	makeStaticData()->AddAlias("g + 12C","g+12C");
	makeStaticData()->AddDecay(-1, "g + 12C -> (g + n) + 11C (quasi-free)",
				   "g + 12C","g + n,11C", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gn_in_12C@g + 12C_to_g + n_11C",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("12C,grandparent,target");
	pmodel->Add("11C,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//(p spectator)
	makeStaticData()->AddDecay(-1, "g + 12C -> (g + p) + 11B (quasi-free)",
				   "g + 12C","g + p,11B", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gp_in_12C@g + 12C_to_g + p_11B",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("12C,grandparent,target");
	pmodel->Add("11B,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//40Ca (n spectator)
	//*****************
	makeStaticData()->AddParticle(670001,"g + 40Ca",37.214694); 
	makeStaticData()->AddAlias("g + 40Ca","g+40Ca");
	makeStaticData()->AddDecay(-1, "g + 40Ca -> (g + n) + 39Ca (quasi-free)",
				   "g + 40Ca","g + n,39Ca", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gn_in_40Ca@g + 40Ca_to_g + n_39Ca",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("40Ca,grandparent,target");
	pmodel->Add("39Ca,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//(p spectator)
	makeStaticData()->AddDecay(-1, "g + 40Ca -> (g + p) + 39K (quasi-free)",
				   "g + 40Ca","g + p,39K", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gp_in_40Ca@g + 40Ca_to_g + p_39K",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("40Ca,grandparent,target");
	pmodel->Add("39K,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//93Nb (n spectator)
	//*****************
	makeStaticData()->AddParticle(754001,"g + 93Nb",86.520784); 
	makeStaticData()->AddAlias("g + 93Nb","g+93Nb");
	makeStaticData()->AddDecay(-1, "g + 93Nb -> (g + n) + 92Nb (quasi-free)",
				   "g + 93Nb","g + n,92Nb", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gn_in_93Nb@g + 93Nb_to_g + n_92Nb",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("93Nb,grandparent,target");
	pmodel->Add("92Nb,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//(p spectator)
	makeStaticData()->AddDecay(-1, "g + 93Nb -> (g + p) + 92Zr (quasi-free)",
				   "g + 93Nb","g + p,92Zr", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gp_in_93Nb@g + 93Nb_to_g + p_92Zr",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("93Nb,grandparent,target");
	pmodel->Add("92Zr,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//3He
	makeStaticData()->AddParticle(49001,"g + 3He",2.808391); 
	makeStaticData()->AddAlias("g + 3He","g+3He");
	makeStaticData()->AddDecay(-1, "g + 3He -> (g + p) + d (quasi-free)",
				   "g + 3He","g + p,d", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gp_in_3He@g + 3He_to_g + p_d",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("3He,grandparent,target");
	pmodel->Add("d,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	//4He
	makeStaticData()->AddParticle(47001,"g + 4He",3.727379); 
	makeStaticData()->AddAlias("g + 4He","g+4He");
	makeStaticData()->AddDecay(-1, "g + 4He -> (g + n) + 3He (quasi-free)",
				   "g + 4He","g + n,3He", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("gn_in_4He@g + 4He_to_g + n_3He",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("g,grandparent,beam");                            
	pmodel->Add("4He,grandparent,target");
	pmodel->Add("3He,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("g,granddaughter,p2");
	pdmutil->Add(pmodel);

	

	gamma_active = kTRUE;
	return kTRUE;
    } else if (strcmp (command,"proton") == 0) {

	if (proton_active) return kTRUE;

	pdmutil->AddSubGroup("proton_beam", "Physics with proton-beam", "root");
	pdmutil->SetGroup("proton_beam");

	Double_t p_mass = makeStaticData()->GetParticleMass("p");

	//Li6
	makeStaticData()->AddParticle(602014,"p + 7Li",6.533833 + p_mass);
	makeStaticData()->AddAlias("p + 7Li","p+7Li");
	makeStaticData()->AddDecay(-1, "p + 7Li -> (p + n) + 6Li (quasi-free)",
				   "p + 7Li","p + n,6Li", 1.0 );
	PFermiMomentumGA *pmodel = 
	    new PFermiMomentumGA("pn_in_7Li@p + 7Li_to_p + n_6Li",
				 "Quasi-free particle production <nucleus_fermi>",-1);	
	// Now add all particles
	// Define spectators and final decay products (the granddaughters)
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("7Li,grandparent,target");
	pmodel->Add("6Li,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");	
	pdmutil->Add(pmodel);

	//Li6 (p spectator): 
	makeStaticData()->AddDecay(-1, "p + 7Li -> (p + p) + 6He (quasi-free)",
				   "p + 7Li","p + p,6He", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pp_in_7Li@p + 7Li_to_p + p_6He",
				 "Quasi-free particle production <nucleus_fermi>",-1);	
	// Now add all particles
	// Define spectators and final decay products (the granddaughters)
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("7Li,grandparent,target");
	pmodel->Add("6He,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");	
	pdmutil->Add(pmodel);

	//12C
	makeStaticData()->AddParticle(614014,"p + 12C",11.174862 + p_mass); 
	makeStaticData()->AddAlias("p + 12C","p+12C");
	makeStaticData()->AddDecay(-1, "p + 12C -> (p + n) + 11C (quasi-free)",
				   "p + 12C","p + n,11C", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pn_in_12C@p + 12C_to_p + n_11C",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("12C,grandparent,target");
	pmodel->Add("11C,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);
	//p spectator:
	makeStaticData()->AddDecay(-1, "p + 12C -> (p + p) + 11B (quasi-free)",
				   "p + 12C","p + p,11B", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pp_in_12C@p + 12C_to_p + p_11B",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("12C,grandparent,target");
	pmodel->Add("11B,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);

	//40Ca
	makeStaticData()->AddParticle(670014,"p + 40Ca",37.214694 + p_mass); 
	makeStaticData()->AddAlias("p + 40Ca","p+40Ca");
	makeStaticData()->AddDecay(-1, "p + 40Ca -> (p + n) + 39Ca (quasi-free)",
				   "p + 40Ca","p + n,39Ca", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pn_in_40Ca@p + 40Ca_to_p + n_39Ca",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("40Ca,grandparent,target");
	pmodel->Add("39Ca,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);

	//p spectator:
	makeStaticData()->AddDecay(-1, "p + 40Ca -> (p + p) + 39K (quasi-free)",
				   "p + 40Ca","p + p,39K", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pp_in_40Ca@p + 40Ca_to_p + p_39K",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("40Ca,grandparent,target");
	pmodel->Add("39K,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);

	//93Nb
	makeStaticData()->AddParticle(754014,"p + 93Nb",86.520784 + p_mass); 
	makeStaticData()->AddAlias("p + 93Nb","p+93Nb");
	makeStaticData()->AddDecay(-1, "p + 93Nb -> (p + n) + 92Nb (quasi-free)",
				   "p + 93Nb","p + n,92Nb", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pn_in_93Nb@p + 93Nb_to_p + n_92Nb",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("93Nb,grandparent,target");
	pmodel->Add("92Nb,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);
	//p spectator:
	makeStaticData()->AddDecay(-1, "p + 93Nb -> (p + p) + 92Zr (quasi-free)",
				   "p + 93Nb","p + p,92Zr", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pp_in_93Nb@p + 93Nb_to_p + p_92Zr",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("93Nb,grandparent,target");
	pmodel->Add("92Zr,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);


	//3He
	makeStaticData()->AddParticle(49014,"p + 3He",2.808391 + p_mass); 
	makeStaticData()->AddAlias("p + 3He","p+3He");
	makeStaticData()->AddDecay(-1, "p + 3He -> (p + p) + d (quasi-free)",
				   "p + 3He","p + p,d", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pp_in_3He@p + 3He_to_p + p_d",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("3He,grandparent,target");
	pmodel->Add("d,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("p,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);

	//4He
	makeStaticData()->AddParticle(47014,"p + 4He",3.727379 + p_mass); 
	makeStaticData()->AddAlias("p + 4He","p+4He");
	makeStaticData()->AddDecay(-1, "p + 4He -> (p + n) + 3He (quasi-free)",
				   "p + 4He","p + n,3He", 1.0 );
	pmodel = 
	    new PFermiMomentumGA("pn_in_4He@p + 4He_to_p + n_3He",
				 "Quasi-free particle production <nucleus_fermi>",-1);
	pmodel->Add("q,parent");                                       
	pmodel->Add("p,grandparent,beam");                            
	pmodel->Add("4He,grandparent,target");
	pmodel->Add("3He,daughter,spectator");
	pmodel->Add("q,daughter,composite");
	pmodel->Add("n,granddaughter,participant"); 
	pmodel->Add("p,granddaughter,p2");
	pdmutil->Add(pmodel);


	proton_active = kTRUE;
	return kTRUE;
    } 





    return kFALSE;
}



ClassImp(PNucleusFermiPlugin)



