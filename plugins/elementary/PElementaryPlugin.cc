/////////////////////////////////////////////////////////////////////
//  A plugin to activate the NN composite particles
//  and decays needed for other plugins,
//  and total cross sections for specific reactions
//  
//  Total cross sections for NN->Delta N are added from Teis et al. 
//  (Ref 10)
//
//  For the eta cross sections, the following references have been used:
//
//  [L1] P. Moskal et al.,
//       "Near threshold production of the eta meson via the quasi-free 
//       pn --> pn eta reaction"
//       PRC 79 (2009) 015208, arXiv:0807.0722 [hep-ex]
//  [L2] H. Calen et al.,
//       "Measurement of the quasifree p + n $\to$ p + n + eta reaction near
//       threshold"
//       PRC 58 (1998) 2667.
//  [L3] H. Calen et al.,
//       "Measurement of the quasifree p + n --> d + eta reaction near threshold"
//       PRL 79} (1997) 2642.
//
//                             Author:  I. Froehlich
//                             Written: 21.8.2008
//                           
//////////////////////////////////////////////////////////////////////



#include "PDataBase.h"
#include "PElementaryPlugin.h"
#include "PTCrossWeight.h"


PElementaryPlugin::PElementaryPlugin() {
}


PElementaryPlugin::PElementaryPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {
}

PElementaryPlugin::~PElementaryPlugin() {
}





Bool_t PElementaryPlugin::Activate(void) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    //************** Composite particles ******************//

    //mass is important for the emin in PChannel
    if (!makeStaticData()->IsParticleValid("p + p")) {
	makeStaticData()->AddParticle(14014, "p + p",
				      2.*makeStaticData()->GetParticleMass("p"));
	makeStaticData()->AddAlias("p + p","p+p");
    }
    if (!makeStaticData()->IsParticleValid("p + n")) {
	makeStaticData()->AddParticle(13014, "p + n",
				      makeStaticData()->GetParticleMass("p")+
				      makeStaticData()->GetParticleMass("n"));
	makeStaticData()->AddAlias("p + n","p+n");
    }
    if (!makeStaticData()->IsParticleValid("n + p")) {
	makeStaticData()->AddParticle(14013, "n + p",
				      makeStaticData()->GetParticleMass("p")+
				      makeStaticData()->GetParticleMass("n"));
	makeStaticData()->AddAlias("n + p","n+p");
    }
    if (!makeStaticData()->IsParticleValid("d + p")) {
	makeStaticData()->AddParticle(14045, "d + p",
				      makeStaticData()->GetParticleMass("p")+
				      makeStaticData()->GetParticleMass("d"));
	makeStaticData()->AddAlias("d + p","d+p");
    }
    if (!makeStaticData()->IsParticleValid("p + d")) {
	makeStaticData()->AddParticle(45014, "p + d",
				      makeStaticData()->GetParticleMass("p")+
				      makeStaticData()->GetParticleMass("d"));
	makeStaticData()->AddAlias("p + d","p+d");
    }

    if(!makeStaticData()->IsParticleValid("g + n")) {
	makeStaticData()->AddParticle(13001, "g + n",
				      makeStaticData()->GetParticleMass("n"));
	makeStaticData()->AddAlias("g + n","g+n");
    }
    if(!makeStaticData()->IsParticleValid("g + p")) {
	makeStaticData()->AddParticle(14001, "g + p",
				      makeStaticData()->GetParticleMass("p"));
	makeStaticData()->AddAlias("g + p","g+p");
    }
    if(!makeStaticData()->IsParticleValid("n + g")) {
	makeStaticData()->AddParticle(1013, "n + g",
				      makeStaticData()->GetParticleMass("n"));
	makeStaticData()->AddAlias("n + g","n+g");
    }
    if(!makeStaticData()->IsParticleValid("p + g")) {
	makeStaticData()->AddParticle(1014, "p + g",
				      makeStaticData()->GetParticleMass("p"));
	makeStaticData()->AddAlias("p + g","p+g");
    }

    //************** Composite decays ******************//
    //*                                                *//
    //* Do not change the order here, otherwise
    //* the decays_idx will be re-sorted                
    //**************************************************//

    Int_t ipid[5];
    ipid[0]=makeStaticData()->GetParticleID("p + p");
    ipid[1]=makeStaticData()->GetParticleID("p");
    ipid[2]=makeStaticData()->GetParticleID("D+");

    //Delta(1232) production

    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "p + p -> p + D+", 
				   "p + p","p,D+", 1.0 );
    ipid[1]=makeStaticData()->GetParticleID("n");
    ipid[2]=makeStaticData()->GetParticleID("D++");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "p + p -> n + D++", 
				   "p + p","n,D++", 1.0 );
        
    ipid[0]=makeStaticData()->GetParticleID("p + n");
    ipid[1]=makeStaticData()->GetParticleID("n");
    ipid[2]=makeStaticData()->GetParticleID("D+");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "p + n -> n + D+", 
				   "p + n","n,D+", 1.0 );
    ipid[1]=makeStaticData()->GetParticleID("p");
    ipid[2]=makeStaticData()->GetParticleID("D0");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "p + n -> p + D0", 
				   "p + n","p,D0", 1.0 );
    ipid[0]=makeStaticData()->GetParticleID("n + p");
    ipid[1]=makeStaticData()->GetParticleID("n");
    ipid[2]=makeStaticData()->GetParticleID("D+");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "n + p -> n + D+", 
				   "n + p","n,D+", 1.0 );
    ipid[1]=makeStaticData()->GetParticleID("p");
    ipid[2]=makeStaticData()->GetParticleID("D0");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "n + p -> p + D0", 
				   "n + p","p,D0", 1.0 );
    //This causes the BR the be equally distributed. This is a caveat
    //since the exp_w_mean will be folded in the PChannel.
    //Therefore it acts as an additional factor

    //Eta physics:
    ipid[0]=makeStaticData()->GetParticleID("p + p");
    ipid[1]=makeStaticData()->GetParticleID("p");
    ipid[2]=makeStaticData()->GetParticleID("p");
    ipid[3]=makeStaticData()->GetParticleID("eta");
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "p + p -> p + p + eta", 
				   "p + p","p,p,eta", 1.0 );
    ipid[0]=makeStaticData()->GetParticleID("n + p");
    ipid[1]=makeStaticData()->GetParticleID("n");
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "n + p -> n + p + eta", 
				   "n + p","n,p,eta", 1.0 );
    ipid[0]=makeStaticData()->GetParticleID("p + n");
    if (makeStaticData()->GetDecayKey(ipid, 3) < 0)
	makeStaticData()->AddDecay(-1, "p + n -> n + p + eta", 
				   "p + n","n,p,eta", 1.0 );

    //pn -> d eta
    ipid[0]=makeStaticData()->GetParticleID("p + n");
    ipid[1]=makeStaticData()->GetParticleID("d");
    ipid[2]=makeStaticData()->GetParticleID("eta");

    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "p + n -> d + eta", 
				   "p + n","d,eta", 1.0 );

    ipid[0]=makeStaticData()->GetParticleID("n + p");
    if (makeStaticData()->GetDecayKey(ipid, 2) < 0)
	makeStaticData()->AddDecay(-1, "n + p -> d + eta", 
				   "n + p","d,eta", 1.0 );

// x=4.98953, y=0.156364
// x=4.72892, y=0.209786
// x=4.49976, y=0.254305
// x=4.22567, y=0.343342
// x=3.98303, y=0.441283
// x=3.73141, y=0.610455
// x=3.49326, y=0.815241
// x=3.23714, y=1.17139
// x=2.97204, y=1.71452
// x=2.74737, y=2.41791
// x=2.64852, y=2.81858
// x=2.53169, y=3.3439
// x=2.42834, y=3.73567
// x=2.3879, y=3.85142
// x=2.35645, y=3.88703
// x=2.28456, y=3.70896
// x=2.21716, y=2.78297
// x=2.17223, y=1.74123
// x=2.1228, y=0.601551


   Double_t x[20]={
	2.0,
	2.1228,
	2.17223,
	2.21716,
	2.28456,
	2.35645,
	2.3879,
	2.42834,
	2.53169,
	2.64852,
	2.74737,
	2.97204,
	3.23714,
	3.49326,
	3.73141,
	3.98303,
	4.22567,
	4.49976,
	4.72892,
	4.98953
    };

    Double_t y[20]={
	0.0,
	0.601551,
	1.74123,
	2.78297,
	3.70896,
	3.88703,
	3.85142,
	3.73567,
	3.3439,
	2.81858,
	2.41791,
	1.71452,
	1.17139,
	0.815241,
	0.610455,
	0.441283,
	0.343342,
	0.254305,
	0.209786,
	0.156364
};


    pdmutil->AddSubGroup("tcross", "Total cross sections", "root");
    pdmutil->SetGroup("tcross");

    TGraph *delta_gra = new TGraph(20,x,y);

    PTCrossWeight * delta_cross = 
	new PTCrossWeight("p + p_to_D+_p/tcross",
			  "D+ cross section in the p+p reaction",-1);
    delta_cross->Add("p,grandparent,beam");
    delta_cross->Add("p,grandparent,target");
    delta_cross->Add("q,parent");
    delta_cross->Add("p,daughter");
    delta_cross->Add("D+,daughter");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) / 1000.); //Isospin coefficient
    //factor 1000 to convert in [b]
    pdmutil->Add(delta_cross);

    delta_cross = 
	new PTCrossWeight("p + p_to_D++_n/tcross",
			  "D++ cross section in the p+p reaction",-1);
    delta_cross->Add("p,grandparent,beam");
    delta_cross->Add("p,grandparent,target");
    delta_cross->Add("q,parent");
    delta_cross->Add("n,daughter");
    delta_cross->Add("D++,daughter");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) * 3. / 1000.); //Isospin coefficient
    pdmutil->Add(delta_cross);

    delta_cross = 
	new PTCrossWeight("p + n_to_D0_p/tcross",
			  "D0 cross section in the p+n reaction",-1);
    delta_cross->Add("p,grandparent,beam");
    delta_cross->Add("n,grandparent,target");
    delta_cross->Add("q,parent");
    delta_cross->Add("p,daughter");
    delta_cross->Add("D0,daughter");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) / 1000.); //Isospin coefficient
    pdmutil->Add(delta_cross);
    
    delta_cross = 
	new PTCrossWeight("n + p_to_D0_p/tcross",
			  "D0 cross section in the p+n reaction [clone]",-1);
    delta_cross->Add("dummy,parent");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) / 1000.); //Isospin coefficient
    pdmutil->Add(delta_cross);


    delta_cross = 
	new PTCrossWeight("p + n_to_D+_n/tcross",
			  "D+ cross section in the p+n reaction",-1);
    delta_cross->Add("p,grandparent,beam");
    delta_cross->Add("n,grandparent,target");
    delta_cross->Add("q,parent");
    delta_cross->Add("n,daughter");
    delta_cross->Add("D+,daughter");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) / 1000.); //Isospin coefficient
    pdmutil->Add(delta_cross);
    
    delta_cross = 
	new PTCrossWeight("n + p_to_D+_n/tcross",
			  "D+ cross section in the p+n reaction [clone]",-1);
    delta_cross->Add("dummy,parent");
    delta_cross->SetCrossSection(delta_gra,kTRUE); 
    delta_cross->SetScaling((3./2.) / 1000.); //Isospin coefficient
    pdmutil->Add(delta_cross);

    //pp -> pp eta case
    //TODO: Parameterization>120MeV
   //  x=1.0574, y=0.138652
// x=2.32627, y=0.367414
// x=4.86402, y=0.840432
// x=9.09359, y=1.45904
// x=19.0331, y=2.48683
// x=28.5496, y=3.33745
// x=44.8335, y=7.77583
// x=56.0419, y=10.8264
// x=65.347, y=15.0739
// x=75.498, y=19.4994
// x=91.3589, y=29.2215
// x=103.202, y=43.7908
// x=117.371, y=62.1022

    Double_t x_pp_eta[13]={
	0,1.0574,2.32627,4.86402,9.09359,19.0331,
//	28.5496,
	44.8335,56.0419,
	65.347,75.498,91.3589,103.202,117.371
    };

    Double_t y_pp_eta[13]={
	0,0.138652,0.367414,0.840432,1.45904,2.48683,
//	3.33745,
	7.77583,10.8264,
	15.0739,19.4994,29.2215,43.7908,62.1022
    };

    Double_t thr=makeStaticData()->GetParticleMass("p")*2+
	makeStaticData()->GetParticleMass("eta");
    //Better is to convert to the SECONDARY_MODELS standard
    for (int i=0;i<13;i++) {
	x_pp_eta[i]*=0.001; //in GeV
	x_pp_eta[i]+=thr;
    }
	    
    TGraph *pp_eta_gra = new TGraph(13,x_pp_eta,y_pp_eta);

    PTCrossWeight * pp_eta_cross = 
	new PTCrossWeight("p + p_to_p_p_eta/tcross",
			  "Eta cross section in the p+p reaction (<120 MeV)",-1);
    pp_eta_cross->Add("p,grandparent,beam");
    pp_eta_cross->Add("p,grandparent,target");
    pp_eta_cross->Add("q,parent");
    pp_eta_cross->Add("p,daughter");
    pp_eta_cross->Add("p,daughter");
    pp_eta_cross->Add("eta,daughter");
    pp_eta_cross->SetCrossSection(pp_eta_gra,kFALSE); 
    pp_eta_cross->SetScaling(.000001); //convert from [ub] to [b]
//    pp_eta_cross->SetExcessEnergy(kTRUE,kTRUE);  //Q in MeV
    pdmutil->Add(pp_eta_cross);



    //Now the pn -> pn eta case:
   //  x=2.74923, y=1.02883
// x=7.40176, y=3.79588
// x=12.6887, y=9.51891
// x=25.8004, y=20.2298
// x=35.7399, y=35.1202
// x=45.468, y=51.6716
// x=55.4075, y=68.0823
// x=63.4437, y=91.3697
// x=73.8061, y=124.898
// x=81.6308, y=139.466
// x=90.09, y=164.565
// x=98.5492, y=220.854
// x=108.489, y=229.127

    Double_t x_pn_eta[14]={
	0.,2.74923,7.40176, 12.6887,25.8004, 35.7399, 45.468,55.4075,
	63.4437,73.8061, 81.6308,90.09,   98.5492,108.489};
    Double_t y_pn_eta[14]={ 
	0.,1.02883,3.79588, 9.51891,20.2298, 35.1202,51.6716, 68.0823,
	91.3697,124.898, 139.466,164.565, 220.854,229.127};
    
    thr=makeStaticData()->GetParticleMass("p")+
	makeStaticData()->GetParticleMass("n")+
	makeStaticData()->GetParticleMass("eta");
    //Better is to convert to the SECONDARY_MODELS standard

    for (int i=0;i<14;i++) {
	x_pn_eta[i]*=0.001; //in GeV
	x_pn_eta[i]+=thr;
    }


    TGraph *pn_eta_gra = new TGraph(14,x_pn_eta,y_pn_eta);
    
    PTCrossWeight * pn_eta_cross = 
	new PTCrossWeight("p + n_to_p_n_eta/tcross",
			  "Eta cross section in the p+n reaction (<120 MeV)",-1);
    pn_eta_cross->Add("p,grandparent,beam");
    pn_eta_cross->Add("n,grandparent,target");
    pn_eta_cross->Add("q,parent");
    pn_eta_cross->Add("p,daughter");
    pn_eta_cross->Add("n,daughter");
    pn_eta_cross->Add("eta,daughter");
    pn_eta_cross->SetCrossSection(pn_eta_gra,kFALSE); 
    pn_eta_cross->SetScaling(.000001); //convert from [ub] to [b]
//    pn_eta_cross->SetExcessEnergy(kTRUE,kTRUE);  //Q in MeV
    pdmutil->Add(pn_eta_cross);

    pn_eta_cross = 
	new PTCrossWeight("n + p_to_p_n_eta/tcross",
			  "Eta cross section in the p+n reaction (<120 MeV) [clone]",-1);
    pn_eta_cross->Add("dummy,parent");
    pn_eta_cross->SetCrossSection(pn_eta_gra,kFALSE); 
    pn_eta_cross->SetScaling(.000001); //convert from [ub] to [b]
//    pn_eta_cross->SetExcessEnergy(kTRUE,kTRUE);  //Q in MeV
    //BUGBUG: I want to have an automatic clone here....
    pdmutil->Add(pn_eta_cross);

    //Array from calen et al.

//     x=0.000195955, y=21.1751
// x=0.00803413, y=28.9042
// x=0.0158723, y=39.4545
// x=0.0248862, y=57.7113
// x=0.0350758, y=61.8431
// x=0.0466371, y=85.8883
// x=0.0548672, y=98.6266
// x=0.0654487, y=98.6266
// x=0.0744626, y=92.0373
// x=0.0846523, y=88.9097
// x=0.09445, y=88.9097
// x=0.10268, y=92.0373
// x=0.113262, y=96.9363

    Double_t x_d_eta[14]={0,0.000195955,0.00803413,0.0158723,0.0248862,0.0350758,
			  0.0466371,0.0548672,0.0654487,0.0744626,0.0846523,
			  0.09445,0.10268,0.113262}; //in Q [MeV]
    
    
    Double_t y_d_eta[14]={0,21.1751,28.9042,39.4545,57.7113,61.8431,85.8883,
			  98.6266,98.6266,92.0373,88.9097,88.9097,92.0373,96.9363}; //in [ub]
    thr=makeStaticData()->GetParticleMass("d")+
	makeStaticData()->GetParticleMass("eta");
    //convert in total c.m. energy
    for (int i=0;i<14;i++) {
	x_d_eta[i]+=thr;
    }

    TGraph *d_eta_gra = new TGraph(14,x_d_eta,y_d_eta);
    PTCrossWeight * d_eta_cross = 
	new PTCrossWeight("p + n_to_d_eta/tcross",
			  "Eta cross section in the p+n->d eta reaction (<120 MeV)",-1);
    d_eta_cross->Add("p,grandparent,beam");
    d_eta_cross->Add("n,grandparent,target");
    d_eta_cross->Add("q,parent");
    d_eta_cross->Add("d,daughter");
    d_eta_cross->Add("eta,daughter");
    d_eta_cross->SetCrossSection(d_eta_gra,kFALSE); 
    d_eta_cross->SetScaling(.000001); //convert from [ub] to [b]
    pdmutil->Add(d_eta_cross);

    d_eta_cross = 
	new PTCrossWeight("n + p_to_d_eta/tcross",
			  "Eta cross section in the p+n->d eta reaction (<120 MeV) [clone]",-1);
    d_eta_cross->Add("dummy,parent");
    d_eta_cross->SetCrossSection(d_eta_gra,kFALSE); 
    d_eta_cross->SetScaling(.000001); //convert from [ub] to [b]

    pdmutil->Add(d_eta_cross);

    return kTRUE;

};


Bool_t PElementaryPlugin::ExecCommand(const char * command, Double_t value) {
    return kFALSE;

}

ClassImp(PElementaryPlugin)



