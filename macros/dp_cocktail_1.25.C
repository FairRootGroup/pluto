//TITLE <b>Total cross sections:</b> The dp cocktail at 1.25AGeV
//
//
{

//for explanations what is going on here see pp_cocktail_1.25.C
    
    makeDistributionManager()->Exec("elementary");
    makeDistributionManager()->Exec("brems: kaptari; weighting");

//makeDistributionManager()->Exec("brems: delta");
    makeDistributionManager()->Exec("brems: sum");  
//makeDistributionManager()->Exec("brems: elastic");
    
//makeDistributionManager()->Exec("brems: fsi");  
    
    PTCrossWeight * my_cross = 
	new PTCrossWeight("n + p_to_n_p_pi0_pi0/tcross",
			  "Cross section for double pi0 production",-1);
    my_cross->SetCrossSection(0.1 * 0.001 * 2); //0.1mb times 2 (see below) 
    makeDistributionManager()->Add(my_cross);
    
    makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
    makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");
    
//Once we are done, print the collection just to see if the plugins worked:    
    
    makeDistributionManager()->Print("root");

//Here we have to add all isospin channels

    PReaction *my_reaction = 
	new PReaction("2.50","d","p","n D+ [p pi0 [g dilepton [e+ e-]]] (p)","dp_1.25",1,0,0,0);
    my_reaction->AddReaction("p D0 [n pi0 [g dilepton [e+ e-]]] (p)");
    my_reaction->AddReaction("p n pi0 pi0 [g dilepton [e+ e-]] (p)");

    //Add Delta, if you want that:
    //my_reaction->AddReaction("n D+ [p dilepton [e+ e-]] (p)");
    //my_reaction->AddReaction("p D0 [n dilepton [e+ e-]] (p)");
    
    //Add Bremsstrahlung, if you want that:
    my_reaction->AddReaction("n p dilepton [e+ e-] (p)");
    //keep in mind double counting, if you enabled brems with "delta" or "sum"

//Attach a control histogram
    TH1F * pn_sum = 
	new TH1F ("pn_sum","pn DiLepton mass",100,0.,0.6);
    pn_sum->Sumw2();

    my_reaction->Do(pn_sum,"_x=[dilepton]->M()");

    my_reaction->Print();

    my_reaction->Preheating(100);
    my_reaction->Loop(10000);

    //Subtreshold eta production
    //make an extra reaction to have an independent histo
    PReaction *my_reaction2 = 
	new PReaction("2.50","d","p","n p eta [g dilepton [e+ e-]] (p)",NULL,1,0,0,0);
    my_reaction2->AddReaction("d eta [g dilepton [e+ e-]] (p)");

    my_reaction2->Print();

//Only eta
    TH1F * eta_sum = 
	new TH1F ("eta_sum","pn DiLepton mass (coherent sum)",100,0.,0.6);

    eta_sum->Sumw2();

    my_reaction2->Do(eta_sum,"_x=[dilepton]->M()");

    my_reaction2->Preheating(100);

    my_reaction2->Loop(10000);

    
//It is very crucial to correct the histo for binsize
    pn_sum->Add(eta_sum);
    PUtils::correct(pn_sum);
    PUtils::correct(eta_sum);


    TCanvas *c1 = new TCanvas("ee_mass", "ee invmass",800,800);
    c1->SetLogy(1);
    pn_sum->Draw();
    eta_sum->Draw("same");

}
