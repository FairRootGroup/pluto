//TITLE <b>Total cross sections:</b> The pp cocktail at 1.25GeV
//
//This macros demonstrates how to set up a cocktail with total cross sections
//
{

    
//The first action we have to take is to activate the plugins
//this is done in Pluto via commands. In the same stap we
//can already forward the options to the plugins

//As an first example let us activate the elementary plugin
//This adds a lot of composite particles into the data base
//and some models for total cross sections we need later

    makeDistributionManager()->Exec("elementary");


//Next step is to activate the bremsstrahlung.
//In our case the Kaptari/Kaempfer version
//We want to use weighting, otherwise there is no total
//normalization

    makeDistributionManager()->Exec("brems: kaptari; weighting");

//There are 3 "sub-versions": Please choose one out
//of them:

//makeDistributionManager()->Exec("brems: delta");
makeDistributionManager()->Exec("brems: sum");  
//makeDistributionManager()->Exec("brems: elastic");

//makeDistributionManager()->Exec("brems: fsi");  

//The Delta production is done via
//the Teis model
//for double pi0 production there is model
//let us add a model
//The corresponding decay is added automatically into the data base
    PTCrossWeight * my_cross = 
	new PTCrossWeight("p + p_to_p_p_pi0_pi0/tcross",
			  "Cross section for double pi0 production",-1);
    my_cross->SetCrossSection(0.2 * 0.001 * 2); //0.2mb times 2 (see below) 
    makeDistributionManager()->Add(my_cross);
    
//Activate the Krivoruchenko Dalitz
    makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");

//Enable the generator mod: It is very usefull to enable a flat
//generator to acquire enough statistics for the high mass region

//There is another topic in this context: For broad resonances we
//should disable the static b.r. scaling and use the dGamma/dM as
//a mass-dependent b.r.
//Let us put the threshold to 100 MeV

    makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");
    
//The krivoruchenko mod enables by default also the QED form factor
//If you want the VMD version, uncomment the following line:
//makeDistributionManager()->Enable("vmd");

//Once we are done, print the collection just to see if the plugins worked:    
    
    makeDistributionManager()->Print("root");

//Set up the reaction:
//N.B: if you want to write a file, add a filename at the "NULL" position

    PReaction *my_reaction = 
//	new PReaction("1.25","p","p","p D+ [p pi0 [g dilepton [e+ e-]]]",NULL,1,0,0,0);
	new PReaction("1.25","p","p","p D+ [p pi0 [g dilepton [e+ e-]]]","pp_1.25",1,0,0,0);
    my_reaction->AddReaction("p p pi0 pi0 [g dilepton [e+ e-]]");
    //pi0 decay as defined above. Double pi0 already taken into account
    //-> therefore we leave one pi0 undecayed

    //Add Delta, if you want that:
    my_reaction->AddReaction("p D+ [p dilepton [e+ e-]]");
    //Add Bremsstrahlung, if you want that:
    //my_reaction->AddReaction("p p dilepton [e+ e-]");
    //keep in mind double counting, if you enabled brems with "delta" or "sum"

    my_reaction->Print();

//Attach a control histogram
    TH1F * pp_sum = 
	new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",100,0.,0.6);
    pp_sum->Sumw2();

    my_reaction->Do(pp_sum,"_x=[dilepton]->M()");

    my_reaction->Preheating(100);

    my_reaction->Loop(10000);

    
//It is very crucial to correct the histo for binsize
    PUtils::correct(pp_sum); 


    TCanvas *c1 = new TCanvas("ee_mass", "ee invmass",800,800);
    c1->SetLogy(1);
    pp_sum->Draw();

}
