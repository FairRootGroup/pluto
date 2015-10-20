//TITLE <b>Usage of weights/generators:</B> another way of eta Dalitz simulation

{
    gROOT->Reset();
    
    //Disable old eta Dalitz model and enable it for weighting:
    makeDistributionManager()->GetDistribution("eta_dalitz")->EnableWeighting();
    
    //TF1 object representing the di-lepton statistics:
    TF1 *flat = new TF1("flat","1",0,1);
    
    //The "PInclusiveModel" can be used as a generator:
    PInclusiveModel *dilepton_generator = 
	new PInclusiveModel("flat@eta_to_g_dilepton/generator","Dilepton generator",-1);
    
    //The distribution template:
    dilepton_generator->Add("eta,parent");
    dilepton_generator->Add("g,daughter");
    dilepton_generator->Add("dilepton,daughter,primary");
    dilepton_generator->SetSampleFunction(flat);

    //Enable distribution as a generator
    dilepton_generator->EnableGenerator();
    
    makeDistributionManager()->Add(dilepton_generator);
    
    PReaction *r =new PReaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]", 
				"eta_generator",1,0,0,0);
    
    r->Print();

    //100 dummy events to have a good starting point for the normalization:
    r->Preheating(100);


    r->loop(10000);
 
}

