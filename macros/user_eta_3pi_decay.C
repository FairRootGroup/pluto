//TITLE Demo using the batch syntax to change the Dalitz plot eta -> 3pi

{

    makeDistributionManager()->Disable("eta_hadronic_decay");


    PDalitzDistribution* decay = 
	new PDalitzDistribution("my_hadronic_decay",
				"Eta matrix element for decay into charged pions");
    decay->Add("eta,    parent");
    decay->Add("pi0,    daughter,    primary");
    decay->Add("pi+,    daughter,    s1");
    decay->Add("pi-,    daughter,    s2");
    //A "step function"
    //decay->AddEquation("_f = 1.; m = (_s1 + _primary)->M2() ; echo $m; if m > 0.12; _f = 0.5");
    decay->AddEquation("_f = 1.; m = (_s1 + _primary)->M2(); if m > 0.12; _f = 0.2");
    decay->SetMax(1);
    makeDistributionManager()->Add(decay);
 

    PReaction my_reaction("2.2","p","p","p p eta [pi+ pi- pi0]",NULL,1,0,0,0);
    
    
    TH2F *hf1= new TH2F("hf1","",100,0.06,.2,100,0.06,.2);
    my_reaction.Do(hf1,"_x = ([pi-] + [pi0])->M2() ; _y = ([pi+] + [pi0])->M2()");
    my_reaction.Print();
    my_reaction.Loop(10000);
    hf1->Draw("box");


}
