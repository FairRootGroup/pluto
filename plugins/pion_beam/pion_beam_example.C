//TITLE Example how to use the pion beam plugin

{

    gSystem->Load("../../libPluto.so");
    //Enable the pion beam plugin:
    makeDistributionManager()->Exec("pion_beam");
    
    //Set the term
    PPionBeamAmplitude* ampl_orig=(PPionBeamAmplitude*) makeDistributionManager()->GetDistribution("Pi_minusBeamAmplitude");
    ampl_orig->SetTerm(0);
    //0=coherent sum, 1=rho, 2=omega

    //Add particle smearing if needed
    //PHadesParticleSmearer * smearer = new PHadesParticleSmearer();
    //smearer->SetResolutionFactor(2);

    TH1F *histo1 = new TH1F("histo1","Dilepton Mass",100,0.,0.8);
    histo1->Sumw2();

    //0.832 => sqrt=1.65
    PReaction my_reaction("0.832","pi-","p","n dilepton [e+ e-]");
    //my_reaction.AddBulk(smearer);
    
    my_reaction.Do(histo1,"_x = ([e+]+[e-])->M()");
    my_reaction.Print();   
    my_reaction.Loop(10000);

    PUtils::correct(histo1); 
    histo1->Draw("");

}
