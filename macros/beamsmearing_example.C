//TITLE <b>Usage of Pluto build-in templates:</b> beam smearing

{
    //This macro demonstrates how the beam smearing model
    //of Pluto can b used.
    //
    //In addition, one can play around what happens to the 
    //reconstructed eta mass if the unkown beam momentum is used


    PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
    smear->SetReaction("p + p");

    //Now one has to set the beam parameters
    //There are 3 options:
    //   * standard set (from the old Pluto)
    //     only angular smearing
    //   * A distribution function for angular
    //     smearing (the rho in polar coordinates)
    //   * A momentum smearing function
    //     (stretching, e.g. ==1 means no change)
    // 
    
    //Set "standard" parameters:
    //smear->SetBeamParameters(1., 0. , 0.5);
    //Parameters 
    // 1.) Tilt theta in degree
    // 2.) Tilt phi in degree
    // 3.) Sigma (of gaus) in degree

    //An angular spot
    //Range of theta in unites of degree
    //N.B. the radial distance probability is not included, so
    //if you want to have a flat distribution
    //one has to take this into account by multiplying with x
    //This produces a flat beam spot around around +/- 1deg
    smear->SetAngularSmearing(new TF1("delme","1    *x",0,1.));

    //momentum smearing +/- 10%
    //smear->SetMomentumSmearing(new TF1("delme","1",0.9,1.1));

    makeDistributionManager()->Add(smear);

    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0);

    //This histogram shows the beam profile:
    TH2F * histo1 = new TH2F ("histo1","Px vs. Py of beam",100,-.1,.1,100,-.1,.1);
    my_reaction.Do(histo1,"_x = [p + p]->Px(); _y  = [p + p]->Py();");
    
    //This histogram shows how a wrong assumption about the beam
    //momentum can influence the reconstruction in an exclusive
    //reaction
    TH1F * histo2 = new TH1F ("histo2","Reconstructed eta mass",100,0.3,.0.7);
    my_reaction.Do("wrong_cm = P3E(0.000000,0.000000,2.994728,4.076545);");
    my_reaction.Do(histo2,"_x = wrong_cm - ([p,1] + [p,2])->M();");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(10000);

    histo1->Draw();
    
    new TCanvas();

    histo2->Draw();
    

}
