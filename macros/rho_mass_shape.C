//TITLE Choosing the rho0 -> e+e- partial width model

{

    TCanvas c1;
    c1.SetLogy(1);


    //Remove VDM M**3 scaling
    //Not part of standard pluto -> Add this by hand
    PFixedDecay *pmodel = new PFixedDecay("rho0_novdm_e-_e+",
					  "Rho decay without VDM",-1);
    pmodel->Add("rho0, parent");
    pmodel->Add("e+, daughter");
    pmodel->Add("e-, daughter");
    makeDistributionManager()->Add(pmodel);

    //Print out what we have:
    makeDistributionManager()->Print("decay_models");

    //Choosing the model using the distribution manager
    makeDistributionManager()->Enable("rho_picutoff_e-_e+");

    //Write the primary models to the data base:
    makeDistributionManager()->LinkDB();

    //Read the primary model for the rho meson:
    PChannelModel * rho = makeDynamicData()->GetParticleModel(makeStaticData()->GetParticleID("rho0"));

    //Print rho0 properties:
    makeStaticData()->PrintParticle("rho0");

    //Set the decay index of the primary model:
    rho->SetParameter(1,85);
    rho->SetDidx(85);
    
    //Save rho model as ROOT will not evaluate the GetRandom again
    PChannelModel * rho2 = rho->Clone();
    PChannelModel * rho3 = rho->Clone();

    TH1F h1("rho0","rho0",100,0.,1.2);
    rho->GetRandom();
    for (int i=0;i<10000;i++) h1->Fill(rho->GetRandom());
 
    h1->Draw();
    
    //Change it!
    makeDistributionManager()->Enable("rho0_ee_e-_e+");
    makeDistributionManager()->LinkDB();

    TH1F h2("rho0_2","rho0_2",100,0.,1.2);
    for (int i=0;i<15000;i++) h2->Fill(rho2->GetRandom());
    h2.SetLineStyle(7);
    h2.Draw("same");

    makeDistributionManager()->Enable("rho0_novdm_e-_e+");
    makeDistributionManager()->LinkDB();

    TH1F h3("rho0_3","rho0_3",100,0.,1.2);
    for (int i=0;i<10000;i++) h3->Fill(rho3->GetRandom());
    h3.SetLineStyle(2);
    h3.Draw("same");

    //scale all histograms that they match in the limit of higher masses
    h2.Scale(h1.GetBinContent(90)/h2.GetBinContent(90));
    h3.Scale(h1.GetBinContent(90)/h3.GetBinContent(90));

}
