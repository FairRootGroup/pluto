//TITLE <b>Virtual detectors and batch packs</b> This macro makes a demo filter root file

{
    //Let us create a "virtual detector"

    //First, we create dummy root files
    //For a real detector, they can be derived e.g.
    //from a full GEANT simulation and re-opened
    //here from a ROOT file

    //Our dummy detector should have an acceptance
    //for electrons in forward directions
    //between 18 - 95 degrees. The efficiency
    //of electrons should be 80%,
    //in a momentum range below 100 MeV only 50,
    //and below 50 MeV no acceptance

    //Protons should have a flat acceptance of 90%
    //using a forward detector (a plane) in 
    //a distance of 2m downstream with 3x3m

    //Other particles are rejected

    //Remember, this is NOT HADES

    //This filter should come in 2 modi operandi:
    //1. For inclusive measurements (e+e-)
    //2. For exclusive measurements (ppe+e-)

    
    //Make a pseudo-Monto-Carlo simulation:
    //Generate particles between: 
    //mom = 0 ... 3000 MeV   (x)
    //cos(theta) = -1 ... 1  (y)
    //phi = 0 ... 360 deg    (z)

    //On purpose, we do not use the Pluto convention
    //to demonstrate the convertion inside the command list

    //In addition, we apply a 1% momentum resolution...
    //...just fixed, to keep it easy. Of course it is possible
    //to use matrices as well
    //This is controlled via "filter_smear_factor"

    //Create the efficiency matrices:
    TH3F * gen = 
	new TH3F ("gen","Generated events",30,0,1500,20,-1,1,20,0,360);
    TH3F * eff_elec = 
	new TH3F ("eff_elec","Efficiency matrix electrons",30,0,1500,20,-1,1,20,0,360);    
    TH3F * eff_protons = 
	new TH3F ("eff_protons","Efficiency matrix protons",30,0,1500,20,-1,1,20,0,360);

    //Start the Monte-Carlo
    int nev = 1000000;
    for (int i=0; i<nev ; i++) {
	Double_t mom = PUtils::sampleFlat()*3000;
	Double_t costheta = PUtils::sampleFlat()*2-1;
	Double_t phi = PUtils::sampleFlat()*360;
	
	gen->Fill(mom,costheta,phi);

	Double_t theta = acos(costheta)*180/TMath::Pi();

	//electrons
	if ((theta > 18) && (theta<95)) {
	    Double_t rnd=PUtils::sampleFlat();
	    if ((mom>100) && (rnd<0.8))
		eff_elec->Fill(mom,costheta,phi);
	    else if ((mom>50) && (rnd<0.5))
		eff_elec->Fill(mom,costheta,phi);	    
	}

	//protons
	Double_t x = tan(acos(costheta)) * 2 * cos(phi*TMath::Pi()/180);
	Double_t y = tan(acos(costheta)) * 2 * sin(phi*TMath::Pi()/180);
	
	if ((fabs(x) < 1.5) && (fabs(y) < 1.5)) {
	    Double_t rnd=PUtils::sampleFlat();
	    if (rnd<0.90)
		eff_protons->Fill(mom,costheta,phi);
	}


    }

    eff_elec->Divide(gen);
    eff_protons->Divide(gen);


    //Now we create the filter file
    TFile *f = new TFile("pluto_demo_filter.root","RECREATE");

    //Commands with the keyword "_main" are called when the filter is attached
    PCommandList *p = new PCommandList("_main", "echo **** This is our demo filter, v1");
    
    //A little bit of hello information:
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Usage:      This filter works in 2 modi:");
    p->AddCommand("echo Usage:      1.) inclusive:");
    p->AddCommand("echo Usage:          default");
    p->AddCommand("echo Usage:      2.) exclusive measurement:");
    p->AddCommand("echo Usage:          Add in your macro after filter attachment:");
    p->AddCommand("echo Usage:          makeDistributionManager()->Startup(\"_filter_exclusive=1\");");
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Usage:      Momentum smearing is applied, you can change:");
    p->AddCommand("echo Usage:          makeDistributionManager()->Startup(\"_filter_smear_factor=....\");");
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Debug mode:     makeDistributionManager()->Startup(\"_filter_debug=1\");");
    p->AddCommand("echo ********************************************************************************");
    //Generate variables
    p->AddCommand("_filter_debug=0");
    p->AddCommand("_filter_exclusive=0");
    p->AddCommand("_filter_smear_factor=1.");

    p->Write();

    //Commands with the keyword "_startup" are called when event loop is started
    //Here, it is used to add another warning
    PCommandList *s = 
	new PCommandList("_startup", 
			 "echo ********************************************************************************");
    
    s->AddCommand("echo This event loop uses the demo filter ");
    s->AddCommand("echo ********************************************************************************");
    s->Write();

    //Commands with the keyword "_loop" are called within event loop
    //Labels should be unique. Avoid common names
    PCommandList *l = new PCommandList("_loop");

    //The first thing we have to do is to calculate for each particle under test the axes values

    l->AddCommand("if (_filter_debug); echo **** New event called"); //do not forget the ";" after the if-statement
    l->AddCommand("_start: removed=0");                   //labels with colon
    l->AddCommand("_x = [e+]->P() * 1000.; ");  //Conversion of units
    l->AddCommand("_y = [e+]->CosTheta(); ");
    l->AddCommand("_z = (([e+]->Phi()) *180/TMath::Pi()); if (_z<0); _z = _z + 360");
    l->AddCommand(eff_elec,"eff_ep = Eval(); rnd =  sampleFlat(); if (rnd>eff_ep); [e+]->SetInActive(); removed=1");
    l->AddCommand("mom = [e+]->P(); newmom = sampleGaus(mom,0.01*mom*_filter_smear_factor);[e+]->SetMom(newmom)");
    l->AddCommand("if (_filter_debug && removed); echo e+: $_x, $_y, $_z eff: $eff_ep <removed>");
    l->AddCommand("if (_filter_debug && !removed); echo e+: $_x, $_y, $_z eff: $eff_ep");

    l->AddCommand("_x = [e-]->P() * 1000.; removed=0");
    l->AddCommand("_y = [e-]->CosTheta(); ");
    l->AddCommand("_z = (([e-]->Phi()) *180/TMath::Pi()); if (_z<0); _z = _z + 360;");
    l->AddCommand(eff_elec,"eff_em = Eval(); rnd = sampleFlat(); if (rnd>eff_em); [e-]->SetInActive(); removed=1");
    l->AddCommand("mom = [e-]->P(); newmom = sampleGaus(mom,0.01*mom*_filter_smear_factor);[e-]->SetMom(newmom)");
    l->AddCommand("if (_filter_debug && removed); echo e-: $_x, $_y, $_z eff: $eff_em <removed>");
    l->AddCommand("if (_filter_debug && !removed); echo e-: $_x, $_y, $_z eff: $eff_em");


    l->AddCommand("formore e+; goto _start");

    l->AddCommand("if (!_filter_exclusive); goto _end;");

    l->AddCommand("_start_proton_loop: removed=0");

    l->AddCommand("_x = [p]->P() * 1000.; ");
    l->AddCommand("_y = [p]->CosTheta(); ");
    l->AddCommand("_z = (([p]->Phi()) *180/TMath::Pi()); if (_z<0); _z = _z + 360;");
    //l->AddCommand("[p]->Print();");
    l->AddCommand(eff_protons,"eff_p = Eval(); rnd = sampleFlat(); if (rnd>eff_p); [p]->SetInActive(); removed=1");
    l->AddCommand("mom = [p]->P(); newmom = sampleGaus(mom,0.01*mom*_filter_smear_factor);[p]->SetMom(newmom)");
    // l->AddCommand("echo $mom , $newmom");
    //l->AddCommand("[p]->Print();");
    l->AddCommand("if (_filter_debug && removed); echo p: $_x, $_y, $_z eff: $eff_p <removed>");
    l->AddCommand("if (_filter_debug && !removed); echo p: $_x, $_y, $_z eff: $eff_p");

    l->AddCommand("formore p; goto _start_proton_loop");

    l->AddCommand("_end:");
    
    

    l->Write();
    eff_elec->Write();
    eff_protons->Write();

    f->Write();
 

}
