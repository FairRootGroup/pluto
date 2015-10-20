//TITLE <b>Cocktails for 3.5GeV:</b> the pp reaction

{
    gROOT->Reset();

    PDecayChannel *c;
    PDecayManager *pdm = new PDecayManager;

    PParticle *beam   = new PParticle(14,3.5);      // beam  14-proton
    PParticle *target = new PParticle(14);          // target
    PParticle *s      = new PParticle(*beam + *target);


    c = new PDecayChannel;                    //define reaction channels
    c->AddChannel(0.61,"p","D+");
    c->AddChannel(0.03,"p","p","eta");                
    c->AddChannel(0.011,"p","p","w");
    c->AddChannel(0.02,"p","p","rho0");

    c_d= new PDecayChannel;
    c_d->AddChannel(1.0,"e+","e-");
    pdm->AddChannel("dilepton",c_d);

    c_pi0= new PDecayChannel;                      //pi0
    c_pi0->AddChannel(.012,"g","dilepton");
    pdm->AddChannel("pi0",c_pi0);

    c_delta= new PDecayChannel;
    c_delta->AddChannel(6.0e-5,"p","dilepton");    // Delta
    c_delta->AddChannel(0.666,"p","pi0");    // add decay mode with weight + products
    pdm->AddChannel("D+",c_delta);

    c_eta= new PDecayChannel;
    c_eta->AddChannel(6.0e-3,"g","dilepton");    // eta
    pdm->AddChannel("eta",c_eta);

    c_omega= new PDecayChannel;
    c_omega->AddChannel(7.0e-5,"e+","e-");    // omega
    c_omega->AddChannel(5.9e-4,"pi0","dilepton");
    pdm->AddChannel("w",c_omega);

    c_rho= new PDecayChannel;                  //rho0
    c_rho->AddChannel(4.5e-5,"e+","e-");
    pdm->AddChannel("rho0",c_rho);


    pdm->InitReaction(s,c);               // cocktail production in p + p

    Int_t n = pdm->loop(10000,1,"pp35_cocktail",1,0,0,1,1); // make events + vertices
    cout << "Events processed: " << n << endl;

    //Draw the spectrum:
    //data.Draw("M()","ID() == 51 || ID()==52 || ID()==41");

}

