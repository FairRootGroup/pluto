{
gROOT->Reset();
PDecayChannel *c;
PDecayManager *pdm = new PDecayManager; 
pdm->SetVerbose(1);                   // comment out to skip details
pdm->SetDefault("eta'");              // enable the desired channels
pdm->SetDefault("w");
pdm->SetDefault("rho0");
pdm->SetDefault("eta");
pdm->SetDefault("pi0");
pdm->SetDefault("dilepton");
pdm->SetDefault("dimuon");
pdm->SetDefault("g");

PParticle *p=new PParticle("p",2.62); // proton beam, Tlab (GeV)
PParticle *d=new PParticle("d");      // deuteron target
PParticle *s=new PParticle(*p + *d);
c = new PDecayChannel;                // define the decay channel
c->AddChannel(1.0, "p", "d", "eta'"); // add decay mode with weight + products

pdm->InitReaction(s,c);               // eta' production in p + d
//pdm->setMaxFileSize(100000);
Int_t n = pdm->loop(10000,0,"etap",0,0,0,1);    // make events + vertices
cout << "Events processed: " << n << endl;
}



