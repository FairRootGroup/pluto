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

PParticle *p=new PParticle("p",2.62); // proton beam, Tlab (GeV)
PParticle *d=new PParticle("d");      // deuteron target
PParticle *s=new PParticle(*p + *d);
c = new PDecayChannel;                // define the decay channel
c->AddChannel(0.1, "p", "d", "eta'"); // add decay mode with weight + products
c->AddChannel(0.9, "p", "d", "w");    // add decay mode with weight + products

pdm->InitReaction(s,c);               // eta' production in p + d
Int_t n = pdm->loop(10000,0,"etap+omeg",0,0,0,1); // make events + vertices
cout << "Events processed: " << n << endl;
}



