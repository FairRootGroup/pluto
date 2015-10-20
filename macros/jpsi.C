{   //test J/Psi decays
gROOT->Reset();
PDecayChannel *c;
PDecayManager *pdm = new PDecayManager; 
pdm->SetVerbose(1);                    // comment out to skip details
pdm->SetDefault("J/Psi");
pdm->SetDefault("pi0");
pdm->SetDefault("dilepton");

PParticle *p1=new PParticle("p",25.0);  // proton beam, Tlab (GeV)
PParticle *p2=new PParticle("p");       // proton target
PParticle *s=new PParticle(*p1 + *p2);
c = new PDecayChannel;                 // define the decay channel
c->AddChannel(1.0, "p", "p", "J/Psi"); // add decay mode with weight + products

pdm->InitReaction(s,c);                // J/Psi production in p + p
Int_t n = pdm->loop(10000,0,"jpsi",0,0,0,0);    // make events
cout << "Events processed: " << n << endl;
}



