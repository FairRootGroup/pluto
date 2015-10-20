{
gROOT->Reset();
PDecayChannel *c;
PDecayManager *pdm = new PDecayManager; 
pdm->SetVerbose(1);                   // comment out to skip details
pdm->SetDefault("pi0");
pdm->SetDefault("dilepton");
pdm->SetDefault("pi+");

PParticle *p=new PParticle("p",2.62); // proton beam, Tlab (GeV)
PParticle *d=new PParticle("d");      // deuteron target
PParticle *s=new PParticle(*p + *d);
c = new PDecayChannel;                // define the decay channel
c->AddChannel(0.9, "p", "d", "pi0");  // add decay mode with weight + products
c->AddChannel(0.1, "p", "d", "pi+");  // add decay mode with weight + products

pdm->InitReaction(s,c);               // pi0 production in p + d
Int_t n = pdm->loop(1000,1,"weight",0,0,1,0);    // make events + vertices
cout << "Events processed: " << n << endl;
}



