//TITLE Omega Dalitz with decay manager in the p+d reaction

{


    PDecayChannel *c = new PDecayChannel;
    PDecayManager *pdm = new PDecayManager;
    pdm->SetVerbose(1);          // Print really useful info
    
    PParticle *p = new PParticle("p",3.5);  // proton beam
    PParticle *d = new PParticle("d");      // deuteron target
    PParticle *s = new PParticle(*p + *d);  // composite quasiparticle
    
    pdm->SetDefault("w");        // include omega decay modes
    pdm->SetDefault("pi0");      // include pi0 decay modes
    pdm->SetDefault("dilepton"); // e+e- production
    
    c = new PDecayChannel;    
    c->AddChannel(0.1,"p","d","w");        // coherent reaction
    c->AddChannel(0.9,"p","p","pi0","n");  // breakup reaction (spec. last)
    pdm->InitReaction(s,c);                // initialize the reaction
    
    pdm->loop(10000,0,"pdomega",1,0,0,1,0);   
    //arguments: num_events, weighting, filename, f0-f3, random_flag

    //data.Draw("M()","(ID()==52) * W()");
}
