{ // test decay manager with 1AGeV d+p collisions
  // R.H. 10/8/2000

  PDecayManager *pdm = new PDecayManager;
  pdm->SetVerbose();
//  pdm->SetDefault("D0"); // use default Delta0 decays (not implemented!)
  pdm->SetDefault("dilepton");

  PDecayChannel *d1 = new PDecayChannel;
  d1->AddChannel(1.0,"dilepton","n");
  pdm->AddChannel("D0",d1);
 
  PDecayChannel *d2 = new PDecayChannel; 
  d1->AddChannel(1.0,"dilepton","p");
  pdm->AddChannel("D+",d2);

  PDecayChannel *c = new PDecayChannel; // set up decay of d+p -> Delta -> e+e-
  c->AddChannel(0.5*0.5,"p","D0","p");
  c->AddChannel(0.5*0.5,"D+","n","p");
  c->AddChannel(0.5*0.25,"p","D+","n");
  c->AddChannel(0.5*0.25,"D+","p","n");

  PParticle *d = new PParticle("d",2.0);
  PParticle *p = new PParticle("p");
  PParticle *s = new PParticle(*d+*p);

  pdm->InitReaction(s,c);
  pdm->Print();

  pdm->loop(10000,0,"d+p",1,0,0,0);  

  
}




