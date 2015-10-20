{ // test decay manager with thermal pi0 source in 2AGeV Ca+Ca collisions
  // R.H. 17/8/2000

    //Make funny use of the decay manager for a thermal macro

    PDecayManager *pdm = new PDecayManager;
    pdm->SetVerbose();
    //  pdm->SetDefault("pi0");       // use default pi0 decays
    //  pdm->SetDefault("dilepton");  // use default dilepton decays

    TH1F h1("Npart","Npart",100,0.,100.);
    TH1F h2("b","Impact Parameter",100,0.,20.);
    PFireball *source1 = new PFireball("pi0",2.,0.07); // pure Boltzmann
    PFireball *source2 = new PFireball("pi0",2.,0.07); // pure Boltzmann
    source2->setRandomB(40.,40.,0.1);                  // Ca + Ca, P(pi0) = 10%
    for(Int_t i=0;i<100000;i++) {
	float b = source2->sampleB();
	Int_t n = source2->sampleNProd(b);
	h1->Fill(n);
	h2->Fill(b);
    }
    Float_t sum = 0.;
    for(Int_t j=2;j<=5;j++) sum += h1->GetBinContent(j);
    h1 = 1./sum * h1;    // normalize multiplicity bins >0 to 1
    //  for(Int_t j=1;j=<100;j++) printf("%d: %f\n",j,h1->GetBinContent(j)); 

    char* products[] = {
	"pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0",
	"pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0",
	"pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0",
	"pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0",
	"pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0","pi0"  };

    PDecayChannel *c = new PDecayChannel;  // set up decay of thermal source
    for(Int_t k=2; k<=5 && h1->GetBinContent(k)>0.01;k++) {
	printf("Mult=%d  Prob=%f\n",k-1,h1->GetBinContent(k));
	c->AddChannel(h1->GetBinContent(k),k-1,products);
    }

    pdm->InitReaction(source1,c);
    pdm->Print();

    pdm->loop(10000,0,"fire",0,0,0,0);
}






