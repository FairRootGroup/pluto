{ // test thermal b sampling
  // R.H. 13/4/2000
  //
  TH1F h1("Npart","Npart",100,0.,100.);
  TH1F h2("b","Impact Parameter",100,0.,20.);
  PFireball source("pi0",0.1,0);
  source->setRandomB(200.,200.,0.03);
  for(Int_t i=0;i<100000;i++) {
     float b = source->sampleB();
     Int_t n = source->sampleNProd(b);
     h1->Fill(n);
     h2->Fill(b);
  }
  h1 = 1./100000. * h1;   // normalize to 1
  for(Int_t j=0;j<30;j++) printf("%d: %f\n",j,h1->GetBinContent(j+1));
}
