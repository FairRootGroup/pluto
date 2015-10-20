//root -l 3He_breakup_1500_MeV_in_GeV.root analyze_g3He.C

{

  gROOT->LoadMacro("KineticRecoilEnergy.C");

  TFile *f = (TFile*) gFile;

  TTree *Reaction = (TTree*)gDirectory->Get("data");

  TClonesArray *evt=new TClonesArray("PParticle", 11);
  Reaction->SetBranchAddress("Particles", &evt);

  PParticle *par[11];

  const double r2d=180./3.14159265358979323846; 

  TH1F *hMomAf       = new TH1F("hMomAf"     ,"Momentum final Nucleus"     , 1000,0.,1.);
  TH1F *hMomAfRec    = new TH1F("hMomAfRec"  ,"Fermi Momentum Reconstructed",1000,0.,1.);

  TH1F *hThetaEta    = new TH1F("hThetaEta","Polar angle theta of the Eta",180,0.,180.);
  TH1F *hThetaPi0    = new TH1F("hThetaPi0","Polar angle theta of the Pi0",180,0.,180.);
  TH1F *hThetaPro    = new TH1F("hThetaPro","Polar angle theta of the Proton",180,0.,180.);
  TH1F *hPhiEta      = new TH1F("hPhiEta","Azimuthal angle phi of the Eta",360,-180.,180.);
  TH1F *hPhiPi0      = new TH1F("hPhiPi0","Azimuthal angle phi of the Pi0",360,-180.,180.);
  TH1F *hPhiPro      = new TH1F("hPhiPro","Azimuthal angle phi of the Proton",360,-180.,180.);
 
  TCanvas *can = new TCanvas("can","g+3He->Pi0+Eta+p+d",200,10,600,400);
  can->Divide(3, 3, .001, .001);

  TLorentzVector NewAfinal;
  TLorentzVector GrandParent;
  TLorentzVector RecProton;
  TLorentzVector Pi0plusEta;
  TLorentzVector NewProton;
  TLorentzVector Eta2g;
  TLorentzVector Pi02g;
  TLorentzVector Ainitial;
  TLorentzVector Afinal;
  TLorentzVector ggg(0,0,1.5,1.5);
  
  TVector3 MomAf;
  TVector3 NewAfinalMom;

  double massP;
  double massAi;
  double massAf;
  double BeamEnergy;

  Int_t nentries = Reaction->GetEntries();
  for (Int_t i=0; i<nentries;i++) 
  {
      Reaction->GetEntry(i);
      par[0] = (PParticle*)evt->At(0);
      par[1] = (PParticle*)evt->At(1);
      par[2] = (PParticle*)evt->At(2);
      par[3] = (PParticle*)evt->At(3);
      par[4] = (PParticle*)evt->At(4);
      par[5] = (PParticle*)evt->At(5);
      par[6] = (PParticle*)evt->At(6);
      par[7] = (PParticle*)evt->At(7);
      par[8] = (PParticle*)evt->At(8);
      par[9] = (PParticle*)evt->At(9);
      par[10] = (PParticle*)evt->At(10);

      Eta2g         = par[9]->Vect4() + par[10]->Vect4(); 
      Pi02g         = par[7]->Vect4() + par[8]->Vect4();
      Pi0plusEta    = Eta2g + Pi02g; 
      RecProton     = par[6]->Vect4();

      Ainitial      = par[0]->Vect4() - ggg;
      Afinal        = par[1]->Vect4();

      GrandParent   = par[0]->Vect4();
      massAi        = Ainitial.M();
      massAf        = Afinal.M();
      massP         = par[6]->M();

      MomAf         = par[1]->Vect();

      BeamEnergy    = ggg.E();
      NewProton     = doCalEnergy(BeamEnergy,Pi0plusEta,RecProton,massAi,massP,massAf);

      NewAfinal     = GrandParent - Eta2g - Pi02g - NewProton;
      NewAfinalMom  = NewAfinal.Vect();

      hMomAfRec     -> Fill( NewAfinalMom.Mag() );
      hMomAf        -> Fill( MomAf.Mag() );

      hThetaEta->Fill(Eta2g.Theta()*r2d);
      hThetaPi0->Fill(Pi02g.Theta()*r2d);
      hThetaPro->Fill(par[6]->Theta()*r2d);
      hPhiEta->Fill(Eta2g.Phi()*r2d);
      hPhiPi0->Fill(Pi02g.Phi()*r2d);
      hPhiPro->Fill(par[6]->Phi()*r2d);

  }

  double px,py,pz;
  TH1F *h1 = new TH1F("h1","Distribution of Fermi Momentum inside Helium-3",1000,0.,1.);
  makeDistributionManager()->Exec("nucleus_fermi:gamma");
  makeDistributionManager()->LinkDB();
  PFermiMomentumGA* f1 = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gp_in_3He"));   
  for(int i=0;i<10000;i++) h1->Fill(f1->GetRandomFermiMomentum(px,py,pz));
 
  can->cd(1);
  hThetaEta->SetLineColor(3);
  hThetaEta->Draw();
  can->cd(2);
  hThetaPi0->SetLineColor(3);
  hThetaPi0->Draw();
  can->cd(3);
  hThetaPro->SetLineColor(3);
  hThetaPro->Draw();
  can->cd(4);
  hPhiEta->SetLineColor(4);
  hPhiEta->Draw();
  can->cd(5);
  hPhiPi0->SetLineColor(4);
  hPhiPi0->Draw();
  can->cd(6);
  hPhiPro->SetLineColor(4);
  hPhiPro->Draw();

  can->cd(7);
  hMomAf->SetLineColor(2);
  hMomAf->Draw();

  can->cd(8);
  hMomAfRec->SetLineColor(2);
  hMomAfRec->Draw();

  can->cd(9);
  h1->SetLineColor(2);
  h1->Draw();

  can->Print("g3He_BreakUp.eps");
  can->Print("g3He_BreakUp.root");

  can->Modified();

}
