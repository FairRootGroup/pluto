//TITLE Analyze the polar angle as created with useAngularDistribution

{
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("angular_distribution.root");
if (!f) {  f= new TFile("angular_distribution.root"); }
TTree *Reaction = (TTree*)gDirectory->Get("data");

TClonesArray *evt=new TClonesArray("PParticle",20);
Reaction->SetBranchAddress("Particles",&evt);
const double r2d=180./3.14159265358979323846; 
PParticle *par[3];
TH1F *hf1=new TH1F("hf1","#Theta_{lab} (deg.)",20,-1.,1.);
//TH1F *hf1=new TH1F("hf1","#Theta_{lab} (deg.)",20,0,1.);
hf1->SetMinimum(0);
Int_t nentries = Reaction->GetEntries();
for (Int_t i=0; i<nentries;i++) {
  Reaction->GetEntry(i);
  par[0] = (PParticle*)evt.At(0);
  par[1] = (PParticle*)evt.At(1);
  par[2] = (PParticle*)evt.At(2);
  par[1]->Boost(-par[0]->BoostVector());
  par[2]->Boost(-par[0]->BoostVector());
//  hf1->Fill(fabs(cos(par[1]->Theta())));
  hf1->Fill((cos(par[1]->Theta())));
} 
hf1->Draw("e1"); 
}

