{
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_delta_dalitz.root");
if (!f) f = new TFile("pp_delta_dalitz.root");
TTree *Reaction = (TTree*)gDirectory->Get("data");
TClonesArray *evt=new TClonesArray("PParticle",7);
Reaction->SetBranchAddress("Particles",&evt);
PParticle *par[7];
TArrayI *Status=new TArrayI(3);
Reaction->SetBranchAddress("Filters",Status->GetArray());
const double r2d=180./3.14159265358979323846; 
TVector3 v1,v2;

TH1F *hf1=new TH1F("hf1","#Delta mass",200,1.05,1.4);
TH1F *hf2=new TH1F("hf2","#Delta lab prod angle (deg)",200,0.,45.);
TH1F *hf3=new TH1F("hf3","Dilepton mass",200,0.,.37);
TH1F *hf4=new TH1F("hf4","accepted e^{+} and e^{-} scattering angle (deg)",200,10.,90.);
TH1F *hf5=new TH1F("hf5","accepted dilepton rapidity",200,0,2);
TH1F *hf6=new TH1F("hf6","accepted e^{+}e^{-} opening angle",200,0.,160.);
c1 = new TCanvas("c1","p + p -> p + Delta -> p + p + gamma* -> p + p + e- + e+",200,10,900,900);
c1->Divide(3,2);
c1->SetFillColor(10);
c1_3->SetLogy();
c1_6->SetLogy();
hf1->SetFillColor(10);
hf2->SetFillColor(10);
hf3->SetFillColor(10);
hf4->SetFillColor(10);
hf5->SetFillColor(10);
hf6->SetFillColor(10);
Int_t nentries = Reaction->GetEntries();
for (Int_t i=0; i<nentries;i++) {
  Reaction->GetEntry(i);
  par[2]=(PParticle*)evt.At(2);
  par[4]=(PParticle*)evt.At(4);
  par[5]=(PParticle*)evt.At(5);
  par[6]=(PParticle*)evt.At(6);
  hf1->Fill(par[2]->M());
  hf2->Fill(r2d*(par[2]->Theta()));
  hf3->Fill(par[4]->M());
  if ( (*Status)[0] == 0) {
    hf4->Fill(r2d*(par[5]->Theta()));
    hf4->Fill(r2d*(par[6]->Theta()));
    hf5->Fill(par[4]->Rapidity());
    v1=par[5]->Vect();
    v2=par[6]->Vect();
    hf6->Fill(r2d*(acos(v1.Dot(v2)/(v1.Mag()*v2.Mag()))));
  }
}
c1_1->cd();
c1_1->SetFillColor(10);
hf1->Draw();
c1_2->cd();
c1_2->SetFillColor(10);
hf2->Draw();
c1_3->cd();
c1_3->SetFillColor(10);
hf3->Draw();
c1_4->cd();
c1_4->SetFillColor(10);
hf4->Draw();
c1_5->cd();
c1_5->SetFillColor(10);
hf5->Draw();
c1_6->cd();
c1_6->SetFillColor(10);
hf6->Draw();
}
