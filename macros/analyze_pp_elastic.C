{

    //Analyze the output of pp_elastic.C

    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_elastic.root");
    if (!f) {  f= new TFile("pp_elastic.root"); }
    TTree *Reaction = (TTree*)gDirectory->Get("data");
    TClonesArray *evt=new TClonesArray("PParticle",3);
    Reaction->SetBranchAddress("Particles",&evt);
    const double r2d=180./3.14159265358979323846; 
    PParticle *par[3];
    TH1F *hf1=new TH1F("hf1","#Theta_{lab} (deg.)",20,0,0.9);
    TCanvas *can = new TCanvas("canvas","pp elastic",200,10,600,400);
    can->SetLogy(1);
    Int_t nentries = Reaction->GetEntries();
    for (Int_t i=0; i<nentries;i++) {
	Reaction->GetEntry(i);
	par[0] = (PParticle*)evt.At(0);
	par[1] = (PParticle*)evt.At(1);
	par[2] = (PParticle*)evt.At(2);
	par[1]->Boost(-par[0]->BoostVector());
	par[2]->Boost(-par[0]->BoostVector());
	hf1->Fill(fabs(cos(par[1]->Theta())));
    } 
    hf1->Draw(); 


}

