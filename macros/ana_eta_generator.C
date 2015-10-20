//TITLE Compare eta Dalitz using weights with parse_eta.C

{

    //First, one has to run parse_eta.C and use_eta_generator.C
    //If one uses the output of hello_world.C one can
    //see very nicely the effect of a detector acceptance

    gStyle->SetOptLogy(1);


    TFile * f = new TFile("eta_dalitz.root");
    TTree *Reaction = (TTree*)gDirectory->Get("data");
    TClonesArray *evt=new TClonesArray("PParticle",7);
    Reaction->SetBranchAddress("Particles",&evt);
    PParticle *par[7];

    Int_t nentries = Reaction->GetEntries();
    TH1F *mass=new TH1F("mass","Dilepton mass",100,0,1.);
    mass->Sumw2();


    PParticle *dil,*d;

    for (Int_t i=0; i<nentries;i++) {
    
	Reaction->GetEntry(i);
    
	for (int j=0;j<7;j++) {
	    par[j]=(PParticle*)evt.At(j);
	    if (par[j]->is("dilepton")) dil    =(PParticle*)evt.At(j);
	    if (par[j]->is("D+")) d    =(PParticle*)evt.At(j);
	}

	mass->Fill(dil->M(),(Stat_t) dil->W()*mass->GetNbinsX());
    
    }


    mass->SetMinimum(0.000000001);
    mass->Draw("e1");
    mass->SetMarkerColor(2);
    mass->SetLineColor(2);


    TFile *f2 = new TFile("eta_generator.root");
    Reaction = (TTree*)gDirectory->Get("data");

    Reaction->SetBranchAddress("Particles",&evt);

    nentries = Reaction->GetEntries();
    TH1F *mass2=new TH1F("mass2","Dilepton mass",100,0,1.);
    mass2->Sumw2();

    for (Int_t i=0; i<nentries;i++) {
    
	Reaction->GetEntry(i);
    
	for (int j=0;j<7;j++) {
	    par[j]=(PParticle*)evt.At(j);
	    if (par[j]->is("dilepton")) dil    =(PParticle*)evt.At(j);
	    if (par[j]->is("D+")) d    =(PParticle*)evt.At(j);
	}

	mass2->Fill(dil->M(),(Stat_t) dil->W()*mass2->GetNbinsX());
    
    }

    mass2->SetXTitle("M [GeV/c^{2}]");
    mass2->SetYTitle("counts * BR");
    mass2->Draw("");
    
    mass->SetMarkerStyle(4); //open circles
    mass->Draw("samee1");
    mass2->SetMarkerColor(4);
    mass2->SetLineColor(4);


}
