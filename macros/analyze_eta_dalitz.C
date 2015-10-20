//TITLE Analyze the eta Dalitz output file

{
    // This macro uses the output of 
    // * complete_eta.C
    // * parse_eta.C
    // * hello_world.C (can have an filter)

    TFile *f =  new TFile("eta_dalitz.root");
    TTree *Reaction = (TTree*)gDirectory->Get("data");
    TClonesArray *evt=new TClonesArray("PParticle",100);
    Reaction->SetBranchAddress("Particles",&evt);
    PParticle *ep,*em;



    TH1F *hf1= new TH1F("hf2","Di-Lepton mass",100,0.,0.6);
    Int_t nentries = Reaction->GetEntries();
    if (nentries>1000) nentries = 1000; //limit number of events
    
    for (Int_t i=0; i<nentries;i++) {
	ep = NULL;
	em = NULL;
	Reaction->GetEntry(i);

	for (int j=0;j<evt->GetEntriesFast();j++) {
	    PParticle *current =(PParticle*)evt.At(j);
	    //current->Print();
	    if (current->Is("e+")) ep=current;
	    if (current->Is("e-")) em=current;
	}
	
	if (ep && em) {
	    //particles found

	    //It is very important to parse the PParticle
	    //to a TLorentzVector, because the "+" operator
	    //is reserved for a reaction (i.e. a compound
	    //particle is created
	    
	    TLorentzVector  dilepton =
		(*(TLorentzVector *) ep) + (*(TLorentzVector *) em);
	    hf1->Fill(dilepton->M());

	}

    }

    hf1->Draw();
}
