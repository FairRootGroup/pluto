//TITLE The pp->ppeta reaction with N* contribution from DISTO
//This macro generates 4 different models for this reaction

{
    


    Bool_t eta_helicity = kTRUE;
    Bool_t pp_model = kTRUE;
    
    Int_t ascii = 1;
    
    // beam
    
    PParticle p1("p",2.2);
    //PParticle p1("p",2.85);
    PParticle p2("p");

    PParticle q=p1+p2;


    PDecayChannel * c = new PDecayChannel;
    PDecayManager * pdm = new PDecayManager;
    pdm->SetVerbose(1);
    PDistributionManager * dim = pdm->GetDistributionManager();

    //primary meson production NONRESONANT
    c->AddChannel(0.42,"p","p","eta");


    //VIA N*(1535)
    c->AddChannel(0.58,"p","NS11+");

    //IMPORTANT: If N* is decaying into more particles, ratios have to scaled!!!

    //decay of the N*
    PDecayChannel * nstar_decay = new PDecayChannel;
    nstar_decay->AddChannel(1.0,"eta","p");
    pdm->AddChannel("NS11+",nstar_decay);

    //decay of the eta
    PDecayChannel * eta_dalitz_decay = new PDecayChannel;
    eta_dalitz_decay->AddChannel(1,"g","dilepton");
    pdm->AddChannel("eta",eta_dalitz_decay);

    //decay of the virtual photon:
    PDecayChannel * eta_dilepton_decay = new PDecayChannel;
    eta_dilepton_decay->AddChannel(1.0,"e+","e-");

    pdm->AddChannel("dilepton",eta_dilepton_decay);

    //dim->Disable("eta_physics");

    TString filename("eta_dalitz_a");
    if (!eta_helicity && !pp_model) {
	dim->Disable("eta_dilepton_helicity");
	dim->Disable("pp_eta_pp_align");
	dim->Disable("pp_ns_eta_pp_align");
    } else if (eta_helicity && !pp_model) {
	filename = TString ("eta_dalitz_b");
	dim->Disable("pp_eta_pp_align");
	dim->Disable("pp_ns_eta_pp_align");
    } else if (!eta_helicity && pp_model) {
	filename = TString ("eta_dalitz_c");
	dim->Disable("eta_dilepton_helicity");
    } else {
	filename = TString ("eta_dalitz_d");
    }

    pdm->InitReaction(&q,c);

    pdm->loop(50000,0,filename.Data(),0,0,0,ascii,1);
  
}


