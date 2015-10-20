//TITLE Use the long example to produce eta Dalitz events

{
    Double_t beam_energy=2.2;
    PParticle p1("p",beam_energy);
    PParticle p2("p");
    
    PParticle q=p1+p2;       //construct the beam particle
    
    // eta production
    PParticle p3("p");
    PParticle p4("p");
    PParticle eta("eta");
    PParticle *eta_part[]={&q,&eta,&p4,&p3};
    PChannel eta_prod(eta_part,3,1);
    
    // eta dalitz decay
    PParticle di_eta("dilepton");
    PParticle g_eta("g");
    
    PParticle *dalitz_part_eta[]={&eta,&di_eta,&g_eta};
    PChannel dalitz_decay_eta(dalitz_part_eta,2,1);
    
    // decay of the eta dilepton   
    PParticle em_eta("e-");
    PParticle ep_eta("e+");
    PParticle *dileptons_eta[]={&di_eta,&em_eta,&ep_eta};
    PChannel dilepton_decay_eta(dileptons_eta,2,1);
    
    PChannel *c[ ]={&eta_prod,&dalitz_decay_eta,&dilepton_decay_eta};
    PReaction r(c,"eta_dalitz",3,0,0,0,1);
    r.Print();
    
    r.loop(10000);
}
