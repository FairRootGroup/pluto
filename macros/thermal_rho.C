//TITLE <b>Thermal macros (PFireball):</b> Test thermal source of rho's
// (R.H. 21/3/2006)

{ 

    //Add user-defined model without VDM M**3 scaling
    PFixedDecay *pmodel = new PFixedDecay("rho0_novdm_e-_e+", "Rho decay without VDM",-1);
    pmodel->Add("rho0, parent");
    pmodel->Add("e+, daughter");
    pmodel->Add("e-, daughter");
    makeDistributionManager()->Add(pmodel); //This enables the model

    //Enable this for removing the pi-cutoff:
    //makeDistributionManager()->Enable("rho0_ee_e-_e+");

    //Re-enabling VDM support:
    makeDistributionManager()->Enable("rho_picutoff_e-_e+");

    //Finally check what we have done:
    makeDistributionManager()->Print("decay_models");

    //PFireball *source=new PFireball("rho0", 0., 0.02, 0.0, 1, 0., 0., 0.0, 0.0, 0.0);
    PFireball *source=new PFireball("rho0",1.0,0.055,0.0,1.0,0.0,0.,0.,0.,0.); //1AGeV temp.
    source->setTrueThermal(kFALSE);

    PParticle *rho = new PParticle("rho0");
    PParticle *s1[] = {source,rho};
    PChannel  *c1 = new PChannel(s1,1,1);

    PParticle *ep = new PParticle("e+");
    PParticle *em = new PParticle("e-");
    PParticle *s2[] = {rho,ep,em};
    PChannel  *c2 = new PChannel(s2,2,1);

    PChannel  *cc[] = {c1,c2};
    PReaction *r = new PReaction(cc,"thermal_rho",2,1,1);
    r->Print();
    r->SetHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(10000);

    //data.Draw("M()","ID()==41","");

}
