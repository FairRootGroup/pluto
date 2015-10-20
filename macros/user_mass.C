//TITLE Play with the mass shape of a dummy resonance

Int_t decay_index, pid_resonance;

double f_my_function(double *x, double * par)
{
    Double_t mass = x[0];

    //this is a dummy box:
    if ((mass < 1.2) || (mass > 1.5)) return 0;

    //add a smooth trailing edge - if wanted
    //(if not return a 1)
    Double_t width = makeDynamicData()->GetDecayPartialWidth(mass,decay_index);
    return width;
}

user_mass() {

    //Add a user-defined particle with mass = 1.25 GeV
    pid_resonance = makeStaticData()->AddParticle(-1,"MyResonance", 1.25);
    //Set the width (important to enable the m1-decay model for the production)
    makeStaticData()->SetParticleTotalWidth("MyResonance",0.25);
    //It must be a hadron:
    makeStaticData()->SetParticleBaryon("MyResonance",1);

    //Add some user-defined decays:
    decay_index = makeStaticData()->AddDecay("MyResonance --> p + pi-", "MyResonance", 
					     "p,pi-", .5 );
    makeStaticData()->AddDecay("MyResonance --> n + pi0", "MyResonance", 
			       "n,pi0", .5 );
    makeStaticData()->AddDecay("MyResonance --> n + pi0 + pi0", "MyResonance", 
			       "n,pi0,pi0", 1 );
    //N.B.: Branching ratios are re-normalized

    listParticle("MyResonance"); //List the data base content
    
    //Add a user-defined mass distribution:
    PMassSampling * model = new PMassSampling("mymodel@MyResonance","My model",-1);
    model->SetSamplingFunction(new TF1("mass_function",f_my_function,0.,10.,1));
    makeDistributionManager()->Add(model);

    //Create the reaction and a control histogram
    TH1F *histo = new TH1F("histo", "MyResonance mass", 100, 1., 1.6);
    PReaction r("1.25", "p", "p", "p MyResonance [p pi-]");
    //(Try what happens if you increase the beam energy!)
    r.Do(histo,"_x = [MyResonance]->M()");
    r.Print();
    r.Loop(10000);
    
    histo->Draw();

}
