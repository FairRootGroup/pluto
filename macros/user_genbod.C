//TITLE Add a user-defined mass distribution to the genbod decay

{
    
    //Add a bunch of dummy particles
    makeStaticData()->AddParticle(-1,"A", 1.2);
    makeStaticData()->AddParticle(-1,"a", 0.1);
    makeStaticData()->AddParticle(-1,"b", 0.2);
    makeStaticData()->AddParticle(-1,"c", 0.3);


    makeStaticData()->AddDecay(-1,"A -> a b c ", "A", "a,b,c",1.);

    PGenBod *d = new PGenBod("my_genbod@A_genbod_a_b_c/genbod","User Genbod",-1);
    d->Add("A,  parent");
    d->Add("a,  daughter,corr1");
    d->Add("b,  daughter,corr2"); // --> these 2 particles are formed to have the mass distribution as below
    d->Add("c,  daughter");
    makeDistributionManager()->Add(d);

    PFunction *f = new PFunction("my_corr@A_genbod_a_b_c/correlation",
				 "User mass correlation",-1);
    //f->SetFunction(new TF1("f","x*x",
    // makeStaticData()->GetParticleMass("a") +
    // 	makeStaticData()->GetParticleMass("b"),
    // 	makeStaticData()->GetParticleMass("A")));

    f->AddEquation("_f = _x * _x");
    f->SetRange(makeStaticData()->GetParticleMass("a") +
		makeStaticData()->GetParticleMass("b"),
		makeStaticData()->GetParticleMass("A"));


    makeDistributionManager()->Add(f);

    TH1F *hf1= new TH1F("hf1","",100,0.1,1.2);
    PReaction my_reaction("_T1=3.5","p","p","p A [a b c]");
    my_reaction.Print();
    my_reaction.Do(hf1,"_x = ([a] + [b])->M() ; ");
    my_reaction.Loop(100000);
    
    hf1->Draw("");
    
    


}
