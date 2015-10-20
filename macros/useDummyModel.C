//TITLE Use the dummy channel model

{

    //Add a bunch of dummy particles
    makeStaticData()->AddParticle(-1,"A", 1.2);
    makeStaticData()->SetParticleTotalWidth("A",0.3);
    makeStaticData()->SetParticleBaryon("A",1);
    makeStaticData()->AddParticle(-1,"b", 0.5);
    makeStaticData()->AddParticle(-1,"c", 0.3);


    makeStaticData()->AddDecay(-1,"A -> b + c", "A", "b,c", 1.);
    makeStaticData()->AddDecay(-1,"A -> b + b", "A", "b,b", 1.);

    //Load the model into the root session
    //Here, only compiled code works, the ".L" does not work
    gSystem->CompileMacro( "./PDummyModel.C");

    //Construct the model from the class
    //Use the parser
    PDummyModel *newmodel = new PDummyModel("A_dummy_b_c", "My private model",-1);

    //make it known to the Pluto world:
    makeDistributionManager()->Add(newmodel);
    makeDistributionManager()->Print("user");

    //We look to the decay A -> b+b
    PReaction my_reaction(3.13,"p","p","p A [b b]", "mydecay",1,0,0,0);
    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(20000);

    //data.Draw("M()","ID()==70","e1");
}
