//TITLE <PPlutoBulkDecay> Demonstrates the bulk decay option

{
    PReaction my_reaction(6,"p","p","p NS11+","n1535_sample_bulk",1,0,0,0);

    PPlutoBulkDecay *pl = new PPlutoBulkDecay();
    pl->SetRecursiveMode(1);  //Let also the products decay
    pl->SetTauMax(0.001);     //maxTau in ns
    my_reaction.AddBulk(pl);

    my_reaction.loop(100000);

}
