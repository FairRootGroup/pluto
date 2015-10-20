//TITLE Convert an TNtuple (from the example above) to a PParticle and make a bulk decay


{
    //Using the ntuple from:
    // * hello_world.C
    // * batch_write_ntuple.C

    TFile *f = new TFile("ntuple.root");

    //Define an "empty" reaction
    PReaction my_reaction("output");
    
    my_reaction.Input(ntuple);

    //reconstruct the eta:
    my_reaction.Do("myeta = P3M(eta_px,eta_py,eta_pz,0.54745)");

    //Add eta to particle stream:
    my_reaction.Do("[eta,+] = myeta");
 


    //Bulk decay of eta
    PPlutoBulkDecay *pl = new PPlutoBulkDecay();
    pl->SetRecursiveMode(1);  //Let also the products decay
    pl->SetTauMax(0.001);     //maxTau in ns
    my_reaction.AddBulk(pl);

    //This is for debugging the decay chain:
    //my_reaction->Do("echo *******new event");
    //my_reaction->Do("foreach(*); id = [*]->ID(); echo PID: $id");


    my_reaction.Print();

    cout << my_reaction.Loop() << " events recovered" << endl;


}
