//TITLE Re-read an TNtuple (from the example above) and make a histogram

{

    //Reads ntuple written by
    // * hello_world.C
    // * batch_write_ntuple.C


    TFile *f = new TFile("ntuple.root");
    TH1F * histo = new TH1F ("histo","cos theta of eta",100,-1.,1.);

    //Define an "empty" reaction
    PReaction my_reaction;
    
    my_reaction.Input(ntuple);

    //reconstruct the eta:
    my_reaction.Do("myeta = P3M(eta_px,eta_py,eta_pz,0.54745)");

    //NB: The following numbers can be obtained via:
    //PParticle p("p"); PParticle p2("p",3.5); q=p+p2 ; q.Print()
    my_reaction.Do("cm = P3E(0.000000,0.000000,4.337961,5.376545) ; myeta->Boost(cm);");
    my_reaction.Do(histo,"if (myeta->M() > 0); _x= myeta->CosTheta();");
    
    my_reaction.Print();

    cout << my_reaction.Loop() << " events recovered" << endl;
    
    histo->Draw();

}
