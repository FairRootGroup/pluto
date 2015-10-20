//TITLE <b>Hello World</b>

{

    //This macro combines several features of Pluto:
    // 1.) Setting up a reaction 
    // 2.) Adding a filter
    // 3.) Writing an additional NTuple
    // 4.) Making an online histogram

    //Create a file for the NTuple:
    TFile *f = new TFile("ntuple.root","RECREATE");
    //Create an NTuple with several variables
    TNtuple *ntuple = new TNtuple("ntuple","data from Pluto events","eta_px:eta_py:eta_pz:opang");
    
    //Create a control histo
    TH1F * histo1 = new TH1F ("histo1","dilepton mass with opening angle < 9deg",100,0.0,0.7);

    //Define the reaction: pp -> pp eta @ 3.5 kinetic beam energy
    PReaction my_reaction("3.5","p","p","p p eta [g dilepton [e+ e-]]","eta_dalitz");

    //Adding a filter
    //It is very simple: all variables starting with "#" are recognized as an file event filter
    my_reaction.Do("theta_ep = ([e+]->Theta() * 180.)/TMath::Pi()");
    my_reaction.Do("theta_em = ([e-]->Theta() * 180.)/TMath::Pi()");
    my_reaction.Do("#acc_filter = 1; if theta_ep<18 || theta_ep>85 || theta_em<18 || theta_em>85; #acc_filter = 0");
    
    //Writing variables to an NTuple
    my_reaction.Do("eta_px = [eta]->Px() ; eta_py = [eta]->Py() ; eta_pz = [eta]->Pz();");
    my_reaction.Do("opang = [e+]->Angle([e-])");
    my_reaction.Output(ntuple);

    //An additional filter on opening angle...
    my_reaction.Do("#opang_filter = 0; if opang > (9./180.)*TMath::Pi();  #opang_filter = 1");

    //Some control histo
    my_reaction.Do(histo1,"if opang > (9./180.)*TMath::Pi(); _x = ([e+] + [e-])->M()");

    cout << my_reaction.Loop(100000) << " events recorded" << endl;
    
    histo1->Draw();
    f->Write();

}
