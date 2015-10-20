//TITLE Writing TNtuple(s)

{

    TFile *f = new TFile("ntuple.root","RECREATE");
    TNtuple *ntuple = new TNtuple("ntuple","data from Pluto events","eta_px:eta_py:eta_pz:eta_m");
    
    TFile *f2 = new TFile("ntuple2.root","RECREATE");
    TNtuple *ntuple2 = new TNtuple("ntuple2","data from Pluto events","cos_theta_eta");

    //Define the reaction
    PReaction my_reaction("3.5","p","p","p p eta [pi+ pi- pi0]");
   
    
    my_reaction.Output(ntuple,"eta_px = [eta]->Px() ; eta_py = [eta]->Py() ; eta_pz = [eta]->Pz(); eta_m = [eta]->M()");
    my_reaction.Output(ntuple2,"myeta = [eta]; myeta->Boost([p + p]); cos_theta_eta = myeta->CosTheta()");

    my_reaction.Print();

    cout << my_reaction.Loop(100000) << " events recorded" << endl;

    f->Write();
    f2->Write();

}
