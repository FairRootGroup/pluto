//TITLE pp elastic scattering with different beam tilts

{

    PData::SetWeightVersion(0); //Disable 1/N_ev * BR weighting

    TH1F * histo2 = new TH1F ("histo2","cos theta of pp",20,-1.,1.);

    PReaction my_reaction("_T1=0.4; _T2=0.4; ","p","p","p p");
    //PReaction my_reaction("_T1=0.4; _T2=0.4; _theta1=90*TMath::DegToRad(); _theta2=-90*TMath::DegToRad(); ","p","p","p p");

    my_reaction.Do(histo2,"p1=[p,1]; p1->Boost([p + p]); _x=cos(p1->Theta())");

    my_reaction.Print();

    my_reaction.Loop(5000);

    histo2->Draw("e1");
}
