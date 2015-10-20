//TITLE <b>The event loop interface:</b> Make histograms in one single line without ana macro with the PBatch script language

{

    PData::SetWeightVersion(0); //Disable 1/N_ev * BR weighting

    //First define our histograms
    TH2F * histo = new TH2F ("histo","Dalitz",30,1.,4.,30,1.,4.);
    TH1F * histo1 = new TH1F ("histo1","p/pi0 missing mass",100,0.1,2.0);
    TH1F * histo2 = new TH1F ("histo2","cos theta of pp",20,-1.,1.);
    TH1F * histo3 = new TH1F ("histo3","cos theta of D+ decay",20,-1.,1.);

    TH2F * histo_p1 = new TH2F ("histo_p1","px vs. vy of the pp pair",30,-1.,1.,30,-1.,1.);

    //If you want another decay angle of the D, uncomment this line:
    ((PScatterDistribution *)makeDistributionManager()->GetDistribution("pp_delta_waves1"))->SetAngleFunction(new TF1("delme","(1+3*x*x)*.25",-1,1));
 
    //No D+ angular distribution?
    //makeDistributionManager()->Disable("pp_delta+_angle");

    //Define the reaction as usual
    PReaction my_reaction(3.13,"p","p","p D+ [p pi0]", "delme",1,0,0,0);

    //Combine the masses of p,1 and pi0, p,2, pi0 to the Dalitz plot
    my_reaction.Do(histo,"_x = ([p,1] + [pi0])->M2() ; _y = ([p,2] + [pi0])->M2() ");

    //Missing mass of the p2 and pi0 pair:
    my_reaction.Do(histo1,"miss= [p + p]- ( [p,2]+ [pi0] );_x=miss->M()");

    //pp aligment angle: Important: copy the proton before boosting it, otherwise you
    //will have the boosted particle on tape
    my_reaction.Do(histo2,"p1=[p,1]; p1->Boost([p + p]); _x=cos(p1->Theta())");

    //The decay angle of the D+:
    my_reaction.Do(histo3,"_pi0=[pi0]; _d=[D+]; _pi0->Rot(_d); _d->Rot(_d); _pi0->Boost(_d); _x= cos(_pi0->Theta())");

    //Momentum of a pair
    //All methods of a PParticle of the type ->XXX(void) can be used
    //Syntax can be ->XXX() or .XXX()  (as you like)
    my_reaction.Do(histo_p1,"pp_pair= [p,1] + [p,2];_x=pp_pair.Px(); _y=pp_pair->Py()");

    my_reaction.Loop(50000);


    TCanvas c1;

    c1.Divide(2,2);

    c1.cd(1); 
    gStyle->SetPalette(1,0);
    histo->Draw("colz");
    c1.cd(2); 
    histo2->Draw("e1");
    c1.cd(3); 
    histo1->Draw("e1");
    c1.cd(4); 
    histo3->Draw("e1");

    TCanvas c2;

    histo_p1->Draw("colz");

}
