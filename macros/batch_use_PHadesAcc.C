{

    gSystem->CompileMacro( "PHadesAcc.C");
    
    TH1F * histo1 = new TH1F ("histo1","opang",50,0,3.14);
    histo1->Sumw2();
    TH1F * histo2 = new TH1F ("histo2","opang",50,0,3.14);
    histo2->Sumw2();
    
    PReaction my_reaction(3.13,"p","p","p p pi0 [dilepton [e+ e-] g]");    
    my_reaction.AddBulk(new PHadesAcc());
   
    //Use prepared variable to supress output:
    my_reaction.Do("#filter = _hadacc * _opang");
    //Use prepared variable to select histogram content:
    my_reaction.Do(histo1,"if(_opang &&  _hadacc);_ep=[e+]; _em=[e-]; _x=_ep->Angle(_em)");

    cout << my_reaction.Loop(10000) << endl;
        
    histo1->Draw("e1");

    makeDistributionManager()->Disable("helicity_angles");

    PReaction my_reaction2(3.13,"p","p","p p pi0 [dilepton [e+ e-] g]");    
    my_reaction2.AddBulk(new PHadesAcc());
   
    //Use prepared variable to supress output:
    my_reaction2.Do("#filter = _hadacc * _opang");
    //Use prepared variable to select histogram content:
    my_reaction2.Do(histo2,"if(_opang &&  _hadacc);_ep=[e+]; _em=[e-]; _x=_ep->Angle(_em)");

    my_reaction2.Loop(10000);
        
    histo2->Draw("e1same");


}
