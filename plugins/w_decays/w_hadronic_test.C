//TITLE Test macro for omega -> pi+pi-pi0

{

    gStyle->SetPalette(1);
    gStyle->SetOptStat(1);
    gStyle->SetOptTitle(0);

    
    PReaction my_reaction("3.5","p","p","p p w [pi+ pi- pi0]");
    TH2F *histo2 = new TH2F ("histo2","DalitzOmega",100,0.,0.5,100,0.,0.5);
    
    my_reaction.Do(histo2, "_x = ([pi0] + [pi-])->M2() ; _y = ([pi+] + [pi0])->M2()");
    
    my_reaction.Print();
    my_reaction.Loop(100000);
    histo2->Draw("colz");
}
