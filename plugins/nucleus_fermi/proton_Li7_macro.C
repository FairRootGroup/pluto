//TITLE Example for a p + 7Li reaction

{

    makeDistributionManager()->Exec("nucleus_fermi:proton");

    PReaction my_reaction("3.5","p","7Li","p (n) phi [K+ K-] (6Li)","delme");
    //both participant and spectator must be in brackets "( )"
    //the last particle is the spectator

    my_reaction.Print();

    //check momentum distribution of K+
    TH1F * histo1 = new TH1F ("histo1","K+ momentum",100,0,3);
    my_reaction.Do(histo1,"_x = [K+]->P()");

    //mass of phi
    TH1F * histo2 = new TH1F ("histo2","phi mass",100,0.9,1.2);
    my_reaction.Do(histo2,"_x = [phi]->M()");

    my_reaction.Loop(1000);

    histo2->Draw();


    //same at much lower energy
    //slightly below threshold, this is really time consuming
    PReaction my_reaction2("2.5","p","7Li","p (n) phi [K+ K-] (6Li)","delme");
    TH1F * histo2a = new TH1F ("histo2a","phi mass",100,0.9,1.2);
    my_reaction2.Do(histo2a,"_x = [phi]->M()");
    my_reaction2.Loop(1000);
    histo2a->Draw("same");

}
