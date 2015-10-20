//TITLE <b>Setting up reactions:</b> Use the one-liner parser to produce eta Dalitz events

{
    PUtils::SetSeed(123); //this is to have a fixed SEED. By default, the systime is used....

    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0);
//If first number in quotation marks, it is the beam energy
//If not, it is the momentum 
    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(100000);
}
