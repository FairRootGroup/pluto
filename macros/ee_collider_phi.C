//TITLE e+e- -> phi

{

    //create phi in limits:
    PReaction my_reaction("_T1=0.51  ; _T2=0.51; _theta1=2.*TMath::DegToRad(); _theta2=2.*TMath::DegToRad(); _phi=20.*TMath::DegToRad();","e+","e-","phi", "test",1,0,0,1);

    //or not within limits:
    //PReaction my_reaction("_T1=0.7  ; _T2=0.7; _theta1=2.*TMath::DegToRad(); _theta2=2.*TMath::DegToRad(); _phi=20.*TMath::DegToRad();","e+","e-","phi", "test",1,0,0,1);
    

    my_reaction.Do("[e+ + e-]->Print(); [phi]->Print();");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(1);
}
