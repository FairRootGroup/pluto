//TITLE <b>Angular distribution interface:</b> Add user-defined angular distributions

double f_my_function(double *x, double * par)
{
    // (1+x*x)/2
    return (1+x[0]*x[0])/2;
    // It is important that the test fuction does not exceed 1 (but comes close to it)
}


void useAngularDistribution(void) {

    //define the angular function
    //it should not exceed 1!!
    //This can be a parsed TF1 object:
    TF1 *angles1=new TF1("angles1","(1+x*x)/2",-1,1);
    //or a hand-made function:
    TF1 *angles2=new TF1("angles2",f_my_function,-1,1,1);
    //or based on a TGraph:
    const Int_t n = 3;
    Double_t x[n], y[n];
    for (Int_t i=0;i<n;i++) {
	x[i] = i*(1./Double_t(n-1));
	y[i] = (1+x[i]*x[i])/2;
    }
    //another idea:
    y[0]=1.;y[1]=0.;y[2]=1.;

    angles3 = new TGraph(n,x,y);




    PAngularDistribution *ang = new PAngularDistribution("my_angle","My angular dist");
    ang->Add("q,parent,reference");
    ang->Add("p,daughter");
    ang->Add("NP11+,daughter,primary");

    //uncomment here one version
    //ang->SetAngleFunction(angles1);
    //ang->SetAngleFunction(angles2);
    ang->SetAngleFunction(angles3,kTRUE,kTRUE);
    //ang->SetAngleFunction(angles3,kFALSE,kTRUE);
    //options for TGraph version: useSpline, useSymmetry
    ang->Print();

    makeDistributionManager()->Add(ang);

    PReaction my_reaction (3,"p","p","p NP11+ [p pi0]","angular_distribution",1,0,0,0);


    makeDistributionManager()->Print();
    my_reaction.Print();
    ang->Draw();
    my_reaction.Loop(10000);
}
