//TITLE Two-step process (phi in medium -> K K -> K K_scatter)

double f_my_function(double *x, double * par)
{
    if (x[0] > 0.8) return 1;
    return 0;
}

fermi_rescattering() {

    gErrorIgnoreLevel=kWarning;

    makeDistributionManager()->Exec("nucleus_fermi:proton");
    Double_t fragment_mass = makeStaticData()->GetParticleMass("6Li") -
        makeStaticData()->GetParticleMass("n");
    makeStaticData()->AddParticle(999, "fragment",fragment_mass);

    //Now we have to do some preparation:
    PFermiDistributions* model = 
        new PFermiDistributions("6Li_fermi@6Li/fermi","Fermi momentum of 6Li",-1);
    model->SetRange(0,2.);
    makeDistributionManager()->Add(model);

    //This is the real reaction:
    PParticle * beam  =new PParticle("p",3.5);
    PParticle * target=new PParticle("7Li");
    PParticle * s= new PParticle(*beam+*target);

    //Quasi-free sub-reaction:
    PParticle * beam2  =new PParticle("p");
    PParticle * target2=new PParticle("n");
    PParticle * s2= new PParticle(*beam2+*target2);

    PParticle * spectator=new PParticle("6Li");
    PParticle * cc1[]={s,s2,spectator};

    //The outgoing products of the p-n scattering:
    PParticle * p1  =new PParticle("p");
    PParticle * p2  =new PParticle("n");
    PParticle * phi =new PParticle("phi");

    PParticle * cc2[]={s2,p1,p2,phi};

    //decay of the phi:
    PParticle * kp  =new PParticle("K+");
    PParticle * km  =new PParticle("K-");

    PParticle * cc3[]={phi,kp,km};

    //scatter the K+ again, this is like treating it as a "beam particle"
    //the "real reaction" in this case is K+ + 6Li
    //Use the "Scatter" method to keep the original pointers
    PParticle * rescatter= new PParticle("dummy");
    rescatter->Scatter(kp,spectator);

    //Next quasi-free sub-reaction:
    PParticle * beam3  =new PParticle("K+");
    PParticle * target3=new PParticle("n");
    PParticle * s3= new PParticle(*beam3+*target3);

    PParticle * spectator3=new PParticle("fragment");
    PParticle * cc4[]={rescatter,s3,spectator3};

   //The outgoing products of the final re-scattering:
    PParticle * final_kp  =new PParticle("K+");
    PParticle * final_n  =new PParticle("n");

    PParticle * cc5[]={s3,final_kp,final_n};

    //after this point our new particles and decays are pre-defined

    PFermiMomentumGA *pmodel = 
        new PFermiMomentumGA("K+n_in_6Li@K+ + 6Li_to_K+ + n_fragment",
                             "Quasi-free particle production",-1);      
    // Now add all particles
    // Define spectators and final decay products (the granddaughters)
    pmodel->Add("q,parent");                                       
    pmodel->Add("K+,grandparent,beam");                            
    pmodel->Add("6Li,grandparent,target");
    pmodel->Add("fragment,daughter,spectator");
    pmodel->Add("q,daughter,composite");
    pmodel->Add("n,granddaughter,participant"); 
    pmodel->Add("K+,granddaughter,p2"); 
    makeDistributionManager()->Add(pmodel);


    //Add a dummy angular distribution
    //In this case it is just a step function, it is without any physical meaning!
    PAngularDistribution *ang = new PAngularDistribution("my_angle","My angular dist");
    ang->Add("q,parent,reference");
    ang->Add("n,daughter");
    ang->Add("K+,daughter,primary");
    ang->Add("n,grandparent,base_reference");
    TF1 *angles2=new TF1("angles2",f_my_function,-1,1,1);
    ang->SetAngleFunction(angles2);
    makeDistributionManager()->Add(ang);


    PChannel *c1=new PChannel(cc1,2);
    PChannel *c2=new PChannel(cc2,3);
    PChannel *c3=new PChannel(cc3,2);
    PChannel *c4=new PChannel(cc4,2);
    PChannel *c5=new PChannel(cc5,2);

    PChannel *cc[]={c1,c2,c3,c4,c5};


    //Some online histograms: compare the undistorted phi shape with
    //the reconstructed one
    PReaction *r=new PReaction(cc,"fermi_rescattering",5);
    TH1F * histo1 = new TH1F ("histo1","phi mass",100,0.9,1.2);
    r->Do(histo1,"_x = [phi]->M()");
    TH1F * histo2 = new TH1F ("histo2","phi mass (scatter)",100,0.9,1.2);
    r->Do(histo2,"_x = ([K+,2] + [K-])->M()");

    r->Print();


    r->loop(10000);

    histo1->Draw();
    histo2->Draw("same");

}
