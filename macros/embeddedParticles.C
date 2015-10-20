//TITLE <PEmbeddedParticles> Adds embedded particles in the eta Dalitz decay

{

    PReaction my_reaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]","eta_dalitz_embedded",1,0,0,0);

    //Construct the bulk container:
    PEmbeddedParticles * embedded = new PEmbeddedParticles();

    //Add an e+ which we emit at a single point:
    PParticle * e_plus = new PParticle("e+",1.,2.,3.);  
    //Just add the particle to the container:
    embedded->AddParticle(e_plus);

    //We can also add a "white" dilepton, which we emit in a small cone:
    PParticle * dilepton = new PParticle("dilepton");
    embedded->AddParticle(dilepton);
    embedded->SetSampling(0, 1.,   //pmin and pmax in lab frame 
			  TMath::Pi()/1000., //opening angle
			  TMath::Pi()/2.,    //Theta of pointing vect.
			  TMath::Pi()/2.,    //Phi of pointing vect.
			  0.2, 1.5  //Mass sampling (optional)
			  ); 

    //Add our container to the reaction:
    my_reaction.AddBulk(embedded);

    my_reaction.loop(100000);
}
