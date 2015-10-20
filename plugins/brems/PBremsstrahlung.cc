/////////////////////////////////////////////////////////////////////
//
// Parameterization of the virtual pn/pp Bremsstrahlung
// Ref. 19 and Ref. 20
// 
//
//                                  Author:  Wuestenfeld/Dohrmann
/////////////////////////////////////////////////////////////////////


#include "PBremsstrahlung.h"
#include <cmath>

Double_t PBremsstrahlung::Pol1(Double_t z, Double_t par1, Double_t par2) const
{
    Double_t x=z;
    return par1+par2*x;
}

Double_t PBremsstrahlung::Pol2(Double_t z, Double_t par1, Double_t par2, Double_t par3) const
{
    Double_t x=z;
    return par1+par2*x+par3*x*x;
}

#if 0
Double_t PBremsstrahlung::Pol3(Double_t z, Double_t par1, Double_t par2, Double_t par3, Double_t par4)
{
    Double_t x=z;
    return par1+par2*x+par3*x*x+par4*x*x*x;
}
#endif

Double_t PBremsstrahlung::DLines(Double_t x, Double_t par0, Double_t par1,Double_t par2,Double_t par3) const
{
    Double_t m,n;
    
    if(x<1.9)
    {
	m=(par1-par0)/(0.65);
	n=par0-m*1.25;
	return m*x+n;
	
    }
    else
    {
	if(x<2.2)
	{
	    m=(par2-par1)/(0.3);
	    n=par1-m*1.9;
	    return m*x+n;
	    
	}
	else
	{
	    m=(par3-par2)/(1.3);
	    n=par2-m*2.2;
	    return m*x+n;
	    
	}
	
    }
}

Double_t PBremsstrahlung::AGauss(  Double_t x,   Double_t *par) const
{
    Double_t result;
    
    if(par[0]==0) return 0;
    
    if(x-par[2]>0)
    {
	result = par[0]*exp(-0.5*pow((x-par[2])/(TMath::Max(1.e-10,((1.0 + par[3])*par[1]))),2));
    }
    else
    {
	result = par[0]*exp(-0.5*pow((x-par[2])/(TMath::Max(1.e-10,((1.0 - par[3])*par[1]))),2));
    }
    //cout << x << "," << par[0] << "," << par[1] << "," << par[2] << "," << par[3] << "," <<result << endl;
    //cout << (x-par[2]) << "," << ((1.0 + par[3])*par[1]) << "," <<  exp(-0.5*pow((x-par[2])/(TMath::Max(1.e-10,((1.0 - par[3])*par[1]))),2)) << endl;
    return result;
}


PDistribution* PBremsstrahlung::Clone(const char*) const
{
    //clone the object
    return new PBremsstrahlung((const PBremsstrahlung &)* this);
};

PBremsstrahlung::PBremsstrahlung(const Char_t *id, const Char_t *de, Int_t key) : 
    PChannelModel(id, de,key) {
    //Constructor

    if (is_channel<0)
	Fatal("PBremsstrahlung","The model (%s) should be bound to CHANNELS only",de);

    //try to figure out if the model is pp or pn
    //in addition make some consistency checks
    Int_t tid[11];
    tid[0]=10; 

    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info

    if (tid[0]!=3)
    {
	Fatal("PBremsstrahlung","Only 3 body decay");
    }
    dilepton_position = -1;
    neutron_position = -1;
    for (int i=1;i<=3; i++) {
	if (makeStaticData()->IsParticle(tid[i],"dilepton")) dilepton_position = i;
	if (makeStaticData()->IsParticle(tid[i],"n")) neutron_position = i;
    }

    if (dilepton_position == -1)
    {
	Fatal("PBremsstrahlung","No dilepton found");
    }
    threshold = 2*makeStaticData()->GetParticleMass("p");
    if (neutron_position != -1)
    {
	threshold = makeStaticData()->GetParticleMass("p")+ makeStaticData()->GetParticleMass("n");
    }

    p2_energy =2.2;
    bin=NULL;

    //...just in the case that they are needed
    mp = makeStaticData()->GetParticleMass("p");
    mn = makeStaticData()->GetParticleMass("n");
    model=0;
    author=0;
    graph=NULL;
    graph2d=NULL;
}

void PBremsstrahlung::SetMode(Char_t mode) {
    // model [0..4] where 0 => bare gNN Bremsstrahlung [default]
    //                    1 =>      gDN Delta contribution included
    //                    2 =>      cSum coherent sum of bare an gDN Bremsstrahlung
    //                    3 =>      N1520 contribution (Shyam/Mosel only)
    // The following modes are untested (IF, 26.9.2008)
    //                    4 =>      FSI final state interactions included
    //                    5 =>      VMD Vector Meson Dominance contribution
    //                    6 =>      VMD and FSI contribution included
    model=mode;
}


Bool_t PBremsstrahlung::Init(void)
{
    //looking for the dilepton. This is mandatory
    dilepton = GetParticle("dilepton");
    if (!dilepton) {
	Warning("Init","Dilepton not found");
	return kFALSE;
    }

    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    //get in/outgoing nucleons...
    p1= GetParticle("p1");
    p2= GetParticle("p2");
    p3= GetParticle("p3");
    p4= GetParticle("p4");


 

    return kTRUE;
}

Double_t PBremsstrahlung:: EvalPar(const double* x, const double* par)
{
    return Eval(x[0]);
}


Double_t PBremsstrahlung::EvalSM(Double_t x, Double_t y, Double_t z, Double_t t) const
{
        Double_t kinLim;
        Double_t M = x;
        Double_t b0, b1, b2, b3;

        if(neutron_position > -1)
        {
                // p + n case

                kinLim = TMath::Sqrt(pow((mn + mp),2)+2*mp* p2_energy)-(mn+mp);
                if(M>kinLim)// FD this is somewhat clumsy, but ...
                {
                        return 0;
                }
                if(model == gNN)
                {
                        b0 = Pol2(p2_energy, 24.4393, -7.58498, 0.819536);
                        b1 = Pol2(p2_energy, -33.4414, 14.8302, -1.75497);
                        b2 = Pol2(p2_energy, 4.13174, -5.99696, 1.95849);
                        b3 = Pol2(p2_energy, 0.942944, 1.5077, -0.207335);
                }
                else
                {
                        if(model == gDN)
                        {
                                b0 = Pol2(p2_energy, 5.40183, 4.9095, -1.06483);
                                b1 = Pol2(p2_energy, 27.9443, -25.8368, 4.04411);
                                b2 = Pol2(p2_energy, 0.859852, -1.2709, 0.437889);
                                b3 = Pol2(p2_energy, 1.3053, 0.788249, -0.109513);
                        }
                        else
                        {
                                if(model == gN1520)
                                {
                                        b0 = Pol2(p2_energy, 41.7807, -22.0668, 2.91033);
                                        b1 = Pol2(p2_energy, -41.2204, 22.8689, -3.01134);
                                        b2 = Pol2(p2_energy, 6.69874, -9.94155, 3.39074);
                                        b3 = Pol2(p2_energy, 0.795684, 1.37388, -0.172581);
                                        //cout << b0 << " " << b1 << " " << b2 << " " << b3 << endl;
                                }
                                else
                                {
                                        if(model == cSum)
                                        {
                                                b0 = Pol2(p2_energy,7.65485, 3.43015, -0.811125);
                                                b1 = Pol2(p2_energy, 14.7932, -15.5442, 2.47197);
                                                b2 = Pol2(p2_energy, 0.379128, -0.559714, 0.193903);
                                                b3 = Pol2(p2_energy, 1.60404, 0.522244, -0.0678947);
                                        }
                                        else
                                        {
                                                Error("EvalSM","Invalid Bremsstrahlung model!");
                                                return 0;
                                        }
                                }
                        }
                }
                //cout << x << endl;
                return (1./1000000.)*pow(fabs(M - kinLim),b3)/fabs((b2*exp(b0*M+b1*M*M)));
        }
        else


        {
                // p + p case
                kinLim = TMath::Sqrt(pow((mp + mp),2)+2*mp* p2_energy)-(mp+mp);
                if(M>kinLim)// FD this is somewhat clumsy again, but ...
                {
                        Error("EvalSM","M>kinLim");
                        return 0;
                }
                if(model == gNN)
                {
                        b0 = Pol2(p2_energy, 27.4811, -10.6609, 1.25601);
                        b1 = Pol2(p2_energy, -38.7574, 19.8781, -2.52863);
                        b2 = Pol2(p2_energy, 9.29954, -13.6681, 4.60585);
                        b3 = Pol2(p2_energy, 1.85104, 0.518799, -0.0646308);
                }
                else
                {
                        if(model == gDN)
                        {
                                b0 = Pol2(p2_energy, 2.13375, 7.31236, -1.44497);
                                b1 = Pol2(p2_energy, 49.7263, -40.9453, 6.24139);
                                b2 = Pol2(p2_energy, 1.93576, -2.96657, 1.02646);
                                b3 = Pol2(p2_energy, 0.621623, 1.22859, -0.160213);

                        }
                        else
                        {
                                if(model == gN1520)
                                {
                                        b0 = Pol2(p2_energy, 43.3482, -24.036, 3.24966);
                                        b1 = Pol2(p2_energy, -43.0502, 26.0555, -3.58266);
                                        b2 = Pol2(p2_energy, 8.37618, -12.8982, 4.7091);
                                        b3 = Pol2(p2_energy, 1.47027, 0.557303, -0.0442349);
                                }
                                else
                                {
                                        if(model == cSum)
                                        {
                                                b0 = Pol2(p2_energy, 6.2699, 3.50518, -0.808666);
                                                b1 = Pol2(p2_energy, 30.9533, -25.3514, 3.81197);
                                                b2 = Pol2(p2_energy, 2.07252, -3.04563, 1.02282);
                                                b3 = Pol2(p2_energy, 2.69517, -0.512592, 0.119903);
                                        }
                                        else
                                        {
                                                Error("EvalSM","Invalid Bremsstrahlung model!");
                                                return 0;
                                        }
                                }
                        }
                }
                return (1./1000000.)*pow(fabs(M - kinLim),b3)/fabs((b2*exp(b0*M+b1*M*M)));
        }
}


Double_t PBremsstrahlung::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{

    if (graph)   return graph->Eval(x);
	
    if (graph2d) return graph2d->Interpolate(x,p2_energy);

    if (author == BREMS_SHYAM_MOSEL) return EvalSM(x, y ,z ,t );

    Double_t kinLim;
    Double_t M = x;
    Double_t b0, b1, b2, b3;
    Double_t parA[4];
    Double_t parB[4];

    if(neutron_position > -1)
    {
	// p + n case
	//cout << p2_energy << endl;
	kinLim = TMath::Sqrt(pow((mn + mp),2)+2*mp* p2_energy)-(mn+mp);
	if(M>kinLim)// FD this is somewhat clumsy, but ...
	{
	    return 0;
	}
	if(model == gNN)
	{
	    b0 = Pol2( p2_energy,      16.9411, -5.93276, 0.7807);
	    b1 = Pol2( p2_energy,      -33.19, 18.8141, -2.91396);
	    b2 = Pol2( p2_energy,0.0654566, -0.102335, 0.0592914);
	    b3 = Pol2( p2_energy,  1.99472, 0.447733, -0.0727639);
	}
	else
	{
	    if(model == gDN)
	    {
		b0 = Pol2( p2_energy,   20.240726,-7.9911127,0.94649249);
		b1 = Pol2( p2_energy,   -27.103907,14.986975,-2.3138011);
		b2 = Pol2( p2_energy,0.057434153,-0.12646526,0.07317695);
		b3 = Pol2( p2_energy,   4.2753172,-1.5160941,0.33231404);
	    }
	    else
	    {
		if(model == cSum)
		{
		    b0 = Pol2( p2_energy,     12.370295,-2.379364,0.045619611);
		    b1 = Pol2( p2_energy,    -10.650187,2.7560029,-0.16295315);
		    b2 = Pol2( p2_energy,0.041419424,-0.069788419,0.038073331);
		    b3 = Pol2( p2_energy,    1.9655238,0.43345186,-0.04576372);
		}
		else
		{
		    if(model == FSI)
		    {
						
		    }
		    else
		    {
			if(model == VMD)
			{
			}
			else
			{
			    if(model == cFsiVmd)
			    {
				b0 = DLines(p2_energy,       11.73, 9.873, 6.253,7.787);
				b1 = DLines(p2_energy,    -7.589, -8.468, -30.1,-11.75);
				b2 = DLines(p2_energy,0.01877,0.03573, 0.007521,0.2954);
				b3 = DLines(p2_energy,         1.595,2.317,11.43,7.124);

				if(p2_energy > 1.9)
				{
				    parA[0] = Pol2(p2_energy, -0.65,0.09,0.116);
				    parA[1] = Pol1(p2_energy,-0.008866,0.00908);
				    parA[2] = Pol1(p2_energy,  0.7693,0.004265);
				    parA[3] = Pol1(p2_energy,   0.5262,-0.2392);
									
				    parB[0] = Pol2(p2_energy,-0.7279,0.4763,-0.04907);
				    parB[1] = Pol1(p2_energy,        -0.1483,0.09447);
				    parB[2] = Pol1(p2_energy,        0.7329,0.009148);
				    parB[3] = Pol1(p2_energy,        0.08683,-0.2189);
				}
				else
				{
				    for(Int_t i=0;i<4;i++)
				    {
					parA[i] = 0;
					parB[i] = 0;
				    }
				}
				return (1./1000000.)*pow(fabs(M - kinLim),b3)/(b2*exp(b0*M+b1*M*M)) + AGauss(p2_energy,parA) + AGauss(p2_energy,parB);
			    }
			    else
			    {
				Error("Eval","Invalid Bremsstrahlung model!");
				return 0;
			    }
			}
		    }
		    //return pow(fabs(M - kinLim),b3)/(b2*exp(b0*M+b1*M*M)) + AGauss(p2_energy,parA) + AGauss(p2_energy,parB);
		}
	    }
	}
	return (1./1000000.)*pow(fabs(M - kinLim),b3)/(b2*exp(b0*M+b1*M*M));
	//FD probably could omit 'fabs' function here 
    }
    else
    {
	// p + p case

	kinLim = TMath::Sqrt(pow((mp + mp),2)+2*mp* p2_energy)-(mp+mp);
	if(M>kinLim)// FD this is somewhat clumsy again, but ...
	{
	    return 0;
	}

	//FD 24/01/08 Take either this block (note return function)
	//FD do NOT delete this set of parameters
	/*       b0 = Pol2(Teff,33.6382,-15.3112,2.10402);
		 b1 = Pol2(Teff,-48.9918,28.9945,-4.43253);
		 b2 = Pol3(Teff,-0.700971,1.71958,-1.39966,0.389064);
		 b3 = Pol2(Teff,2.16083,0.835047,-0.216302);
		 return pow(fabs(M *M - kinLim * kinLim),b3)/(b2*exp(b0*M+b1*M*M));
		 //FD cannot omit 'fabs' function there
		 */

	//FD 24/01/08 Or that block (note different return function)
	if(model == gNN)
	{
	    b0 = Pol2( p2_energy,     28.572,-15.518,2.41939);
	    b1 = Pol2( p2_energy,  -61.4211,43.4459,-7.57627);
	    b2 = Pol2( p2_energy, 0.226395,-0.510502,0.35198);
	    b3 = Pol2( p2_energy,    4.2694,-1.5697,0.338622);
	}
	else
	{
	    if(model == gDN)
	    {
		b0 = Pol2( p2_energy,   17.5256, -5.47771, 0.389949);
		b1 = Pol2( p2_energy,    -31.254, 18.2722, -2.80298);
		b2 = Pol2( p2_energy,0.0928247, -0.188658, 0.102528);
		b3 = Pol2( p2_energy,   4.63276, -1.19151, 0.210233);
	    }
	    else
	    {
		if(model == cSum)
		{
		    b0 = Pol2( p2_energy,     16.1204,-6.72895,0.927332);
		    b1 = Pol2( p2_energy,    -13.0257,4.86514,-0.614537);
		    b2 = Pol2( p2_energy,0.0421482,-0.0873232,0.0551647);
		    b3 = Pol2( p2_energy,   2.48228,0.668239,-0.0493283);
		}
		else
		{
		    if(model == FSI)
		    {
						
		    }
		    else
		    {
			if(model == VMD)
			{
			}
			else
			{
			    if(model == cFsiVmd)
			    {
			    }
			    else
			    {
				Error("Eval","Invalid Bremsstrahlung model!");
				return 0;
			    }
			}
		    }
		}
		//return pow(fabs(M - kinLim),b3)/(b2*exp(b0*M+b1*M*M))  + AGauss(p2_energy,parA) + AGauss(p2_energy,parB);
	    }
	}
	return (1./1000000.)*pow(fabs(M - kinLim),b3)/(b2*exp(b0*M+b1*M*M));
	//FD probably could omit 'fabs' function here
    }

}

Double_t PBremsstrahlung::GetWeight(void) {
    //boost particle 2 such to be in particle 1 rest frame

    Double_t m1=p1->M();
    Double_t m2=p2->M();

    Double_t mq=parent->M();

    //m2 should be the target:
    p2_energy = (mq*mq-m1*m1-m2*m2-2*m1*m2)/(2*m2);

    if (p1->Is("n")) //neutron at rest
	p2_energy = (mq*mq-m1*m1-m2*m2-2*m1*m2)/(2*m1);

    Double_t val=Eval(dilepton->M());
    //cout << p2_energy << ":" << dilepton->M() << ":" << val << endl;
    return val;

    
}

Bool_t PBremsstrahlung::SampleMass(void)
{

    //boost particle 2 such to be in particle 1 rest frame
  
    PParticle tmp(p2);  //particle under investigation. Make better a copy
    tmp.Boost(-p1->BoostVector());
    p2_energy = tmp.KE();

    if (!bin) {
	bin = new PAdaptiveMesh(0,3,0.5);
	
	bin->SetTF1((TF1*)this);
//bin->Divide(5,3);
	bin->Divide(3,3);
	bin->Print();
    }
   

    //dilepton->SetM(this->GetRandom());

    dilepton->SetM(bin->GetRandom());


    return kTRUE;
}


ClassImp(PBremsstrahlung)
