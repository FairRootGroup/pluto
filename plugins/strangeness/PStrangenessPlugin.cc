/////////////////////////////////////////////////////////////////////
//Plugin to add the strangeness particles Lambda(1405) and Sigma(1385)
//
//Pids are 70-73
//
//This plugin is automatically enabled. If you don't like it, call
//makeDistributionManager()->Disable("strangeness");
//in the very beginning of your macro
//
//////////////////////////////////////////////////////////////////////

#include "PStrangenessPlugin.h"




PStrangenessPlugin::PStrangenessPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {
    is_initialized =0;
    //    RequiresPlugin("elementary");
}

Bool_t PStrangenessPlugin::Activate(void) {

    return kTRUE;
}

PStrangenessPlugin::~PStrangenessPlugin() {
}


Bool_t PStrangenessPlugin::ExecCommand(const char * command, Double_t value) {

    if (strcmp (command,"init") == 0) {
	if (!is_initialized) {
	    is_initialized=1;

	    Double_t *system_alloc_verbosity = makeStaticData()->GetBatchValue("_system_alloc_verbosity");
	    Double_t old = (*system_alloc_verbosity);
	    (*system_alloc_verbosity)=0;

	    if (makeStaticData()->AddParticle(70,"Lambda1405", 1.406) > 0) {
		makeStaticData()->AddAlias("Lambda1405","Lambda(1405)");
		
		makeStaticData()->SetParticleTotalWidth("Lambda1405",0.05);
		makeStaticData()->SetParticleBaryon("Lambda1405",1);
		
		
		makeStaticData()->AddDecay("Lambda(1405 )0 --> Sigma+ + pi-", "Lambda1405", 
					   "Sigma+,pi-", .3333 );
		makeStaticData()->AddDecay("Lambda(1405)0 --> Sigma- + pi+", "Lambda1405", 
					   "Sigma-,pi+", .3333 );
		makeStaticData()->AddDecay("Lambda(1405)0 --> Sigma0 + pi0", "Lambda1405", 
					   "Sigma0,pi0", .3333 );
	    
	    
		makeStaticData()->AddParticle(71,"Sigma1385+", 1.3828);
		makeStaticData()->AddAlias("Sigma1385+","Sigma(1385)+");
		makeStaticData()->SetParticleTotalWidth("Sigma1385+",0.0358);
		makeStaticData()->SetParticleBaryon("Sigma1385+",1);
	    
		makeStaticData()->AddDecay("Sigma(1385)+ --> Lambda + pi+", "Sigma1385+",
					   "Lambda,pi+", .8815 );    
		makeStaticData()->AddDecay("Sigma(1385)+ --> Sigma+ + pi0", "Sigma1385+",
					   "Sigma+, pi0", .05925);
		makeStaticData()->AddDecay("Sigma(1385)+ --> Sigma0 + pi+", "Sigma1385+",
					   "Sigma0, pi+", .05925);
	    
	    
	    
		makeStaticData()->AddParticle(72,"Sigma1385-", 1.3872);
		makeStaticData()->AddAlias("Sigma1385-","Sigma(1385)-");
		makeStaticData()->SetParticleTotalWidth("Sigma1385-",0.0394);
		makeStaticData()->SetParticleBaryon("Sigma1385-",1);
	    
		makeStaticData()->AddDecay("Sigma(1385)- --> Lambda + pi-", "Sigma1385-",
					   "Lambda,pi-", .8815 );    
		makeStaticData()->AddDecay("Sigma(1385)- --> Sigma- + pi0", "Sigma1385-",
					   "Sigma-, pi0", .05925);
		makeStaticData()->AddDecay("Sigma(1385)- --> Sigma0 + pi-", "Sigma1385-",
					   "Sigma0, pi-", .05925);
	


		makeStaticData()->AddParticle(73,"Sigma13850", 1.3837);
		makeStaticData()->AddAlias("Sigma13850","Sigma(1385)0");
		makeStaticData()->SetParticleTotalWidth("Sigma13850",0.0365);
		makeStaticData()->SetParticleBaryon("Sigma13850",1);
	
		makeStaticData()->AddDecay("Sigma(1385)0 --> Lambda + pi0", "Sigma13850",
					   "Lambda,pi0", .87000 );  
		makeStaticData()->AddDecay("Sigma(1385)0 --> Lambda + g", "Sigma13850",
					   "Lambda,g", .01300 );   
		makeStaticData()->AddDecay("Sigma(1385)0 --> Sigma+ + pi-", "Sigma13850",
					   "Sigma+, pi-", .05850);
		makeStaticData()->AddDecay("Sigma(1385)0 --> Sigma- + pi+", "Sigma13850",
					   "Sigma-, pi+", .05850);
	
	    } else {
		Error("ExecCommand","PIDs blocked, plugin disabled");
	    }

	    (*system_alloc_verbosity)=old; 
	    	
	}
	return kTRUE;
    }

    
    return kFALSE;
}



ClassImp(PStrangenessPlugin)



    
