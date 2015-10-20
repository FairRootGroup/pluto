/////////////////////////////////////////////////////////////////////
//  PDistributionManager Class implementation file
//
//  PDistributionManager keeps information about all distributions
//  Users can add user-defined distributions
// 
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>

#include "PDistributionManager.h"

#include "Plugins.h"
#include "TFile.h"

PDistributionManager *gDM=0;          // global pointer 

PDistributionManager& fDistributionManager()
 {
   static PDistributionManager* ans = new PDistributionManager();
   return *ans;
 }

PDistributionManager * makeDistributionManager()
 {
   return &fDistributionManager();
 }



PDistributionManager::PDistributionManager() {

    collect_pointer=0;
    makeStaticData();

  
    

    pdmutil = makeDistributionManagerUtil();
    
    pdmutil->no_warning=kTRUE;
	
    //some bugfixes .... TODO: Move this to extra plugin
    //makeStaticData()->AddDecay("eta' Dalitz", "eta'","dilepton,g",0.0009);

    std_models = new PStdModels();
    std_models->Add(pdmutil);

    pdmutil->AddGroup("plugins", "Plugins");
    SetGroup("plugins");

    batch=NULL;
    loop_start=NULL;



#include "Plugins.cc"

    pdmutil->no_warning=kFALSE; //workaround in Add()

    from_pdecaymanager =0;

    SetGroup("user");

    Info("PDistributionManager()","(%s)", PRINT_CTOR);

}

int PDistributionManager::Attach(PChannel * ch) {
    //Attach the PChannel "ch" to the distribution manager
    //This enables all distributions in the PChannel
    //and fills the PChannel with life


    //before Attaching, make the data base link
    pdmutil->LinkDB(); 



    if (!makeStaticData()->GetFreezeOut()) {
        pdmutil->no_warning=kTRUE;
        //remove warning, otherwise the Add() will cause warning
	//Add(makeStdModels()->GetModels());
	std_models->Add(pdmutil);

        pdmutil->no_warning=kFALSE;
	Info("Attach","Re-iteration of std plugin done");
    }

    return pdmutil->Attach(ch);
}


Bool_t PDistributionManager::AddPlugin(PDistributionCollection *plugin) {
  if (collect_pointer == PDISTRIBUTIONMANAGER_MAX_COLLECT) {
      Warning("AddPlugin","PDISTRIBUTIONMANAGER_MAX_COLLECT reached");
      return kFALSE;
  }
  collect[collect_pointer]=plugin;
  collect_pointer++;

  plugin->SetEnable(0);
  if (makeDistributionManagerUtil()->Add(plugin)) {
      return kFALSE;
  }

  return kTRUE;

}

Bool_t PDistributionManager::Activate(const char * name) {
    //is model existing?
    if (!GetDistribution(name)) {
	Warning("Activate","plugin %s does not exist",name);
	return kFALSE;
    }
    //first check if model is already activated
    if (GetDistribution(name)->GetEnable()) return kTRUE;
    
    for (int i=0;i<collect_pointer;i++) {
	if (strcmp(name,collect[i]->GetIdentifier()) == 0) {
	    //first we check the dependencies
	    Int_t pointer =0;
	    const char * depname=NULL;
	    while ((depname=collect[i]->GetDependency(&pointer))!=NULL) {
		Activate(depname);
	    }
	    if (collect[i]->Activate()) {
		Info("Exec","Plugin <%s> activated",name);
		collect[i]->SetEnable(1);
		SetGroup("user");
		return kTRUE;
	    }
	    Warning("Activate","plugin %s could not be activated",name);
	    return kFALSE;
	}
    }
    Warning("Activate","plugin %s does not exist",name);
    return kFALSE;

}

#if 1
Bool_t PDistributionManager::ExecAll(const char * command) {
    //sends "command" to all enabled plugins

    for (int i=0; i<collect_pointer; i++) {
	if (collect[i] -> GetEnable()) {
	    collect[i] -> ExecCommand(command,0);
	}
    }
    return kTRUE;
}
#endif

Bool_t PDistributionManager::Exec(const char * command) {
    //Executes a command like "model: command"
    //If "model" is a plugin, it is activated first

    //split at ":"
    char *array2[200];
    Int_t array2_s=200; 
    PUtils::Tokenize(command, ":", array2, &array2_s);
    
    if (array2_s > 2) {
	Warning("Exec","Syntax error: Too many :'s");
	return kFALSE;
    } else {
	//First check for a distribution
	PDistribution * dist = GetDistribution(array2[0]);
	if (!dist) {
	    Warning("Exec","model/plugin %s does not exist",array2[0]);
	    return kFALSE;
	}

	if (Activate(array2[0])) {
	    if (array2_s == 2) {
		dist = GetDistribution(array2[0]);
		if (dist->Exec(array2[1])) {
		    SetGroup("user");
		    return kTRUE;
		}
		SetGroup("user");
		return kFALSE;
	    } else  { //Consider also "empty" executions
		dist = GetDistribution(array2[0]);
		if (dist->Exec("")) {
		    SetGroup("user");
		    return kTRUE;
		}
		SetGroup("user");
		return kFALSE;
	    }
	    //	    return kTRUE;
	}
    }

    return kFALSE;

}

void PDistributionManager::PluginInfo(const char * info) {
    
    Info("PDistributionManager",info);

}


Bool_t PDistributionManager::Startup(char * command) {
    if (!batch) 
	batch = new PBatch;
    batch->AddCommand(command);
    return kTRUE;
};

Bool_t PDistributionManager::Startup(void) {
    if (!batch) return kFALSE;
    return batch->Execute();
}


Bool_t PDistributionManager::Unpack(const char * filename) {
    //Unpack the command list
    //The "_main" object is executed immediately


    //BUGBUG: Avoid that the same pack is used 2times
    TFile * current = gROOT->GetFile();
    new TFile(filename);

    //Unpack "_main"
    PCommandList*main = (PCommandList*)gROOT->FindObject("_main");

    if (main) {
	PBatch *priv = new PBatch();
	int level=0;
	char *cmd;
	while (main->GetCommand(&cmd,level)) {
	    level++;
	    if (cmd) {
		//cout << cmd << endl;
		priv->AddCommand(cmd);
	    }
	    if (level>MAX_COMMAND_POINTER) return kFALSE; //breakup
	}
	priv->Execute();
    }

    //Unpack "_startup"
    PCommandList*startup = (PCommandList*)gROOT->FindObject("_startup");

    if (startup) {
	int level=0;
	char *cmd;
	while (startup->GetCommand(&cmd,level)) {
	    level++;
	    if (cmd) {
		Startup(cmd);
	    }
	    if (level>MAX_COMMAND_POINTER) return kFALSE; //breakup
	}
    }

    //Unpack first projector for event loop
    PCommandList*loop_start_cmd = 
	(PCommandList*)gROOT->FindObject("_loop");

    if (loop_start_cmd) {
	if (loop_start) {
	    Error("Unpack","Loop command list already defined");
	    //TODO: stack?
	} else {
	    loop_start = new PProjector();
	    loop_start->SetPriority(FILTER_PRIORITY);
	    int level=0;
	    char *cmd;
	    TObject *obj=NULL;
	    while (loop_start_cmd->GetCommand(&cmd,level,&obj)) {
		level++;
		if (cmd) {
		    //cout << cmd << endl;
		    if (!obj)
			loop_start->AddCommand(cmd);
		    else if (!strncmp(obj->ClassName(),"TH1",3)) {
			Info("Unpack","Recovered TH1 <%s>",obj->GetName());
			loop_start->AddHistogram((TH1*)obj,cmd,0); //nofill flag set
		    } else if (!strncmp(obj->ClassName(),"TH2",3)) {
			Info("Unpack","Recovered TH2 <%s>",obj->GetName());
			loop_start->AddHistogram((TH2*)obj,cmd,0); //nofill flag set
		    } else if (!strncmp(obj->ClassName(),"TH3",3)) {
			Info("Unpack","Recovered TH3 <%s>",obj->GetName());
			loop_start->AddHistogram((TH3*)obj,cmd,0); //nofill flag set
		    } else {
			Warning("Unpack","Object <%s> is of type <%s>. Unsupported and ignored",
				obj->GetName(), obj->ClassName());
			loop_start->AddCommand(cmd);
		    }

		}
		if (level>MAX_COMMAND_POINTER) return kFALSE; //break
	    }
	    
	}

    }


    if (current)
	current->cd();

    return kTRUE;

}

ClassImp(PDistributionManager)
