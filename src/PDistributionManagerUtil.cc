/////////////////////////////////////////////////////////////////////
//  PDistributionManagerUtil Class implementation file
//
//  PDistributionManagerUtil is the private implementation of 
//  the PDistributionManager
// 
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>

#include "PDistributionManagerUtil.h"

PDistributionManagerUtil& fDistributionManagerUtil()
 {
   static PDistributionManagerUtil* ans = new PDistributionManagerUtil();
   return *ans;
 }

PDistributionManagerUtil * makeDistributionManagerUtil()
 {
   return &fDistributionManagerUtil();
 }

PDistributionManagerUtil::PDistributionManagerUtil() {


    makeStaticData();
    position=0;
    verbosity = 1;

    group_position=0;
    corr_position=0;
    AddGroup("none", "No Group");
    AddGroup("user", "User defined distributions");


    for (int i=0;i<MAX_DISTRIBUTION_LIST;i++)
	alt_distribution[i]=-1;
    linkdb_done=0;


    SetGroup("user");

}

void PDistributionManagerUtil::AlternativeTo(const char * a, const char *b) {
    //Makes the distributions "a" and "b" strongly alternative
    //This means only one of these 2 (or more) can be enabled at the same time

    Int_t apos=-1, bpos=-1;

    for (int i=0; i<position; i++) {
	if (strcmp(distribution[i]->GetIdentifier(), a)==0) {
	    apos=i;
	}
	if (strcmp(distribution[i]->GetIdentifier(), b)==0) {
	    bpos=i;
	}
    }

    if (apos==-1) {
	Warning("AlternativeTo","First argument %s not found",a);
	return;
    }
    if (bpos==-1) {
	Warning("AlternativeTo","Second argument %s not found",b);
	return;
    }
    if (alt_distribution[bpos] == -1) {
	alt_distribution[bpos] = apos;
	alt_distribution[apos] = bpos;
    } else if (alt_distribution[apos] == -1) {
	alt_distribution[bpos] = apos;
	alt_distribution[apos] = bpos;
    } else {
	Int_t bold=alt_distribution[bpos];	
	alt_distribution[bpos] = apos;
	alt_distribution[apos] = bold;
    }
    int mylist=alt_distribution[bpos];
    while (mylist != bpos) {

	if (distribution[mylist]->GetEnable()) //BUGBUG: Warning needed
	    distribution[mylist]->SetEnable(0);
	mylist = alt_distribution[mylist];
	if (mylist == -1) mylist=bpos;
    }
}

int PDistributionManagerUtil::Add(PDistribution * dist) {
    //Adds a distribution to the list 
    //return value -1 on failure

    Int_t num_channels;
    
    if (makeDynamicData()->GetPChannels(&num_channels) && (no_warning==kFALSE)) {
        Warning("Add","You add a distribution (%s) after PChannels have been created",dist->GetIdentifier());
        Warning("Add","Please add the PDistribution in the very beginning");
    }

    if (position == MAX_DISTRIBUTION_LIST) {
	Warning("Add","MAX_DISTRIBUTION_LIST reached");
	return -1;
    }

    //check if my identifier already exists...

    for (int i=0; i<position; i++) {
	if (strcmp(distribution[i]->GetIdentifier(), dist->GetIdentifier())==0) {
	    Warning("Add","Identifier %s already exists",distribution[i]->GetIdentifier());
	    return -1;
	}

    }

    distribution[position]=dist;
    alt_distribution[position]=-1;
    AddCorrelation(current_group,position);    
    position++;

    //If key exists already, make it alternative
    for (int i=0; i<position-1; i++) {
	if ((distribution[i]->GetKey() == dist->GetKey()) && (dist->GetKey()>=0))
	    AlternativeTo(distribution[i]->GetIdentifier(), dist->GetIdentifier());
	}

    return 0;

}

int  PDistributionManagerUtil::Add(PDistribution * dist, const Char_t *gr) {
    
    Int_t save=current_group;
    if (GetGroup(gr)<0) return -1;
    current_group = GetGroup(gr);

    if (Add(dist) < 0) {
	current_group = save;
	return -1;
    }
    current_group = save;
    return 0; 

}


int PDistributionManagerUtil::Add(TObjArray * arr) {
    //Adds a complete list of distributions

    for (int pat = 0; pat < arr->GetEntriesFast(); pat++) {
	const char * group = ((PDistribution *) (*arr)[pat])->GetGroupID();
	if (group) SetGroup(group);
        Bool_t saved_warning=no_warning;
        no_warning=kTRUE;
	if (Add((PDistribution *) (*arr)[pat]) == -1) {
	    Warning("Add(TObjArray *)","Adding %s failed ",((PDistribution *) (*arr)[pat])->GetDescription());
	}
        no_warning=saved_warning;
    }
    return 0;
}

int PDistributionManagerUtil::Attach(PChannel * ch) {
    //Attach the PChannel "ch" to the distribution manager
    //This enables all distributions in the PChannel
    //and fills the PChannel with life


    if ( ch->GetQuasi() ) {
	Attach(ch->GetQuasi()); //Attach also dummy quasi-free scattering
    }
    
    for (int i=0; i<position; i++) {
	if (distribution[i] -> GetEnable()) {
	    if (ch->SetDistribution(distribution[i]) == 0) {
//		distribution[i] = (PDistribution*) distribution[i]->Clone();
		distribution[i]->Reset();
	    }
	}
    }
    return 0;
}


void PDistributionManagerUtil::DisableAlts(int id) {
    if (alt_distribution[id] != -1 ) {
	int mylist=alt_distribution[id];
	while (mylist != id) {
	    if (distribution[mylist]->GetEnable()) distribution[mylist]->SetEnable(0);
	    mylist = alt_distribution[mylist];
	    if (mylist == -1) mylist=id;
	}
    }
}


Int_t PDistributionManagerUtil::PrintGroup(Int_t group_id, Int_t width,Int_t indent,
				       const Char_t* name,
				       Int_t *num_enabled_mods,Int_t *num_total_mods,
				       Int_t *num_subs, Int_t *will_print) const {
    //Print a group (group_id) with all sub-groups
    //if width=0 make a "dry" test
    //width>0 :column width
    //returns the new column width (if width=0)
    //updates num_enabled_mods and num_subs

    //First, we loop over the groups and try to find
    //the sub-groups

    Int_t return_width=0;
    Int_t dummy_will_print=0,local_will_print=0;

    //In a first step we find out how many valid subgroups and models
    //we have

    Int_t local_num_enabled_mods=0,local_num_total_mods=0,local_num_subs=0;

    for (int gr=0; gr<group_position; gr++) {
	if (group_corr[gr] == group_id) {
	    Int_t local_width = PrintGroup(gr,0,indent+2,name,&local_num_enabled_mods,
					   &local_num_total_mods, &local_num_subs, 
					   &dummy_will_print);
	    if (local_width>return_width) return_width=local_width;
	    if (dummy_will_print) local_will_print=1;
	}
    }    

    for (int corr=0; corr<corr_position; corr++) {
	if (corr_gr[corr] == group_id) { //has found correlation
	    local_num_total_mods++;
	    if (distribution[corr_dis[corr]]->GetEnable())
		local_num_enabled_mods++;
	    Bool_t matched=kFALSE;
	    if (name)
		if (strcmp(name,group_identifier[group_id]) ==0) matched=kTRUE;
	    if (group_expanded[group_id] ||
		matched || local_will_print) { //Print List
		Int_t local_width = strlen(distribution[corr_dis[corr]]->GetIdentifier())+indent;
		if (local_width>return_width) return_width=local_width;
		if (alt_distribution[corr_dis[corr]] != -1 ) {

		    //have to take into account correlations
		    int mylist=alt_distribution[corr_dis[corr]];
		    while (mylist != corr_dis[corr]) {
			local_width = strlen(distribution[mylist]->GetIdentifier())+indent;
			if (local_width>return_width) return_width=local_width;
			mylist = alt_distribution[mylist];
			if (mylist == -1) mylist=corr_dis[corr];
		    }
		}

	    }
	}
    }

    //for the return statistics:
    if (local_num_total_mods && (width==0)) {
	*num_enabled_mods += local_num_enabled_mods;
	*num_total_mods += local_num_total_mods;
	*num_subs         += local_num_subs+1;
	Int_t local_width = strlen(group_identifier[group_id])+indent;
	if (local_width>return_width) return_width=local_width;

	if (name)
	    if ((strcmp(name,group_identifier[group_id]) ==0)  || (group_expanded[group_id]) 
		|| local_will_print) {
		*will_print =1;
		return return_width;
	    }
	*will_print =0;
	return local_width+2;
    }



    if (!local_num_total_mods) return 0; //nothing to print
    
    
    
    //Print group header
    cout << "    ";
    for (Int_t  k=0;k<indent;k++) cout << " ";
    cout << group_identifier[group_id];

    if ((width-indent-strlen(group_identifier[group_id]))>0) {
	for (UInt_t  k=0;k<((width-indent-strlen(group_identifier[group_id])));k++)
	    cout << " ";
    }
    cout << group_description[group_id] << " ";
    cout << ": "<<local_num_enabled_mods << " objects (of " << local_num_total_mods << ")";
    if (local_num_subs)
	cout << ", subgroups: " << local_num_subs << endl;
    else cout << endl;
    
    if (name != NULL)

	if ((strcmp(name,group_identifier[group_id]) !=0)  && (group_expanded[group_id] ==0) 
	    && (local_will_print ==0)) {
	    return 0; //Unexpanded group
	}

    //first the sub-groups
    for (int gr=0; gr<group_position; gr++) {
	if (group_corr[gr] == group_id) {
	    PrintGroup(gr,width,indent+2,name,&local_num_enabled_mods,&local_num_total_mods, 
		       &local_num_subs, 
		       &dummy_will_print);
	}
    }    

    //Now we loop over the distributions

    for (int corr=0; corr<corr_position; corr++) {
	if (corr_gr[corr] == group_id) { //has found correlation
	    if (distribution[corr_dis[corr]]->GetEnable()) (*num_enabled_mods)++;

	    Bool_t matched=kFALSE;
	    if (name)
		if (strcmp(name,group_identifier[group_id]) ==0) matched=kTRUE;
	    if (group_expanded[group_id] ||
		matched) { //Print List
		
		if (alt_distribution[corr_dis[corr]] == -1 ) {
		    if (distribution[corr_dis[corr]]->GetEnable())
			cout << "[X] ";
		    else
			cout << "[ ] ";
		} else {
		    if (distribution[corr_dis[corr]]->GetEnable())
			cout << "(X) ";
		    else
			cout << "( ) ";
		}
		for (Int_t  k=0;k<indent;k++) cout << " ";
		cout << distribution[corr_dis[corr]]->GetIdentifier();

		if ((width-indent-strlen(distribution[corr_dis[corr]]->GetIdentifier()))>0) {
		    for (UInt_t k=0;k<((width-indent-strlen(distribution[corr_dis[corr]]->GetIdentifier())));k++)
			cout << " ";
		}
//		cout << (width-indent-strlen(distribution[corr_dis[corr]]->GetIdentifier())) << endl;
		cout << distribution[corr_dis[corr]]->GetDescription() << endl;;
		if (verbosity == 2)  { //Print out parameters etc... for the given Model
		    distribution[corr_dis[corr]]->Print();
		}		
		if (alt_distribution[corr_dis[corr]] != -1 ) {
		    int mylist=alt_distribution[corr_dis[corr]];
		    while (mylist != corr_dis[corr]) {
			if (distribution[mylist]->GetEnable()) {
			    cout << " (X)";
			    for (Int_t  k=0;k<indent;k++) cout << " ";
			    cout << distribution[mylist]->GetIdentifier();

			    if ((width-indent-strlen(distribution[mylist]->GetIdentifier()))>0) {
				for (UInt_t k=0;k<((width-indent-strlen(distribution[mylist]->GetIdentifier())));k++) cout << " ";
				cout << distribution[mylist]->GetDescription() << endl;
			    }
			}
			else {
			    cout << " ( )";
			    for (Int_t  k=0;k<indent;k++) cout << " ";
			    cout << distribution[mylist]->GetIdentifier();

			    if ((width-indent-strlen(distribution[mylist]->GetIdentifier()))>0) {
				for (UInt_t k=0;k<((width-indent-strlen(distribution[mylist]->GetIdentifier())));k++) 
				    cout << " ";
			    }
			    cout << distribution[mylist]->GetDescription() << endl;
			}
			mylist = alt_distribution[mylist];
			if (mylist == -1) mylist=corr_dis[corr];
		    }
		}
		
	    } //end print list
	}
	
    }

    return 0;

}


void PDistributionManagerUtil::Print(const Option_t* delme) const {

    UInt_t max_id_length =0;

    for (int i=0; i<group_position; i++) {
	if (strlen(group_identifier[i])>max_id_length)
	    max_id_length=strlen(group_identifier[i]);
    }
    for (int i=0; i<position; i++) {
	if (strlen(distribution[i]->GetIdentifier())>max_id_length)
	    max_id_length=strlen(distribution[i]->GetIdentifier());
//	cout << alt_distribution[i] << endl;
    }

    cout << "--------------------" << endl;
    cout << "PDistributionManager" << endl;
    cout << "--------------------" << endl;

    Int_t max_width=0;
    Int_t a,b,c,d;
    for (int gr=0; gr<group_position; gr++) {
	//do not print if I'm a sub-group
	if (group_corr[gr] <0 ) {
	    Int_t local_width=PrintGroup(gr, 0,0,
					 (char*)delme,
					 &a,&b,&c,&d);
	    if (local_width>max_width) max_width=local_width;
	    //cout << "max_width:" << max_width << endl;
	}
    }

    for (int gr=0; gr<group_position; gr++) {
//	cout << max_width << endl;
	//return;
	if (max_width && (group_corr[gr] <0 ))
	    PrintGroup(gr,max_width+2 ,0,
		   (char*)delme,
		   &a,&b,&c,&d);
    }
    
//     if (!group_expanded[gr] && (print == 0) && (sum>0)) { //Print Summary
// 	print=1;
// 	cout << endl << "    ";
// 	cout << group_identifier[gr];
// 	for (UInt_t  k=0;k<((max_id_length-strlen(group_identifier[gr])+2));k++)
// 	    cout << " ";
// 	cout << group_description[gr] << ": " << enable << " enabled (out of " << sum << ")" << endl;
//     }
// } 
//end loop over groups

    cout << "--------------------" << endl;

    for (int i=0; i<position; i++) {
	if (delme)
	    if (strcmp(distribution[i]->GetIdentifier(),delme)==0)
		distribution[i]->Print();
    }

}


Bool_t PDistributionManagerUtil::Disable(const Char_t * id) {
    //Disable the distribution "id"
    //If "id" is a group, disable all distributions
    //which are members

    Bool_t retval = kFALSE;

    for (int i=0; i<position; i++) {
	if (strcmp(id, distribution[i]->GetIdentifier()) == 0) {
	    retval = kTRUE;
	    if (distribution[i]->GetEnable())
		distribution[i]->SetEnable(0);
	    
	}
    }
    //look for groups
    if (!retval) {
	Int_t mygr=GetGroup(id,0);
	if (mygr>=0) {
	    for (int corr=0; corr<corr_position; corr++) {
		if (strcmp(id, group_identifier[corr_gr[corr]]) == 0) {
		    retval = kTRUE;
		    distribution[corr_dis[corr]]->SetEnable(0);
		}
	    }
	    for (int i=0; i<group_position; i++) {
		if (group_corr[i]==mygr ) {
		    //i has id as parent
		    if (Disable(group_identifier[i])) 
			retval = kTRUE;
		}
	    }
	}
    }
    
    if (!retval) Warning("Disable","%s not found",id);
    return retval;
}

Bool_t PDistributionManagerUtil::Enable(const Char_t * id) {
    //Enable the distribution "id"
    //If "id" is a group, enable all distributions
    //which are members
    //If the distribution is part of an alternatice chain, disable all 
    //other distributions

    Bool_t retval = kFALSE;

    for (int i=0; i<position; i++) {
	if  (strcmp(id, distribution[i]->GetIdentifier()) == 0) {
	    retval = kTRUE;
	    if ((distribution[i]->GetEnable()==0)) {
		distribution[i]->SetEnable(1);
		DisableAlts(i);
	    }
	}
    }
    //look for groups
    if (!retval) {
	Int_t mygr=GetGroup(id,0);
	if (mygr>=0) {
	    for (int corr=0; corr<corr_position; corr++) {
		if (strcmp(id, group_identifier[corr_gr[corr]]) == 0) {
		    retval = kTRUE;
		    distribution[corr_dis[corr]]->SetEnable(1);
		    DisableAlts(corr_dis[corr]);
		}
	    }
	    for (int i=0; i<group_position; i++) {
		if (group_corr[i]==mygr ) {
		    //i has id as parent
		    if (Enable(group_identifier[i])) 
			retval = kTRUE;
		}
	    }
	}
    }


    if (!retval) Warning("Enable","%s not found",id);
    return retval;
}


PDistribution * PDistributionManagerUtil::GetDistribution(const Char_t * id) {
    for (int i=0; i<position; i++) {
//	if (distribution[i]->GetEnable() && (strcmp(id, distribution[i]->GetIdentifier()) == 0)) {
	if (strcmp(id, distribution[i]->GetIdentifier()) == 0) {
	    return distribution[i];
	}
    }
    return NULL;
}

void PDistributionManagerUtil::LinkDB(void) {
    //Links the list of distributions to the PDataBase
    //This makes coupled-channel calculations possible

    for (int i=0; i<position; i++) 
	if (distribution[i]->GetKey()>=0 && distribution[i]->GetEnable()) {

 	    if (makeDynamicData()->GetDecayModelByKey(distribution[i]->GetKey()) && (linkdb_done==0)) {
		//Print warning in first loop
		Warning("LinkDB","'%s' would overwrite'%s' ",distribution[i]->GetDescription(),
			makeDynamicData()->GetDecayModelByKey(distribution[i]->GetKey())->GetDescription());
 	    } else {

	      //	      if (!distribution[i]->LinkDBdone()) {
		makeDynamicData()->SetDecayModelByKey(distribution[i]->GetKey(), 
						      (PChannelModel *)distribution[i] );
		//		distribution[i]->LinkDBdone(1);
		//}
	    }
	}
    
    for (int i=0; i<position; i++) 
	if (distribution[i]->GetKey()>=0 && distribution[i]->GetEnable()) {
	    distribution[i]->FreezeOut();
	}
    
    linkdb_done=1;
    
}

//BUGBUG: Check if 2 keys are used

int PDistributionManagerUtil::AddGroup(const Char_t *id, const Char_t *de) {
    if ((GetGroup(id,0)) >= 0)
	return -1;

    if (group_position == MAX_GROUP_LIST) {
	Warning("AddGroup","MAX_GROUP_LIST reached");
	return -1;
    }
    group_identifier[group_position]=id;
    group_description[group_position]=de;
    group_expanded[group_position]=0;
    group_corr[group_position]=-1;
    group_position++;
    return 0;
};

int PDistributionManagerUtil::AddSubGroup(const Char_t *id,const  Char_t *de, const Char_t *parent_group) {
    if (AddGroup(id, de) <0) return -1;
    Int_t mygr=GetGroup(id);
    Int_t par=GetGroup(parent_group);
    if (par<0) return -1;
    group_corr[mygr]=par;
    
     return 0;
}

Int_t PDistributionManagerUtil::AddCorrelation(Int_t gr, Int_t dis) {
    if (corr_position == MAX_CORR_LIST) {
	Warning("AddCorrelation","MAX_CORR_LIST reached");
	return -1;
    }
    corr_gr[corr_position]=gr;
    corr_dis[corr_position]=dis;
    corr_position++;
    return 0;
}

int PDistributionManagerUtil::GetGroup(const Char_t *id, Int_t warning) {
     for (int i=0; i<group_position; i++) {
	 if (strcmp(id, group_identifier[i]) == 0)
	     return i;
     }
     if (warning) Warning("GetGroup","id %s not found",id);
     return -1;
};

void PDistributionManagerUtil::ExpandGroup(const char *id, Int_t ex) {
    Int_t pos=GetGroup(id);
    if (pos<0) return;
    group_expanded[pos]=ex;
}

void PDistributionManagerUtil::SetGroup(const Char_t *id) {
    current_group = GetGroup(id);
};

ClassImp(PDistributionManagerUtil)
