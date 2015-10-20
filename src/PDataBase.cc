////////////////////////////////////////////////////////
//  Particle ID, properties and decay data base
//
//  
//
//                    Author: I. Froehlich
//                    Written: 11.04.2007
//                    Revised: 
//
////////////////////////////////////////////////////////


#include "PDataBase.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "PMesh.h"

PDataBase& fDataBase()
{
    static PDataBase* ans = new PDataBase();
    return *ans;
}

PDataBase * makeDataBase()
{
    return &fDataBase();
}


PDataBase::PDataBase() {
    param_double_pointer=0;
    param_string_pointer=0;
    param_int_pointer=0;
    param_tobj_pointer=0;
    lastkey=0;

    //Init arrays to NULL
    for (int i=0; i<PDATABASE_MAX_LINES; i++) {
	for (int j=0; j<PDATABASE_MAX_STRING_PARAM; j++)
	    strings[i][j]=NULL;
	for (int j=0; j<PDATABASE_MAX_DOUBLE_PARAM; j++)
	    doubles[i][j]=NULL;
	for (int j=0; j<PDATABASE_MAX_INT_PARAM; j++)
	    ints[i][j]=NULL;
	for (int j=0; j<PDATABASE_MAX_TOBJ_PARAM; j++)
	    tobjs[i][j]=NULL;
    }

    for (int j=0; j<PDATABASE_MAX_INT_PARAM; j++) {
	param_int_key[PDATABASE_MAX_INT_PARAM]=NULL;
	param_int_key_max[PDATABASE_MAX_INT_PARAM]=0;
    }

    MakeParamString("name", "Database name"); 
//this is absolute minimum and required for testing for existence
    Info("PDataBase()","(%s)", PRINT_CTOR);
}


void    PDataBase::Print(const Option_t* delme) const {
    cout << param_int_pointer << " INT's booked (out of " << PDATABASE_MAX_INT_PARAM <<")" << endl;
    for (int i=0;i<param_int_pointer;i++)
	cout << "INT:" << param_int_name[i] << ":" << param_int_descr[i] << endl;

    cout << param_double_pointer << " DOUBLE's booked (out of " << PDATABASE_MAX_DOUBLE_PARAM <<")" << endl;
    for (int i=0;i<param_double_pointer;i++)
	cout << "DOUBLE:" << param_double_name[i] << ":" << param_double_descr[i] << endl;

    cout << param_string_pointer << " STRING's booked (out of " << PDATABASE_MAX_STRING_PARAM <<")" << endl;
    for (int i=0;i<param_string_pointer;i++)
	cout << "STRING:" << param_string_name[i] << ":" << param_string_descr[i] << endl;

    cout << param_tobj_pointer << " TOBJ's booked (out of " << PDATABASE_MAX_TOBJ_PARAM <<")" << endl;
    for (int i=0;i<param_tobj_pointer;i++)
	cout << "TOBJ:" << param_tobj_name[i] << ":" << param_tobj_descr[i] << endl;
    
}

void    PDataBase::SetFastKey(Int_t pkey, Int_t maxkey) {
    //set if param(pkey) provides a unique fast key
    //this is the case for e.g. "pid"
    //set the maxkey to max "pid" value
    param_int_key[pkey] = new Int_t[maxkey] ;
    param_int_key_max[pkey]=maxkey;
    for (int i=0;i<maxkey;i++) param_int_key[pkey][i]=-1;
}

Bool_t PDataBase::CheckEntry(Int_t key){
    //returns a kTRUE if a line with the "key" is already initialized
    if ((key >= lastkey) || (key<-1)) return kFALSE;
//    Int_t pp=getParamString("name");
    Int_t pp=0;
    if (strings[key][pp] == NULL) return kFALSE;
    return kTRUE;

};

const char* PDataBase::GetName(Int_t key){
    //returns the primary name
    if ((key >= lastkey) || (key<-1)) return NULL;
//    Int_t pp=getParamString("name");
    Int_t pp=0;
    return strings[key][pp];
};

Int_t PDataBase::MakeParamDouble(const char* paramname, const char *descr) {
    //check maximum size
    if (param_double_pointer==PDATABASE_MAX_DOUBLE_PARAM) {
	Fatal("MakeParamDouble","PDATABASE_MAX_DOUBLE_PARAM reached");
    }
    //check if paramname already exists
    if (GetParamDouble(paramname) >= 0) {
	Warning("MakeParamDouble","Value %s already exists",paramname);
	return -1;
    }
    param_double_name[param_double_pointer]=paramname;
    param_double_descr[param_double_pointer]=descr;
    param_double_pointer++;
    return  param_double_pointer-1;
}


Int_t PDataBase ::MakeParamString(const char* paramname, const char *descr) {
    //check maximum size
    if (param_string_pointer==PDATABASE_MAX_STRING_PARAM) {
	Fatal("MakeParamString","PDATABASE_MAX_STRING_PARAM reached");
    }
    //check if paramname already exists
    if (GetParamString(paramname) >= 0) {
	Warning("MakeParamString","Value %s already exists",paramname);
	return -1;
    }
    param_string_name[param_string_pointer]=paramname;
    param_string_descr[param_string_pointer]=descr;
    param_string_pointer++;
    return  param_string_pointer-1;
}

Int_t PDataBase ::MakeParamInt(const char* paramname, const char *descr) {
    //check maximum size
    if (param_int_pointer==PDATABASE_MAX_INT_PARAM) {
	Fatal("MakeParamInt","PDATABASE_MAX_INT_PARAM reached");
    }
    //check if paramname already exists
    if (GetParamInt(paramname) >= 0) {
	Warning("MakeParamInt","Value %s already exists",paramname);
	return -1;
    }
    param_int_name[param_int_pointer]=paramname;
    param_int_descr[param_int_pointer]=descr;
    param_int_pointer++;
    return  param_int_pointer-1;
}

Int_t PDataBase ::MakeParamTObj(const char* paramname, const char *descr) {
    //check maximum size
    if (param_tobj_pointer==PDATABASE_MAX_TOBJ_PARAM) {
	Fatal("MakeParamTObj","PDATABASE_MAX_TOBJ_PARAM reached");
    }
    //check if paramname already exists
    if (GetParamTObj(paramname) >= 0) {
	Warning("MakeParamTObj","Value %s already exists",paramname);
	return -1;
    }
    param_tobj_name[param_tobj_pointer]=paramname;
    param_tobj_descr[param_tobj_pointer]=descr;
    param_tobj_pointer++;
    return  param_tobj_pointer-1;
}

Int_t PDataBase ::GetParamDouble(const char* paramname) {
    for (int i=0; i<param_double_pointer; i++)
	if (strcmp(paramname,param_double_name[i])==0) return i;
    return -1;
}

Int_t PDataBase ::GetParamString(const char* paramname) {
    //BUGBUG: If paramname is a pointer like "d1:name" this does not work->write getDescription
    for (int i=0; i<param_string_pointer; i++)
	if (strcmp(paramname,param_string_name[i])==0) return i;
    return -1;
}

Int_t PDataBase ::GetParamInt(const char* paramname, Int_t length) {
    for (int i=0; i<param_int_pointer; i++) {

	if (length<0) {
	    if (strcmp(paramname,param_int_name[i])==0) return i;
	} else {
	    if (strncmp(paramname,param_int_name[i],length)==0) return i;
	}
    }
    return -1;    
}

Int_t PDataBase ::GetParamTObj(const char* paramname) {
    for (int i=0; i<param_tobj_pointer; i++)
	if (strcmp(paramname,param_tobj_name[i])==0) return i;
    return -1;
}

Int_t  PDataBase ::ConvertParamKey(const char * &newparamname, Int_t key) {
    //checks for the ":" delmiter, sets the pointer to the remaing string
    //evaluate the new key
    //cout << "newparamname init " <<newparamname << endl;
    Int_t *newkey;
    for (unsigned int i=0;i<strlen(newparamname);i++) {
	if (newparamname[i] == ':') {
	    if (!GetParamInt(key,newparamname,&newkey,i)) {
		Warning("ConvertParamKey","Value %s for key %i not existing",newparamname,key);
		return -1;
	    }
	    newparamname = newparamname + i + 1;

	    return *newkey;
	}	
    }
    return -1;
};

char*   PDataBase ::GetDescription(const char* paramname) {
    //Interpretation of pattern
    TString spattern(paramname);
    TObjArray *array = spattern.Tokenize(TString(":"));
    TString bla;
    Int_t done=0;
    for (int pat = 0; pat < array->GetEntriesFast(); pat++) {
	if (done)  bla.Append("->");
	TObjString *patoption = (TObjString *) (*array)[pat];
	char * options=(char*)patoption->GetString().Data();
	if (GetParamString(options)>-1) {
	    bla.Append(param_string_descr[GetParamString(options)]);
	    done=1;
	}
	if (GetParamInt(options)>-1) {
	    bla.Append(param_int_descr[GetParamInt(options)]);
	    done=1;
	}
	if (GetParamDouble(options)>-1) {
	    bla.Append(param_double_descr[GetParamDouble(options)]);
	    done=1;
	}
	if (GetParamTObj(options)>-1) {
	    bla.Append(param_tobj_descr[GetParamTObj(options)]);
	    done=1;
	}
    }
    return (char*) bla.Data();
}

Bool_t  PDataBase ::GetParamDouble (Int_t key, const char*  paramname, Double_t ** result){
    //return kFALSE if 1.) key not existing or 2.) Param not existing or 3.) Param not used for key
    if (!CheckEntry(key)) return kFALSE;
    Int_t newkey=-1;
    const char *newparamname=paramname;

    if ((newkey = ConvertParamKey(newparamname,key))>=0) {
	return GetParamDouble (newkey,newparamname,result);
    }
    paramname=newparamname;
    Int_t pp=GetParamDouble(paramname);
    if (pp<0) return kFALSE;

    if (doubles[key][pp] == NULL) return kFALSE;

    *result = doubles[key][pp];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamDouble (Int_t key, Int_t pkey, Double_t ** result) {
    if (!CheckEntry(key)) return kFALSE;
    if (pkey<0) return kFALSE;
    if (doubles[key][pkey] == NULL) return kFALSE;
    *result = doubles[key][pkey];
    return kTRUE;
};


Bool_t  PDataBase ::GetParamString (Int_t key, const char*  paramname, const char ** result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) return kFALSE;
    Int_t newkey=-1;
    const char *newparamname=paramname;

    if ((newkey = ConvertParamKey(newparamname,key))>=0) {
	return GetParamString (newkey,newparamname,result);
    }
    paramname=newparamname;
    Int_t pp=GetParamString(paramname);
    if (pp<0) return kFALSE;

    if (strings[key][pp] == NULL) return kFALSE;

    *result = strings[key][pp];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamString (Int_t key, Int_t pkey, const char ** result) {
    if (!CheckEntry(key)) return kFALSE;
    if (pkey<0) return kFALSE;
    if (strings[key][pkey] == NULL) return kFALSE;
    *result = strings[key][pkey];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamInt (Int_t key, const char*  paramname, Int_t ** result, Int_t length){
    //return kFALSE if 1.) key not existing or 2.) Param not existing or 3.) Param not used for key
    if (!CheckEntry(key)) return kFALSE;

    Int_t newkey=-1;
    const char *newparamname=paramname;
    if (length<0) {
	if ((newkey = ConvertParamKey(newparamname,key))>=0) {

	    return GetParamInt (newkey,newparamname,result);
	}
    }
    paramname=newparamname;
    Int_t pp=GetParamInt(paramname,length);
    if (pp<0) return kFALSE;
    
    if (ints[key][pp] == NULL) return kFALSE;

    *result = ints[key][pp];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamInt (Int_t key, Int_t pkey, Int_t ** result) {
    if (!CheckEntry(key)) return kFALSE;

    if (pkey<0) return kFALSE;
    if (ints[key][pkey] == NULL) return kFALSE;
    *result = ints[key][pkey];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamTObj (Int_t key, const char*  paramname, TObject ** result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) return kFALSE;
    Int_t newkey=-1;
    const char *newparamname=paramname;

    if ((newkey = ConvertParamKey(newparamname,key))>=0) {
	return GetParamTObj (newkey,newparamname,result);
    }
    paramname=newparamname;
    Int_t pp=GetParamTObj(paramname);
    if (pp<0) return kFALSE;

    if (tobjs[key][pp] == NULL) return kFALSE;

    *result = tobjs[key][pp];
    return kTRUE;
};

Bool_t  PDataBase ::GetParamTObj (Int_t key, Int_t pkey, TObject ** result) {
    if (!CheckEntry(key)) return kFALSE;
    if (pkey<0) return kFALSE;
    if (tobjs[key][pkey] == NULL) return kFALSE;
    *result = tobjs[key][pkey];
    return kTRUE;
};



Bool_t  PDataBase ::SetParamDouble (Int_t key, const char*  paramname, Double_t * result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) {
	cout << "PDataBase::SetParamDouble: key " << key << " not existing" << endl;
	return kFALSE;
    }
    
    Int_t pp=GetParamDouble(paramname);
    if (pp<0) {
	cout << "PDataBase::SetParamDouble: paramname " << paramname << " not existing" << endl;
	return kFALSE;
    }

    doubles[key][pp] = result;
    return kTRUE;
};

Bool_t  PDataBase ::SetParamString (Int_t key, const char*  paramname, char * result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) {
	Error("SetParamString","key %i not existing",key);
	return kFALSE;
    }
    
    Int_t pp=GetParamString(paramname);
    if (pp<0) {
	Error("SetParamString","paramname %s not existing",paramname);
	return kFALSE;
    }

    strings[key][pp]= result;
    return kTRUE;
};

Bool_t  PDataBase ::SetParamInt (Int_t key, const char*  paramname, Int_t * result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) {
	Error("SetParamInt","key %i not existing",key);
	return kFALSE;
    }
    
    Int_t pp=GetParamInt(paramname);
    if (pp<0) {
	Error("SetParamInt","paramname %s not existing",paramname);
	return kFALSE;
    }

    ints[key][pp] = result;
    return kTRUE;
};

Bool_t  PDataBase ::SetParamTObj (Int_t key, const char*  paramname, TObject * result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) {
	Error("SetParamTObj","key %i not existing",key);
	return kFALSE;
    }
    
    Int_t pp=GetParamTObj(paramname);
    if (pp<0) {
	Error("SetParamTObj","paramname %s not existing",paramname);
	return kFALSE;
    }

    tobjs[key][pp] = result;
    return kTRUE;
};

Bool_t  PDataBase ::SetParamTObj (Int_t key, Int_t pp, TObject * result){
    //return kFALSE if 1.) KEY not existing or 2.) Param not existing or 3.) Param not used for KEY
    if (!CheckEntry(key)) {
	Error("SetParamTObj","key %i not existing",key);
	return kFALSE;
    }
    
    tobjs[key][pp] = result;
    return kTRUE;
};

Int_t PDataBase ::GetEntry(const char* name){
    for (int i=0; i<lastkey; i++) {
	if (strings[i][0])
	    if (strcmp(name,strings[i][0])==0) return i;
    }
    return -1;	
}; 


Int_t   PDataBase ::GetEntryInt(const char* paramname, Int_t value){
    Int_t * result;
    for (int i=0; i<lastkey; i++) {
	if (GetParamInt (i, paramname,&result)) {
	    if (*result == value)
		return i;
	}
    }
    return -1;
}; 

Int_t   PDataBase ::GetEntryInt(Int_t pkey, Int_t value){
//return the index key if param(pkey) is matching the value
//otherwise -1
    Int_t * result;

    //use fast checking is available
    if (param_int_key[pkey]) {
	if ((value>=0) && (value<param_int_key_max[pkey])) {
	    if (param_int_key[pkey][value]!=-1) {
		return param_int_key[pkey][value];
	    } else {
		for (int i=0; i<lastkey; i++) {
		    if (GetParamInt (i, pkey,&result)) {
			if (*result == value) {
			    //cout << "init " << pkey << ":" << value << endl;
			    param_int_key[pkey][value] = i;
			    return i;
			}
		    }
		}
	    }
	}
    }
    //slow version: loop over entries

    for (int i=0; i<lastkey; i++) {
	if (GetParamInt (i, pkey,&result)) {
	    if (*result == value)
		return i;
	}
    }
    return -1;
}; 

Bool_t  PDataBase ::AddEntry(Int_t key, const char * name){
    if (CheckEntry(key)) {
	Warning("AddEntry","Key %i already existing",key);
	return kFALSE;
    }
    if (GetEntry(name) >=0) {
	Warning("AddEntry","An entry with name %s already exists",name);
	CRASH;
	return kFALSE;
    }
    Int_t pp=GetParamString("name");
    strings[key][pp] = name;
    return kTRUE;
};

Int_t   PDataBase ::AddEntry(const char * name) {
    if (AddEntry(lastkey, name)) {
	lastkey++;
	return lastkey-1;
    }
    return -1;
};

Int_t   PDataBase::AddListEntry(const char * name, const char * count, 
				const char * link, const char * newname){
    //This is used to add linked-lists the the entry "name"
    //The "count" (which should be an int param) holds the number of links
    //"link" is used for:
    //1.) setting the key to the last entry of the lists (updated on the fly)
    //2.) setting the key from the last list entry to the last-but-one
    //The first list entry points back to the original "name", thus forming a ring
    //"newname" is the name of the new list entry

    Int_t *i_count;
    Int_t *i_link;
    Int_t *i_newlink;
    Int_t key=GetEntry(name);
    Int_t targetkey=GetEntry(newname);
    if (key<0) {
	Warning("AddListEntry","Unable to get entry %s",name);
	//Fatal("AddListEntry","Unable to get entry %s",name);
	return -1;
    }

    //test if linked list already initialized
    if (!GetParamInt (key,count, &i_count)) {
	//initialize
	//first add the new entry before modifying the old one
	if (targetkey < 0) targetkey = AddEntry(newname);
	if (targetkey < 0) {
	    Warning("AddListEntry","Unable to name entry %s",newname);
	    return -1;
	}
	//if this worked, set the link from the new to the old entry
	//TODO: More consitency checks

	i_count = new Int_t(1);
	SetParamInt (key, count,i_count);
	i_link  = new Int_t(targetkey);
	SetParamInt (key, link,i_link);
	i_newlink= new Int_t(key);
	SetParamInt (targetkey, link,i_newlink);
    } else {
	//first add the new entry before modifying the old one
	if (targetkey < 0) targetkey = AddEntry(newname);
	if (targetkey < 0) {
	    Warning("AddListEntry","Unable to name entry %s",newname);
	    return -1;
	}
	//increment target count
	GetParamInt (key, count,&i_count);
	(*i_count)++;

	//Find the last entry in the chain
	Int_t listkey=-1,lastkey;
	while (makeDataBase()->MakeListIterator(key, count,link , &listkey)) {
	    lastkey = listkey;
	}

	//copy old pointer to new entry
	//GetParamInt (key, link,&i_link);
	GetParamInt (lastkey, link,&i_link);
	i_newlink= new Int_t(*i_link);
	SetParamInt (targetkey, link,i_newlink);
	//set new entry
	i_link  = new Int_t(targetkey);
	SetParamInt (lastkey, link,i_link);
    }
    return targetkey;
};

Bool_t  PDataBase ::MakeListIterator(Int_t key, const char * count, 
				     const char * link, Int_t * listkey) {
    //get the list entries which belongs to "key" and is described by "counts" and "link"
    //(both has to be defined in a proper way)
    //It is very important that on the 1st call *listkey is -1
    //on the value kTRUE, *listkey contains the key link to the list entry
    //on the value kFALSE, the iteration has been finished (or not started due to an error)
    Int_t *i_count,*loc_listkey_p;
    
    if (key==-1) return kFALSE;
    if (*listkey==-1) {//first call: check list header entry
	if (count) {
	    if (!GetParamInt (key, count,&i_count)) {
		Warning("MakeListIterator","count %s not found",count);
		return kFALSE;
	    }
	}
	if (!GetParamInt (key, link,&loc_listkey_p)) {
	    Warning("MakeListIterator","link %s not found",link);
	    return kFALSE;
	}
    } else {
	if (!GetParamInt (*listkey, link,&loc_listkey_p)) return kFALSE;
    }
    //now copy to external
    *listkey=*loc_listkey_p;
    if (*listkey==key) {
	*listkey=-1;
	return kFALSE;
    }
    return kTRUE;
}

Bool_t  PDataBase ::MakeListIterator(Int_t key, Int_t count, Int_t link, Int_t * listkey) {
    //get the list entries which belongs to "key" and is described by "counts" and "link"
    //(both has to be defined in a proper way)
    //It is very important that on the 1st call *listkey is -1
    //on the value kTRUE, *listkey contains the key link to the list entry
    //on the value kFALSE, the iteration has been finished (or not started due to an error)
    Int_t *i_count,*loc_listkey_p;
    
    if (key==-1) return kFALSE;
    if (*listkey==-1) {//first call: check list header entry
	if (count>=0) {
	    if (!GetParamInt (key, count,&i_count)) {
		Warning("MakeListIterator","count %i not found",count);
		return kFALSE;
	    }
	}
	if (!GetParamInt (key, link,&loc_listkey_p)) {
	    Warning("MakeListIterator","link %i not found",link);
	    return kFALSE;
	}
    } else {
	if (!GetParamInt (*listkey, link,&loc_listkey_p)) return kFALSE;
    }
    //now copy to external
    *listkey=*loc_listkey_p;
    if (*listkey==key) {
	*listkey=-1;
	return kFALSE;
    }
    return kTRUE;
}



Bool_t   PDataBase ::ListEntries(Int_t key,Int_t option, const char * pattern) {
    //key=line (or -1 for all)
    //option=0 : line break =1: no line break
    //pattern like "mass,width"
    
    Int_t start=0;
    Int_t end=PDATABASE_MAX_LINES-1;
    Double_t *result;
    const char *result2;
    Int_t *result3;
    TObject *result4;
    Int_t sz[PDATABASE_MAX_LINES][PDATABASE_MAX_DOUBLE_PARAM+PDATABASE_MAX_STRING_PARAM]; 
    Int_t max_sz[PDATABASE_MAX_DOUBLE_PARAM+PDATABASE_MAX_STRING_PARAM]; 
    Int_t valid_key[PDATABASE_MAX_LINES];

    for (int i=0; i<(PDATABASE_MAX_DOUBLE_PARAM+PDATABASE_MAX_STRING_PARAM); i++) {
	max_sz[i]=0;
	for (int j=0; j<PDATABASE_MAX_LINES; j++) {
	    sz[j][i]=0;
	    valid_key[j]=0;
	}
    }


    //Interpretation of pattern
    TString spattern(pattern);
    TObjArray *array = spattern.Tokenize(TString(","));

    if (key>=0) {start=end=key;}
    for (int run=0;run<2; run++) {
	//run0: check size etc.
	//run1: print info
	for (int i=start; i<=end; i++) {
	    //TODO: check if at least one param is there
	    if (run && valid_key[i]>0) { if (option) {
		cout << i << ":";
		if (i<10) cout << " ";
		if (i<100) cout << " ";
		if (i<1000) cout << " ";
	    }
	    else
		cout << "Database key=" << i << endl;

	    }
	    for (int pat = 0; pat < array->GetEntriesFast(); pat++) {
		if (pat>=(PDATABASE_MAX_DOUBLE_PARAM+
			  PDATABASE_MAX_STRING_PARAM+
			  PDATABASE_MAX_INT_PARAM+
			  PDATABASE_MAX_TOBJ_PARAM)) {
		    //Too many string->Something is wrong
		    cout << "PDataBase::listEntries: Too many pattern strings" << endl;
		    return kFALSE;
		}
		TObjString *patoption = (TObjString *) (*array)[pat];
		char * options=(char*)patoption->GetString().Data();
		Bool_t checkline=kTRUE;
		Bool_t invert=kFALSE;
		if (options[0] == '*') {
		    //check paramteter, but not used for line selection
		    options++;
		    checkline=kFALSE;
		}
		if (options[0] == '!') { 
		    options++;
		    invert=kTRUE;
		}

		if (GetParamDouble (i, options, &result)) {
		    if (run) {
			if (valid_key[i]>0) {
			   
			    if (option) cout << options <<  "=";
			    else  cout << GetDescription(options) <<  "=";
			    printf("%f",*result);
			    if (option) cout << " "; else cout << " \n";
			}
		    } else {
			char bla[1000]; //I dont know a better way to get the length
			// (if somebody has an idea -> help yourself)
			sprintf(bla,"%f",*result);
			sz[i][pat]=strlen(bla);
			if (sz[i][pat] >  max_sz[pat]) max_sz[pat]=sz[i][pat];
			if (checkline && !invert) valid_key[i]++;
			else if (invert && checkline) valid_key[i]=-999;
		    }
		} 

		else if (GetParamInt (i, options,&result3)) {
		    if (run) {
			if (valid_key[i]>0) {
			    if (option) cout << options <<  "=";
			    else  cout << GetDescription(options) <<  "=";
			    printf("%i",*result3);
			    if (option) cout << " "; else cout << " \n";
			}
		    } else {
			char bla[1000]; //I dont know a better way to get the length
			// (if somebody has an idea -> help yourself)
			sprintf(bla,"%i",*result3);
			sz[i][pat]=strlen(bla);
			if (sz[i][pat] >  max_sz[pat]) max_sz[pat]=sz[i][pat];
			if (checkline && !invert) valid_key[i]++;
			else if (invert && checkline) valid_key[i]=-999;
		    }
		} 

		else if (GetParamString (i, options,&result2)) {
		    if (run) {
			if (valid_key[i]>0) {
			    if (option) cout << options <<  "=" << result2;
			    else  cout << GetDescription(options) <<  "="<< result2;
			    if (option) cout << " "; 
			    else cout << " \n";
			}
		    } else {
			sz[i][pat]=strlen(result2);
			if (sz[i][pat] >  max_sz[pat]) max_sz[pat]=sz[i][pat];
			if (checkline && !invert) valid_key[i]++;
			else if (invert && checkline) valid_key[i]=-999;
		    }
		} 

		else if (GetParamTObj (i, options,&result4)) {
		    if (run) {
			if (valid_key[i]>0) {
			   
			    if (option) cout << options <<  "=<yes>";
			    else  cout << GetDescription(options) <<  "=<yes>";
			    if (option) cout << " "; else cout << " \n";
			}
		    } else {
			sz[i][pat]=4; //for the yes
			if (sz[i][pat] >  max_sz[pat]) max_sz[pat]=sz[i][pat];
			if (checkline && !invert) valid_key[i]++;
			else if (invert && checkline) valid_key[i]=-999;
		    }
		} 
		 
//		 else if (invert && checkline) valid_key[i]++;
		 

		if (option && run && valid_key[i]>0) {
		    //fill missing gaps
		    int ga=sz[i][pat];
		    //account for missing option string
		    if (!ga) ga-=strlen(options)+2;
		    for (int x=ga; x<max_sz[pat];x++) {
			cout << " ";
		    }
		}
	    }
	     
	    if (run && valid_key[i]>0) cout << "\n"; 
	}
    }
    delete array;
    return kTRUE;
};



void PDataBase ::Performance(void) {
    TStopwatch timer;
    //way1: loop over pids, and get the name
    timer.Start();
    Int_t mypid=0, pid_param=0,*i_result;
    for (int r=0;r<10000;r++)
	for (int i=0;i<50;i++) {
	    //get the pid from name
	    int k=0;
	    for (k=0; (k<50) && strcmp(strings[i][0],strings[k][0]); k++) {
	    }
	    mypid+=k;
//	    cout << strings[i][0] << " is " << k << endl;
	}
    cout << timer.RealTime() << endl;
    cout << mypid << endl;
    mypid=0;

    timer.Start();
    //now the new way
    for (int r=0;r<10000;r++)
	for (int i=0;i<50;i++) {
	    GetFastParamInt ("pid",&pid_param);
	    if (! GetParamInt (GetEntry(strings[i][0]),pid_param, &i_result)) {
//	    if (! getFastParamInt (pid_param,i, &i_result)) {
		cout << "PStaticData::ID: " << i << " not found" << endl;	
	    } 
	    mypid+=*i_result;
//	    cout << strings[i][0] << " is " << *i_result << endl;
	}
    cout << timer.RealTime() << endl;
    cout << mypid << endl;

}


ClassImp(PDataBase)
