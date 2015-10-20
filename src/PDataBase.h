// Author: I. Froehlich
// Written: 11.04.2007
// Revised: 
// PDataBase
// Replacement for the particle data base in PData

#ifndef _PDATABASE_H_
#define _PDATABASE_H_

#include "TROOT.h"
#include <iostream>
#include "Pdefines.h"

#define PDATABASE_MAX_DOUBLE_PARAM 20
#define PDATABASE_MAX_STRING_PARAM 20
#define PDATABASE_MAX_INT_PARAM 150
#define PDATABASE_MAX_TOBJ_PARAM 10

#define PDATABASE_MAX_LINES 1000

using namespace std;

class PDataBase;
PDataBase* makeDataBase();
PDataBase& fDataBase();

class PDataBase : public TObject {

 private:

    //Entries properties
    const char*   param_double_name[PDATABASE_MAX_DOUBLE_PARAM]; //the parameter name matching the ID
    const char*   param_string_name[PDATABASE_MAX_STRING_PARAM]; 
    const char*   param_int_name[PDATABASE_MAX_INT_PARAM];       
    const char*   param_tobj_name[PDATABASE_MAX_INT_PARAM]; 
    const char*   param_double_descr[PDATABASE_MAX_DOUBLE_PARAM]; //the parameter long description
    const char*   param_string_descr[PDATABASE_MAX_STRING_PARAM]; 
    const char*   param_int_descr[PDATABASE_MAX_INT_PARAM];       
    const char*   param_tobj_descr[PDATABASE_MAX_INT_PARAM]; 
    Int_t   param_double_pointer;
    Int_t   param_string_pointer;
    Int_t   param_int_pointer;
    Int_t   param_tobj_pointer;

    //TODO: Params need description

    Int_t  *param_int_key[PDATABASE_MAX_INT_PARAM];
    Int_t   param_int_key_max[PDATABASE_MAX_INT_PARAM];

    //The entries array itself
    //for each entry, we add the pointer array to our params
    const char* strings[PDATABASE_MAX_LINES][PDATABASE_MAX_STRING_PARAM];
    Double_t*   doubles[PDATABASE_MAX_LINES][PDATABASE_MAX_DOUBLE_PARAM];
    Int_t*      ints[PDATABASE_MAX_LINES][PDATABASE_MAX_INT_PARAM];
    TObject*    tobjs[PDATABASE_MAX_LINES][PDATABASE_MAX_TOBJ_PARAM];

    Int_t lastkey;

    Bool_t CheckEntry(Int_t key);  //check if key with unique name already existing
    Int_t  ConvertParamKey(const char * &newparamname,Int_t key);

 public:

    //constructor
    PDataBase();

    void Performance(void);

    void    SetFastKey(Int_t pkey, Int_t maxkey);

    //return value of the makeParam's: the param_id of -1 on failure
    Int_t   MakeParamDouble (const char* paramname, 
			     const char *descr);//add "paramname" to double database
    Int_t   MakeParamString(const char* paramname, 
			    const char *descr); //add "paramname" to string database
    Int_t   MakeParamInt(const char* paramname, 
			 const char *descr);    //add "paramname" to int database
    Int_t   MakeParamTObj(const char* paramname, 
			  const char *descr);   //add "paramname" to object database

    Int_t   GetParamDouble (const char*  paramname); //get "paramname" from double database
    Int_t   GetParamString(const char* paramname);   //get "paramname" from string database
    Int_t   GetParamInt(const char* paramname, Int_t length=-1);  //get "paramname" from string database
    Int_t   GetParamTObj(const char* paramname);

    char*   GetDescription(const char* paramname);

    void    GetFastParamInt(const char* paramname, Int_t *pkey)
	{
	    if ((*pkey)<0) *pkey=GetParamInt(paramname);
	};  
    void    GetFastParamString(const char* paramname, Int_t *pkey)
	{
	    if ((*pkey)<0) *pkey=GetParamString(paramname);
	}; 
    void    GetFastParamDouble(const char* paramname, Int_t *pkey)
	{
	    if ((*pkey)<0) *pkey=GetParamDouble(paramname);
	}; 
    void    GetFastParamTObj(const char* paramname, Int_t *pkey)
	{
	    if ((*pkey)<0) *pkey=GetParamTObj(paramname);
	}; 
    
    using TObject::GetName;
    const char* GetName(Int_t key);
    //getting particle/decay properties from data base 
    //return kFALSE if 1.) key not existing or 2.) Param not existing or 3.) Param not used for entry
    Bool_t  GetParamDouble (Int_t key, const char*  paramname, Double_t ** result);
    Bool_t  GetParamString (Int_t key, const char*  paramname, const char ** result);
    Bool_t  GetParamInt  (Int_t key, const char*  paramname, Int_t ** result, Int_t length=-1);
    Bool_t  GetParamTObj (Int_t key, const char*  paramname, TObject ** result);

    //same as above but with the pkey as obtained by the getFastParam's
    Bool_t  GetParamInt (Int_t key, Int_t pkey, Int_t ** result);
    Bool_t  GetParamString (Int_t key, Int_t pkey, const char ** result);
    Bool_t  GetParamDouble (Int_t key, Int_t pkey, Double_t ** result);
    Bool_t  GetParamTObj   (Int_t key, Int_t pkey, TObject ** result);

    //setting value instead of pointer
    Bool_t  GetParamInt (Int_t key, Int_t pkey, Int_t * result) {
	Int_t *p;
	Bool_t b=GetParamInt (key, pkey, &p);
	if (b && p)
	    *result=*p;
	else 
	    *result=0;
	return b;
    };

    //getting param by primary name
    Bool_t  GetParamInt (const char * name, const char*  paramname, Int_t ** result){
	return GetParamInt (GetEntry(name), paramname, result);
	};
    Bool_t  GetParamDouble (const char * name, const char*  paramname, Double_t ** result){
	return GetParamDouble (GetEntry(name), paramname, result);
	};
    Bool_t  GetParamString (const char * name, const char*  paramname, const char ** result){
	return GetParamString (GetEntry(name), paramname, result);
	};
    Bool_t  GetParamTObj (const char * name, const char*  paramname, TObject ** result){
	return GetParamTObj (GetEntry(name), paramname, result);
	};

    //Getting param by matching integer value
    Bool_t  GetParamInt (const char*  paramname1, Int_t value1, const char * paramname2, Int_t ** result){	
	return GetParamInt (GetEntryInt(paramname1,value1), paramname2, result);
    };
    Bool_t  GetParamDouble (const char*  paramname1, Int_t value1, const char * paramname2, Double_t ** result){	
	return GetParamDouble (GetEntryInt(paramname1,value1), paramname2, result);
    };
    Bool_t  GetParamString (const char*  paramname1, Int_t value1, const char * paramname2, const char ** result){	
	return GetParamString (GetEntryInt(paramname1,value1), paramname2, result);
    };
    Bool_t  GetParamTObj (const char*  paramname1, Int_t value1, const char * paramname2, TObject ** result){	
	return GetParamTObj (GetEntryInt(paramname1,value1), paramname2, result);
    };

    //same as above when the pkey are known...
    Bool_t  GetParamInt (Int_t pkey1, Int_t value1, Int_t pkey2, Int_t ** result){
//	cout << pkey1 << ":" << value1<< ":" << pkey2<< ":" << endl;
	return GetParamInt (GetEntryInt(pkey1,value1), pkey2, result);
    };
    Bool_t  GetParamDouble (Int_t pkey1, Int_t value1, Int_t pkey2, Double_t ** result){	
	return GetParamDouble (GetEntryInt(pkey1,value1), pkey2, result);
    };
    Bool_t  GetParamString (Int_t pkey1, Int_t value1, Int_t pkey2, const char ** result){	
	return GetParamString (GetEntryInt(pkey1,value1), pkey2, result);
    };
    Bool_t  GetParamTObj (Int_t pkey1, Int_t value1, Int_t pkey2, TObject ** result){	
	return GetParamTObj (GetEntryInt(pkey1,value1), pkey2, result);
    };


    Bool_t  SetParamDouble (Int_t key, const char*  paramname, Double_t * result);
    Bool_t  SetParamString (Int_t key, const char*  paramname, char * result);   
    Bool_t  SetParamInt (Int_t key, const char*  paramname, Int_t * result);   
    Bool_t  SetParamTObj (Int_t key, const char*  paramname, TObject * result);   

    //faster:
    Bool_t  SetParamTObj (Int_t key, Int_t pp, TObject * result);

    Bool_t  SetParamDouble (const char * name, const char*  paramname, Double_t result) {
	Int_t key = GetEntry(name);
	if (key<0) return kFALSE;
	return SetParamDouble (key, paramname, new Double_t(result));
    }
    Bool_t  SetParamInt (const char * name, const char*  paramname, Int_t result) {
	Int_t key = GetEntry(name);
	if (key<0) return kFALSE;
	return SetParamInt (key, paramname, new Int_t(result));
    }
    Bool_t  SetParamString (const char * name, const char*  paramname, char * result) {
	Int_t key = GetEntry(name);
	if (key<0) return kFALSE;
	return SetParamString (key, paramname, result);
    };   
    Bool_t  SetParamTObj (const char * name, const char*  paramname, TObject * result) {
	Int_t key = GetEntry(name);
	if (key<0) return kFALSE;
	return SetParamTObj (key, paramname, result);
    }; 

    //Dealing with the entries
    Bool_t  AddEntry(Int_t key, const char * name);
    Int_t   AddEntry(const char * name);
    Int_t   GetEntry(const char* name); 
    Int_t   GetEntryInt(const char* paramname, Int_t value); 
    Int_t   GetEntryInt(Int_t pkey, Int_t value); 

    //listings
    Int_t   AddListEntry(const char * name, const char * count, const char * link, const char * newname);
    Bool_t  MakeListIterator(Int_t key, const char * count, const char * link, Int_t * listkey);
    Bool_t  MakeListIterator(Int_t key, Int_t count, Int_t link, Int_t * listkey);

    Bool_t  ListEntries(Int_t key=-1,Int_t option=0, const char * pattern=NULL);

    void Print(const Option_t* delme=NULL) const ;


    ClassDef(PDataBase,0)  //Relational data base for pluto
};

#endif
