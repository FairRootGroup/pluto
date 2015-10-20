// Author: I. Froehlich
// Written: 6.10.2009


#ifndef _PCOMMANDLIST_H_
#define _PCOMMANDLIST_H_

#include "TObjString.h"


class PCommandList : public TObjString {

 private:

    TObjString fCmd;  //command

    TObjString fNext;        //Next command 
    PCommandList *fNextP; //! Pointer to next command
    Int_t fLevel;         //! Level of nested objects

    TObjString fTool;  //command

 public:

    //constructor
    PCommandList();
    PCommandList(TString name);
    PCommandList(TString name, TString cmd);
    PCommandList(TObject * obj, TString name, TString cmd);
    ~PCommandList();
 
    const char*	GetCmd() const {
	return fCmd.GetString().Data();
    };

    Bool_t AddCommand(const char * cmd, const char * basename = NULL, int level = 0);
    Bool_t AddCommand(TObject * obj, const char * cmd, const char * basename = NULL, int level = 0);
    Bool_t GetCommand(char ** cmd, int level = 0,TObject ** obj = NULL);


    using TObject::Write; 
    Int_t Write(const char *name, Int_t option, Int_t bufsize); 
    void Print(const Option_t* delme= NULL) const;

    ClassDef(PCommandList,1)  //Tool class to make lists of commands
};

#endif
