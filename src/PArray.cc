////////////////////////////////////////////////////////
//  Small tool class to convert text files into TGraph 
//  objects
//  
//
//                    Author: I. Froehlich
//                    Written: 1.08.2007
//                    Revised: 
//
////////////////////////////////////////////////////////


#include "PArray.h"
#include <iostream>

using namespace std;

PArray::PArray(Int_t dimension) {
    dim = dimension;
    fp  = NULL;
    real_size_1d =0;
    scaling=1.;

    for (int i=0;i<PARRAY_MAX_COLUMNS;i++)
	vals_1d[i]=NULL;
}

PArray::~PArray() {
};

Bool_t PArray::OpenFile(const char * filename, Int_t syntax, Int_t num_columns, Double_t y) {

    if ((dim == 1) && vals_1d[0]) {
	Warning("OpenFile","Cannot read more then one file for 1dim arrays");
	return kFALSE;
    }

    fp = fopen(filename,"r");   // open file

    Double_t val[PARRAY_MAX_COLUMNS];
    Int_t    row=0;

    if (!vals_1d[0]) {
	for (int i=0;i<num_columns;i++) {
	    vals_1d[i]=new TArrayD(PARRAY_GRANULARITY);
	}
	if (dim == 2)
	    vals_2d=new TArrayD(PARRAY_GRANULARITY);
    }

    if (fp==NULL) {
	Error("OpenFile","Cannot open file %s",filename);
	return kFALSE;
    }
    
    Int_t ret = 0;
    while (ret != EOF) {
	if (num_columns == 5) {
	    ret = fscanf(fp,"%le %le %le %le %le", &val[0], &val[1], &val[2], &val[3], &val[4]);
	}

	for (int i=1;i<num_columns;i++)
	    val[i] *= scaling;

	if (ret != EOF) {

	    if (dim ==1) {
		if (vals_1d[0]->GetSize() > row) {
		    for (int i=0;i<num_columns;i++)
			vals_1d[i]->Set(vals_1d[0]->GetSize() + PARRAY_GRANULARITY);
		}
		for (int i=0;i<num_columns;i++)
		    vals_1d[i]->AddAt(val[i],row);

	    }

	    if (dim ==2) {
		if (vals_1d[0]->GetSize() > (real_size_1d+row)) {
		    for (int i=0;i<num_columns;i++) {
			vals_1d[i]->Set(vals_1d[0]->GetSize() + PARRAY_GRANULARITY);
		    }
		    vals_2d->Set(vals_2d->GetSize() + PARRAY_GRANULARITY);
		}
		for (int i=0;i<num_columns;i++) {
		    vals_1d[i]->AddAt(val[i],real_size_1d+row);
		}
		vals_2d->AddAt(y,real_size_1d+row);
	    }
	    row++;
	}
    }
    real_size_1d+=row;

    Info("OpenFile","%d lines read from file %s",row, filename);
    
    fclose(fp);
    return kTRUE;

}


TGraph* PArray::GetTGraph(Int_t xcol, Int_t ycol) {
    //NB: TGraph is not killed, it must be deleted in the app.
    
    if (dim == 1) {
	if (!vals_1d[0]) {
	    Error ("GetTGraph","No (1dim) file opened");
	    return NULL;
	}
	TGraph * delme = new TGraph(real_size_1d,vals_1d[xcol]->GetArray(),vals_1d[ycol]->GetArray() );
	return delme;

    }
    return NULL;
}

TGraph2D* PArray::GetTGraph2D(Int_t xcol, Int_t ycol) {
    //NB: TGraph is not killed, it must be deleted in the app.
    
    if (dim == 2) {
	if (!vals_1d[0]) {
	    Error ("GetTGraph","No (2dim) files opened");
	    return NULL;
	}

	TGraph2D * delme = new TGraph2D(real_size_1d,vals_1d[xcol]->GetArray(),
					vals_2d->GetArray(),vals_1d[ycol]->GetArray() );
// 	for (int i=0;i<real_size_1d;i++)
// 	    cout << (vals_1d[xcol]->GetArray())[i] << ":"
// 		 << (vals_2d->GetArray())[i] << ":"
// 		 << (vals_1d[ycol]->GetArray())[i] << ":"
// 		 << endl;
	return delme;

    }
    return NULL;
}

ClassImp(PArray)
