#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP verify_ancestral_C(SEXP RinMatrix){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

SEXP ancestral   = allocVector(REALSXP,J);
double *anc      = REAL(ancestral);
PROTECT(ancestral);

// Init ancestral
for(int j=0; j< J; j++){
REAL(ancestral)[j]=2; // default polymorph -> not included
}

int nullen;
int einsen;
int alle;

for (int i = 0; i < J; i++){

nullen  =0;
einsen  =0;
alle    =0;
 
 //count 0,1
 for (int j = 0; j < I; j++){
   if(Rval[j +I*i]==0){nullen++;}
   if(Rval[j +I*i]==1){einsen++;}
 }
 
 alle = nullen + einsen;  
 if(einsen==alle){
    REAL(ancestral)[i] = 1; // 1 is ancestral
 }
 if(nullen==alle){
    REAL(ancestral)[i] = 0; // 0 is ancestral
 }
 
}

UNPROTECT(1);

return ancestral;

}



