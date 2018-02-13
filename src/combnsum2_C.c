#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP combnsum2_C(SEXP RinMatrix, SEXP RinMatrix2){

// APPROX SNP SEARCH

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;
SEXP Rvalue2;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

int value1;
int value2;
int prod  = 0;
int summe = 0; 

Rvalue           = coerceVector(RinMatrix, INTSXP);
int *Rval        = INTEGER(Rvalue);
Rvalue2          = coerceVector(RinMatrix2, INTSXP);
int *Rval2       = INTEGER(Rvalue2);



PROTECT(ret = allocVector(INTSXP,1));

// Init ret
INTEGER(ret)[0]=0; // diversity


//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){

   value1  = Rval[i];

 for (int j = 0; j < J; j++){

   if(i!=j){
    value2 = Rval2[j];
    prod   = value1 * value2;
    summe  = summe + prod;
   }

 }
}

INTEGER(ret)[0] = summe;


UNPROTECT(1);

return ret;

}

