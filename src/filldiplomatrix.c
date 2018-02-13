#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP filldiplomatrix(SEXP HapMatrix, SEXP DiplMatrix){

SEXP null = R_NilValue;
PROTECT(null);
int I;
int J;
SEXP Rdim;
Rdim = getAttrib(HapMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

double *hap  = REAL(HapMatrix);
double *dipl = REAL(DiplMatrix);

double diplvalue;
double hapvalue1;
double hapvalue2;
int inter;
int inter2;
int count = 0;


for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){

 //convert 
 diplvalue =  hap[j +I*i];
 inter     = (int)(diplvalue/10);
 hapvalue1 = (double)(inter);
 inter2    = (int)diplvalue;
 hapvalue2 = (double)(inter2%10); 

 //fill diplo Matrix
 dipl[ 2*j       +  2*I*i]  = hapvalue1;
 dipl[ 2*j+1     +  2*I*i]  = hapvalue2;

 }
}
 

UNPROTECT(1);
return null;

}

