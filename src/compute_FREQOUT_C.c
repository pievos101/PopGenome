
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP compute_FREQOUT_C(SEXP RinMatrix){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen           
// BUT !
// last row is ancestral
J    = INTEGER(Rdim)[1]; // Spalten

double value1;
double value2;

Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

SEXP freq    = allocVector(REALSXP,3);
PROTECT(freq);
double *fr   = REAL(freq);
SEXP sfreq   = allocMatrix(REALSXP,I,J);
PROTECT(sfreq);
double *sfr  = REAL(sfreq);


// Init freq
for(int i=0; i< 3; i++){
REAL(freq)[i]=0; 
}

// Init Matrix sfreq
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
  REAL(sfreq)[j +I*i]=0;
 }
}

int nullen;
int einsen;
int alle;

for (int i = 0; i < J; i++){

nullen  = 0;
einsen  = 0;
alle    = 0;
 
 //count 0,1
 for (int j = 0; j < I-1; j++){
   if(Rval[j +I*i]==0){nullen++;}
   if(Rval[j +I*i]==1){einsen++;}
 }
 
 alle = nullen + einsen;  

 if(Rval[(I-1)+I*i]==0){ // letzte Reihe
     REAL(sfreq)[einsen+I*i] = 1.0; // Anzahl der Einsen in der Spalte
 }

 if(Rval[(I-1)+I*i]==1){ // letzte Reihe
     REAL(sfreq)[nullen+I*i] = 1.0; // Anzahl der Nullen
 }
 
}

UNPROTECT(2);
return sfreq;

}





