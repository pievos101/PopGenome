#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>


SEXP C_get_sfreqh_C(SEXP RRuniquematrix, SEXP RRmatrix)
{

SEXP freq = R_NilValue; // frequenz of each haplotype

SEXP Rdim;
SEXP Runiquematrix;
SEXP Rmatrix;

int I1;
int J1;

Rdim  = getAttrib(RRuniquematrix, R_DimSymbol);
I1    = INTEGER(Rdim)[0]; // Reihen 
J1    = INTEGER(Rdim)[1]; // Spalten

int I2;
int J2;

Rdim  = getAttrib(RRmatrix, R_DimSymbol);
I2    = INTEGER(Rdim)[0]; // Reihen 
J2    = INTEGER(Rdim)[1]; // Spalten



PROTECT(freq     = Rf_allocVector(INTSXP,I1));


Runiquematrix           = coerceVector(RRuniquematrix, INTSXP);
int *R1                 = INTEGER(Runiquematrix);

Rmatrix                 = coerceVector(RRmatrix, INTSXP);
int *R2                 = INTEGER(Rmatrix);

int count;
int equal;

for (int j = 0; j < I1; j++){

 count = 0; 

 for (int j2 = 0; j2 < I2; j2++){
 
  equal = 1 ;

  for (int i = 0; i < J1; i++){ // compare each element of the two sequences
    
     if(R1[j + I1*i] != R2[j2 + I2*i]){
       equal = 0;
       break;
     }

     //Rval[j +I*i]  
  }

 if(equal){count ++;}

 }

INTEGER(freq)[j] = count ; 

}

UNPROTECT(1);

return(freq);

}

