#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP Ccompare(SEXP Rvector1, SEXP Rvector2)
{

SEXP val = R_NilValue;
PROTECT(val= Rf_allocVector(INTSXP,1));
INTEGER(val)[0] = 1;

// Rvector1 = coerceVector(Rvector1, INTSXP);
// Rvector2 = coerceVector(Rvector2, INTSXP);

int size;
const double *vec1 = REAL(Rvector1);
const double *vec2 = REAL(Rvector2);
size            = length(Rvector1);

for(int i=0; i < size; i++) {

	if(vec1[i]!=vec2[i]){          
	   INTEGER(val)[0] = 0;
           break;
	}
}

UNPROTECT(1);
return val;


}

