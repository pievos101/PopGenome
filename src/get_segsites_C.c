#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP get_segsites_C(SEXP RinMatrix){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

double value1;
double value2;

Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

PROTECT(ret = allocVector(LGLSXP,J));

// Init ret

for(int i=0; i< J; i++){
 LOGICAL(ret)[i]=FALSE; // default monomorph
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){
   value1 = Rval[I*i];
 for (int j = 0; j < I; j++){
   value2 = Rval[j +I*i];
  
   //printf("%f",value1);
   if(value1!=value2){
        //printf("Treffer");
	LOGICAL(ret)[i]=TRUE; // polymorphic
        break;      
   } 
   
 }
}



UNPROTECT(1);

return ret;

}




/*SEXP Ccompare(SEXP Rvector1, SEXP Rvector2)
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
*/
