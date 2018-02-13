#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP polyC(SEXP RinMatrix){

// APPROX SNP SEARCH

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
double value3;
int    treffer;
int    ispoly;
Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

PROTECT(ret = allocVector(INTSXP,J));

// Init ret

for(int i=0; i< J; i++){
INTEGER(ret)[i]=0; // default monomorph
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){

   value1  = Rval[I*i];
   treffer = 0;
   ispoly  = 0;

 for (int j = 0; j < I; j++){

   value2 = Rval[j + I*i];
  
   //printf("%f",value1);
  
   if(value2==6){//gaps
      INTEGER(ret)[i]=2;  
      break;
   }

   if(value2==5){//unknown
     INTEGER(ret)[i]=3;
     break;
   }
  
   if((value1!=value2)){  
      
        if((treffer) && (value2!=value3)){
          INTEGER(ret)[i]=4;// polyallelic
          ispoly=1;
        }
        value3 = value2;
        if(!ispoly){
	  INTEGER(ret)[i]=1; // polymorphic/biallelic !
        }
        treffer=1;                
   } 
   
 }
}

UNPROTECT(1);

return ret;

}

// include unknown !
SEXP polyCinclude(SEXP RinMatrix){

// APPROX SNP SEARCH

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
double value3;
int    treffer;
int    ispoly;
Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

PROTECT(ret = allocVector(INTSXP,J));

// Init ret

for(int i=0; i< J; i++){
INTEGER(ret)[i]=0; // default monomorph
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){

   value1  = Rval[I*i];
   treffer = 0;
   ispoly  = 0;

 for (int j = 0; j < I; j++){

   value2 = Rval[j + I*i];
  
   //printf("%f",value1);
  
   if(value2==6){//gaps
      INTEGER(ret)[i]=2;  
      break;
   }

   if(value2==5){//unknown
     // INTEGER(ret)[i]=3;
     continue;
   }
  
   if((value1!=value2)){  

	if(value1==5){
         value1=value2;
         continue;
    	}
      
        if((treffer) && (value2!=value3)){
          INTEGER(ret)[i]=4;// polyallelic
          ispoly=1;
        }
        value3 = value2;
        if(!ispoly){
	  INTEGER(ret)[i]=1; // polymorphic/biallelic !
        }
        treffer=1;                
   } 
   
 }
}

UNPROTECT(1);

return ret;

}





