#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP fittingGFFC(SEXP RinMatrix, SEXP positions){

//printf("Geilo");

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;
SEXP Rpositions;



Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

//printf("Böö");

//Rlang             = coerceVector(lang, INTSXP);
//long int *ll      = INTEGER(Rlang);    

int ll = length(positions);

Rpositions            = coerceVector(positions, INTSXP);
int *Rpos             = INTEGER(Rpositions);

Rvalue                = coerceVector(RinMatrix, INTSXP);
int *Rval             = INTEGER(Rvalue);

PROTECT(ret = allocMatrix(INTSXP,I,J));

// Init Matrix
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
   INTEGER(ret)[j +I*i]=-1;
 }
}


int value1;
int value2;
int ende;
int treffer;
int kk;
int PosVal;
int start = 0;

//printf("%d",ll);

//printf("Hallo");

for (int j = 0; j < I; j++){

 treffer= 0;
 value1 =  Rval[j +I*0];
 value2 =  Rval[j +I*1];

   //  printf("%d",value1);
   //  printf("\n");
   //  printf("%d",value2);
   //  printf("\n");
   //  printf("%d",Rpos[1]);

// Searching if position exists
 
 for( int q = start; q < ll; q++ ){
    
    if((Rpos[q]>=value1) && (Rpos[q]<=value2)){

      if(treffer){
       INTEGER(ret)[j + I*1] = q+1;
       continue;
      }

     start=q;  // hahah .. geht nicht weil gff nicht sortiert :( // FIXED ! 
     treffer = 1;
     INTEGER(ret)[j + I*0] = q+1;
     INTEGER(ret)[j + I*1] = q+1;
      
    }else{
 
     if(treffer){break;}     
     continue;

    }


 } 
}

UNPROTECT(1);

return ret;

}


