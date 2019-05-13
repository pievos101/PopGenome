#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP my_match_C(SEXP pos, SEXP Rvector)
{

PROTECT(pos);

  pos      = coerceVector(pos, INTSXP);
  PROTECT(pos);
  Rvector  = coerceVector(Rvector, INTSXP);

int size1;
int size2;
int *vec1 = INTEGER(pos);
int *vec2 = INTEGER(Rvector);

size1       = length(pos);
size2       = length(Rvector);

int treffer = 0;

int start = 0;

for(int i=0; i < size1; i++) {

  for (int j=start; j < size2; j++){

	if(vec1[i]==vec2[j]){
           treffer = 1;          
	   INTEGER(pos)[i] = j+1;
           start = j-1;
           break;
	}

  }

 if(treffer==0){
   INTEGER(pos)[i] = -1; 
 }
 treffer=0;

}

UNPROTECT(2);
return pos;


}

// returns the id in Rvector of bigger value than pos

SEXP whichbigger_C(SEXP pos, SEXP Rvector)
{

PROTECT(pos);

  pos      = coerceVector(pos, INTSXP);
  PROTECT(pos);
  Rvector  = coerceVector(Rvector, INTSXP);

int size1;
int size2;
int *vec1 = INTEGER(pos);
int *vec2 = INTEGER(Rvector);

size1       = length(pos);
size2       = length(Rvector);

int treffer = 0;

for(int j=0; j < size2; j++) {

	if(vec1[0]<=vec2[j]){          
	   INTEGER(pos)[0] = j + 1;
           treffer = 1;
           break;
	}
}

if(treffer==0){
  INTEGER(pos)[0] = 1;
}


UNPROTECT(2);
return pos;

}




