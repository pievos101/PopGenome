#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP find_windowC(SEXP Rpositions, SEXP Rstart, SEXP Rend, SEXP RMERKEN)
{

SEXP val        = R_NilValue;
PROTECT(val     = Rf_allocVector(INTSXP,2));
INTEGER(val)[0] = 0;
INTEGER(val)[1] = 0;

 Rstart       = coerceVector(Rstart, INTSXP);
 PROTECT(Rstart);
 Rend         = coerceVector(Rend, INTSXP);
 PROTECT(Rend);
 Rpositions   = coerceVector(Rpositions, INTSXP);
 PROTECT(Rpositions);
 RMERKEN      = coerceVector(RMERKEN, INTSXP);
 PROTECT(RMERKEN);

int raus=0;
int size;
int *start  = INTEGER(Rstart);
int *end    = INTEGER(Rend);
int *pos    = INTEGER(Rpositions);
int *merken = INTEGER(RMERKEN);

size                = length(Rpositions);

int begin = merken[0] - 1;

if(start[0]>pos[size-1]){
    UNPROTECT(5);
 return R_NilValue;
}
if(end[0]<pos[0]){
   UNPROTECT(5);
 return R_NilValue;
}


for(int i=begin; i < size; i++) {

	//if(vec1[i]!=vec2[i]){          
	//   INTEGER(val)[0] = 0;
        //   break;
	//}

//printf("%f",start[0]);
//printf("%f",end[0]);
//printf("%f",pos[i]);

   if(start[0]<=pos[i]){

       if(start[0]<pos[i] && end[0]<pos[i]){
          UNPROTECT(5);
          return R_NilValue;
       }

       INTEGER(val)[0]=i+1;
       int startX = i;
       

       for(int xx = startX; xx < size; xx++){

          if(end[0]<=pos[xx]){

             if(end[0]==pos[xx]){
                INTEGER(val)[1]=xx+1;
                //raus = 1;
                //break;
                goto mark;
             }

            INTEGER(val)[1]=xx;
            //raus = 1;
            //break;
               goto mark;
          }
         
       }

    // rechte Seite ist über der bial-Grenze
    if(!raus){
    INTEGER(val)[1]=size;
    break;
    }
          
   }

}

mark:
UNPROTECT(5);
return val;


}

