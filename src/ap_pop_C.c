#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP ap_pop_C(SEXP RinMatrix){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

//double value1;
//double value2;

Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

// pData = (int*) calloc (i,sizeof(int));

PROTECT(ret  = allocMatrix(REALSXP, 6+I+1, J));

// SEXP freq    = allocVector(REALSXP,3);
// double *fr   = REAL(freq);
// double *freq ;  
// freq  = (double*)malloc(3);

//SEXP sfreq   = allocVector(REALSXP,I+1);
//double *sfr  = REAL(sfreq);

double * sfreq;  
sfreq = (double*) calloc((I+1),sizeof(double));


double S       ;
double thetaS  ;
double thetaT  ;
double thetaFL ;
double thetaSA ;
double thetaTA ;

// Init freq
//for(int i=0; i< 3; i++){
//REAL(freq)[i]=0;
//freq[i]=0;
//}

// Init sfreq
for(int i=0; i< I+1; i++){
//REAL(sfreq)[i]=0; // default monomorph
sfreq[i]=0;

}

// Init Matrix
for (int i = 0; i < J; i++){
 for (int j = 0; j < I+6+1; j++){
  REAL(ret)[j +I*i] = 0;
 }
}

double an ;
int nullen;
int einsen;
int alle;

for (int i = 0; i < J; i++){


// Init sfreq
for(int i=0; i< I+1; i++){
//REAL(sfreq)[i]=0; // default monomorph
 sfreq[i]=0;
}
  
S       = 0;
thetaS  = 0;
thetaT  = 0;
thetaFL = 0;
thetaSA = 0;
thetaTA = 0;
nullen  = 0;
einsen  = 0;
alle    = 0;
 
 //count 0,1
 for (int j = 0; j < I; j++){
   if(Rval[j +I*i]==0){nullen++;}
   if(Rval[j +I*i]==1){einsen++;}
 }
 
 alle = nullen + einsen;  


 // only when include.unknown
 //if(alle>0){

   //REAL(sfreq)[einsen] = 1.0; // Anzahl der Einsen in der Spalte
     sfreq[einsen] = 1.0;

   if((einsen>0) && (einsen < alle)){
    S = 1.0;
   } // S !!!

 // }
 // freq[1]=alle, freq[2]=nullen, freq[3]=einsen

 // only when include unknown

  if((einsen > 0) && (einsen<alle)){

          // calculate an
          //an     = sum(1/1:(freq[1]-1)) // muss gemacht werden
          an=0.0;
          for(int xx=1; xx<alle; xx++){       
             an = an +  1.0/(double)xx;
	  }	
          
          // end of calculate an
          thetaS = 1.0/an;
         

         // thetaT= (freq[2]*freq[3])/((freq[1]*(freq[1]-1))/2)
          thetaT = (double)(nullen*einsen) / (double) ( (double)(alle*(alle-1))/2.0 );

         

         // thetaFL = 1*(freq[1]-1)/freq[1]
          if((einsen==1) || (nullen==1)){
              thetaFL = 1* (double)(alle-1)/(double)alle;
          }

          if((einsen> 1) && (nullen >1)){
             // thetaSA = 1/(an - freq[1]/(freq[1]-1))
             thetaSA = 1.0/(double) (an -  (double)alle/(double)(alle-1));
           
             
             
             // thetaTA <- freq[3]*freq[2]/(freq[1]*(freq[1]-1)/2) * (freq[1]-1)/(freq[1]-3)
             thetaTA = (double)(einsen*nullen) / (double)( (double)(alle*(alle-1))/2.0)       * (double)(alle-1)/(double)(alle-3);

             
          }
   }

 //Fill the matrix

 REAL(ret)[0 +(7+I)*i]= (double)S;
 REAL(ret)[1 +(7+I)*i]= (double)thetaS;
 REAL(ret)[2 +(7+I)*i]= (double)thetaT;
 REAL(ret)[3 +(7+I)*i]= (double)thetaFL;
 REAL(ret)[4 +(7+I)*i]= (double)thetaSA;
 REAL(ret)[5 +(7+I)*i]= (double)thetaTA;

// fill the rest sfreq
for(int yy=0; yy<I+1; yy++){
 REAL(ret)[(6+yy) +(7+I)*i]= (double)sfreq[yy]; //(double) REAL(sfreq)[yy];  
}
// End of filling sfeq


 //  matrixerg[1] <- S
 //  matrixerg[2] <- thetaS
 //  matrixerg[3] <- thetaT
 //  matrixerg[4] <- thetaFL
 //  matrixerg[5] <- thetaSA
 //  matrixerg[6] <- thetaTA

//matrixerg2 <- c(matrixerg,sfreq)
  
   //if(value2==6){break;}
   //if(value2==5){break;}  
   //if(value1!=value2){
   //	INTEGER(ret)[i]=1; // polymorphic
   //     break;      
   //} 
   
}

free(sfreq);

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
