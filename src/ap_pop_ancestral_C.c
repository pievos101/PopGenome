
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP ap_pop_ancestral_C(SEXP RinMatrix){

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

PROTECT(ret  = allocMatrix(REALSXP, 8+I, J));
SEXP freq    = allocVector(REALSXP,3);
PROTECT(freq);
double *fr   = REAL(freq);
SEXP sfreq   = allocVector(REALSXP,I);
PROTECT(sfreq);
double *sfr  = REAL(sfreq);


double S       ;
double thetaS  ;
double thetaT  ;
double thetaFL ;
double thetaSA ;
double thetaTA ;
double thetaFW ;
double thetaL  ;

// Init freq
for(int i=0; i< 3; i++){
REAL(freq)[i]=0; 
}
// Init sfreq
for(int i=0; i< I; i++){
REAL(sfreq)[i]=0; // default monomorph
}

// Init Matrix
// FIXME J I
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
  REAL(ret)[j +I*i]=0;
 }
}


double an ;
int nullen;
int einsen;
int alle;

for (int i = 0; i < J; i++){


 // Init sfreq
  for(int ff=0; ff< I; ff++){
  REAL(sfreq)[ff]=0; // default monomorph
  }
  
S       = 0;
thetaS  = 0;
thetaT  = 0;
thetaFL = 0;
thetaSA = 0;
thetaTA = 0;
thetaFW = 0;
thetaL  = 0;
nullen  = 0;
einsen  = 0;
alle    = 0;
 
 //count 0,1
 for (int j = 0; j < I-1; j++){
   if(Rval[j +I*i]==0){nullen++;}
   if(Rval[j +I*i]==1){einsen++;}
 }
 
 alle = nullen + einsen;  

 //printf("%d",nullen);

 if(Rval[(I-1)+I*i]==0){ // letzte Reihe
     REAL(sfreq)[einsen] = 1.0; // Anzahl der Einsen in der Spalte
         if((einsen>0) && (einsen < alle)){
            S = 1.0;
         } // S !!!
  }

 if(Rval[(I-1)+I*i]==1){ // letzte Reihe
     REAL(sfreq)[nullen] = 1.0; // Anzahl der Nullen
         if((nullen>0) && (nullen < alle)){
            S = 1.0;
         } // S !!!
    // change mutation
    int fr       = nullen;
    nullen       = einsen;
    einsen       = fr;
 }
 
// printf("%d",einsen);

if((einsen > 0) && (einsen<alle)){

          // calculate an
          //an     = sum(1/1:(freq[1]-1)) // muss gemacht werden
          an=0;
          for(int xx=1; xx<alle; xx++){       
             an = an +  1.0/(double)xx;
	  }	
          
          // end of calculate an
          thetaS = 1.0/an;
         

         // thetaT= (freq[2]*freq[3])/((freq[1]*(freq[1]-1))/2)
          thetaT = (double)(nullen*einsen) / (double) ( (double)(alle*(alle-1))/2.0 );

      if(Rval[(I-1)+I*i]==1 || Rval[(I-1)+I*i] ==0){   //ancestral information

          if(einsen==1){
              thetaFL = 1;
          }

          // thetaFW <-(freq[3]*freq[3])/(freq[1]*(freq[1]-1)/2)
          thetaFW = (double)(einsen*einsen)/(double) ( (double)(alle*(alle-1))/2.0 );
          // thetaL <- freq[3]/(freq[1]-1)
          thetaL  = (double)(einsen)/(double)(alle-1);
 
          if(einsen> 1){
             thetaSA = 1.0/(double) (an - 1);          
             // -2 mit outgroup
             thetaTA = (double)(einsen*nullen) / (double)( (double)(alle*(alle-1))/2.0) * (double)(alle-1)/(double)(alle-2);
          }
      }

  }

 //Fill the matrix

 REAL(ret)[0 +(8+I)*i]= (double)S;
 REAL(ret)[1 +(8+I)*i]= (double)thetaS;
 REAL(ret)[2 +(8+I)*i]= (double)thetaT;
 REAL(ret)[3 +(8+I)*i]= (double)thetaFL;
 REAL(ret)[4 +(8+I)*i]= (double)thetaSA;
 REAL(ret)[5 +(8+I)*i]= (double)thetaTA;
 REAL(ret)[6 +(8+I)*i]= (double)thetaFW;
 REAL(ret)[7 +(8+I)*i]= (double)thetaL;

// fill the rest sfreq
for(int yy=0; yy<I; yy++){
 REAL(ret)[(8+yy) +(8+I)*i]= (double) REAL(sfreq)[yy];  
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

UNPROTECT(3);

return ret;

}





