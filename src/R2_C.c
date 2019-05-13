#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP R2_C(SEXP RinMatrix, SEXP EINSEN, SEXP NULLEN){

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
double *einsen   = REAL(EINSEN);
double *nullen   = REAL(NULLEN);

double  ones;
double  zeros;
double     id;
double     id2;
double     count;
double valid_comp;

double  freqsite1;
double  freqsite2;
double  freqsite1_low;
double  freqsite2_low;
double  site_length = I;
double  d_raw;
double  r2;
int     fill;

PROTECT(ret = allocVector(REALSXP,(J*(J-1))/2));

// Init ret

for(int i=0; i < (J*(J-1))/2; i++){
REAL(ret)[i]=0; // default monomorph
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}


fill = 0;
for (int m=0; m < J-1; m++){
 
 // sites 1
    ones  = einsen[m];
    zeros = nullen[m];
    site_length = ones + zeros;

//    zeros = site_length - ones;

    if(ones>=zeros){
       freqsite1 = ones/site_length;
       id = 1;
    }else{
       freqsite1 = zeros/site_length;
       id = 0;
    } 
 // End of sites 1
 
 for (int i = m+1; i < J; i++){

 // sites 2
    ones  = einsen[i];
    zeros = nullen[i];
    site_length = ones + zeros;

//    zeros = I - ones;

    if(ones>=zeros){
       freqsite2 = ones/site_length;
       id2 = 1;
    }else{
       freqsite2 = zeros/site_length;
       id2 = 0;
    } 
 // End of sites 2
    
   count = 0;
   valid_comp = 0;

   for (int j = 0; j < I; j++){
   
    value1 = Rval[j +I*m];
    value2 = Rval[j +I*i];
	
    //count valid comparisons
    if((value1==0 && value2==0) || (value1==1 && value2==0) || (value1==0 && value2==1) || (value1==1 && value2==1)){
    valid_comp ++;
    }    

     if(value1 == id && value2 == id2){
       count ++;
     }
    
   }

   //d_raw  = (double)count/(double)site_length - freqsite1*freqsite2;
   if(valid_comp==0){continue;}
   d_raw  = (double)count/(double)valid_comp - freqsite1*freqsite2;

   freqsite1_low = 1-freqsite1;
   freqsite2_low = 1-freqsite2;
   r2 = (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low);
   REAL(ret)[fill] = r2;
   fill ++;

   // printf("%f",value1);
   // if(value2==6){break;}
   // if(value2==5){break;}  
   // if(value1!=value2){
        //printf("Treffer");
   //	INTEGER(ret)[i]=1; // polymorphic
   //     break;       
 }
}

//

UNPROTECT(1);

return ret;

}


SEXP R2_between_C(SEXP RRbial1,SEXP RReinsen1, SEXP RRbial2, SEXP RReinsen2){

int I1; // n.rows of bial1
int J1; // n.cols of bial1
int I2; // n.rows of bial2
int J2; // n.cols of bial2

// Get dimensions
SEXP Rdim;
Rdim  = getAttrib(RRbial1, R_DimSymbol);
I1    = INTEGER(Rdim)[0]; // Reihen 
J1    = INTEGER(Rdim)[1]; // Spalten
Rdim = getAttrib(RRbial2, R_DimSymbol);
I2    = INTEGER(Rdim)[0]; // Reihen 
J2    = INTEGER(Rdim)[1]; // Spalten


SEXP Rbial1;
Rbial1            = coerceVector(RRbial1, REALSXP);
double *bial1     = REAL(Rbial1);
double *einsen1   = REAL(RReinsen1);

SEXP Rbial2;
Rbial2            = coerceVector(RRbial2, REALSXP);
double *bial2     = REAL(Rbial2);
double *einsen2   = REAL(RReinsen2);


SEXP R2;
R2 = allocVector(REALSXP,J1*J2);
PROTECT(R2);
SEXP M;
M = allocMatrix(INTSXP,J1*J2,4); // for the fisher exact test
PROTECT(M);

// Init R2
for(int i=0; i< J1*J2; i++){
REAL(R2)[i]=0; 
}

// Init M
for (int i = 0; i < 4; i++){
 for (int j = 0; j < J1*J2; j++){
   INTEGER(M)[j +(J1*J2)*i]=0;
 }
}

// Algorithm 
double  ones;
double  zeros;
double  id;
double  id2;
double  count;
double  value1;
double  value2;
double  freqsite1;
double  freqsite2;
double  freqsite1_low;
double  freqsite2_low;
double  site_length = I1;
double  d_raw;
double  r2;
int     fill;

int m00=0;
int m01=0;
int m10=0;
int m11=0;


fill = 0;
for (int m=0; m < J1; m++){
 
 // sites 1 of bial 1
    ones  = einsen1[m];
    zeros = site_length - ones;
    if(ones>=zeros){
       freqsite1 = ones/site_length;
       id = 1;
    }else{
       freqsite1 = zeros/site_length;
       id = 0;
    } 
 
 for (int i = 0; i < J2; i++){

 // sites 2 of bial 2
    ones  = einsen2[i];
    zeros = site_length - ones;
    if(ones>=zeros){
       freqsite2 = ones/site_length;
       id2 = 1;
    }else{
       freqsite2 = zeros/site_length;
       id2 = 0;
    } 

   count = 0;
   m00=0;
   m01=0;
   m10=0;
   m11=0;

   for (int j = 0; j < I1; j++){
   
    value1 = bial1[j +I1*m];
    value2 = bial2[j +I1*i];
     if(value1==id && value2 ==id2){
       count ++;
     } 
     // fill M Matrix
     if(value1==0 && value2 == 0){
	m00++;
     }
     if(value1==0 && value2 == 1){
	m01++;
     }
     if(value1==1 && value2 == 0){
	m10++;
     }
     if(value1==1 && value2 == 1){
	m11++;
     }

   }
 
  // Fill Matrix M
  INTEGER(M)[fill+(J1*J2)*0] = m00;
  INTEGER(M)[fill+(J1*J2)*1] = m01;
  INTEGER(M)[fill+(J1*J2)*2] = m10;
  INTEGER(M)[fill+(J1*J2)*3] = m11;

   d_raw  = (double)count/(double)site_length - freqsite1*freqsite2;
   freqsite1_low = 1-freqsite1;
   freqsite2_low = 1-freqsite2;
   r2 = (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low);
   REAL(R2)[fill] = r2;
   fill ++;
 
 }

}

SEXP list = R_NilValue;

   // Creating a list with 3 vector elements:
   PROTECT(list = allocVector(VECSXP, 2)); 
     // attaching myint vector to list:
   SET_VECTOR_ELT(list, 0, R2); 
     // attaching mydouble vector to list:
   SET_VECTOR_ELT(list, 1, M); 
     // and attaching the vector names:
     //setAttrib(list, R_NamesSymbol, list_names); 
   //SET_VECTOR_ELT(list, 2, SUBST);
 
UNPROTECT(3);
   
return list;

}


SEXP R2_C_plus(SEXP RinMatrix, SEXP EINSEN, SEXP NULLEN, SEXP bialsites){

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
double *einsen   = REAL(EINSEN);
double *nullen   = REAL(NULLEN);
double *bialpos  = REAL(bialsites);

double  ones;
double  zeros;
double     id;
double     id2;
double     count;
double valid_comp;

double  freqsite1;
double  freqsite2;
double  freqsite1_low;
double  freqsite2_low;
double  site_length = I;
double  d_raw;
double  r2;
int     fill;
int m00=0;
int m01=0;
int m10=0;
int m11=0;

ret = allocVector(REALSXP,(J*(J-1))/2);
PROTECT(ret);
SEXP M;
M = allocMatrix(INTSXP,(J*(J-1))/2,4); // for the fisher exact test
PROTECT(M);
SEXP Dist;
Dist = allocVector(REALSXP,(J*(J-1))/2);
PROTECT(Dist);
SEXP SNP1;
SNP1 = allocVector(REALSXP,(J*(J-1))/2);
PROTECT(SNP1);
SEXP SNP2;
SNP2 = allocVector(REALSXP,(J*(J-1))/2);
PROTECT(SNP2);

// Init Dist
for(int i=0; i < (J*(J-1))/2; i++){
REAL(Dist)[i]=0; 
}
// Init SNP1
for(int i=0; i < (J*(J-1))/2; i++){
REAL(SNP1)[i]=0; 
}
// Init SNP2
for(int i=0; i < (J*(J-1))/2; i++){
REAL(SNP2)[i]=0; 
}
// Init ret
for(int i=0; i < (J*(J-1))/2; i++){
REAL(ret)[i]=0; 
}

// Init M
for (int i = 0; i < 4; i++){
 for (int j = 0; j < ((J*(J-1))/2); j++){
   INTEGER(M)[j +((J*(J-1))/2)*i]=0;
 }
}

fill = 0;
for (int m=0; m < J-1; m++){
 
 // sites 1
    ones  = einsen[m];
    //zeros = site_length - ones;
    zeros = nullen[m];
    site_length = ones + zeros;
	
    if(ones>=zeros){
       freqsite1 = ones/site_length;
       id = 1;
    }else{
       freqsite1 = zeros/site_length;
       id = 0;
    } 
 // End of sites 1
 
 for (int i = m+1; i < J; i++){

 // sites 2
    ones  = einsen[i];
    zeros = nullen[i];
    site_length = ones + zeros;

    //zeros = I - ones;
    if(ones>=zeros){
       freqsite2 = ones/site_length;
       id2 = 1;
    }else{
       freqsite2 = zeros/site_length;
       id2 = 0;
    } 
 // End of sites 2
    
   count = 0;
   valid_comp = 0;
   m00=0;
   m01=0;
   m10=0;
   m11=0;

   for (int j = 0; j < I; j++){
   
    value1 = Rval[j +I*m];
    value2 = Rval[j +I*i];
    
     if(value1==id && value2 ==id2){
       count ++;
     }
   	  // fill M Matrix
     if(value1==0 && value2 == 0){
        valid_comp ++;
	m00++;
     }
     if(value1==0 && value2 == 1){
	m01++;
        valid_comp ++;
     }
     if(value1==1 && value2 == 0){
	m10++;
	valid_comp ++;
     }
     if(value1==1 && value2 == 1){
	m11++;
        valid_comp ++;
     }
   }


 
  // Fill Matrix M
  INTEGER(M)[fill+((J*(J-1))/2)*0] = m00;
  INTEGER(M)[fill+((J*(J-1))/2)*1] = m01;
  INTEGER(M)[fill+((J*(J-1))/2)*2] = m10;
  INTEGER(M)[fill+((J*(J-1))/2)*3] = m11;
  REAL(Dist)[fill] = bialpos[i] - bialpos[m];
  REAL(SNP1)[fill] = bialpos[i];
  REAL(SNP2)[fill] = bialpos[m];

   if(valid_comp==0){continue;}
   //d_raw  = (double)count/(double)site_length - freqsite1*freqsite2;
     d_raw  = (double)count/(double)valid_comp - freqsite1*freqsite2;
	
   freqsite1_low = 1-freqsite1;
   freqsite2_low = 1-freqsite2;
   r2 = (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low);
   REAL(ret)[fill] = r2;
   fill ++;
  
 }
}

SEXP list = R_NilValue;

   // Creating a list with 3 vector elements:
   PROTECT(list = allocVector(VECSXP, 5)); 
     // attaching myint vector to list:
   SET_VECTOR_ELT(list, 0, ret); 
     // attaching mydouble vector to list:
   SET_VECTOR_ELT(list, 1, M); 
   SET_VECTOR_ELT(list, 2, Dist);
   SET_VECTOR_ELT(list, 3, SNP1);
   SET_VECTOR_ELT(list, 4, SNP2);

     // and attaching the vector names:
     //setAttrib(list, R_NamesSymbol, list_names); 
   //SET_VECTOR_ELT(list, 2, SUBST);
 
//UNPROTECT(1);
UNPROTECT(6);

return list;

}

SEXP count_congruent(SEXP RinMatrix){

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
Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);

PROTECT(ret = allocVector(INTSXP,J-1));

// Init ret
for(int i=0; i< J-1; i++){
INTEGER(ret)[i]=1; // congruent
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J-1; i++){
 for (int j = 0; j < I; j++){

   value1 = Rval[j + I*i];     // site 1
   value2 = Rval[j + I*(i+1)]; // site 2
  
   //printf("%f",value1);
   //printf(":");
   //printf("%f",value2);
   //printf("\n");
   	
   // one of the values are Nan
   if((value1 != 0 && value1!=1) || (value2 != 0 && value2!=1)){continue;}  

   if((value1!=value2)){  
   INTEGER(ret)[i] = 0;
   break;                   
   } 
   
 }

//printf("new site pair");
//printf("\n");

}

UNPROTECT(1);

return ret;

}







// Stuff
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
