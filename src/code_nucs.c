#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP code_nucs(SEXP RinMatrix){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten


Rvalue                        = coerceVector(RinMatrix, STRSXP);
PROTECT(Rvalue);

//currentfilename = STRING_ELT(filenames,i);
//filename = (char*)CHAR(currentfilename);

//const char *Rval              = CHAR(Rvalue);

PROTECT(ret = allocMatrix(INTSXP,I,J));

// Init ret
// Init Matrix

/*
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
   INTEGER(ret)[j +I*i]=5;
 }
}
*/

char	nucleotide_mapping[] = {

	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 5	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 16	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,6,5,5,				// 32	: !"#$%&'()*+’-./
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 48	:5123456789:;<=>?
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 64	:@ABCDEFGHIJKLMNO
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 80	:PQRSTUVWXYZ[\]^_
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 96	:`abcdefghijklmno
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 112	:pqrstuvwxyz{|}~

	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 128	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 144	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 160	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 176	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 192	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 208	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 224	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5					// 240	: - mostly nonprintable - 

};


SEXP val;
char *value;

for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){

   val   =  STRING_ELT(Rvalue,j+I*i);  
   value =  (char*)CHAR(val);

   INTEGER(ret)[j+I*i] = nucleotide_mapping[(int)*value];

  /* if(*value=='A'||*value=='a'){INTEGER(ret)[j +I*i]=4;continue;}
   if(*value=='C'||*value=='c'){INTEGER(ret)[j +I*i]=2;continue;}
   if(*value=='G'||*value=='g'){INTEGER(ret)[j +I*i]=3;continue;}
   if(*value=='U'||*value=='u'||*value=='T'||*value=='t'){INTEGER(ret)[j +I*i]=1;continue;}
   //if(*value=='N'||*value=='n'||*value=='?'){INTEGER(ret)[j +I*i]= 5;continue;}
   if(*value=='-'){INTEGER(ret)[j +I*i]= 6;continue;}
   
 */


 }
}

UNPROTECT(2);

return ret;

}


