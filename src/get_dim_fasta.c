#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP get_dim_fasta (SEXP RRfilename) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

  char *file;
  int  ch;
  int  n_ind=0;

  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);


  FILE *fp;
  fp = fopen( file ,"r");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }


//Rprintf("Start count individuals \n");

  // check rows of ff matrix
  // ##################################
  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='>'){
      n_ind ++;     
    }

  }

   INTEGER(ret)[0] = n_ind;

//Rprintf("Individuals checked \n");

  // ##################################


// SAVING NAMES
// ######################################

 rewind(fp);

//Rprintf("Start get individual names \n");

 SEXP row_names;

 PROTECT(row_names = allocVector(STRSXP,n_ind));
 
  char str[100]={""};
  char temp[2]={""};
       temp[1]='\0';

  int count;
  int count_ind=-1;

  while((ch = fgetc( fp )) != EOF) {
    
    if(ch=='>'){
         count_ind ++;
         count = 0;
         while(1){
            ch = fgetc( fp );              
           if(ch=='\n'|| ch=='\t' || ch=='\r'){            
              break;
           }else{
             temp[0] = ch;
             strcat(str,temp);                  
           }
         }              
	 SET_STRING_ELT(row_names,count_ind,mkChar(str));
         strcpy(str,""); 	            
    }
 }

//Rprintf("Individual names checked \n");

//##################################################
// check columns of ff matrix
// ############################################
  rewind(fp); 

//Rprintf("Start count Nucleotides \n");

  int n_nucs  = 0;
 
   // Get first nucleotide   
   while(1){
     ch = fgetc( fp );
     //printf("%c", ch);
     if(ch=='\n'|| ch=='\t'){
        break;
     }
   }
   
   // Count nucleotides
   while(1){
   ch = fgetc( fp );

     if(ch=='>' || ch==EOF){
      break;
     }

     if((ch!='\n') && (ch=='A' || ch=='a'|| ch=='C' || ch=='c'|| ch=='G'|| ch=='g' || ch=='U'|| ch=='u'|| ch=='T' || ch=='t' || ch=='N'||ch=='n'||ch=='?' || ch=='-')){      
        n_nucs ++;
     }
    
  // printf("%c", ch);
   } 

//Rprintf("Nucleotide number checked \n");

  INTEGER(ret)[1] = n_nucs; 
  fclose(fp);
//###################################################

SEXP list;
PROTECT(list = allocVector(VECSXP, 2));
SET_VECTOR_ELT(list, 0, ret);    
SET_VECTOR_ELT(list, 1, row_names); 

UNPROTECT(3);
return list;

}
SEXP FASTA_open(SEXP name, SEXP mode, SEXP fun)

{
static SEXP FILE_type_tag;
    FILE *f = fopen(CHAR(STRING_ELT(name, 0)), CHAR(STRING_ELT(mode, 0)));
    if (f == NULL){
	Rprintf("Cannot open file");
        return R_NilValue;
    }else{ 
	SEXP val = R_MakeExternalPtr(f,FILE_type_tag, R_NilValue);
        //R_RegisterCFinalizer(val, (R_CFinalizer_t) FASTA_close);  
	//R_RegisterCFinalizerEx(res, vcff_finalize, Rboolean_TRUE) 
        return val;
    }
}

SEXP FASTA_getNextIndividual(SEXP s, SEXP RRn_nucs){

int ch;
SEXP ret;
SEXP Rn_nucs;
Rn_nucs     = coerceVector(RRn_nucs, INTSXP);
int *n_nucs = INTEGER(Rn_nucs);
PROTECT(ret = allocVector(INTSXP,n_nucs[0]));

FILE *f;
f = R_ExternalPtrAddr(s);

// Skip to next individual
while(1){
ch = fgetc(f);
if(ch=='>' || ch==EOF){break;}
}
while(1){
ch = fgetc(f);
if(ch=='\n'||ch=='\r'||ch=='\t'||ch==EOF){break;}
}
// -----------------------

int xx = 0;
while((xx< n_nucs[0]) && (ch!=EOF)){
 ch = fgetc(f);
 //printf("%c", ch);
 if(ch=='A'||ch=='a'){INTEGER(ret)[xx]=4;xx++;continue;}
 if(ch=='C'||ch=='c'){INTEGER(ret)[xx]=2;xx++;continue;}
 if(ch=='G'||ch=='g'){INTEGER(ret)[xx]=3;xx++;continue;}
 if(ch=='U'||ch=='u'||ch=='T'||ch=='t'){INTEGER(ret)[xx]=1;xx++;continue;}
 if(ch=='N'||ch=='n'||ch=='?'){INTEGER(ret)[xx]= 5;xx++;continue;}
 if(ch=='-'){INTEGER(ret)[xx]= 6;xx++;continue;}
}

UNPROTECT(1);
return ret;
}// End of Function 


SEXP FASTA_close(SEXP s)
{
    FILE *f;
   // CHECK_FILE_STREAM(s);
    f = R_ExternalPtrAddr(s);
    if (f != NULL) {
        fclose(f);
        R_ClearExternalPtr(s);
    }
    return R_NilValue;
}

