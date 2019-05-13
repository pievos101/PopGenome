#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP get_ind_fasta (SEXP RRfilename,SEXP RRind, SEXP RRn_nucs) {

  SEXP ret    = R_NilValue;
  SEXP Rind;
  SEXP Rn_nucs;

  Rind        = coerceVector(RRind, INTSXP);
  PROTECT(Rind);
  Rn_nucs     = coerceVector(RRn_nucs, INTSXP);
 
  int *n_nucs = INTEGER(Rn_nucs);
  int *ind    = INTEGER(Rind);
  
  PROTECT(ret = allocVector(INTSXP,n_nucs[0]));

// Init 

for(int x = 0; x < n_nucs[0]; x++){
    INTEGER(ret)[x]=6;
} 

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
    UNPROTECT(2);
    return R_NilValue;
  }

  // check rows of ff matrix
  // ##################################
  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='>'){

      n_ind ++;

      if(n_ind==ind[0]){

        // Get first nucleotide

          while(1){
            ch = fgetc( fp );
           //printf("%c", ch);
            if(ch=='\n'|| ch=='\t'){
               break;
            }
          }           
 
         for(int xx=0; xx < n_nucs[0]; xx++){
                 ch = fgetc( fp );
                 if(ch=='\n'){ch=fgetc( fp );}
                 if(ch=='A'||ch=='a'){INTEGER(ret)[xx]=4;continue;}
  		 if(ch=='C'||ch=='c'){INTEGER(ret)[xx]=2;continue;}
   		 if(ch=='G'||ch=='g'){INTEGER(ret)[xx]=3;continue;}
   		 if(ch=='U'||ch=='u'||ch=='T'||ch=='t'){INTEGER(ret)[xx]=1;continue;}
   		 if(ch=='N'||ch=='n'||ch=='?'){INTEGER(ret)[xx]= 5;continue;}
   		 if(ch=='-'){INTEGER(ret)[xx]= 6;continue;}
         }
       break;
      }
    }

  }

fclose(fp);

// ##################################

UNPROTECT(2);
return ret;

}

