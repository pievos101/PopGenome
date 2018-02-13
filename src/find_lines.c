#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP find_lines_GFF_Human2 (SEXP RRfilename, SEXP RRchromosome) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

// init
for(int xx =0; xx<2;xx++){
INTEGER(ret)[xx] = 0;
}

  //char ch;
  int ch;
  char *file;
  int line_count = 0;
  //int chr2;
  int raus=0;
  char *chr;
  int FOUND  = 0; 
  

  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);

  SEXP Rchromosome;

  Rchromosome   = STRING_ELT(RRchromosome,0);
  chr           = (char*)CHAR(Rchromosome);
  

  // printf("%c", *file);
 
  FILE *fp   = fopen( file ,"rt");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }

char ident[1000] = {""};
char temp[2] = {""};
//Rprintf("Identifier not found");

  while(1) {
    
    ch      = fgetc(fp);

    if(ch==EOF){break;} 

    if(ch=='\n'){
	
      line_count ++;
        
             while(1){	                 
	      ch      = fgetc(fp);
              //printf("%c", ch);

	      //
              if(ch=='#'){goto skip;}
 	      //
	             
	      if(ch==EOF && FOUND==0){ Rprintf("Identifier not found");goto end;}
	      if(ch==EOF && FOUND==1){ INTEGER(ret)[1]=line_count-1;goto end;}	
              //if(ch==EOF){ Rprintf("Identifier not found");goto end;}
	      if(ch=='\t'){break;}	
              temp[0] = ch;
	      strcat(ident,temp);
             }// end of while 

	     //printf("%c", *ident);
	     //break;		

	     // compare identifiers
             if(FOUND==0){	
             	if(strcmp(ident,chr)==0){ // identical  
              	INTEGER(ret)[0] = line_count;
              	FOUND = 1;
             	}
	     }else{	
             	if(strcmp(ident,chr)!=0){ // not identical  
                //printf("%c", *ident);
              	INTEGER(ret)[1] = line_count-1;
              	break;
             	}
	     }
     
            skip: // skip line
            strcpy(ident, ""); // delete array 
	    
    } // one line 



 } // End of while            

  //printf("%d",line_count);
  //printf("%d",chr2);

  end:
  fclose(fp);

  UNPROTECT(1);
 
  return ret;

}



SEXP find_lines_GFF (SEXP RRfilename, SEXP RRchromosome) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

  //char ch;
  int ch;
  char *file;
  int line_count = 1;
  int chr2;
  int raus=0;
  char *chr;
 
  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);

  SEXP Rchromosome;
  Rchromosome  = STRING_ELT(RRchromosome,0);
  chr          = (char*)CHAR(Rchromosome);



  // printf("%c", *file);
 
  FILE *fp   = fopen( file ,"rt");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }

  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='\n'){

      line_count ++;
       
             ch = fgetc(fp);

             if(ch=='c' || ch=='C'){  
              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);
             }

             if(ch==*chr){
               if(line_count==2){
               INTEGER(ret)[0] = 1; 
               break; 
               }
               INTEGER(ret)[0] = line_count;
               break;
             }        
    }
  }

 while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='\n'){

      line_count ++;
        
         

             ch = fgetc(fp);

             if(ch=='c'||ch=='C'){  
              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);
             }		

             if(ch!=*chr){
               INTEGER(ret)[1] = line_count - 1;
               break;
             }

               
    }
 }



  //printf("%d",line_count);
  //printf("%d",chr2);

  fclose(fp);

  UNPROTECT(1);
 
  return ret;

}

SEXP find_lines_SNP (SEXP RRfilename, SEXP RRchromosome) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

  //char ch;
  int ch;
  char *file;
  int line_count = 1;
  int chr2;
  int raus=0;
  char *chr;
 
  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);

  SEXP Rchromosome;
  Rchromosome  = STRING_ELT(RRchromosome,0);
  chr          = (char*)CHAR(Rchromosome);



  // printf("%c", *file);
 
  FILE *fp   = fopen( file ,"rt");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }

  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='\n'){

      line_count ++;
                         

             while((ch = fgetc(fp))!='\t');

             ch = fgetc(fp);
             
             if(ch=='c' || ch=='C'){  
              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);
             }

             if(ch==*chr){
               if(line_count==2){
               INTEGER(ret)[0] = 1; 
               break; 
               }
               INTEGER(ret)[0] = line_count;
               break;
             }

               
    }
  }

 while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='\n'){

        line_count ++;
             
            while(1){
              
             ch = fgetc(fp);

             if(ch=='\t'){break;}
             if(ch==EOF){
                INTEGER(ret)[1] = line_count - 1;
                goto MARKE;
             } 

            }
   
           //while( (ch = fgetc(fp)) !='\t');

           ch = fgetc(fp);
             
             if(ch=='c'||ch=='C'){  

              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);

             }		

             if(ch!=*chr){

               INTEGER(ret)[1] = line_count - 1;
               break;

             }

               
    }
 }

MARKE:

  //printf("%d",line_count);
  //printf("%d",chr2);

  fclose(fp);

  UNPROTECT(1);
 
  return ret;

}

SEXP find_lines_GFF_Human (SEXP RRfilename, SEXP RRchromosome) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

  //char ch;
  int ch;
  char *file;
  int line_count = 1;
  //int chr2;
  int raus=0;
  char *chr;
  char *chr2;
  int FOUND  = 0; 
  

  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);

  SEXP Rchromosome;
  SEXP Rchromosome2;

  Rchromosome   = STRING_ELT(RRchromosome,0);
  Rchromosome2  = STRING_ELT(RRchromosome,1);

  chr           = (char*)CHAR(Rchromosome);
  chr2          = (char*)CHAR(Rchromosome2);


  // printf("%c", *file);
 
  FILE *fp   = fopen( file ,"rt");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }

  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", *chr2);

    if(ch=='\n'){

      line_count ++;
                         
             ch = fgetc(fp);

             if(ch=='c'||ch=='C'){  
              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);
             }

             if(ch==*chr){

               if(line_count==2 && *chr2=='z'){
               INTEGER(ret)[0] = 1;
               FOUND = 1; 
               break; 
               }
               // check second position of chrXX
               if(*chr2=='z'){
                INTEGER(ret)[0] = line_count;
                FOUND = 1;
                break;
               }else{
                ch = fgetc(fp);
                //printf("%c", ch);
                if(ch==*chr2){
	        INTEGER(ret)[0] = line_count;
                FOUND=1;
                break;  	
                }
               }

             }

               
    }
  }

 // If the chromosome id was not found
 
 if(FOUND==0) {
    UNPROTECT(1);
    return R_NilValue;
 }



 while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='\n'){

      line_count ++;
              
             ch = fgetc(fp);

             if(ch=='c'||ch=='C'){  
              ch = fgetc(fp);
              ch = fgetc(fp);
              ch = fgetc(fp);
             }		

             if(ch!=*chr){
               INTEGER(ret)[1] = line_count - 1;
               break;
             }
            
             if(ch==*chr){               
                if(*chr2=='z'){
                  continue;
                }
                ch = fgetc(fp);
                if(ch!=*chr2){
                  INTEGER(ret)[1] = line_count - 1;
                  break;
                }
             }
               
    }
 }



  //printf("%d",line_count);
  //printf("%d",chr2);

  fclose(fp);

  UNPROTECT(1);
 
  return ret;

}


