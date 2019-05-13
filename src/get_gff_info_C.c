#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP get_gff_info_C (SEXP RRstart, SEXP RRend, SEXP RRfilename, SEXP RRposition) {

  SEXP info = R_NilValue;
  
  SEXP Rposition;
  Rposition   = coerceVector(RRposition, INTSXP);
  PROTECT(Rposition);
  SEXP Rstart;
  Rstart      = coerceVector(RRstart, INTSXP);
  PROTECT(Rstart);
  SEXP Rend;
  Rend        = coerceVector(RRend, INTSXP);  
  PROTECT(Rend);
  
  int  *cstart = INTEGER(Rstart);
  int  *cend   = INTEGER(Rend);
  int  *pos    = INTEGER(Rposition);
  char *file;

  int  ch;

  
  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);


  FILE *fp;
  fp = fopen( file ,"r");

PROTECT(info = allocVector(STRSXP,1));

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(4);
    return R_NilValue;
  }

// PROTECT(info = allocVector(STRSXP,1));

//   PROTECT(info = allocVector(INTSXP,2));
//   INTEGER(info)[0] = -1;
//   INTEGER(info)[1] = -1;

  
  int  count    = 0;
  int  count2   = 0;
  int  intcount = 0;
  char str[2000];
  char temp[10];

  int poss1;
  int poss2;

// Jump to the right chromosome
int jumpii = 1;
while(jumpii<cstart[0]){ 
  ch = fgetc( fp );
  if(ch=='\n'){
     jumpii++;
  }
}
// ----------------------------

int line_count = 0;

  while((ch = fgetc( fp )) != EOF) {
    
    if(ch=='\n'){
        
       line_count ++;
       // Wenn schon nächstes Chromosom break;   
       if((cstart[0]+line_count)>cend[0]){break;}  
       
       //-------------------------------------

       //jump to the start-position
         while(1){
            ch = fgetc( fp );
            if(ch==EOF){goto SCHLUSS;}              
           if(ch=='\t'){
             count ++;
           }
          if(count==3){count=0;break;}
         }
       // ------------------------     
 
              // save start pos as integer
              while((ch=fgetc(fp))!='\t'){
                // printf("%c", ch);
                 temp[0] = (char)ch;
                 strcat(str,temp);
                 //strcpy(str,""); 
                 //intcount++;   
              }
           //printf("-----------");
           //printf("%s",str);
            

            //UNPROTECT(1);
            //fclose(fp);
            //return(info);

            poss1 = atoi(str);
            //printf("%s",str);
            //printf("/");
            strcpy(str,"");

             // save end position
             while((ch=fgetc(fp))!='\t'){
                 // printf("%c", ch);
                 temp[0] = (char)ch;
                 strcat(str,temp);
                 //strcpy(str,""); 
                 //intcount++;   
              }

           
           poss2 = atoi(str);
           //printf("%s",str);
           //goto SCHLUSS; 
           strcpy(str,"");

          
            if((poss1 <= pos[0]) && (poss2 >= pos[0])){
               // hier infos zusammensammeln  
      	                    //jump to the info-position
     			    while(1){
            			ch = fgetc( fp );
            			// if(ch==EOF){goto SCHLUSS;}              
           			if(ch=='\t'){
             		            count2 ++;
                	        }
          			if(count2==3){count2=0;break;}
                              }
                            // ------------------------
                            // Get the informations
		            while((ch=fgetc(fp))!='\n'){
                 		// printf("%c", ch);
                 	        temp[0] = (char)ch;
                 		strcat(str,temp);
                 		//strcpy(str,""); 
                 		//intcount++;   
              	            }
                            // --------------------------- 
  
                            SET_STRING_ELT(info,0,mkChar(str));

	                    //INTEGER(info)[1] = poss2;
                            //INTEGER(info)[0] = poss1;

                            goto SCHLUSS;
                            //UNPROTECT(1);
               		    //fclose(fp);
               		    //return(info);

            }
  
           //printf("-----------");
           //printf("%s",str);
           
    } // END if Zeile zuende 	
	
 } // while bis EOF

SCHLUSS:

fclose(fp);
UNPROTECT(4);
return(info);

}//END OF FUNCTION


//           }else{
//             temp[count] = ch;
//             count ++;
//             strcat(str,temp);
//             SET_STRING_ELT(row_names,count_ind,mkChar(str));
//             strcpy(str, "");
//           }
//         }
//        
//    }

//  }
//##################################################


  // Count nucleotides
  
/*   while(1){
   ch = fgetc( fp );

     if(ch=='>' || ch==EOF){
      break;
     }
 
     if(ch!='\n'){      
        n_nucs ++;
     }
    
  // printf("%c", ch);
   } 

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
*/

