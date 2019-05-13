//#define _FILE_OFFSET_BITS		64

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <string.h>

//#define			PLATFORM_BITS		32

// Get individual names
SEXP getIndividuals(char *VCFbuffer, int ploidy){

SEXP row_names = R_NilValue;

char ind[1000];
int  x = 0;
int  n_ind=1; 
int  skip=0;
int y;

// Count individuals
while(1){
	if(VCFbuffer[x]=='\n'){
		if(VCFbuffer[x+1]=='#' && VCFbuffer[x+2]!='#'){
              
			while(1){
					if(VCFbuffer[x]=='\t'){
						skip++;
 						x++;
					}
					if(skip==9){
						y=x;
						// count individuals
						while(VCFbuffer[x]!='\n'){
							if(VCFbuffer[x]=='\t'){n_ind++;}
						x++;
                                                }
						goto end;					
					}
				x++;
                                } 
		//CHAR(ret)[0] = 'a';
                //CHAR(ret)[1] = 'b';
				                
		//printf("%c",VCFbuffer[x+1]);
                //UNPROTECT(1);                
		//return(ret);
		}
	x++;
	}
x++;
}
end:


// save individual names:
PROTECT(row_names = allocVector(STRSXP,ploidy*n_ind));
int count = 0;
char str[100]={""};
char temp[2]={""};
char aChar[2]={""};
char str2[100]={""};
char dott[2] ={"."};

while(1){

	if(VCFbuffer[y]!='\t' && VCFbuffer[y]!='\n'){
	     if(VCFbuffer[y]!='\r'){temp[0] = VCFbuffer[y];}
             strcat(str,temp);	
	     y++;	
	}else{
	     SET_STRING_ELT(row_names,count,mkChar(str));
	     // polyploid	
	     for(int xx=2; xx<=ploidy;xx++){
	     strcpy(str2, str);	
	     count++;	
	     aChar[0] = '0' + (xx);
	     strcat(str2,dott);
             strcat(str2,aChar);
	     SET_STRING_ELT(row_names,count,mkChar(str2));
	     strcpy(str2, "");	
             }
	     if(VCFbuffer[y]=='\n'){break;}	
       	     strcpy(str, "");
             count++;
             y++;	     	
	}
}
//printf("%d",n_ind);
UNPROTECT(1);
return(row_names);
} // end of function


// FUNCTION: Count SNPs ---------------------------------------------------------
int countSNPs(char *VCFbuffer, long int size){

 long int x=0;
 long int count=-1;

// Skip to first SNP
while(1){
	if(VCFbuffer[x]=='\n'){
		if(VCFbuffer[x+1]=='#' && VCFbuffer[x+2]!='#'){
                break;
		}
	}
x++;
}
//count SNPs
for(long int i = x+1; i <= size;i++){

//printf("%c",VCFbuffer[i]);
     if(VCFbuffer[i]=='\n'){
     count++;
     }
  }	
  //if(VCFbuffer[size-1]!='\n'){count++;}

return(count);
}//----------------------------------------------------- End of Function


// FUNCTION: Check Ploidy -----------------------------------------------
int checkPloidy(char *VCFbuffer, long int size){

long int count = 1;
// Check Ploidy from the end of the file
// Jump to the first TAB
while(VCFbuffer[size]!='\t'){
size--;
}

while(VCFbuffer[size]!=':' && VCFbuffer[size]!='\n'){
	if(VCFbuffer[size]=='/'||VCFbuffer[size]=='|'){
        count++;
	}
size++;
}
return(count);
}
// END OF FUNCTION // ------------------------------------------------------------------


SEXP myReadVCFC (SEXP RRfilename) {


// Nucleotide Mapping 
static char	nucleotide_mapping[] = {
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 5	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 16	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,6,5,5,				// 32	: !"#$%&'()*+’-./
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 48	:5123456789:;<=>?
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 64	:@ABCDEFGHIJKLMNO
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 80	:PQRSTUVWXYZ[\]^_
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 96	:`abcdefghijklmno
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 112	:pqrstuvwxyz{|}~
	//FIXME : U = T => 1 ?
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 128	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 144	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 160	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 176	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 192	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 208	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 224	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5					// 240	: - mostly nonprintable - 
};



  SEXP list = R_NilValue;
  //PROTECT(ret = allocVector(INTSXP,2));

// Get the size of the file
  char *file;
  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);
  FILE *fp;
  fp = fopen( file ,"rt");
  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    //UNPROTECT(1);
    return R_NilValue;
  }
//  printf("0");

  fseeko(fp,0,SEEK_END);
  /*unsigned int size*/ off_t size = ftello(fp);
// printf("1");
  rewind(fp);

//printf("size is: %d",size);
//printf("2");

// Read file into buffer
char *VCFbuffer;
size_t result;
VCFbuffer = calloc( 1, size+1 );
if (VCFbuffer == NULL) {Rprintf("Memory error !!! ");Rprintf("\n");
//UNPROTECT(1); 
return R_NilValue;
}

//printf("3");

result=fread(VCFbuffer , size, 1, fp);//FIXME warnings ...

//printf("4");

fclose(fp);

//
// Count SNPs
//
long int n_snps = countSNPs(VCFbuffer,size);
//printf("n.biallelic.sites: %d",n_snps);
//printf("\n");

SEXP positions = allocVector(INTSXP,n_snps);
PROTECT(positions);

//printf("5");

//printf("N.SNPS=: %d",n_snps);
//
// Check Ploidy
//
int ploidy = checkPloidy(VCFbuffer,size);
//printf("Ploidy=: %d",ploidy);

//printf("6");

//
// Get individuals
//
SEXP row_names;
row_names = getIndividuals(VCFbuffer, ploidy);
PROTECT(row_names);

//printf("7");
//
// Init Matrix
//
SEXP matrix;
int I = length(row_names);
int J = n_snps;
PROTECT(matrix = allocMatrix(INTSXP,I,J));
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
   INTEGER(matrix)[j +I*i]=6;
 }
}
//
// Fill Matrix ------------------------------------------
//
//printf("8");

//skip to first individual
int nucmap[10];
long int x = 0;
while(1){
	if(VCFbuffer[x]=='\n' && VCFbuffer[x+1]!='#'){
        break;
        }
x++;
}
// iterate over snps
int tabcount = 0;
long int count;
char str[100]={""};
char temp[2]={""};

//printf("9");

for(long int i=0; i<J;i++){

tabcount = 0;
x = x+1; // skpip the End of line \n

		// SNP line
		while(VCFbuffer[x]!='\n'){
			if(VCFbuffer[x]=='\t'){
                        tabcount++;
                        }
			if(tabcount==1){ // Save SNP position 
			x++;
				while(VCFbuffer[x]!='\t'){
				temp[0] = VCFbuffer[x];
             			strcat(str,temp);
				x++;
				}
				tabcount++;

			INTEGER(positions)[i]= atoi(str);
			strcpy(str, "");
			}
			if(tabcount==3){ // REF
 			x = x+1;
			//printf("REF:%c",VCFbuffer[x]); 
                        nucmap[0] = nucleotide_mapping[(int)VCFbuffer[x]];
			// check if REF is an insertion, if yes, set to unknknown 5
				if(VCFbuffer[x+1]!='\t'){
					while(VCFbuffer[x]!='\t'){x++;}
					x--; // x is one position before \t
				nucmap[0] = 5 ;
				}	   
                        }
 			if(tabcount==4){ // ALT
			x = x+1;
			count = 1;
				while(VCFbuffer[x]!='\t'){
                                	if(VCFbuffer[x]!=','){
                                        	// if insertion skip to next
						if(VCFbuffer[x+1] != '\t' && VCFbuffer[x+1] != ','){
					 	 nucmap[count] = 5;
						 count ++;
 						 while(VCFbuffer[x]!=',' && VCFbuffer[x]!='\t'){x++;} 
						 x--;	                                        	
						}else{
						//printf("ALT:%c",VCFbuffer[x]);
					 	 nucmap[count] = nucleotide_mapping[(int)VCFbuffer[x]];
					 	 count++;
						}	
					}
				x++;					
				}
				tabcount++;
			}
                        if(tabcount==9){ // Individuals
			//printf("*%c",VCFbuffer[x]);
			//printf("---");break;
			 for(long int j=0; j < I; j++){
				while(VCFbuffer[x]!='\t' && VCFbuffer[x]!='\n'){x++;} //skip Tabs
				if(VCFbuffer[x]=='\n'){break;}
				x= x+1;
					//fill matrix one column
					//printf("%c",VCFbuffer[x]);
					for(int p=0; p< ploidy ;p++){
					 if(VCFbuffer[x + 2*p]=='.'){
					 INTEGER(matrix)[(ploidy*j+p) +I*i] = 5; //unknown position
					 }else{	
					 INTEGER(matrix)[(ploidy*j+p) +I*i] = nucmap[VCFbuffer[x + 2*p]-'0'];
 					 } //printf("--------------%d",val);
					}
			 }
			}
		if(VCFbuffer[x]=='\n'){break;}
		x++; // increase Buffer incr x
		}// End of while 

}

//for (int x=0; x<3;x++){
//printf("/%d",nucmap[x]);
//}

   // Creating a list with 3 vector elements:
   PROTECT(list = allocVector(VECSXP, 3)); 
   SET_VECTOR_ELT(list, 0, matrix); 
   SET_VECTOR_ELT(list, 1, positions); 
   SET_VECTOR_ELT(list, 2, row_names); 
   UNPROTECT(4);


free(VCFbuffer);
//printf("Test: %c",CHAR(test)[1]);

return(list);
}

SEXP pimpMatrix(SEXP RinMatrix, SEXP RinMatrix2){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;   // input matrix
SEXP Rvalue2;  // output matrix

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

//Rvalue        = coerceVector(RinMatrix, INTSXP);
int    *Rval    = INTEGER(RinMatrix);
//Rvalue2       = coerceVector(RinMatrix2, INTSXP);
int    *Rval2   = INTEGER(RinMatrix2);

PROTECT(ret = allocVector(INTSXP,1));
INTEGER(ret)[0] = 0;

int value;
//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){

  // value = (Rval[j + I*i]);
  // printf("%f",value);

   Rval2[2*j + 2*I*i]     = (int)(Rval[j + I*i]/10);
   Rval2[(2*j+1) + 2*I*i] = Rval[j + I*i] % 10;

  
  
 }
}

INTEGER(ret)[0] = 1;
UNPROTECT(1);
return ret;

}


