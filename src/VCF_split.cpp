/*
 *
 *  Created on: 10.02.2015
 *      Author: bastian
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace std;

// Get Contig of the current line
string getContig(string  l) {

	char first = l[0];

	if(first=='#'){
		//printf("1");
	 	return "#";
	}

	int i=0;
	string contig;
    char currentChar;

	while(1){
		// printf("uff");
		// cout << l[2] << endl;
		currentChar = l[i];
		// cout << currentChar << endl;
		if(currentChar=='\t'){
			//cout << l[i+1] << endl;
			//printf("Treffer");
			break;
		}
		i++;

	contig = contig + currentChar;
	}
	return contig;
}


//void splitData(const char *inputFile){
extern "C" SEXP split_VCF_scaffolds(SEXP RRfilename, SEXP RRfolder){

	SEXP    ret = R_NilValue;
	PROTECT(ret = allocVector(INTSXP,1));
	INTEGER(ret)[0] = 1;	

	SEXP Rfilename;
	Rfilename   = STRING_ELT(RRfilename,0);
        char *inputFile;
	inputFile   = (char*)CHAR(Rfilename);

	//OutputFolder
	SEXP Rfolder;
	Rfolder     = STRING_ELT(RRfolder,0);
        char *outputFolder;
	outputFolder   = (char*)CHAR(Rfolder);
        string s_outputFolder = outputFolder;
 
        //mkdir(outputFolder);

	// input file
	 //   cout << inputFile << endl;
	    ifstream iFile(inputFile);
	    if(iFile){
	    	Rprintf("Success open Input file \n");
	    }else{
	    	Rprintf("Error::opening Input file");
                UNPROTECT(1);
                return(ret);
//	    	exit(1);
	    }

	    // output file

	    ofstream myfile;

	    //myfile.close(); // Closing file

	    string line;
	    string s;
	    string header;
	    string ssave;

	    // Get Header information from VCF
	    while(getline(iFile, line)) {
	            s  = getContig(line);
	            //cout << line << endl;
	            if(s=="#"){
	            header = header +  line + "\n";
	            }else{break;}
	    }

	    // Get the very first Contig and open new file
	    ssave = getContig(line);
	    string outfilebase = inputFile;
	    string outfilefirst;
	    outfilefirst =  s_outputFolder + ssave; //+ ".vcf";//+ outfilebase + ssave;
	    //cout << ssave << endl;
	    //cout << outfilefirst << endl;


	    const char * c = outfilefirst.c_str();
	    myfile.open(c);

	    if(myfile){
	         	Rprintf("Success open Output file \n");
	    }else{
	           	Rprintf("Error::opening output file");
			UNPROTECT(1);
			return(ret);
	           	//exit(1);
	    }
	    // -----------------------------------------------

	    string outfile;
	    int count = 1;
	    myfile << header << endl;
	    //myfile << "\n";
	    myfile << line  << endl;
	    //myfile << "\n";
	    /* Iterate over SNPs within the Contigs */
	    while(getline(iFile, line)) {
	        s  = getContig(line);
	        if(s==ssave){
	        //cout << line << endl;
	        myfile << line << endl;
	        //myfile << "\n"; //endl;
	        }else{
	        	myfile.close();
	        	ssave = s;
	        	//open new file and write infos
	        	outfile =   s_outputFolder + ssave; //+ ".vcf";//outfilebase + ssave ; 
	        	//cout << outfile << endl;
	        	myfile.open (outfile.c_str());
	        	myfile << header << endl;
	        	//myfile << "\n";
	        	myfile << line<< endl;
	        	//myfile << "\n";
	        	count ++;
	        	//if(count==2){exit(1);}
	        }
	    }
	    iFile.close();
	   // cout << "Finished" << endl;

UNPROTECT(1);
return(ret);
}



