/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		Whopgen - Reading data into matrices
**
**
**

	TODO:



**
*/

//*
//*			INCLUDES
//*

#include	"whopgen_common.h"




	//

using namespace std;


//*
//*			DEFINES
//*



//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*



//*
//*			DATA
//*



SEXP	nucleotide_A	=	R_NilValue;
SEXP	nucleotide_C	=	R_NilValue;
SEXP	nucleotide_G	=	R_NilValue;
SEXP	nucleotide_T	=	R_NilValue;
SEXP	nucleotide_N	=	R_NilValue;	//for missing fields


//*
//*			EXTERNS
//*



//*
//*			CODE
//*

/*

		ids <- c("T" ,"t",	"U","u",	"C","c",	"G","g",	"A","a",	"N","n","?",	"-")
		nuks <- c(1,1,		1,1,		2,2,    	3,3,		4,4,    	 5,5,5,			6)
*/

//!	For quick mapping of nukleotide character codes into numeric constants for the biallelic matrix
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

/*


bld(); v= .Call("VCF_open","/media/data715/ALL.chr1.phase1.projectConsensus.genotypes.vcf.gz"); .Call("VCF_selectSamples",v,c("HG00096","HG00098")); ndat = matrix(nrow=2,ncol=2,rep("A",4) )


.Call("VCF_exportAsFasta",v,ndat,0,0,0)

*/


/*!	Reads a coding matrix for all biallelic SNPs of the given samples
**
**	- vcff : a VCFhandle
**
**
*/
EXPORT	SEXP	VCF_readIntoCodeMatrix( SEXP vcfptr, SEXP mat )
{
	// statistics
	int	nonbialcols=0,
		bialcols=0;

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 1 is not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("VCF_readIntoCodeMatrix :: VCF does not appear to have a FORMAT field!\n");
		return RBool::False();
	}
	
	//
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		return RBool::False();
	}
	
	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		return RBool::False();
	}
	
	//	enough rows in matrix for all samples ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("VCF_readIntoCodeMatrix :: No samples selected!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("VCF_readIntoCodeMatrix :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return RBool::False();
	}

	//
	//
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
	{
		Rprintf("WhopGenome::VCF_readIntoCodeMatrix : WARNING : matrix has no column names vector! Cannot set SNP positions in matrix!\n");
		return RBool::False();
	}

	//-
	//-
	//-		parse VCF and fill biallelic matrix
	//-
	//-
	//-

	unsigned int	ncol = m.numCols();
	int*		ptr = m.getIntPtr();
	
	//
	char		*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//when too little data exists
	
	unsigned int	column_stepsize = nrow;
	
	int			snppos=-1;//unused
	SEXP			minus1_char = mkChar("-1");
        PROTECT(minus1_char);
	
//	df1("ncol=%d, nrow=%d, wanted=%d\n",ncol,nrow,f->num_wanted_samples);

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		const char* refptr;
		const char* altptr;
		bool	bLineParsed = f->parseNextLine();
		//
		while( bLineParsed )
		{

			//
			//	- make sure its a bi-allelic SNP line
			//
			refptr = (char*)f->getFieldPtr( REF );
			altptr = (char*)f->getFieldPtr( ALT );
			if( refptr && refptr[1] == '\t' )	
			{
				if( altptr && altptr[1] == '\t' )
				{
					break;
				}
				//else
				//{
				//	//df1("Not a biallelic SNP! (REF=%9s)\n",refptr);
				//}
			}
			//else
			//{
			//	//df1("Not a SNP! (REF=%2s)\n",refptr);
			//}

			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		//	if could not read another line from VCT, exit with error
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------

		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( fieldptr )
		{
			snppos = atoi( fieldptr );
			if( snppos == 0 )
			{
				df1("VCF_readIntoCodeMatrix :: SNPpos=%d\n",snppos);
			}
		}

		//-
		//
		//	Identify the GT-colon-field used in the per-individual fields
		//
		//-

		//
		//	- decompose FORMAT field
		//		OPT: skip buffer-copying it
		//		? what to do if missing ??
		//	- find GT in FORMAT
		//		! break with error if not found!
		//	- memorise field-pos of GT
		//
		fieldptr = (char*)f->getFieldPtr( FORMAT );
		int GTidx=0;
		int i=0;
		for( ; fieldptr[i]!=0 && fieldptr[i]!='\t';i++ )
		{
			if( (fieldptr[i]=='G') && (fieldptr[i+1]=='T') && (fieldptr[i+2]==':'||fieldptr[i+2]=='\t'||fieldptr[i+2]==0) )
				break;
			if( fieldptr[i]==':' )
				GTidx++;
		}

		//
		if( fieldptr[i]==0 || fieldptr[i]=='\t' )
		{
			df0("VCF_readIntoCodeMatrix :: NO GT FIELD DEFINED!\n");
			UNPROTECT(1);
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------


		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-

		bool	bHadRef = false,
				bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			
			//
			if( fieldptr == 0 )
			{
				Rprintf("VCF_readIntoCodeMatrix ::  Problem with reading sample's data!\n");
				Rprintf("	debug info : per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				UNPROTECT(1);
				return RBool::False();
			}

			//	- find the GT subfield (skip ':' until found)
			//
			while( GTidx > 0 )
			{
				if( fieldptr[0] == ':' )
					GTidx--;
				fieldptr++;
			}
			
			//
			//	- parse the GT subfield ( regexp: [0-9]+[/|][0-9]+ )
			//
			int left_allele = fieldptr[0]  - '0';
			int right_allele = fieldptr[2] - '0';
			//assert( left_allele >= 0 && left_allele <= 9 )
			//assert( right_allele >= 0 && right_allele <= 9 )
			
			//
			if( (fieldptr[1] != '|' && fieldptr[1] != '/') || (fieldptr[3] != '\t' && fieldptr[3] != ':') )
			{
				df0("VCF_readIntoCodeMatrix :: Malformed GT field!\n");
				UNPROTECT(1);
				return RBool::False();
			}

			// any chromosome has alt allele -> result is alt allele
			//
			///char	actual_nucleotide;
			if( ((left_allele==1)||(right_allele==1)) )
			{
				bHadAlt=true;
				ptr[per_row] = nucleotide_mapping[(unsigned)*altptr];
			}
			else
			{
				bHadRef=true;
				ptr[per_row] = nucleotide_mapping[(unsigned)*refptr];
			}

			//
		}
		
		//
		//	clear any unused matrix rows
		//
		if( (bHadRef&&bHadAlt) )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = -2;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(posbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{			
			//
			per_column--;
			nonbialcols++;
		}

	}//...for( each column in the matrix )
	
	//	Reset/clear unused columns
	//
	for( unsigned int rcol = per_column; rcol < ncol; rcol ++ )
	{
		//
		for( unsigned int per_row=0; per_row < nrow ; per_row ++ )
		{
			ptr[per_row] = -2;
		}
		ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != colnamvec )
		{
			SET_STRING_ELT( colnamvec, rcol, minus1_char );
		}
	}
	
	//	print some statistics
	//
	df1("VCF_readIntoCodeMatrix ::\n\t%d nonbial columns\n",nonbialcols);
	df1("\t%d bial columns\n",bialcols);
	df1("\t%d total columns\n",bialcols+nonbialcols);
	UNPROTECT(1); //minus1_char
	//
	return RBool::True();
}//...

EXPORT	SEXP	VCF_readIntoCodeMatrixdiploid( SEXP vcfptr, SEXP mat )
{
	// statistics
	int	nonbialcols=0,
		bialcols=0;

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 1 is not a VCFhandle EXTPTR!\n");
		
		return RBool::False();
	}

	//
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("VCF_readIntoCodeMatrix :: VCF does not appear to have a FORMAT field!\n");
		
		return RBool::False();
	}
	
	//
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		
		return RBool::False();
	}
	
	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		
		return RBool::False();
	}
	
	//	enough rows in matrix for all samples ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("VCF_readIntoCodeMatrix :: No samples selected!\n");
		
		return RBool::False();
	}

	//
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("VCF_readIntoCodeMatrix :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		
		return RBool::False();
	}

	//
	//
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
	{
		Rprintf("WhopGenome::VCF_readIntoCodeMatrix : WARNING : matrix has no column names vector! Cannot set SNP positions in matrix!\n");
		
		return RBool::False();
	}

	//-
	//-
	//-		parse VCF and fill biallelic matrix
	//-
	//-
	//-


	unsigned int	ncol = m.numCols();
	int*		 	ptr = m.getIntPtr();
	
	//
	char		*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	
	unsigned int	column_stepsize = nrow;
	
	int				snppos=-1;//unused
	SEXP			minus1_char = mkChar("-1");
        PROTECT(minus1_char);
	
//	df1("ncol=%d, nrow=%d, wanted=%d\n",ncol,nrow,f->num_wanted_samples);

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		const char* refptr;
		const char* altptr;
		bool	bLineParsed = f->parseNextLine();
		//
		while( bLineParsed )
		{

			//
			//	- make sure its a bi-allelic SNP line
			//
			refptr = (char*)f->getFieldPtr( REF );
			altptr = (char*)f->getFieldPtr( ALT );


			if( refptr && refptr[1] == '\t' )	
			{
				if( altptr && altptr[1] == '\t' )
				{
					break;
				}
				//else
				//{
				//	//df1("Not a biallelic SNP! (REF=%9s)\n",refptr);
				//}
			}
			//else
			//{
			//	//df1("Not a SNP! (REF=%2s)\n",refptr);
			//}

			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		//	if could not read another line from VCT, exit with error
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------

		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( fieldptr )
		{
			snppos = atoi( fieldptr );
			if( snppos == 0 )
			{
				df1("VCF_readIntoCodeMatrix :: SNPpos=%d\n",snppos);
			}
		}

		//-
		//
		//	Identify the GT-colon-field used in the per-individual fields
		//
		//-

		//
		//	- decompose FORMAT field
		//		OPT: skip buffer-copying it
		//		? what to do if missing ??
		//	- find GT in FORMAT
		//		! break with error if not found!
		//	- memorise field-pos of GT
		//
		fieldptr = (char*)f->getFieldPtr( FORMAT );
		int GTidx=0;
		int i=0;
		for( ; fieldptr[i]!=0 && fieldptr[i]!='\t';i++ )
		{
			if( (fieldptr[i]=='G') && (fieldptr[i+1]=='T') && (fieldptr[i+2]==':'||fieldptr[i+2]=='\t'||fieldptr[i+2]==0) )
				break;
			if( fieldptr[i]==':' )
				GTidx++;
		}

		//
		if( fieldptr[i]==0 || fieldptr[i]=='\t' )
		{
			df0("VCF_readIntoCodeMatrix :: NO GT FIELD DEFINED!\n");
			UNPROTECT(1);
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------


		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-

		bool	bHadRef = false,
				bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			
			//
			if( fieldptr == 0 )
			{
				Rprintf("VCF_readIntoCodeMatrix ::  Problem with reading sample's data!\n");
				Rprintf("	debug info : per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				UNPROTECT(1);
				return RBool::False();
			}

			//	- find the GT subfield (skip ':' until found)
			//
			while( GTidx > 0 )
			{
				if( fieldptr[0] == ':' )
					GTidx--;
				fieldptr++;
			}
			
			//
			//	- parse the GT subfield ( regexp: [0-9]+[/|][0-9]+ )
			//
			int left_allele = fieldptr[0] - '0';
			int right_allele = fieldptr[2] - '0';
			
						
			
			//assert( left_allele >= 0 && left_allele <= 9 )
			//assert( right_allele >= 0 && right_allele <= 9 )
			
			//
			if( (fieldptr[1] != '|' && fieldptr[1] != '/') || (fieldptr[3] != '\t' && fieldptr[3] != ':') )
			{
				df0("VCF_readIntoCodeMatrix :: Malformed GT field!\n");
				UNPROTECT(1);
				return RBool::False();
			}

			// any chromosome has alt allele -> result is alt allele
			//
			///char	actual_nucleotide;
			// diploid
			if( ((left_allele==1)||(right_allele==1)) )
			{
				bHadAlt=true;
				//ptr[per_row] = nucleotide_mapping[(unsigned)*altptr];
				if((left_allele==1) && (right_allele==0)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*altptr]+nucleotide_mapping[(unsigned)*refptr];
				}
				if((left_allele==0) && (right_allele==1)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*refptr]+nucleotide_mapping[(unsigned)*altptr];
				}
				if((left_allele==1) && (right_allele==1)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*altptr]+nucleotide_mapping[(unsigned)*altptr];
				}
					

			}
			else
			{
				bHadRef=true;
				//ptr[per_row] = nucleotide_mapping[(unsigned)*refptr];
				ptr[per_row]   = 10*nucleotide_mapping[(unsigned)*refptr]+nucleotide_mapping[(unsigned)*refptr];
			}

			//
		}
		
		//
		//	clear any unused matrix rows
		//
		if( (bHadRef&&bHadAlt) )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = 77;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(posbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{			
			//
			per_column--;
			nonbialcols++;
		}

	}//...for( each column in the matrix )
	
	//	Reset/clear unused columns
	//
	for( unsigned int rcol = per_column; rcol < ncol; rcol ++ )
	{
		//
		for( unsigned int per_row=0; per_row < nrow ; per_row ++ )
		{
			ptr[per_row] = 77;
		}
		ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != colnamvec )
		{
			SET_STRING_ELT( colnamvec, rcol, minus1_char );
		}
	}
	
	//	print some statistics
	//
	df1("VCF_readIntoCodeMatrix ::\n\t%d nonbial columns\n",nonbialcols);
	df1("\t%d bial columns\n",bialcols);
	df1("\t%d total columns\n",bialcols+nonbialcols);
	UNPROTECT(1);
	//
	return RBool::True();
}//...








#define	DBG(s)	Rprintf(#s);


/*!	
**
**
**	- vcfptr:	a VCFhandle
**	- mat	:	a matrix to receive the "A","C","G","T" entries
**
**
*/
EXPORT	SEXP	VCF_readIntoNucleotideMatrix( SEXP vcfptr, SEXP mat )
{
	//
	//
	int	nonbialcols=0,bialcols=0;

	//
	if( nucleotide_A == R_NilValue )		nucleotide_A = Rf_mkChar("A");
	if( nucleotide_C == R_NilValue )		nucleotide_C = Rf_mkChar("C");
	if( nucleotide_G == R_NilValue )		nucleotide_G = Rf_mkChar("G");
	if( nucleotide_T == R_NilValue )		nucleotide_T = Rf_mkChar("T");
	if( nucleotide_N == R_NilValue )		nucleotide_N = Rf_mkChar("N");

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df1("WhopGenome::VCF_readIntoNucleotideMatrix : Parameter 1 not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : VCF does not have a FORMAT field!\n");
		return RBool::False();
	}
	
	//
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : Parameter 2 not a (string) matrix!\n");
		return RBool::False();
	}

	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != STRSXP )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : Parameter not a string matrix!\n");
		return RBool::False();
	}

	//	enough rows in matrix for all samples ?
	//
	if( f->num_wanted_samples < 1 )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : No samples selected!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return RBool::False();
	}
	
	//-
	//-
	//-		parse VCF and fill biallelic matrix
	//-
	//-
	//-
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
	{
		df0("WhopGenome::VCF_readIntoNucleotideMatrix : WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");
		return RBool::False();
	}


	unsigned int	ncol = m.numCols();
	char			*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists	
	int				snppos=-1;
	bool			bHadRef = false,	//these two bools are used to find out whether the matrix' columns
					bHadAlt = false;	//	contain both alleles (in biallelic SNP lines from the VCF
										//		the set of selected samples may all have the same allele )

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{

		//
		//	Get a valid biallelic SNP line from the file
		//

		const char*	refptr;
		const char*	altptr;
		bool			bLineParsed = f->parseNextLine();
		//
		while( bLineParsed )
		{
			
			//
			//	- make sure its a bi-allelic SNP line
			//
			refptr = (char*)f->getFieldPtr( REF );
			altptr = (char*)f->getFieldPtr( ALT );
			if( refptr && refptr[1] == '\t' )	
			{
				if( altptr && altptr[1] == '\t' )
				{
					//@TODO : check for non-standard characters as ref,alt (besides A,C,T,G,N)
					break;	//both ref and alt allele are single characters (=> biallelic SNP)
				}
				else
				{
					//df0("Not a biallelic SNP! (REF=%9s)\n",refptr);
				}
			}
			else
			{
				//df0("Not a SNP! (REF=%2s)\n",refptr);
			}
			
			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		//	if could not read another line from VCT, exit with error
		//
		if( bLineParsed == false )
		{
			df0("WhopGenome::VCF_readIntoNucleotideMatrix : No more lines!\n");
			break;
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------
		
		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( fieldptr )
			snppos = atoi( fieldptr );
		//@TODO : what to do if getFieldPtr( POS ) fails ? assume entire line invalid!
		//Rprintf("POS(%d)\n",snppos);

		//-
		//
		//	Identify the GT-colon-field used in the per-individual fields
		//NOTE: VCF actually requires the GT field to be the first, but this code allows deviations from the standard
		//NOTE:		-> maybe just check for FORMAT =~ /^GT\:/  and get genotype from first sample field?
		//-

		//
		//	- decompose FORMAT field
		//		OPT: skip buffer-copying it
		//		? what to do if missing ??
		//	- find GT in FORMAT
		//		! break with error if not found!
		//	- memorise field-pos of GT
		//
		fieldptr = (char*)f->getFieldPtr( FORMAT );
		int GTidx=0;
		int i=0;
		for( ; fieldptr[i]!=0 && fieldptr[i]!='\t';i++ )
		{
			if( (fieldptr[i]=='G') && (fieldptr[i+1]=='T') && (fieldptr[i+2]==':'||fieldptr[i+2]=='\t'||fieldptr[i+2]==0) )
				break;
			if( fieldptr[i]==':' )
				GTidx++;
		}

		//
		if( fieldptr[i]==0 || fieldptr[i]=='\t' )
		{
			df0("WhopGenome::VCF_readIntoNucleotideMatrix : NO GT FIELD DEFINED!\n");
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------

		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-
		
		bHadRef = false;
		bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			
			//
			if( fieldptr == 0 )
			{
				Rprintf("WhopGenome::VCF_readIntoNucleotideMatrix : FAILED TO ACCESS SAMPLE (FIELD %d)\n", f->wanted_samples[ per_row ] );
				Rprintf("	debug info: per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("		baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("		numparsedfields=%d\n",f->numParsedFields());
				return RBool::False();
			}

			//	- find the GT subfield (skip ':' until found)
			//
			while( GTidx > 0 )
			{
				if( fieldptr[0] == ':' )
					GTidx--;
				fieldptr++;
			}
			
			//
			//	- parse the GT subfield ( regexp: [0-9]+[/|][0-9]+ )
			//
			

			//
			int left_allele = fieldptr[0] - '0';
			int right_allele = fieldptr[2] - '0';
			if( (fieldptr[1] != '|' && fieldptr[1] != '/') || (fieldptr[3] != '\t' && fieldptr[3] != ':') )
			{
				df0("WhopGenome::VCF_readIntoNucleotideMatrix : Malformed GT field!\n");
				return RBool::False();
			}

			// either chromosome has alt allele -> matrix stores "result is alt allele"
			//
			char	actual_nucleotide;
			if( ((left_allele==1)||(right_allele==1)) )
			{
				bHadAlt=true;	//only columns that had alternate and reference alleles are kept in the matrix; memorize we got an alternate allele
				actual_nucleotide = *altptr;
			}
			else
			{
				bHadRef=true;	//only columns that had alternate and reference alleles are kept in the matrix; memorize we got a reference allele
				actual_nucleotide = *refptr;
			}
			
			//	store nucleotide code
			//
			if( actual_nucleotide == 'A' )
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_A );
			else if( actual_nucleotide == 'C' )
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_C );
			else if( actual_nucleotide == 'G' )
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_G );
			else if( actual_nucleotide == 'T' )
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_T );
			else if( actual_nucleotide == 'N' )
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_N );
			else
				df1("WhopGenome::VCF_readIntoNucleotideMatrix : Unexpected NUCLEOTIDE : [%c] (%d)!\n",actual_nucleotide,actual_nucleotide);
			
			//
		}

		//
		//	clear any unused matrix rows if this column contains both alleles
		//
		if( (bHadRef&&bHadAlt) )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				SET_STRING_ELT( mat, per_row + per_column * nrow, nucleotide_N );
			}
			
			bialcols++;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(posbuffer) );
			}

		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{			
			per_column--;
			nonbialcols++;
		}//..else..if( had both alleles for sample-set )

	}//...for( each column in the matrix )
	
	//
	//
	SEXP m1 = mkChar("-1");
        PROTECT(m1);
	SEXP nix = mkChar("-");
        PROTECT(nix);
	for( unsigned int rcol = per_column; rcol < ncol; rcol ++ )
	{
		//
		for( unsigned int per_row=0; per_row < nrow ; per_row ++ )
		{
			SET_STRING_ELT( mat, per_row + rcol * nrow, nix );
		}//..for all unused rows in this unused column
			
		//	Clear column name
		//
		if( R_NilValue != colnamvec )
		{
			SET_STRING_ELT( colnamvec, rcol, m1 );
		}

	}//...for all unused columns
	
	//	print some statistics
	//
	df1("%d nonbial columns\n",nonbialcols);
	df1("%d bial columns\n",bialcols);
	df1("%d total columns\n",bialcols+nonbialcols);
	UNPROTECT(2);
	//
	return RBool::True();
}











		//-----------------------------------------------------------------------------		Biallelic
		
		
		
		
		
		
		
		
		
//!
inline bool isBiallelic( const char * _REF, const char* _ALT )
{
	if( _REF[1] != '\t' )	return false;
	if( _ALT[1] != '\t' )	return false;
	return true;
}



/*
w = ld("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")
*/

#define		DBI		if( true )
//#define		DBI		if( false )
	
//-
//-
//-		parse VCF and fill biallelic matrix
//-
//-
//-
bool	read_bial( bool bFilterActivated, vcff * f, RMatrix &m )
{

	int	nonbialcols=0,bialcols=0;
	
	
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();

	unsigned int	nrow = m.numRows();
	unsigned int	ncol = m.numCols();
	int*		 	ptr = m.getIntPtr();
	if( ptr == 0 )
	{
		Rprintf("WhopGenome::getBial : ERROR : Could not get access to the matrix in form of an int*!\n");
		return false;
	}
	
	//
	char			*fieldptr=0;
	unsigned int	per_column = 0;		//vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	unsigned int	column_stepsize = nrow;
	int				snppos=-1;
	SEXP			minus1_char = mkChar("-1");
	PROTECT(minus1_char);
	
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
		Rprintf("WhopGenome::getBial : WARNING : matrix has no column names vector! not setting SNP positions in matrix!\n");

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		const char* refptr;
		const char* altptr;
		bool	bLineParsed = f->parseNextLine();
		
		//	read lines until a SNP line is found (or the chromosomal region's end is reached)
		//
		while( bLineParsed )
		{
			//
			//
			/*
				- make sure its a bi-allelic SNP line
				* 
				* REFerence is a single nucleotide in case of a SNP
				* ALTernative is also a single nucleotide in case of a SNP
			*/
			refptr = (char*)f->getFieldPtr( REF );	//REFerence
			altptr = (char*)f->getFieldPtr( ALT );	//ALTernative
			
			//second char must be a \tab to ensure there is only a single REFerence nucleotide
			//
			if( refptr && refptr[1] == '\t' )
			{
				//second char must be a \tab to ensure there is only a single ALTernative nucleotide
				//
				if( altptr && altptr[1] == '\t' )
				{
					break;
				}
				else
				{
					//ALTernative is not just a single nucleotide
					//df1("Not a biallelic SNP! (REF=%9s)\n",refptr);
				}
			}
			else
			{
				//REFerence is not just a single nucleotide
				//df1("Not a SNP! (REF=%2s)\n",refptr);
			}
			
			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		//	if could not read another line from VCT, exit with error
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------
		
		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( fieldptr )
			snppos = atoi( fieldptr );
		else
			snppos=-1;

		//	per VCF format definition, genotype field must always be first in sample columns !
		//
		int GTidx=0;

#ifdef	ALLOW_NONSTANDARD_VCF_GT
		//-
		//
		//	Identify the GT-colon-field used in the per-individual fields
		//
		//-

		/*
			- decompose FORMAT-column for this row
				OPTIM: skip buffer-copying it
				? what to do if missing ??
			- find GT in FORMAT
				! break with error if not found!
			- memorise field-pos of GT
		*/
		fieldptr = (char*)f->getFieldPtr( FORMAT );
		int i=0;
		for( ; fieldptr[i]!=0 && fieldptr[i]!='\t';i++ )
		{
			// if we found the GT field inside the format specification, break out of the loop
			if( (fieldptr[i]=='G') && (fieldptr[i+1]=='T') && (fieldptr[i+2]==':'||fieldptr[i+2]=='\t'||fieldptr[i+2]==0) )
				break;
			//if we have found another : dividing the FORMAT column's elements, increase GT-index
			if( fieldptr[i]==':' )
				GTidx++;
		}

		//	make sure a GT field was defined in the FORMAT column
		//
		if( fieldptr[i]==0 || fieldptr[i]=='\t' )
		{
			df0("NO GT FIELD DEFINED!\n");
			UNPROTECT(1);
			return false;
		}
#endif
		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------
		
		
		//	check whether filtering needs to be done
		//
		if( bFilterActivated && (false == filterLine( f )) )
		{
			df1("get_bial :: Line (pos %d)has been filtered away\n",snppos);
			continue;
		}
		
		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-		
		bool	bHadRef = false, bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			
			//
			if( fieldptr == 0 )
			{
				Rprintf("vcf_get_bial :: ERROR when trying to get sample %d (matrix row %d) in file!\n",f->wanted_samples[per_row],per_row);
				Rprintf("	per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				UNPROTECT(1);
				return false;
			}

			//	- find the GT subfield (skip ':' until found)
			//
			while( GTidx > 0 )
			{
				if( fieldptr[0] == ':' )
					GTidx--;
				fieldptr++;
			}
			
			/*
				- parse the GT subfield ( regexp: [0-9]+[/|][0-9]+ )
			*/

			//
			//
			char* fcopy = fieldptr;
//			Rprintf("%d: %c,%c,%c //",snppos,fcopy[0],fcopy[1],fcopy[2]);
			int left_allele = fcopy[0]-'0';//parseDecInt( fcopy );
			
			//make sure there are at least two alleles, regardless of whether phased or unphased
			//
			if( fcopy[1] != '/' && fcopy[1] != '|' )
			{
				Rprintf("ERROR : unexpected character '%c' in Genotype field at position %d\n",fcopy[1],snppos);
				df0("	=> Syntax error in GT field (%s)!\n",fieldptr);
				UNPROTECT(1);
				return false;
			}
			
			//get second allele
			//
			int right_allele = fcopy[2]-'0';//parseDecInt( fcopy );
			
			//make sure we're handling a truely diploid species,
			//	i.e. exactly two alleles, no more | or / which divide the alleles in the genotype field)
			//
				/*\t separates samples*/
				/*0-byte ends last sample in the line (no \n ?? ) */
				/*: separates genotype from other per-sample data*/
			if( fcopy[3] != '\t' && fcopy[3] != 0 && fcopy[3] != ':' )
			{
				Rprintf("Syntax error in GT field (%s)!\n",fieldptr);
				UNPROTECT(1);
				return false;
			}

//			Rprintf("left=%d, right=%d\n",left_allele,right_allele);

			//
			//
			if( (left_allele==1)||(right_allele==1) )
			{
				bHadAlt = true;	//memorize that we got samples with alternate alleles
				ptr[per_row] = 1;
			}
			else
			{
				bHadRef = true;	//memorize that we got samples with reference alleles
				ptr[per_row] = 0;
			}
			
			//
		}//...for( all rows == all samples )
		
		//
		//	clear any unused matrix rows of this column, if it had both alleles
		//
		if( (bHadAlt && bHadRef) )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = -2;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(posbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{
			//df1("getbial : position %d : no bi-allelic position for the selected samples!\n",snppos);

			if( R_NilValue != colnamvec )
			{
				//char posbuffer[256];
				//snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, minus1_char );
			}
			
			per_column--;	//re-use the matrix column for the next SNP
			nonbialcols++;	//count how many non-biallelic columns we found
		}
		
		//
	}//...for( each column in the matrix )
	
	//? check whether the end of the region has been reached?
	//

	//
	//
	df1("getbial:\n	%d nonbial columns\n",nonbialcols);
	df1("	%d bial columns\n",bialcols);
	df1("	%d total columns\n",bialcols+nonbialcols);
	UNPROTECT(1);
	//
	return true;
}



/*! Reads 
**
**	- 	Overwrites the given matrix <mat>
**
**
*/
EXPORT SEXP VCF_getBial( SEXP vcfptr, SEXP mat )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getBial :: Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}
	
	df1("A\n");
	
	//	any samples selected ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("VCF_getBial :: No samples selected!\n");
		return R_NilValue;
	}
	
	df1("B\n");

	//	without a FORMAT field, it is not clear how to discover the genotype subfield per sample!
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("VCF_getBial :: VCF does not have a FORMAT field!\n");
		return R_NilValue;
	}
	
	df1("C\n");
	
	//	is mat a matrix?
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("VCF_getBial :: Second parameter is not an integer matrix!\n");
		return R_NilValue;
	}
	
	df1("D\n");
	
	//	biallelic matrices are integer-typed
	//
	//
	if( m.getType() != INTSXP )
	{
		Rprintf("VCF_getBial :: Second parameter not an integer matrix!\n");
		return R_NilValue;
	}
	
	df1("E\n");
	
	
	//	enough rows to hold all samples?
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("VCF_getBial :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return R_NilValue;
	}
	
	//
	//
	//
	bool succeeded = read_bial( false, f, m );

	if( ! succeeded )
		return RBool::False();
	return RBool::True();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT	SEXP	VCF_readIntoCodeMatrixdiploid2( SEXP vcfptr, SEXP mat )
{
	// statistics
	int	nonbialcols=0,
		bialcols=0;

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 1 is not a VCFhandle EXTPTR!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	samplefieldindex = f->getFirstSampleFieldIndex();
	if( samplefieldindex <= FORMAT )
	{
		Rprintf("VCF_readIntoCodeMatrix :: VCF does not appear to have a FORMAT field!\n");
		return RBool::False();
	}
	
	//
	//
	RMatrix m(mat);
	if( false == m.isValid() )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		return RBool::False();
	}
	
	//	biallelic matrices are integer-typed
	//		save 4 bytes per entry over standard 'double' datatype
	//
	if( m.getType() != INTSXP )
	{
		Rprintf("VCF_readIntoCodeMatrix :: Parameter 2 not an integer matrix!\n");
		return RBool::False();
	}
	
	//	enough rows in matrix for all samples ?
	//
	if( f->num_wanted_samples < 1 )
	{
		Rprintf("VCF_readIntoCodeMatrix :: No samples selected!\n");
		return RBool::False();
	}

	//
	//
	unsigned int	nrow = m.numRows();
	if( f->num_wanted_samples > (unsigned)nrow )
	{
		Rprintf("VCF_readIntoCodeMatrix :: %d samples selected but matrix offers only rows for %d samples!\n",f->num_wanted_samples,nrow);
		return RBool::False();
	}

	//
	//
	SEXP colnamvec = m.getColNames();
	if( R_NilValue == colnamvec )
	{
		Rprintf("WhopGenome::VCF_readIntoCodeMatrix : WARNING : matrix has no column names vector! Cannot set SNP positions in matrix!\n");
		return RBool::False();
	}

	//-
	//-
	//-		parse VCF and fill biallelic matrix
	//-
	//-
	//-


	unsigned int	ncol = m.numCols();
	int*		ptr  = m.getIntPtr();
	
	//
	char 		maprefalt[1000] = {""};  //      store REF and ALT values
	char		*fieldptr=0;
	unsigned int	per_column = 0;		//      vars here to find out what to clear
	unsigned int	per_row = 0;		//	when too little data exists
	
	unsigned int	column_stepsize = nrow;
	
	int			snppos      = -1;//unused
	SEXP			minus1_char = mkChar("-1");
        PROTECT(minus1_char);
	
//	df1("ncol=%d, nrow=%d, wanted=%d\n",ncol,nrow,f->num_wanted_samples);

	//
	//
	for( ; per_column < ncol ; per_column ++ )
	{
		//-
		//
		//	Get a valid biallelic SNP line from the file
		//
		//-

		const char* refptr;
		const char* altptr;
		bool	bLineParsed = f->parseNextLine();
		//
		//printf("Line: %d", bLineParsed);
		while( bLineParsed )
		{

			//
			//	- make sure its a bi-allelic SNP line
			//
			refptr = (char*)f->getFieldPtr( REF );
			altptr = (char*)f->getFieldPtr( ALT );

			// fill maprefalt ///////////////////////////////////////////////////
			maprefalt[0] = refptr[0];
			int xyz  = 0;
                        int xyz2 = 0;
			while(altptr[xyz]!='\t'){
			
				if(altptr[xyz]!=','){
				maprefalt[xyz2+1] = altptr[xyz];
				xyz2++;
				}
			xyz++;
			}
			//prints
			/* 
			printf(": %c",maprefalt[0]);
			printf(": %c",maprefalt[1]);
			printf(": %c",maprefalt[2]);
			printf(": %c",maprefalt[3]);
			printf("**");
			printf(": %c",altptr[0]);
			printf(": %c",altptr[1]);
			printf(": %c",altptr[2]);
			printf(": %c",altptr[3]);
			*/
			//end prints


			////////////////////////////////////////////////////////////////////
			//FIXME Take all lines ! Also multi allelic
			//printf("REF: %c", refptr[0]);
			//printf("ALT: %c", altptr[0]);
			//printf("\n");			

			if( refptr && refptr[1] == '\t' )	
			{
				//if( altptr && altptr[1] == '\t' )

				//{
					break;
				//}


				//else
				//{
				//	//df1("Not a biallelic SNP! (REF=%9s)\n",refptr);
				//}
			}
			//else
			//{
			//	//df1("Not a SNP! (REF=%2s)\n",refptr);
			//}

			
			//	try next line
			//
			bLineParsed = f->parseNextLine();

		}//while( could read another line from VCF )
		
		//	if could not read another line from VCT, exit with error
		//
		if( bLineParsed == false )
		{
			//df1("No more lines!\n");
			break;
		}

		//--------------
		//
		//	POSTCONDITION : got a line from VCF and it is a valid SNP
		//
		//--------------

		//
		fieldptr = (char*)f->getFieldPtr( POS );
		if( fieldptr )
		{
			snppos = atoi( fieldptr );
			if( snppos == 0 )
			{
				df1("VCF_readIntoCodeMatrix :: SNPpos=%d\n",snppos);
			}
		}

		//-
		//
		//	Identify the GT-colon-field used in the per-individual fields
		//
		//-

		//
		//	- decompose FORMAT field
		//		OPT: skip buffer-copying it
		//		? what to do if missing ??
		//	- find GT in FORMAT
		//		! break with error if not found!
		//	- memorise field-pos of GT
		//
		fieldptr = (char*)f->getFieldPtr( FORMAT );
		int GTidx=0;
		int i=0;
		for( ; fieldptr[i]!=0 && fieldptr[i]!='\t';i++ )
		{
			if( (fieldptr[i]=='G') && (fieldptr[i+1]=='T') && (fieldptr[i+2]==':'||fieldptr[i+2]=='\t'||fieldptr[i+2]==0) )
				break;
			if( fieldptr[i]==':' )
				GTidx++;
		}

		//
		if( fieldptr[i]==0 || fieldptr[i]=='\t' )
		{
			df0("VCF_readIntoCodeMatrix :: NO GT FIELD DEFINED!\n");
			UNPROTECT(1);
			return RBool::False();
		}

		//--------------
		//
		//	POSTCONDITIONS :
		//		- got a line from VCF and it is a valid SNP
		//		- successfully identified the GT subfield
		//
		//--------------


		//-
		//
		//	for each selected individual, get the SNP genotype information
		//		and store a 0 or 1 in the biallelic matrix
		//
		//-

		bool		bHadRef = false,
				bHadAlt = false;
		
		//
		for( per_row = 0; per_row < f->num_wanted_samples ; per_row ++ )
		{
			
			//	- get field of sample
			//
			fieldptr = (char*)f->getFieldPtr( f->wanted_samples[per_row] );
			
			//
			if( fieldptr == 0 )
			{
				Rprintf("VCF_readIntoCodeMatrix ::  Problem with reading sample's data!\n");
				Rprintf("	debug info : per_row =%d\nwanted_sample[per_row]=%d\n",per_row, f->wanted_samples[per_row] );
				Rprintf("	baseindex=%d, field = %d\n",samplefieldindex, (samplefieldindex + f->wanted_samples[per_row]) );
				Rprintf("	numparsedfields=%d\n",f->numParsedFields());
				UNPROTECT(1);
				return RBool::False();
			}

			//	- find the GT subfield (skip ':' until found)
			//
			while( GTidx > 0 )
			{
				if( fieldptr[0] == ':' )
					GTidx--;
				fieldptr++;
			}
			
			//
			//	- parse the GT subfield ( regexp: [0-9]+[/|][0-9]+ )
			//

			int left_allele;
			if(fieldptr[0]=='.'){
			left_allele = -1; // unknown
			}else{
			left_allele = fieldptr[0] - '0';
                        }
			
			int right_allele;
			if(fieldptr[2]=='.'){
			right_allele = -1; // unknown
			}else{
			right_allele = fieldptr[2] - '0';
			}
			
			
			//printf("\n");
			//printf("left: %c", fieldptr[0]);			
			//printf("right: %c", fieldptr[2]);
			//printf("\n");

			//assert( left_allele >= 0 && left_allele <= 9 )
			//assert( right_allele >= 0 && right_allele <= 9 )
			
			//
			if( (fieldptr[1] != '|' && fieldptr[1] != '/')) //|| (fieldptr[3] != '\t' && fieldptr[3] != ':') )
			{
				df0("VCF_readIntoCodeMatrix :: Malformed GT field!\n");
				UNPROTECT(1);
				return RBool::False();
			}

			// any chromosome has alt allele -> result is alt allele
			//
			///char	actual_nucleotide;
			// diploid

			/*
			if( ((left_allele==1)||(right_allele==1)) )
			{
				bHadAlt=true;
				//ptr[per_row] = nucleotide_mapping[(unsigned)*altptr];
				if((left_allele==1) && (right_allele==0)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*altptr]+nucleotide_mapping[(unsigned)*refptr];
				}
				if((left_allele==0) && (right_allele==1)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*refptr]+nucleotide_mapping[(unsigned)*altptr];
				}
				if((left_allele==1) && (right_allele==1)){
				ptr[per_row] =	10*nucleotide_mapping[(unsigned)*altptr]+nucleotide_mapping[(unsigned)*altptr];
				}
					

			}
			else
			{
				bHadRef=true;
				//ptr[per_row] = nucleotide_mapping[(unsigned)*refptr];
				ptr[per_row]   = 10*nucleotide_mapping[(unsigned)*refptr]+nucleotide_mapping[(unsigned)*refptr];
			}
			*/

			///////////////////////////////////////////////////////////////////////////////////////////////////////
			bHadRef=true;
			bHadAlt=true;

			int first;
			if(left_allele==-1){
                        first = 50; 
                        }else{
			first = 10*nucleotide_mapping[(int)maprefalt[left_allele]];
			}
			
			int second;
			if(right_allele==-1){
			second = 5;
			}else{
			second = nucleotide_mapping[(int)maprefalt[right_allele]];
			}				

			//printf("first %d",first);
			//printf("second %d",second);
			//printf("\n");

			ptr[per_row] = first + second;

			// if(fieldptr[3]==EOF){break;}
			// printf("Last:%c", fieldptr[3]);
			//////////////////////////////////////////////////////////////////////////////////////////////////////
				
		}
		
		//
		//	clear any unused matrix rows
		//
		if( (bHadRef&&bHadAlt) )
		{
			
			//
			for( ; per_row < nrow ; per_row ++ )
			{
				ptr[per_row] = 77;
			}
			
			//
			ptr += column_stepsize;
			
			bialcols++;
			
			
			//	Set Column Names = SNP Positions
			//
			if( R_NilValue != colnamvec )
			{
				char posbuffer[256];
				snprintf(posbuffer,sizeof(posbuffer)-2,"%d",snppos);
				SET_STRING_ELT( colnamvec, per_column, mkChar(posbuffer) );
			}
		
			//FIXME any rows beyond the ones we need : zero them out here
			//	or as a separate step after this double for-loop
		}
		//
		//	if the last VCF-line's individuals did not provide both alleles
		//		reuse current column and try the next line
		//
		else
		{			
			//
			per_column--;
			nonbialcols++;
		}

	}//...for( each column in the matrix )
	
	//	Reset/clear unused columns
	//
	for( unsigned int rcol = per_column; rcol < ncol; rcol ++ )
	{
		//
		for( unsigned int per_row=0; per_row < nrow ; per_row ++ )
		{
			ptr[per_row] = 77;
		}
		ptr += column_stepsize;
			
		//	Clear column name
		//
		if( R_NilValue != colnamvec )
		{
			SET_STRING_ELT( colnamvec, rcol, minus1_char );
		}
	}
	
	//	print some statistics
	//
	df1("VCF_readIntoCodeMatrix ::\n\t%d nonbial columns\n",nonbialcols);
	df1("\t%d bial columns\n",bialcols);
	df1("\t%d total columns\n",bialcols+nonbialcols);
	UNPROTECT(1);
	//
	return RBool::True();
}//...



