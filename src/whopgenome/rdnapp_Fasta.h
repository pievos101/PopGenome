/*
**
**
**			FASTA
**
**
**
*/

#ifdef SUPPORT_FASTA

//#define	FDBG
#define	FDBG	ONDBG


/*!	Checks if the given file is a FASTA file and scans it to determine the alignment matrix dimensions we are going to need
**
**
**
*/
bool	determineAlignmentDimensionsFasta( unsigned int &numnucleotides, unsigned int& numsamples )
{
	unsigned char c=0;

	//	read first N bytes of the file
	//
	if( 0 == readFileBlock( 0 ) )
	{
		Rprintf("(!!) determineAlignmentDimensionsFasta : couldnt read from file!\n" );
		return false;
	}
		
	//	identify as Fasta?
	//
	if( currentChar() != '>' )
	{
//		Rprintf("(!!) determineAlignmentDimensionsFasta : file does not begin with a '>'\n" );
		return false;
	}

	//	
	//
	numsamples=0;
	numnucleotides=0;
	
	//	1.	find end of >BLABLA line in memory buffer
	//
	//
	FDBG Rprintf("skiprestofline...");
	while( (c=nextChar()) != '\n' && (c!='\r') )
		;
	FDBG Rprintf("last char=%c=%d\n",c,c );

	//	2.	determine length of first sequence
	//
	FDBG Rprintf("determine a) bytes-per-sequence b) num-nucleotides-per-sequence\n");
	int	bytespersequence = 0;
	while( (c=nextChar()) != '>' && c != 0 )
	{
		if( isNucleotide(c) )
		{
			numnucleotides++;
		}
		bytespersequence++;
	}
	FDBG Rprintf("%d bytes, %d nucs\n",bytespersequence,numnucleotides);
	
	if( numnucleotides == 0 )
	{
		FDBG Rprintf("(!!) determineAlignmentDimensionsFasta : 1st sequence has no nucleotides!\n");
		return false;
	}
	
	//	3.	determine number of sequences
	//		- find end of >BLABLA line
	//		- jump 'bytespersequence' bytes ahead
	//		- expect to find a > after one or more whitespaces
	//		- if > could not be found, return false
	//

	//count number of sequences
	//	check consistency of numnucleotides per sequence
	numsamples=1;
	unsigned int thissample_numnucleotides=0;
	while( ! atEof )
	{	
		numsamples++;
		thissample_numnucleotides=0;
		
		//	skip rest of > line
		//
		if( c == '>' )
		{
			while( (c=nextChar()) && (c != '\r') && (c != '\n') && (c != 0) )
			;
		}

		//	read sequence data until the next >
		//
		while( (c=nextChar()) != '>'  && c != 0 )
		{
			
			//only count valid nucleotides
			//
			if( isIUPACNucleotide(c) )
			{
				thissample_numnucleotides++;
			}
			else
			{
				if( c != '\n' && c != '\t' && c != ' ' && c != '\r' )
				{
					Rprintf("FASTA: NOT COUNTING '%c'(%d) as NUCLEOTIDE! REPLACED BY 'n'! \n",c,c);
					//return false;
					thissample_numnucleotides++;
				}
			}
			bytespersequence++;
		}
		if( thissample_numnucleotides != numnucleotides )
		{
			Rprintf("\tERROR : Sample #%d has %d instead of %d nucleotides!\n" ,numsamples,thissample_numnucleotides,numnucleotides);
			return false;
		}
	}

	FDBG Rprintf("at end of file - numsamples= %d numsamples\n");

	FDBG Rprintf("\tbytespersequence is %d  __ numnucleotides is %d\n",bytespersequence,numnucleotides);
	FDBG Rprintf("blkidx %d blkidx\n");
	FDBG Rprintf("fileidx %d fileidx\n");
	
	//
	return true;
}
	

//!
SEXP	processAlignmentFasta( void )
{
	SEXP	resmat;
	unsigned int 	numnucleotides=0;
	unsigned int		numsamples=0;
	
	//
	if( false == determineAlignmentDimensionsFasta( numnucleotides, numsamples ) )
		return R_NilValue;

	FDBG Rprintf("File has %d nucleotides per sample and  %d samples\n",numnucleotides,numsamples);

	//
//	rownames.reset();

	//	start copying matrix from the start
	//
	if( 0 == readFileBlock(0)  )
		return R_NilValue;
		
	//	allocate result matrix
	//
	RMatrix resmat_inst;
	if( false == resmat_inst.alloc(INTSXP,numsamples,  numnucleotides) )
	{
		Rprintf("(!!) FASTA error : Matrix alloc failure!\n	numsamples=%d\nnumnucleotides=%d\n",numsamples,numnucleotides);
		return R_NilValue;
	}
	int	*mat = resmat_inst.getIntPtr();//allocIntRMatrix( numsamples,  numnucleotides, &resmat );
	if( mat == 0 )
	{
		Rprintf("(!!) FASTA error : Matrix alloc failure!\n	numsamples=%d\nnumnucleotides=%d\n",numsamples,numnucleotides);
		resmat_inst.dealloc();
		return R_NilValue;
	}
	resmat = resmat_inst.get();

	//	loop over file from beginning
	//
	char			lastc='\n';
	char			mode='>';
	unsigned int	i=0;		//FIXME : long long int for matrices with more than 1 GEl
	unsigned int	cursample=0;
	int				*tempmat= mat;

	char			c = currentChar();
	while( ! Eof() )
	{
		//
		//
		switch( mode )
		{
			//	find next >
			//
			case '>':
				if( c == '>' )
				{
					if( lastc != '\n' )
					{
						break;
					}
					mode = '\n';
				}
				break;
			//
			//	scan sample-name from >-line
			//
			case '\n':
				if( c=='\r' || c == '\n' )
				{
					if( c == '\r' )
						c = nextChar();
					rownames.serialize('\0');
					i=0;
					mode = 'A';
					break;
				}
				rownames.serialize(c);
				break;
			//
			//	translate nucleotide sequence into biallelic-codes and copy into R matrix
			//
			case 'A':
				
				//
				//
				if( c == '>' )
				{
					cursample++;
					if( cursample > numsamples )
					{
						Rprintf("(!!) FASTA : more samples in sequence than expected: %d but expected only %d\n",cursample,numsamples);
						resmat_inst.dealloc();
						return R_NilValue;
					}
					tempmat = &mat[cursample];
					mode = '\n';
				}
				//
				//
				else if( c != '\n' && c != ' ' && c != '\t' && c != '\r' )
				{
					tempmat[i*numsamples] = (int)(nucleotide_mapping[(unsigned)c]);
					i++;
					if( i > numnucleotides )
					{
						Rprintf("CHAR=%c\n",c);
						Rprintf("(!!) FASTA : more nucleotides in sequence than expected: %d but expected only %d\n",i,numnucleotides);
						resmat_inst.dealloc();
						return R_NilValue;
					}
				}
				break;
			default:
				break;
		}//...switch( mode )
			

		lastc=c;
		c=nextChar();

	}

	//
	return resmat;
	//
}



#endif	//#ifdef SUPPORT_FASTA

