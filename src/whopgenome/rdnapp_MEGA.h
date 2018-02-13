/*
**
**
**			MEGA
**
**
**
--------------------------------------------------
#mega
TITLE: Written by EMBOSS 31/12/05

#HSFAU                ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
#HSFAU1               ctaccattttccctctcgattctatatgtacactcgggacaagttctcct

#HSFAU                gttcagtcaaaaaaaaaa
#HSFAU1               tcaggtaagaatggggccttggctggatccgaagggcttgtagcaggttg

#HSFAU
#HSFAU1               gctgcggggtcagaaggcgcggggggaaccgaagaacggggcctgctccg

#HSFAU
#HSFAU1               gaggtgtgcttctcgg

--------------------------------------------------

- one or more alignment blocks (i.e. a block of lines beginning with "#" and last line is empty)
- inside an alignment block, sequences have same length, but each alignment block can have its own sequence length


*/


//*
//*
//*

#include <ctype.h>

//*
//*
//*

#define		IS_LC_BASE( ch )	( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' )
#define		IS_UC_BASE( ch )	( ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' )
#define		IS_BASE( ch )		( IS_LC_BASE( ch ) || IS_UC_BASE( ch ) )


#ifdef SUPPORT_MEGA




/*!	Checks if the given file is a MEGA file and scans it to determine the alignment matrix dimensions we are going to need
**
*/
bool	determineAlignmentDimensionsMEGA( unsigned int &numnucleotides, unsigned int& numsamples )
{
 try
 {
	//
	int		sta=0;

	 /*
#mega
TITLE: Written by EMBOSS 31/12/05	oder !Title
zB !Format	DataType=Protein

#HSFAU                ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
#HSFAU1               ctaccattttccctctcgattctatatgtacactcgggacaagttctcct

#HSFAU                gttcagtcaaaaaaaaaa
#HSFAU1               tcaggtaagaatggggccttggctggatccgaagggcttgtagcaggttg

#HSFAU
#HSFAU1               gctgcggggtcagaaggcgcggggggaaccgaagaacggggcctgctccg

#HSFAU
#HSFAU1               gaggtgtgcttctcgg
--------------------------------------------------------------------
		#mega als erste zeile
		finde erste zeile mit sequenz (#blabla  sequenzsequenzsequenzzz)
			erste Zeile nach #mega, die mit # anfängt
		schleife:
			zeile einlesen
			eof?
				schleife ende
			1ste?
				laenge fuer block festlegen
			oder leerzeile?
				blockende, zeilenzähler zurücksetzen
			andernfalls
				name noch zu speichern?
					namen merken
				sequenz einlesen [laenge merken + ob proteinseq erkennen // nur sequenz einlesen ]
				laenge der ersten sequenz des blocks == laenge der andern? und >0 ? -->

	 */

	char 	c;
	int		curseqlen=0;
	//
	int		sta=0;

	//
	//
	numsamples=0;
	numnucleotides=0;

	//	parse and validate header
	//
	if( false == readFileBlock( 0 ) )
	{
		Rprintf( "Load-MEGA : couldnt read fileblock at 0\n");
		return false;
	}

	//	parse and validate header
	//
	if( false == cmpString("#mega")  )
	{
		ONDBG Rprintf("A\n");
		return false;
	}

	//	skip extraneous lines at the top
	//
	while( (c = skipLine()) != '#' )
	{
		ONDBG Rprintf("line skipped\n");
	}


	//---------	begin of alignmentblocks (interleaved)
	int		lastseqlen = -1,
			curnumseqs=0,
			lastnumseqs=0;

	//------------- FORMAT: -----------------------------------------------------
	//#HSFAU                ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
	//#HSFAU1               ctaccattttccctctcgattctatatgtacactcgggacaagttctcct
	//
	//[next block begins here]
	//---------------------------------------------------------------------------

	//
	//guaranteed to have read the first # after a newline
	//
	while( ! Eof() )
	{

		//
		c = nextChar();

		//
		//
		switch( sta )
		{

			//
			//	begin of sequence (after #)
			//
			case 0:
				if( c != ' ' )
					ONDBG Rprintf("SeqNameChar: %c\n",c);
				else
				{
					curnumseqs++;
					sta = 1;
				}
				break;

			//
			//	empty space after #seqname
			//
			case 1:
				if( c == ' ' )		//go to next state if we found a non-space character
					break;

				curseqlen=0;
				sta = 2;
				//remember this bug: a break here without this else, the first sequence-character was lost!

			//
			//	sequence string
			//
			//		stay in this state until either 1) a newline or 2) end of file appears
			//			expect only nucleotide characters here
			//
			case 2:
				c = tolower( c );

				//
				//	if newline, end of this sequence's fragment in the current block, maybe end of current block
				//
				if( c == '\n' )
				{
					//
					if( lastseqlen > 0 && curseqlen != lastseqlen )
					{
						Rprintf("ERROR!:");Rprintf("Sequence has %d instead of %d bases in this block!\n",curseqlen,lastseqlen);
						return false;
					}

					//
					sta = 3;
					lastseqlen = curseqlen;
					curseqlen=0;
				}
				//
				//	if not a nucleotide, ERROR
				//
				else if( c != 't' && c != 'a' && c != 'c' && c != 'g' )
				{
					Rprintf("ERROR!:");Rprintf("(%c) is not a valid character in a DNA/RNA sequence!\n",c);
					return false;
				}

				//
				//	valid character for a sequence
				//
				ONDBG Rprintf("(%c)",c);
				curseqlen++;

				break;

			//
			//	end of line
			//
			case 3:

				//
				//	if the last sequence fragment line is followed by a # in the next line, the next sequence's fragment begins
				//
				if( c == '#' )
				{
					ONDBG Rprintf("Next sequence begins now...\n");
					ONDBG Rprintf("last len=%d thislen=%d newlen=%d\n", numnucleotides,lastseqlen,numnucleotides+lastseqlen);
					sta = 0;
				}
				//
				//	empty line after sequence fragment -> end of sequence fragments block
				//
				else if( c == '\n' )
				{
					ONDBG Rprintf("End of block found\n");

					//
					if( lastnumseqs > 0 && lastnumseqs != curnumseqs )
					{
						Rprintf("ERROR!:");Rprintf("This block has %d instead of the expected %d number of sequences!\n",curnumseqs,lastnumseqs);
						return false;
					}

					//
					numnucleotides+=lastseqlen;
					numsamples=curnumseqs;

					//
					lastnumseqs=curnumseqs;
					curnumseqs=0;
					curseqlen=lastseqlen=0;
					sta=4;
				}
				break;

			//
			//	skip to next block
			//
			case 4:
				ONDBG Rprintf("SKIP[%c]\n",c);
				while( c != '#' && !Eof())
				{
					c = skipLine();
					ONDBG Rprintf("SKIP[%c]\n",c);
				}
//				if( Eof() )
//					return true;
				sta=0;
				break;

		}//...switch

	}//...while

 }//...try
 catch(const char* err )
 {
	 Rprintf("ERROR!:");Rprintf("(!!) Exception caught during MEGA load: '%s'\n",err);
 }
 catch( ... )
 {
	 Rprintf("ERROR!:");Rprintf("(!!) Unknown exception caught during MEGA load!\n");
 }

 return false;

}









/*!
**
**
**
*/
SEXP	processAlignmentMEGA( void )
{
	//
	RMatrix			resmat_inst;
	SEXP			resmat	=	R_NilValue
						;
	char			c=1
						;
	int				*mat=0
						;
	unsigned int	numnucleotides=0,
					numsamples=0,
					curseqlen=0
						;

	//
	//
	if( false == determineAlignmentDimensionsMEGA( numnucleotides, numsamples ) )
	{
		//Rprintf("Not a proper MEGA format file\n");
		return R_NilValue;
	}

	//
	Rprintf("Detected MEGA format file with %d sequences and %d nucleotides per sequence\n",numsamples,numnucleotides);

	//	read first N bytes of the file
	//
	if( false == readFileBlock( 0 ) )
	{
		Rprintf( "Load-MEGA : couldnt read fileblock at 0\n");
		return false;
	}

	//
	//

	//	parse and validate header
	//
	if( false == cmpString("#mega")  )
	{
		ONDBG Rprintf("A\n");
		return R_NilValue;
	}

	//	skip extraneous lines at the top
	//
	c = skipLine();
	while( c != '#' )
	{
		c = skipLine();
	}

	//
	//	allocate result matrix
	//
	if( false == resmat_inst.alloc(INTSXP,numsamples,  numnucleotides) )
	{
		Rprintf("(!!) Result matrix alloc failure! numsamples = %d, numnucleotides = %d\n", numsamples , numnucleotides);
		return R_NilValue;
	}
	mat = resmat_inst.getIntPtr();
	if( mat == 0 )
	{
		Rprintf("(!!) Result matrix alloc failure! numsamples = %d, numnucleotides = %d\n", numsamples , numnucleotides);
		resmat_inst.dealloc();
		return R_NilValue;
	}
	resmat = resmat_inst.get();


	//---------	begin of alignmentblocks (interleaved)
	unsigned int	lastseqlen = -1,
					curnumseqs=0,
					lastnumseqs=0;
	int				sta=0;

	//------------- FORMAT: -----------------------------------------------------
	//#HSFAU                ttcctctttctcgactccatcttcgcggtagctgggaccgccgttcagtc
	//#HSFAU1               ctaccattttccctctcgattctatatgtacactcgggacaagttctcct
	//
	//[next block begins here]
	//---------------------------------------------------------------------------

	//
	//guaranteed to have read the first # after a newline
	//
	while( ! Eof() )
	{

		//
		c = nextChar();

		//
		//
		switch( sta )
		{

			//
			//	begin of sequence (after #)
			//
			case 0:
				if( c != ' ' )
					ONDBG Rprintf("SeqNameChar: %c\n",c);
				else
				{
					curnumseqs++;
					sta = 1;
				}
				break;

			//
			//	empty space after #seqname
			//
			case 1:
				if( c == ' ' )		//go to next state if we found a non-space character
					break;

				curseqlen=0;
				sta = 2;
				//remember this bug: a break here without this else, the first sequence-character was lost!

			//
			//	sequence string
			//
			//		stay in this state until either 1) a newline or 2) end of file appears
			//			expect only nucleotide characters here
			//
			case 2:
				c = tolower( c );

				//
				//	if newline, end of this sequence's fragment in the current block, maybe end of current block
				//
				if( c == '\n' )
				{
					//
					if( lastseqlen > 0 && curseqlen != lastseqlen )
					{
						Rprintf("Sequence has %d instead of %d bases in this block!\n",curseqlen,lastseqlen);
						resmat_inst.dealloc();
						return R_NilValue;
					}

					//
					sta = 3;
					lastseqlen = curseqlen;
					curseqlen=0;
				}
				//
				//	if not a nucleotide, ERROR
				//
				else if( c != 't' && c != 'a' && c != 'c' && c != 'g' )
				{
					Rprintf("(%c) is not a valid character in a DNA/RNA sequence!\n",c);
					resmat_inst.dealloc();
					return R_NilValue;
				}

				//
				//	valid character for a sequence
				//
				ONDBG Rprintf("(%c)",c);
				curseqlen++;

				break;

			//
			//	end of line
			//
			case 3:

				//
				//	if the last sequence fragment line is followed by a # in the next line, the next sequence's fragment begins
				//
				if( c == '#' )
				{
					ONDBG Rprintf("Next sequence begins now...\n");
					ONDBG Rprintf("last len=%d thislen=%d newlen=%d\n", numnucleotides,lastseqlen,numnucleotides+lastseqlen);
					sta = 0;
				}
				//
				//	empty line after sequence fragment -> end of sequence fragments block
				//
				else if( c == '\n' )
				{
					ONDBG Rprintf("End of block found\n");

					//
					if( lastnumseqs > 0 && lastnumseqs != curnumseqs )
					{
						Rprintf("This block has %d instead of the expected %d number of sequences!\n",curnumseqs,lastnumseqs);
						resmat_inst.dealloc();
						return R_NilValue;
					}

					//
					numnucleotides+=lastseqlen;
					numsamples=curnumseqs;

					//
					lastnumseqs=curnumseqs;
					curnumseqs=0;
					curseqlen=lastseqlen=0;
					sta=4;
				}
				break;

			//
			//	skip to next block
			//
			case 4:
				ONDBG Rprintf("SKIP[%c]\n",c);
				while( c != '#' && !Eof())
				{
					c = skipLine();
					ONDBG Rprintf("SKIP[%c]\n",c);
				}
//				if( Eof() )
//					return resmat;
				sta=0;
				break;

		}//...switch

	}//...while

	return resmat;


//============================================

	//
	//	set sequence names as row names
	//
	setMatrixRownames( resmat , numsamples );

	//
	//	read sequence data and encode it in the numeric format
	//


/*
	//	read first N bytes of the file
	//
	blkidx=0;
	fileidx=0;
	if( 0 == readFileBlock( 0 ) )
	{
		L( "Couldn't read first block of data from file!\n" )
		return false;
	}
*/

	//
	return resmat;
}



#endif //ifdef SUPPORT_MEGA
