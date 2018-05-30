/*
**
**
**			MAF
**
**
**
--------------------------------------------------
##maf version=1 scoring=whatever
# stuff on how it is obtained can be in this line
a score=30778.0
s sequence1     1542 460 + 246127941 ctggagattctta-ttagtgatttgggctggggc-ctggccatgtgtattttttta-aatttccactgatgattttgctgcatggccggtgttgagaatgactgCG-CAAATTTGCCGGATTTCCTTTGCTGTTCCTGCATGTAGTTTAAACGAGATTGCCAGCACCGGGTATCATTCACCAT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTTTCTTTTCGTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCC---TGTGCCAGGGTGCAAGCTGAGC-ACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATC-TTCTGGCCCCTGTTGTCT
s sequence2 27723223 600 - 149950539 cTAGGGAGTCTTAGTCAAAGGTTTGGACCAAGTCCCTGGCCATGCAGATCTTTGTAGAATCTCCACTCGTGACTTTCCTGCATAACCAGAGTTGAGCATCTTTGAGTCAAGTGTGCCAACTTTCTT--------------TGCTGTTTAAATAAGGATGCCAACACCGCATGTCATTAACAGTCTCGTAGGTTGATTGATTTGTTGGCTGGCTCAAAAATGAGAG-TTATTTTTCATTTTGTTTTGAt----------------------------------------------------------------------------------------------------------------tgattttttaagtcttgatctagatagcccagctgggttggagcttactatgtagtttaggttgcgctgcaactctcaatcctccagtctccaatcctcaagtgtccctccaggtctatgctactgtactcagGTAAAAAGGAG-TTTTCTGTCTGCTAATTTGCCACCAGTCATTTC---------------CTATT-ACGTGTGTCTGCTGCCTCCTAGCCCAGGCT-----TGCCCTTCCTCCC--TCTTCTGAGGTGTCATAGGGTCGTGAC--------------------TTACCTGGTTTGGGGGAGTAGTTGGAA------------------GCTGAGTGAGTG--GTGGGGTTTTCTTATGCTAAAGACCTGCGTCCAGTATAGGAAGAGCCATGTGCCTCCACTCTGGCCCTTGTGGTCT
s sequence3 29160419 613 - 187371129 CTGGAGAGTCTTATTTGAAGGGTTGGACCAAGCCACTGGCCATGTAGATCTATTCATAATCACTACTGGTGACTTTCATGTATAACCAGAGTTGAGCATCTTTGAGTCAAATGTGCCAAATTTCCT--------------TGCTGTTTAAATAAGGATGCCAACACTGCATATCATTAACAGTCTTGTAGGTTGATTGATTAGTTGGCTGGCTGGGGAACGGGGGAGTATTTTTCATTTTGTTTTGATTTTTAAGTCATGATCTATATAGCCCAGCTGGGCTGGAGCTTACTATGTAGTTTAGATTAGGCTGCAACTCTCAATCCTCCTGTCTCCACTTTCCAGTGTCACTCCAGGTCTATGCTACTGT----------------------------------------------------------------------------------------------------------------------ACCCAGCTAAAAAGTAGTTTTTCTGCCTGCTAATTTGCCACCAGTCTTTTC---------------CTGTT-ACATGTACCCACTGCCTCCTAGCCTAGGCT-----TGTCCTTCCTCCC--TCTTCTGAAAGGTCACAGGGTCTTGAC--------------------TTACCTGGGTTGGGGGAGGGGTTGGAAGCACACGCTGATTTGGATGCTGAGTGACTG--GTGGAGTTTTCTTATGACAAAGACCTGTGTCCAGGATGGGATAGGCCACACGCTTCC-CTCTGGCCCTTGTGGTCT

a score=1111.0
s sequence1     2040 63 + 246127941 ATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTT
s sequence2 64593042 56 +  93529596 GTTGGAGGGAAGATGAGTGAAGGGATCAATTTCTCTGATGACCTGGGCCGGTAGGT-------
s sequence3 29162556 61 - 187371129 ATTGGAGGGAGGGTGAACAAAGAGATAGACT--TCTGGCAACCTGGGCCAGTAGGTAGTGTCT
s sequence5     3040 63 + 246127941 TTAGGTGGATAGATTAGTGATAGCATCATCTTCTCTCTCAACCTAGTCCAGTAAGTATTGCTT

a score=-999.0
s sequence1      2103 19 + 246127941 GTGCTCATCTCCTTGGCTG----
s sequence2  64593098 19 +  93529596 GTGGTGTCCTCTTTGTCTG----
s sequence3 105237468 19 + 113649943 ----GGCCCTAATTGCTAAGGCA
s sequence4 105237468 19 + 113649943 ----GGCCCTAATGCTATAGGCA
s sequence5 105237468 19 + 113649943 ----GGCCTAATTGCTAAGCGCA

--------------------------------------------------

- one or more alignment blocks (i.e. a block of lines, first beginning with "a score=..." and last line is empty)
- an alignment block does not necessarily have to include all sequences (above: sequence4 missing in the first two alignment blocks, sequence 5 in the first)
- inside an alignment block, sequences have same length, but each alignment block can have its own sequence length
- besides "a " and "s " lines, others ("i ", "q ",..) could appear, too


*/

#ifdef SUPPORT_MAF



/*!	Checks if the given file is a MAF file and scans it to determine the alignment matrix dimensions we are going to need
**
**
**
*/
bool	determineAlignmentDimensionsMAF( unsigned int &numnucleotides, unsigned int& numsamples )
{
	char 	c=0;
	char	seqnamebuf[128];
	unsigned int		si=0;
	int		curseqlen=0,
			oldseqlen=0;
	int		ni=0;
	
	numnucleotides=0;
	numsamples=0;

	//	read first N bytes of the file
	//
	if( 0 == readFileBlock( 0 ) )
	{
		Rprintf( "MAF : determineAlignmentDimensionsMAF error : Couldn't read first block of data from file!\n" );
		return false;
	}
	
	//	parse and validate header
	//
	if( cmpString("track ") )		//track .......
		c = skipLine();
		
	//
	if( false == cmpString("##maf version=")  )	//##maf version
	{
		//Rprintf("(!!) MAF error : file does not start with proper header!\n");
		return false;
	}
	c = skipLine();

	while( (c=nextChar()) == '#' )
	{
		c = skipLine();
	}

	//
	//	scan through the file until EOF
	//
	while( atEof == false )
	{
		curseqlen=0;
		oldseqlen=0;

		//
		//----------------	find next alignment block
		//

		//	find first line beginning with "a " == alignment block begin
		//
		while( c != 'a' && c != 0 )
		{
			c = skipLine();
		}
		if( c == 0 || nextChar() != ' ' )
		{
			if( numsamples == 0 )
			{
				Rprintf("(!!) MAF error : No sequences found!\n[%s]\n",memorybuffer);
				return false;
			}
			break;//eof
		}

		//
		//----------------	find next sequence
		//

		while( (c=skipLine()) != '\n' && c!='\r' && c != 0 )
		{
			if( oldseqlen==0 )
				oldseqlen=curseqlen;
			else
				if( curseqlen != oldseqlen )
				{
					Rprintf("(!!) MAF error : Error in file : sequence length in alignment block mismatch : last=%d <=> %d=current\n",oldseqlen,curseqlen);
					return false;
				}

			//	find first "s " line, count its bytelength including \n, count its sequencelength
			//	count number of "s " lines
			//		skip ahead by bytelength bytes to speed up
			//		skip other "? " lines without counting them
			//		at the first empty line, stop counting
			//
	
			//	repeat this process ("a ", "s ", count seqlengths,"\n") until end of file
			//
			
			// if it is a sequence line...
			//
			if( c == 's' )
			{
				
				// expect a ' ' to follow the 's'
				//
				if( (c=nextChar()) != ' ' )
				{
					ONDBG Rprintf("	MAF Line expected to be a sequence line, but does not begin with 's '\n");
					return false;			
				}
				
				//
				//	get name of sequence
				//
				si=0;
				while( (c=nextChar()) != ' ' && si < sizeof(seqnamebuf) )
				{
					seqnamebuf[si++]=c;
				}
				
				//
				if( si >= sizeof(seqnamebuf) )
				{
					ONDBG Rprintf("	sequence name length exceeds available buffer size %d\n",sizeof(seqnamebuf));
					return false;
				}
				
				//
				seqnamebuf[si]=0;
				
				//	new sequence ?
				if( (ni=rownames.hasname( seqnamebuf )) >= 0)
				{
					ONDBG Rprintf("(ii) Sequence name '%s' already stored in rownames as sequence %d\n",seqnamebuf,ni);
					//return false;
				}
				else
				{
					rownames.serialize(seqnamebuf,si+1);
					numsamples++;
				}
			

				//
				//	skip past the extra fields, e.g. "     2048 63 + 246127941 "[seq begins here]
				//
				bool	hadsign=false;
				while( (c=nextChar()) != 0 )
				{
					ONDBG Rprintf("[%c]",c);
					if( c=='\r' || c=='\n' )
					{
						ONDBG Rprintf("Unexpected, premature EOL CR or LF (%d)\n",c);
						return false;
					}
					if( c == '+' ){ hadsign=true;continue;}
					if( c== '-' ){ if( hadsign )break;else hadsign=true; continue;}
					if( isDecDigit(c) )continue;
					if( c == ' ' ) continue;
					
					//anything else must belong to the sequence, so end the 'skipping extra fields'-part
					break;
				}
				
				//@COND: c == first char of sequence now

				//
				//----------------	count number of nucleotides used
				//
			
				ONDBG Rprintf("'%s':\nSeq=",seqnamebuf);
				curseqlen=0;
				
				// read sequence until End of line
				//
				while( c != '\r' && c != '\n' && c != 0 )
				{
		//			if( c != 'A' && c!='C' && c!='T'&&c!='G'&&c!='a'&&c!='c'&&c!='t'&&c!='g'&&c!='-')
		//				ONDBG Rprintf("\n\n***UNEXPECTED SEQUENCE CHAR: '%c' (%d)\n\n",c,c);
					ONDBG{ Rprintf("%c",c);}
					curseqlen++;
					c=nextChar();
				}//..while( not EOL )

				//
			}//...line began with an 's'
			else// if( c == '\n' )
			{
				numnucleotides += oldseqlen;
				curseqlen=oldseqlen=0;
			}

		}//while inside alignment-block

		ONDBG Rprintf("numnucleotides = %d + %d = %d\n",numnucleotides,oldseqlen,numnucleotides+oldseqlen);
		numnucleotides += oldseqlen;
		oldseqlen=0;

	}//while not eof
	
	//
	/* commented out due to CRAN warnings #FIXME
	ONDBG Rprintf("MAF File containing %d distinct sequences; total alignment length = %d\n",numsamples,numnucleotides);
	numsamples=numsamples;
	*/
	//
	return true;
}










SEXP	processAlignmentMAF( void )
{
	//
	RMatrix			resmat_inst;
	SEXP			resmat	=	R_NilValue
						;
	char			c=1,
					seqnamebuf[128]
						;
	int				*mat=0,
					*matbegin=0,
					*matend=0
						;
	uint			usedsequences[128],
					usedsequenceid=0;		//since MAF does not necessarily provide all sequences per alignment-block, we need to fill the missing parts with -
	unsigned int	numnucleotides=0,
					numsamples=0,
					si=0,
					curseqlen=0,
					oldseqlen=0,
					ni,
					currentalignmentoffset=0	//where the current alignment-block begins
					;
	
	//
	//				
	blkidx=0;
	fileidx=0;	
	if( false == determineAlignmentDimensionsMAF( numnucleotides, numsamples ) )
		return R_NilValue;

	//	allocate matrix
	//
	if( false == resmat_inst.alloc(INTSXP,numsamples,  numnucleotides) )
	{
		Rprintf("(!!) MAF error : Matrix alloc failure!\n	numsamples=%d\nnumnucleotides=%d\n",numsamples,numnucleotides);
		return R_NilValue;
	}
	mat = resmat_inst.getIntPtr();//mat = allocIntRMatrix( numsamples,  numnucleotides, &resmat );
	if( mat == 0 )
	{
		Rprintf("(!!) MAF error : Matrix alloc failure!\n	numsamples=%d\nnumnucleotides=%d\n",numsamples,numnucleotides);
		resmat_inst.dealloc();
		return R_NilValue;
	}
	resmat = resmat_inst.get();
	int *debug_mat_offs_begin = mat;
	matbegin=	mat;
	matend	=	&mat[numsamples*numnucleotides];

	//	read first N bytes of the file
	//
	blkidx=0;
	fileidx=0;
	if( 0 == readFileBlock( 0 ) )
	{
		Rprintf( "(!!) MAF error : Couldn't read first block of data from file!\n" );
		resmat_inst.dealloc();
		return R_NilValue;
	}
	
	//	reset the defined-sequences array - used to track which sequences are included in an alignment-block
	//
	for( unsigned int ii=0;ii<numsamples;usedsequences[ii++]=-1)
		;

	//	skip header lines (starting with a #)
	//
	while( (c = skipLine()) == '#' )
		/* No Other oPeration ;) */;

	//
	//	scan through the file until EOF
	//
	int		global_seqpos=0;
	while( c != 0 )
	{

		//
		//----------------	find next alignment block
		//

		//	find first line beginning with "a " == alignment block begin
		//
		while( c != 'a' && c != 0 )
			c = skipLine();

		if( c == 0 || nextChar() != ' ' )
			break;//found the E.O.F.

		//
		//----------------	find next sequence
		//

		//
		//
		while( (c=skipLine()) != '\n' && c!='\r' && c != 0 )
		{
			//	find first "s " line, count its bytelength including \n, count its sequencelength
			//	count number of "s " lines
			//		skip ahead by bytelength bytes to speed up
			//		skip other "? " lines without counting them
			//		at the first empty line, stop counting
			//
	
			//	repeat this process ("a ", "s ", count seqlengths,"\n") until end of file
			//
			//
			//

			//ONDBG Rprintf("At %d of %d\n",matend - matbegin,numsamples*numnucleotides);
			//ONDBG Rprintf("Alignmentblock : Line begins with '%c' (%d)\n",currentChar(),currentChar());
			
			//
			//
			if( c == 's' )
			{
				
				//
				//----------------	skip whitespace character after /^s/
				//

				if( (c=nextChar()) != ' ' )
				{
					Rprintf("(!!) MAF error : Line expected to be a sequence line, but does not begin with 's '\n");
					return R_NilValue;			
				}
				
				//
				//----------------	get name of sequence
				//

				si=0;
				while( (c=nextChar()) != ' ' && si < sizeof(seqnamebuf) )
				{
					seqnamebuf[si++]=c;
				}
				
				//	
				//
				if( si >= sizeof(seqnamebuf) )
				{
					Rprintf("(!!) MAF error : Sequence name length exceeds available buffer size %d\n",sizeof(seqnamebuf));
					return R_NilValue;
				}
				
				//
				seqnamebuf[si]=0;
				
				//	get sequence index
				//
				if( (ni=rownames.hasname( seqnamebuf )) >= 0 && ni < numsamples )
				{
					ONDBG Rprintf("	Sequence name '%s' is sequence #%d\n",seqnamebuf,ni);
					//TODO: mark sequence ni as contributing in this alignment block
					
				}
				else
				{
					Rprintf("(!!) MAF error : Sequence name '%s' not encountered during initial scan or %d >= %d!\n",seqnamebuf,ni,numsamples);
					//TODO : release allocated matrix
					return R_NilValue;
				}
				
				usedsequences[ni]=usedsequenceid;

				//
				//----------------	skip past the extra fields : '     2048 63 + 246127941 '[seq begins directly after that]
				//
				bool	hadsign=false;
				while( (c=nextChar())!=0 )
				{
					//ONDBG Rprintf("[%c]",c);
					if( isDecDigit(c) )	continue;
					if( c == ' ' )	continue;
					if( c == '+' ){ hadsign=true;continue;}
					if( c== '-' ){ if( hadsign )break; hadsign=true; continue;}
					if( c=='\r' || c=='\n' )
					{
						Rprintf("(!!) MAF error : Unexpected, premature EOL CR or LF (%d)\n",c);
						//TODO : release matrix
						return R_NilValue;
					}
					break;
				}
				
				//
				//----------------	copy mapped nucleotide codes to matrix elements
				//

				//
				//

				mat = &matbegin[ni];

				//	copy-loop
				//
				ONDBG Rprintf("ni=%d,'%s':\nSeq=",ni,seqnamebuf);
				curseqlen=0;
				while( c != '\r' && c != '\n' && c != 0 && c != ' ')
				{
					ONDBG Rprintf("%c[cs=%d]",c,curseqlen);
					ONDBG if( ni == 0 )
						Rprintf("*********>>>>>>	%d\n",(mat - debug_mat_offs_begin)/numsamples );
					mat[0] = (int)nucleotide_mapping[(unsigned)c];
					curseqlen++;
					c=nextChar();
					mat+=numsamples;
				}
				ONDBG Rprintf("[ %d old, %d this]\n",oldseqlen,curseqlen);
				
				//
				//
				if( oldseqlen==0 )
					oldseqlen=curseqlen;
				else
					if( curseqlen != oldseqlen && ni > 0 )
					{
						Rprintf("(!!) MAF error : Error in file : sequence length in alignment block mismatch : last=%d <=> %d=current, (ni=%d) \n",ni,oldseqlen,curseqlen);
						return R_NilValue;//resmat;
					}
				
				//
			}//...if( c == 's' )
			else// if( c == '\n' )
			{
				global_seqpos+=oldseqlen;
				ONDBG Rprintf("NON-/^s/-line	[%d]!\n",global_seqpos);
				break;
			}

		}//while( inside alignment-block )
		ONDBG Rprintf("End of Alignment block /////////////////////////\n");
	
		//
		//
		c = (int)nucleotide_mapping['-'];
		for( unsigned int ii=0; ii < numsamples; ii++)
		{
			if( usedsequences[ii] != usedsequenceid )
			{
				ONDBG Rprintf("Filling up empty space of Sequence #%d (%d nucleotides),global_seqpos=%d\n",ii,oldseqlen,global_seqpos);
				mat = &matbegin[ii];
				ONDBG Rprintf("matoffsstart=%d\n",mat-debug_mat_offs_begin);
				for( int jj=oldseqlen; jj>0; jj-- )
				{
					//Rprintf("jj=%d: %d\n",jj,mat-debug_mat_offs_begin);
					ONDBG Rprintf("jj=%d: %d\n",jj,((mat-matbegin)-ii)/numsamples);
					if( mat > matend )
					{
						Rprintf("(!!) MAF error : fill-unused-positions would write %d beyond matrix end!\n",matend-mat);
						break;
					}
					else
						mat[0]=c;
					mat+=numsamples;
				}
			}
		}
	
		//
		//		
		matbegin += oldseqlen*numsamples;//&matbegin[ numsamples * currentalignmentoffset];
		currentalignmentoffset += oldseqlen;
		curseqlen=oldseqlen=0;
		usedsequenceid++;

		//

	}//...while( not eof )
	
	
	//
	//
	setMatrixRownames( resmat , numsamples );

	//
	return resmat;
}



#endif //ifdef SUPPORT_MAF
