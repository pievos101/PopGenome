/*
**
**
**			Phylip
**
**
**
*/



#ifdef SUPPORT_PHYLIP


/*!	
**
**
**
*/
SEXP	processAlignmentPhylip( void )
{
	//
	RMatrix	resmat_inst;
	SEXP	resmat = R_NilValue
				;
	bool	error=false;

	//
	unsigned int		numnucleotides=0,
						numsamples=0,
						numblockparsednucleotides=0,
						numlastseqblockparsednucleotides=0,
						numnucleotideslefttoparse
							;
	char				c;
	
	//
	//
	blkidx=0;
	fileidx=0;

	//	start copying matrix from the start
	//
	;
	if( (fileidx=readFileBlock( 0 )) == 0 )
	{
		ONDBG Rprintf("(!!) processAlignmentPhylip : updateFileBlock failure!\n\tblkidx = %d, fileidx = %d\n",blkidx,fileidx);
		return R_NilValue;
	}

	//	read "<numsamples> <numnucleotides>\n" line
	//
	skipWhitespaces();
	numsamples = scanUInt();
	skipWhitespaces();
	numnucleotides = scanUInt();
	numnucleotideslefttoparse = numnucleotides;
	
	unsigned int	numelems = numsamples * numnucleotides,
					numprocessedelems = 0;

	//	sanity check numsamples and numnucleotides
	//
	if( numsamples <= 1 || numnucleotides < 1 )
	{
		ONDBG Rprintf("PHYLIP : Not a Phylip file : invalid alignment dimensions: numsamples=%d<=1 , numnucleotides=%d<1\n",numsamples,numnucleotides);
		return R_NilValue;
	}
	
	//	allocate matrix
	//
	if( false == resmat_inst.alloc(INTSXP,numsamples,  numnucleotides) )
	{
		Rprintf("PHYLIP error : Matrix alloc failure!\nnumsamples = %d\nnumnucleotides = %d\n",numsamples,numnucleotides);
		return R_NilValue;
	}
	resmat = resmat_inst.get();
	int *mat = resmat_inst.getIntPtr();
	if( 0 == mat )
	{
		resmat_inst.dealloc();
		Rprintf("PHYLIP error : Matrix alloc failure!\nnumsamples = %d\nnumnucleotides = %d\n",numsamples,numnucleotides);
		return R_NilValue;
	}
	int	*matend =&mat[numelems];
	
	//
	int				*matrowptr;

	//	parse first block of the alignment : ("seqname   AAAAAAAAAAAATTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCC\n" * numsamples)\n
	//
	for( unsigned int samplenum=0; samplenum < numsamples ; samplenum++ )
	{
		//
		c = skipWhitespaces();

		//
		//
		for( int samplenameidx=0; c != ' ';samplenameidx++)
		{
			rownames.serialize((char*)&c,1);
			c=nextChar();
		}
		while( c == ' ' )	c=nextChar();

		//
		rownames.serialize((char)0);

		//
		matrowptr=&mat[samplenum];

		//
		//
		numlastseqblockparsednucleotides=0;
		for( ; c!='\n' && c!='\r' && c != 0 ; c=nextChar() )
		{
			if( matrowptr >= matend )
			{
				error=true;
				Rprintf("parsing error - premature end of file OR parsed beyond specified alignment length\n");
				Rprintf("\tlast char=%d\n",c);
				Rprintf("\tdistance in elements to end of matrix=%d\n",matend - matrowptr);
				resmat_inst.dealloc();
				return R_NilValue;
			}
			//
			*matrowptr = (int)(nucleotide_mapping[ (unsigned)c ]);
			ONDBG numprocessedelems++;
			numlastseqblockparsednucleotides++;
			matrowptr+=numsamples;
		}

		//
	}

	//
	//
	mat = matrowptr-numsamples+1;
	numblockparsednucleotides += numlastseqblockparsednucleotides;
	numnucleotideslefttoparse -= numlastseqblockparsednucleotides;


	//	parse following blocks of the alignment : (AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGG\n * numsamples)\n 
	//
	//
		//FIXME : do not just rely on updateFileBlock to help terminate this loop - especially,
		//		when the whole alignment is stored sequentially instead of interleaved!
	while( (numnucleotideslefttoparse > 0) && (numblockparsednucleotides < numnucleotides) )
	{

		//	process an N-sample block of the alignment
		//
		unsigned int firstsamplenumblocknucleotides = 0;
		for( unsigned int samplenum=0; samplenum < numsamples && c!= 0 ; samplenum++ )
		{
//			Rprintf("***SMP %d***\n",samplenum);
			//
			c = skipWhitespaces();
		
			//
			matrowptr=&mat[samplenum];

			//
			//
			numlastseqblockparsednucleotides=0;
			for( ; c != '\n' && c != '\r' && c != 0 ; c=nextChar() )
			{
				if( matrowptr >= matend )
				{
					error=true;
					Rprintf("parsing error - premature end of file OR parsed beyond specified alignment length\n" );
					Rprintf("\tlast char=%d\n",c);
					Rprintf("\tdistance in elements to end of matrix=%" PRIu64 "\n", matend - matrowptr );
					resmat_inst.dealloc();
					return R_NilValue;
				}
				*matrowptr = (int)(nucleotide_mapping[ (unsigned)c ]);
				ONDBG numprocessedelems++;
				matrowptr+=numsamples;
				numlastseqblockparsednucleotides++;
			}//...for all nucleotides in sample's sequence in the alignment sub-block
			if( samplenum == 0 ) firstsamplenumblocknucleotides = numlastseqblockparsednucleotides;
			else if( numlastseqblockparsednucleotides != firstsamplenumblocknucleotides )
			{
				Rprintf("PHYLIP : ERROR : Line in alignment block shorter for the %dth sample than for the others!\n"
						"	(%d = numlastseqblockparsednucleotides != firstsamplenumblocknucleotides = %d)\n",
						samplenum+1,numlastseqblockparsednucleotides , firstsamplenumblocknucleotides);
			}
		}

		//

		if( c == 0 )
		{
			ONDBG Rprintf("processALignmentPhylip : c=0\n");
			break;
		}
		//
		mat = matrowptr-numsamples+1;
		numblockparsednucleotides += numlastseqblockparsednucleotides;
		numnucleotideslefttoparse -= numlastseqblockparsednucleotides;
//		if( numlastseqblockparsednucleotides != 60 && numnucleotideslefttoparse > 0)
	//	{
		//	Rprintf("PHYLIP : ERROR: numlastseqblockparsednucleotides == %d != 60 while %d nucs left to parse!\n"
			//		"	(i.e. an alignment block has at least one line that is too short)\n\n",numlastseqblockparsednucleotides,numnucleotideslefttoparse);
			//return R_NilValue;
		//}

	}
	
	ONDBG Rprintf("numprocessedelems=%d\n",numprocessedelems);
	ONDBG Rprintf("numnucleotideslefttoparse=%d\n",numnucleotideslefttoparse);

	//
	//
	if( error )
	{
		resmat_inst.dealloc();
		return R_NilValue;
	}

	//
	//
	setMatrixRownames( resmat , numsamples );

	//
	return resmat;

	//
}



#endif	//ifdef SUPPORT_PHYLIP

