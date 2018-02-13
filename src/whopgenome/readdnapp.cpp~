/*
**
**		ReadDNAplusplus
**
**		replacement for read.dna for use in PopGen
**
**
**
**
**
**

	TODO:
	
The input data files are multiple aligned DNA sequence data in a number of interleaved formats as
MAF, MGA, XMFA, PHYLIP, or the HapMap genotype format.
MAF (http://genome.ucsc.edu/goldenpath/help/maf.html).
MGA (http://bibiserv.techfak.uni-bielefeld.de/mga/).
XMFA (http://lagan.stanford.edu/).
PHYLIP (http://evolution.genetics.washington.edu/phylip.html).
HapMap genotype format (http://www.hapmap.org/downloads/encode1.html.en)

**
*/

//*
//*			INCLUDES
//*

#include <stdio.h>
#include <stdlib.h>

	//for printing 64-bit values
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

	//
	//
#include "whop_rsupport.h"

//*
//*			DEFINES
//*

#define			SUPPORT_FASTA
#define			SUPPORT_PHYLIP
#define			SUPPORT_MAF
//#define			SUPPORT_MEGA

//#define			DEBUG

	//
	//
#define			PLATFORM_BITS		32

#	define		_FILE_OFFSET_BITS PLATFORM_BITS

	//	
	//
#ifdef			DEBUG
#	define			ONDBG			if( true )
#else
#	define			ONDBG			if( false )
#endif

	//	maximum number of bytes of a file that is loaded into memory at once
	//
#define			MAX_IMMEDIATE_LOAD_HARD_LIMIT		__readdnapp_internal_maxblockloadsize

	//	error reporting macros
	//

#define DEBUG

	//
	//

typedef		unsigned int			uint;

//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*

class	RGuard {
public:

	//
	RGuard() : unprotectcount(0) {
	}
	
	//
	~RGuard(){
		if( unprotectcount )
			UNPROTECT( unprotectcount );
	}
	
	//!
	inline SEXP	allocVec( int type, int num ){
		SEXP res = R_NilValue;
		PROTECT( res = allocVector( type, num ) );
		unprotectcount++;
		return res;
	}
	
	//!
	inline SEXP	allocMatrix( int type, int numrows, int numcols ){
		SEXP res = R_NilValue;
		PROTECT( res = allocMatrix( type, numrows, numcols ) );
		unprotectcount++;
		return res;
	}
	/*
	PROTECT(dimnamesvec = allocVector(VECSXP, 2));
	PROTECT(rownamesvec = allocVector(STRSXP, num_rows));
	PROTECT(ans = allocMatrix(INTSXP, numrows, numcolumns));
	PROTECT(ans = allocMatrix(INTSXP, numsamples, numnucleotides));
			PROTECT(currentfilename = allocVector(STRSXP, 1));
*/
	int	unprotectcount;
};


//!
class dynstorage			//
{
public:

	//
	dynstorage(unsigned int inblocksize ){
		blocksize=inblocksize;
		currentblock = 0;
		currentblockpointer=0;
		currentblockindex=0;
		top=0;
		enlarge();
	}
	
	//! ONLY CALL IF CURRENT BLOCK IS EXHAUSTED - OTHERWISE REMAINING BYTES ARE INVALID DATA AND WILL BE CONSIDERED VALID!
	void	enlarge( void )
	{
		chunkheader * newptr;
		if( currentblock == 0 || currentblock->next == 0 )
		{

			newptr = (chunkheader*)malloc( blocksize+sizeof(chunkheader)-1 );
											//-1 because chunkheader contains 1 databyte
			if( newptr == 0 )
				throw "dynstorage : failed malloc to enlarge !\n";
			if( currentblock )
				currentblock->next = newptr;
			newptr->next = 0;
		}
		else
			newptr = currentblock->next;
		
		//switch to next block
		//
		currentblock = newptr;
		currentblockpointer = &currentblock->databegin[0];
		currentblockindex = 0;
		
		//
		//
		if( top == 0 )
			top = newptr;

		//
	}//enlarge
	
	//
	struct chunkheader {
		chunkheader		*next;
		char			databegin[1];
	};
	
	//! Returns whether the given string at <seqname> is stored here
	int		hasname( const char * seqname )
	{
		int		numm=0;
		char	cc=1;
		int		i=0;
		int		si=0;
		int		ni=0;

		while( get(i,cc) )
		{
			//
			if( cc )
			{
//				Rprintf("[%c]",cc);
				if( seqname[si++] == cc )
				{
					numm++;
				}
			}
			else
			{
//				Rprintf("(%d <=> %d)\n",si,numm);
				if( seqname[si] == 0 && numm==si )
					return ni;
				si=0;
				numm=0;
				ni++;
			}
		
			i++;
		}

		return -1;
	}

	
	//! Get character number <index> in the datablock
	inline	bool	get( int index, char&c )
	{
		int numchunk = index/blocksize;
		unsigned int	numbyte = index%blocksize;
		chunkheader * p = top;
//		printf("get[%d:c%d/b%d]",index,numchunk,numbyte);
		while( p && numchunk-- )
		{
			p = p->next;
		}
		
		//
		if( p == 0 || (p==currentblock && numbyte >= currentblockindex) )
		{
			c=0;
			return false;
		}
		
		//
		c = p->databegin[numbyte];
		
		//
		return true;
		
		//
	}

	//! Store <numbytes> characters stored at <p>
	inline void		serialize( char * p, int numbytes=1 )
	{
	
		//
		//
		while( numbytes > 0 )
		{

			// determine the maximum number of bytes we can copy
			int numbytescopyable = numbytes>(blocksize-currentblockindex) ? (blocksize-currentblockindex) : numbytes;
			
			//
			memcpy(currentblockpointer,p, numbytescopyable);
			
			//
			currentblockpointer += numbytescopyable;
			currentblockindex += numbytescopyable;
			p += numbytescopyable;
			numbytes -= numbytescopyable;
			
			if( currentblockindex >= blocksize )
				enlarge();
		
		}
	
		
		//
	}//serialize( buffer )
	
	//
	inline	void	reset(void)
	{
		chunkheader * p = top;
		while( p )
		{
			memset( &p->databegin[0], 0, blocksize );
			p = p->next;
		}
		currentblock=top;
		currentblockpointer=&currentblock->databegin[0];
		currentblockindex=0;
	}
	
	
	//!
	inline	void	serialize( char c )
	{
//		printf("ser[%02d:%d,%d] -> '%s'\n",c,currentblockindex,blocksize,&top->databegin[0]);
		if( currentblockindex >= blocksize )
			enlarge();
		*currentblockpointer = c;
		currentblockpointer++;
		*currentblockpointer=0;
		currentblockindex++;
	}//serialize( char )
	
	//
	chunkheader			*top;
	chunkheader			*currentblock;
	char				*currentblockpointer;
	unsigned int		currentblockindex;
	
	unsigned int		blocksize;
};


//*
//*			DATA
//*

int				returncode = 0;

FILE			*filehandle=0;
// long long		filebytelength=0;
// Bastian change

unsigned int            filebytelength=0;

char			*filedatabuffer = 0;
unsigned int	currentfileposition = 0;
int				globalAlignLoadError = 0;
bool			atEof = false;

	//
	//
static	unsigned char		*memorybuffer = 0;
static	unsigned int		memorybufferlength = 0;
static	unsigned int		memorybuffervalidsize = 0;
static	unsigned int		blkidx=0,
							fileidx=0;

	//
	//
unsigned int __readdnapp_internal_maxblockloadsize = ( 1000 * 1024 *1024 );


/*

		ids <- c("T" ,"t",	"U","u",	"C","c",	"G","g",	"A","a",	"N","n","?",	"-")
		nuks <- c(1,1,		1,1,		2,2,    	3,3,		4,4,    	 5,5,5,			6)
*/

//!	For quick mapping of nucleotide character codes into numeric constants for the biallelic matrix
static char	nucleotide_mapping[] = {
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 5	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 16	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,6,6,5,				// 32	: !"#$%&'()*+’-./
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 48	:5123456789:;<=>?
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 64	:@ABCDEFGHIJKLMNO
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 80	:PQRS TUVW XYZ[\]^_
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 96	:`abcdefghijklmno
	5,5,5,5,	1,1,5,5,	5,5,5,5,	5,5,5,5,				// 112	:pqrs tuvw xyz{|}~
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

static char	is_iupac_nucleotide[] = {
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 5	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 16	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,1,1,0,				// 32	: !"#/$%&'_()*+_’-./
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 48	:0123_4567_89:;_<=>?
	0,1,0,1,	1,0,0,1,	1,0,0,1,	0,1,1,0,				// 64	:@ABC_DEFG_HIJK_LMNO
	0,0,1,1,	1,1,1,1,	0,1,0,0,	0,0,0,0,				// 80	:PQRS_TUVW_XYZ[_\]^_
	0,1,0,1,	1,0,0,1,	1,0,0,1,	0,1,1,0,				// 96	:`abc_defg_hijk_lmno
	0,0,1,1,	1,1,1,1,	0,1,0,0,	0,0,0,0,				// 112	:pqrs_tuvw_xyz{_|}~

	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 128	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 144	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 160	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 176	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 192	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 208	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,				// 224	: - mostly nonprintable - 
	0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0					// 240	: - mostly nonprintable - 
};

//!
dynstorage		rownames(2048);

unsigned int	maxloopiters=1000;

//*
//*			EXTERNS
//*



//*
//*			CODE
//*



		//************************
		//************************		File helper functions
		//************************
		//************************
/*!
**
**
*/
void*		loadFile( void )
{
	if( filehandle == 0 )
		return 0;
		
	//
	//
	if( filedatabuffer )
	{
		free( filedatabuffer );
		filedatabuffer=0;
	}

	//
	//
	filedatabuffer = (char*)malloc( filebytelength );
	if( filedatabuffer==0 )
	{
		Rprintf("(!!) Failed to allocate %lld bytes to load file into memory!\n",filebytelength);
		return 0;
	}

	//
	//		
	int numbytesloaded = fread( filedatabuffer, 1, (size_t)filebytelength, filehandle );
	if( numbytesloaded < filebytelength )
	{
		Rprintf("(!!) Only %d bytes of %llu could be read!\n",numbytesloaded,filebytelength);
		free( filedatabuffer );
		return 0;
	}
	//
	
	return filedatabuffer;
}


/*
**
*/
uint	readFileBlock( uint startoffs )
{
	if( filehandle == 0 )
	{
		ONDBG Rprintf("(!!) readdnapp::readFileBlock(%d) : filehandle is Zero!\n",startoffs);
		return 0;
	}
	if( memorybuffer == 0 )
	{
		ONDBG Rprintf("(!!) readdnapp::readFileBlock(%d) : memorybuffer is Zero!\n",startoffs);
		return 0;
	}
	if( currentfileposition != startoffs )
		if( 0 != fseek( filehandle, startoffs, SEEK_SET ) )
			return 0;
	currentfileposition = startoffs;
	
	//
	size_t numreadbytes = fread( memorybuffer , 1, memorybufferlength, filehandle );
	memorybuffervalidsize = numreadbytes;
	currentfileposition += numreadbytes;
	blkidx=0;
	fileidx=startoffs;
	atEof = false;
	return numreadbytes;
}




/*!	Ensures that for the requested position in a file the correct block of data is available in the memory buffer
**	[IN]	curpos = current position in the memory block
**	[IN]	filepos = current in-file start-offset of the data in the memory block
**	[OUT]	returns whether the file block could be read
**	[OUT]	implicitly : updates the in-file start-offset of the data that has been loaded into the memory buffer
**
**	[OUT]	implicitly : updates curpos to 0 = beginning of memory buffer, if successful
*/
inline	bool	updateFileBlock( unsigned int &curpos , unsigned int &filepos )
{
//	L( ufB___ );	L( curpos= ) ;I(curpos);L( validsize= );I( memorybuffervalidsize );L( filepos= ) ;I(filepos);LN;
	if( curpos >= memorybuffervalidsize )
	{
		unsigned int nextpos = 	filepos;// + memorybuffervalidsize;
		//L( nextpos= ) ;I(nextpos);L( \nfilebytelength= );I( filebytelength );LN;
		if( (nextpos+curpos) >= filebytelength )		//we would be beyond the end of the file
		{
//			L( nextpos GE filebytelength );LN;
			curpos=filebytelength;
			filepos=filebytelength;
			atEof=true;
			return false;
		}
			//remember: readFileBlock auto-updates 'memorybuffervalidsize' with the actual number of bytes read
		if( 0 == readFileBlock( nextpos ) )	//did not get any data from the file ?
		{
			ONDBG Rprintf("(!!) readdnapp::updateFileBLock(curpos=%d,filepos=%d) : readFileBlock == 0\n",curpos,filepos );
			atEof=true;
			return false;
		}
		filepos += memorybuffervalidsize;
	}
	atEof=false;
	return true;
}

	//
	//
	//

//!	Returns the filesize of the given file
/*
long long	filesize( FILE * fh )
{
	long long res = 0;
	long long curpos = ftello(fh);
	if( fseeko(fh,0,SEEK_END) == 0 )
		res = ftello(fh);
	fseeko(fh,curpos,SEEK_SET);
	return res;
}
*/

//my change Bastian

unsigned int	filesize( FILE * fh )
{
	unsigned int res = 0;
	unsigned int curpos = ftello(fh);

	if( fseeko(fh,0,SEEK_END) == 0 )
		res = ftello(fh);
	fseeko(fh,curpos,SEEK_SET);
	return res;
}




/*!
*/
inline bool		isDecDigit( char c )
{
	return( (c>='0' && c<='9') );
}

/*!
*/
inline	char	currentChar( void )
{
	return memorybuffer[blkidx];
}


/*!
*/
inline	char	nextChar( void )
{
	blkidx++;
	if( blkidx >= memorybuffervalidsize )
	{
		if( false==updateFileBlock( blkidx, fileidx ) )
		{
			return 0;
		}
		blkidx=0;
	}
	return memorybuffer[blkidx];
}


/*!
*/
inline bool	Eof( void )
{
	return atEof;
}


/*!
**
*/
inline char	skipWhitespaces( void )
{
	char c = currentChar();
	while( (c == ' ') || (c == '\t') || (c == '\n') || ( c== '\r') )
	{
		c=nextChar();
	}
	return c;
}


/*!
**
*/
inline char	skipLine( void )
{
	char c = currentChar();
	while( (c != '\n') && ( c!= 0) )
	{
		c=nextChar();
	}
	c = nextChar();
	return c;
}



/*!
**
*/
inline int scanUInt( void )
{
	int res=0;
	char c = currentChar();
	while( c && (c >= '0') && (c <= '9') )
	{
		res = (res*10)+(c-'0');
		c = nextChar();
	}
	return res;
}

/*!
**
*/
inline bool	cmpString( const char * str )
{
	char c = currentChar();
	while( *str != 0 && c == *str ){
		c=nextChar();
		str++;
	}
	return (*str==0)?true:false;
}




/*!
*/
inline bool	isNucleotide( unsigned char c )
{
	return( c=='A' || c=='a' || c=='C' || c=='c' || c=='G' || c=='g' || c=='T' || c=='t' || c=='-' || c=='N' || c=='n' || c=='U' || c=='u' );
}

inline bool	isIUPACNucleotide( unsigned char c )
{
	return ( is_iupac_nucleotide[c] != 0 );
}

		//************************
		//************************			Matrix helper functions
		//************************
		//************************
		

	//
	//
	//

bool	setMatrixRownames( SEXP mat, int num_rows )
{
	SEXP		dimnamesvec,
				rownamesvec
					;
	
	//	allocate vectors to set rownames
	//
	PROTECT(dimnamesvec = allocVector(VECSXP, 2));
	PROTECT(rownamesvec = allocVector(STRSXP, num_rows));

	//	set rownames as vector elements
	//
	char	rownambuf[512];
	unsigned int		rownambufidx=0;
	int		smplnum=0;
	char	cc=0;
	int		i=0;
	while( rownames.get(i,cc) && num_rows > smplnum )
	{
		//
		if( cc )
		{
			if( rownambufidx < sizeof(rownambuf)-1 )
				rownambuf[rownambufidx++] = cc;
		}
		else
		{
			rownambuf[rownambufidx++] = 0;
			SET_STRING_ELT(rownamesvec,smplnum,mkChar(&rownambuf[0]));
			rownambufidx=0;
			smplnum++;
		}
		
		i++;
	}
	
	//	set dimnames
	//
	SET_VECTOR_ELT(dimnamesvec, 0, rownamesvec);
	setAttrib(mat, R_DimNamesSymbol, dimnamesvec);
	
	UNPROTECT(2);
	
	return true;
	
	//
}

	//
	//
	//

SEXP	copyToMatrix( const char *codemat, int numrows, int numcolumns )
{

	//
	SEXP		ans = R_NilValue
					;
	int*		rans;
	int			i,j;


	//	allocate matrix
	//
	//
	PROTECT(ans = allocMatrix(INTSXP, numrows, numcolumns));

	//	copy data
	//
	rans = INTEGER(ans);
	for(i = 0; i < numrows; i++)
	{
		for(j = 0; j < numcolumns; j++)
			rans[i + numrows*j] = *codemat++;
	}
	
	//
	//
	setMatrixRownames(ans, numrows);

	//
	UNPROTECT(1);
	
	return ans;
}

/*!	Allocates an integer R matrix and returns a pointer to its first element
**
**
int	num_ans_rows  = 0;
int*	allocIntRMatrix( int numrows, int numcolumns , SEXP*sexpstore )
{
	int* rans = 0 ;
	SEXP	ans=R_NilValue;
	
	//
	PROTECT(ans = allocMatrix(INTSXP, numrows, numcolumns));
	if( sexpstore )
		*sexpstore = ans;
	
	//
	if( ans != R_NilValue )
	{
		rans = INTEGER(ans);		//get int-pointer to matrix data so that we can fill right now
	}
	
	//
	num_ans_rows = numrows;
	
	//
	UNPROTECT(1);
	
	//
	return rans;
}
*/



		//************************
		//************************			Format loaders
		//************************
		//************************
		

	//
	//
	//
#include		"rdnapp_Fasta.h"


	//
	//
	//
#include		"rdnapp_Phylip.h"


	//
	//
	//
#include		"rdnapp_MAF.h"


	//
	//
	//
#include		"rdnapp_MEGA.h"



		//************************
		//************************			----------
		//************************
		//************************

	//
	//
	//
	
SEXP	processAlignmentAny( void )
{
	SEXP	res = R_NilValue;
	
	//
//FIXME : processAlignmentXXX should set globalAlignLoadError to > 0 to allow this function to try a different method, to < 0 to cancel altogether

	//	Fasta
	//
#ifdef SUPPORT_FASTA
	//Rprintf("Test for FASTA\n");
	res = processAlignmentFasta();
	if( res!=R_NilValue )
	{
		//Rprintf("FASTA recognized\n");
		return res;
	}
#endif


	//	Phylip
	//
#ifdef SUPPORT_PHYLIP
	//Rprintf("Test for PHYLIP\n");
	if( res==R_NilValue )
	{
		res = processAlignmentPhylip();
		if( res!=R_NilValue )
		{
		//	Rprintf("Phylip recognized\n");
			return res;
		}
	}
#endif


	//	MAF
	//
#ifdef SUPPORT_MAF
	//	'##maf version=1 ... \n'
	//
	//Rprintf("Test for MAF\n");
	if( res==R_NilValue )
	{
		res = processAlignmentMAF();
		if( res!=R_NilValue )
		{
		//	Rprintf("MAF recognized\n");
			return res;
		}
	}

#endif

	//	Mega
	//
#ifdef SUPPORT_MEGA
	//	'#mega\n'
	//
	//Rprintf("Test for MEGA\n");
	if( res==R_NilValue )
	{
		res = processAlignmentMEGA();
		if( res!=R_NilValue )
		{
		//	Rprintf("MEGA recognized\n");
			return res;
		}
	}
#endif


	//
	return res;
}


/*!
**
**
**
**
*/
extern "C"  SEXP rdnapp_SetMaxBlockSize( SEXP newmaxblocksize )
{

	int *nmbs = INTEGER( newmaxblocksize );
	if( nmbs == 0 || *nmbs < 0 )
		return R_NilValue;
	__readdnapp_internal_maxblockloadsize = *nmbs;
	return newmaxblocksize;
}

	//
	//
	//

/*!
**
**
**
**
*/
extern "C"  SEXP readdna(SEXP filenames )
{
 
 	//
 	//
	int			i=0
					;
	char		*filename;
					
	
	SEXP		ans = R_NilValue;
	SEXP		currentfilename
					;

//	//
//	//Test YES/NO-Table for 'is ASCII char a IUPAC nucleotide code?'	
//	//
//	for( int i=32; i < 128; i++ )
//	{
//		Rprintf("%03d='%c' -> %s\n",i,i, isIUPACNucleotide(i)?"YES":"no");
//	}
	
	//
	currentfilename = STRING_ELT(filenames,i);
	filename = (char*)CHAR(currentfilename);

	//	open file
	//
	filehandle = fopen( filename ,"rt");
	if( filehandle )
	{
	
		//	load complete file
		//
		filebytelength = filesize( filehandle );
		ONDBG { Rprintf("File='%s'\nFilesize=%llu\n",filename,filebytelength); }

		//
		//
		if( filebytelength <= 0 )
			return(ans);
			
		//
		atEof = false;
			
		//	pre-allocate a memory buffer for loading (parts of) the file into, but do not load anything yet
		//
		if( memorybuffer == 0 )
		{
			memorybuffer = (unsigned char*)malloc( MAX_IMMEDIATE_LOAD_HARD_LIMIT );
			if( memorybuffer == 0 )
			{
				Rprintf(" ReadDNA++ : failed to allocate %d bytes of memory buffer for file '%s'\n",MAX_IMMEDIATE_LOAD_HARD_LIMIT, filename);
				return ans;
			}
			memorybufferlength = MAX_IMMEDIATE_LOAD_HARD_LIMIT;
		}
		
		//
		//
		rownames.reset();
	
		//
		//
		ans = processAlignmentAny();
		
		if( ans == R_NilValue )
		{
			Rprintf("(!!)	Translation failed in file %s!\n",filename);
		}
		else
		{
			//	set the 'path' attribute
			//
			PROTECT(currentfilename = allocVector(STRSXP, 1));
			SET_STRING_ELT(currentfilename,0,mkChar(filename));
			setAttrib(ans, install("path"), currentfilename);
			UNPROTECT(1);
			
			//
			//
			setMatrixRownames( ans, RMatrix::numRows(ans) );
			
			//
		}

		//
		free( filedatabuffer );
		filedatabuffer = 0;

		//
		fclose( filehandle );

	}//if ( filehandle )
	else
		Rprintf("(!!) Could not open file for reading: '%s'\n",filename);


	//
	return(ans);
	
	//
}




