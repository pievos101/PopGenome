/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**
**
**
**

	TODO:






**
*/

//*
//*			INCLUDES
//*


	//	standard includes
	//
#include	<stdlib.h>				//string comparisons,...      
#include	<sys/time.h>			//gettimeofday

#define		STRICT_R_HEADERS		//to solve problems on Windows with Realloc, ERROR definitions made by R headers

	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>

	//
#include	<string>
#include	<vector>


//*
//*			DEFINES
//*


#define			Rboolean_TRUE	((Rboolean)1)

		//	--	debug mode related defines
		//
		//
#define			DEBUG			0

#if DEBUG
#	define		ONDBG			if( true )
#else
#	define		ONDBG			if( false )
#endif

	//
	//
#define		EXPORT				extern "C"

	//
	//
#define		guard				try{
	
#define		unguards			}catch(const char*err){Rprintf("EXCEPTION CAUGHT:'%s' in file %s line %d\n",err,__FILE__,__LINE__);}catch(...){Rprintf("UNKNOWN EXCEPTION CAUGHT in file %s line %d\n",__FILE__,__LINE__);}



//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*


	//	*******************	*************	******************
	//	*******************	*************	******************
	//	*******************	*************	******************




//*
//*			DATA
//*


//*
//*			CODE
//*

//*
//*			EXTERNS
//*


	//
	//	internal functions related to filtering
	//

/*! Checks wether REF and ALT define a SNP (i.e. 1 single-nucleotide REF allele, up to 3 alt )
**
*/
inline bool	_internal_isSNP( const char * refptr, const char* altptr )
{
	//
	if( refptr[1] != '\t' )	//make sure its not a deletion
		return false;
		
	//	a SNP can have more than one alternate allele, but no multi-nucleotide alleles!
	//		=> test for letter followed by comma followed by letter [...] followed by \t
	
	//---------------------
	//examples:
	//
	//REF=A
	//
	//	and 
	//
	//ALT=C_		or
	//ALT=C,T_		or
	//ALT=C,T,G_
	//---------------------
	// ? can '.' or N alt-alleles happen? (i.e. deletion)
	
	int i=0;
	for( ; altptr[i] != 0 && altptr[i] != '\t'; i+=2 )
	{
		if( false == ((altptr[i] >= 'A' && altptr[i] <='Z') || (altptr[i]>='a' && altptr[i]<='z')) )
			break;
		if( altptr[i+1] == '\t' )
			return true;
		if( altptr[i+1] != ',' )
			break;
	}
	return false;
}


/*! Checks wether REF and ALT define a biallelic SNP (i.e. single-nucleotide REF and ALT alleles)
**
*/
inline bool	_internal_isBiallelicSNP( const char * refptr, const char* altptr )
{
	if( refptr[1] == '\t' && altptr[1] == '\t' )//order of tests might accelerate tests
			return true;
	return false;
}

class vcff;
bool	filterLine( vcff* f );


	//	- utility functions -
	//
EXPORT	SEXP	WhopDebugLevel( void );
EXPORT	SEXP	SetWhopDebugLevel( SEXP lev );

		void	df0( const char *fmt, ...);
#ifndef DEBUG
#	define		df1		//
#else
		void	df1( const char *fmt, ...);
#endif
		void	df(int level, const char *fmt, ...);

//		bool	_internal_isSNP( const char * refptr, const char* altptr );

	//	Open/Close
	//
EXPORT	SEXP	VCF_open( SEXP filename );
EXPORT	SEXP	VCF_close( SEXP vcfptr );
EXPORT	SEXP	VCF_eor( SEXP vcfptr );

	//	Region
	//
EXPORT	SEXP	VCF_setRegion( SEXP vcfptr, SEXP tid, SEXP beg, SEXP end );
EXPORT	SEXP	VCF_getRegion( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionTid( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionBegin( SEXP vcfptr );
EXPORT	SEXP	VCF_getCurrentRegionEnd( SEXP vcfptr );
EXPORT	SEXP	VCF_restartRegion( SEXP vcfptr );

	//	Samples
	//
EXPORT	SEXP	VCF_getSampleNames( SEXP vcfptr );
EXPORT	SEXP	VCF_selectSamples( SEXP vcfptr, SEXP samplesvec );
EXPORT	SEXP	VCF_getSelectedSamples( SEXP vcfptr );


	//	Reading
	//
EXPORT	SEXP	VCF_readLineRaw( SEXP vcfptr, SEXP str );
EXPORT	SEXP	VCF_readLineRawFiltered( SEXP vcfptr, SEXP str );
EXPORT	SEXP	VCF_readLineTSV( SEXP vcfptr );
EXPORT	SEXP	VCF_readLineTSVFiltered( SEXP vcfptr );
EXPORT	SEXP	VCF_readLineDF( SEXP vcfptr );
EXPORT	SEXP	VCF_readLineDFFiltered( SEXP vcfptr );
EXPORT	SEXP	VCF_getBial( SEXP vcfptr, SEXP mat );


EXPORT	SEXP	VCF_readIntoNucleotideMatrix( SEXP vcfptr, SEXP mat );
EXPORT	SEXP	VCF_readIntoCodeMatrix( SEXP vcfptr, SEXP mat );

EXPORT	SEXP	VCF_parseNextLine( SEXP vcfptr );
EXPORT	SEXP	VCF_parseNextSNP( SEXP vcfptr );

		//
EXPORT	SEXP	VCF_getChrom( SEXP vcfptr );
EXPORT	SEXP	VCF_getPos( SEXP vcfptr );
EXPORT	SEXP	VCF_getID( SEXP vcfptr );
EXPORT	SEXP	VCF_getRef( SEXP vcfptr );
EXPORT	SEXP	VCF_getAlt( SEXP vcfptr );
EXPORT	SEXP	VCF_getQual( SEXP vcfptr );
EXPORT	SEXP	VCF_getFilter( SEXP vcfptr );
EXPORT	SEXP	VCF_getInfo( SEXP vcfptr );
EXPORT	SEXP	VCF_getInfoField( SEXP vcfptr, SEXP fieldnam );
EXPORT	SEXP	VCF_getFormat( SEXP vcfptr );
EXPORT	SEXP	VCF_getSample( SEXP vcfptr, SEXP stridx );


	//	Misc. utility
	//

EXPORT	SEXP	VCF_countSNPs( SEXP vcfptr );
EXPORT	SEXP	VCF_countBiallelicSNPs( SEXP vcfptr );
EXPORT	SEXP	VCF_isSNP( SEXP vcfptr );
EXPORT	SEXP	VCF_isInDel( SEXP vcfptr );

	//	Filters
	//
EXPORT	SEXP	VCF_describeFilterConfig( SEXP vcfptr );
EXPORT	SEXP	VCF_addFilter( SEXP vcfptr, SEXP fieldnam, SEXP cmptype, SEXP action, SEXP arg1, SEXP arg2 );
EXPORT	SEXP	VCF_clearFilters( SEXP vcfptr );

	//		Tabix for R
	//
EXPORT	SEXP	tabix_open( SEXP filename );
EXPORT	SEXP	tabix_close( SEXP tabix );
EXPORT	SEXP	tabix_reopen( SEXP tabix );
EXPORT	SEXP	tabix_setRegion( SEXP tabix, SEXP tid, SEXP begin, SEXP endpos );
EXPORT	SEXP	tabix_restartRegion( SEXP tabix );
EXPORT	SEXP	tabix_getRegion( SEXP tabix );
EXPORT	SEXP	tabix_readLine( SEXP tabix );


	//		GTF <nonfunctional>
	//
EXPORT	SEXP	gtftest( SEXP vcfptr );
	
	//
	//
#include	"whop_tools.h"

	//
	//
#include	"whop_rsupport.h"

	//
	//
#include	"whop_tabix.h"

	//
	//
#include	"whop_vcf.h"

