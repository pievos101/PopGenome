/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		main module - dll exports
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
#include	<stdarg.h>			//for df1(), df0() debugging functions


//*
//*			DEFINES
//*



#define		CHKTYPE( nam, type, len )		( is##type(nam) && (length(nam) == 1) )
#define		CHKSTR( nam, len )				CHKTYPE( nam, String, len )
#define		CHKINT( nam, len )				CHKTYPE( nam, Integer, len )


//*
//*			STRUCTS
//*


//*
//*			CLASSES
//*



//*
//*			DATA
//*

int			debug_level	= 0;


SEXP		vcfhandle_attrname_filename = R_NilValue;

//*
//*			EXTERNS
//*



//*
//*			CODE
//*





SEXP _internal_VcfGetAttrFilename( void )
{
	if( vcfhandle_attrname_filename == R_NilValue )
	{
		vcfhandle_attrname_filename = install("VCF.filename");
	}
	return vcfhandle_attrname_filename;
}

	//
	//		debug output control
	//



/*!
*/
EXPORT	SEXP	WhopDebugLevel( void )
{
	Rprintf("WhopGen : Current debug level is %d\n",debug_level);
	return R_NilValue;
}



/*!
*/
EXPORT	SEXP	SetWhopDebugLevel( SEXP lev )
{
	if( false == RNumeric::isInt( lev ) )
		return R_NilValue;

	debug_level = RNumeric::getInt( lev );
	
	return R_NilValue;
}




/*!
*/
void	df(int level, const char *fmt, ...)
{
	if( debug_level < level )
		return;
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf( dfbuf, sizeof(dfbuf),fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}




/*!
*/
void	df0( const char *fmt, ...)
{
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	//vfprintf(stderr, fmt, ap);
	vsnprintf( dfbuf, sizeof(dfbuf), fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}




/*!
*/
void	df1( const char *fmt, ...)
{
	if( debug_level < 1 )
		return;
	char dfbuf[1024];
	va_list ap;
	va_start(ap, fmt);
	//vfprintf(stderr, fmt, ap);
	vsnprintf( dfbuf, sizeof(dfbuf), fmt, ap);
	va_end(ap);
	Rprintf("%s",&dfbuf[0] );
}


		//
		//	finalizers
		//

/*!
**
**
**
*/
static void vcff_finalize(SEXP extPtr)
{
	df1("VCFF FINALIZE!\n");
	vcff* f = (vcff*)R_GetExtPtr( extPtr, "VCFhandle" );
	if( f )
	{
		df1("vcff_finalize : Finalizing VCFhandle!\n");
		VCF_close( extPtr );
		df1("vcff_finalize : Successfully finalized VCFhandle!\n");
	}
	else
	{
		df1("vcff_finalize : Could not finalize potential VCFhandle!\n");
	}

	//
}



		//
		//
		//


/*!	Open a Tabix-indexed VCF file and return a VCFhandle, a whopgen-managed EXTPTR SEXP
**
**
**
*/
EXPORT  SEXP VCF_open( SEXP filename )
{
	//
	if(! CHKSTR(filename,1) )
	{
		df0("VCF_open : filename is not a single string!");
		return R_NilValue;
	}
	
	//-----
	
	//open vcf with filename
	//	if failed: return NULL
	//
	vcff * f = new vcff( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("VCF_open: Could not open file '%s' as tabix-indexed!\n", CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}
	
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("VCF_open : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}
	
	//
	//
	//
#if 0
	int numseqs = f->getNumSequenceNames();
	for( int i=0; i < numseqs; i++ )
	{
		//
		const char * thisseqnam = f->getSequenceName(i);
		Rprintf("#%d=%s\n",i,thisseqnam);
		
		unsigned int	minr=0,maxr=2*1024*1024*1024;
		
		//
		for( unsigned int low = minr; low < maxr; low += 100000 )
		{
			const char * s = f->readNextLine();
			if( s )
			{
				while( *s != '\t' && *s != 0 )
					s++;
				int pos = atoi(s);
				
			}
		}
		
		//
	}
#endif

	//
	df1("(VCF_open) opened file '%s' is a VCF!\n",CHAR(STRING_ELT(filename, 0)));
	
	//
	//
	SEXP res;
	PROTECT(
		res = R_MakeExternalPtr( f, install("VCFhandle"), R_NilValue)
	);
	if( res == R_NilValue )
	{
		df0("VCF_open : could not create external pointer SEXP!\n");
		return res;
	}
	
	//
	R_RegisterCFinalizerEx(res, vcff_finalize, Rboolean_TRUE);
	
	//
	setAttrib( res, _internal_VcfGetAttrFilename(), filename );
	
	UNPROTECT( 1 );
	return res;

}


/*!	Close a VCFhandle, an EXTPTR SEXP of a whopgen-managed Tabix-indexed VCF file 
**
**
*/
EXPORT SEXP	VCF_close( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_close : parameter is not a VCFhandle or nil!\n");
		return RBool::False();
	}
	
	//
	R_ClearExternalPtr(vcfptr);
	
	//
	delete f;
	
	//
	return RBool::True();
}





/*!	Reopens a VCf file if the VCFhandle got stale
**
**
**
**
**
*/
EXPORT SEXP	VCF_reopen( SEXP vcfptr )
{
	//
	if( false == RType::IsExtPtr( vcfptr ) || ( strcasecmp( RExtPtr::getTag(vcfptr), "VCFhandle" ) != 0 ) )
	{
		df1("VCF_reopen : parameter is not an externalptr VCFhandle!\n");
		return RBool::False();
	}
	
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( f != 0 )
	{
		return RBool::True();
	}
	
	//
	SEXP filename = getAttrib( vcfptr, _internal_VcfGetAttrFilename() );
	
	//
	//open vcf with filename
	//	if failed: return False
	//
	f = new vcff( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df0("VCF_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}
	
	//	if loading failed, return with error
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df0("VCF_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}

	R_SetExternalPtrAddr( vcfptr, f );
	
	//
	return RBool::True();
}







/*!	Returns the next header-line, which was preparsed and stored in
**		 the VCFhandle object.
**	Returns NilValue, if last header-line was returned but also resets
**		the counter to begin at the first line again with next invocation.
**
**
*/
EXPORT SEXP	VCF_getHeaderLine( SEXP vcfptr, SEXP num )
{

	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_getHeaderLine : Parameter 1 is not a VCFhandle!\n");
		return R_NilValue;
	}
	
	//
	if( false == isInteger( num ) )
	{
		df0("VCF_getHeaderLine : parameter 2 needs to be an integer!\n");
		return R_NilValue;
	}

	//
	int idx = INTEGER(num)[0];
	const char * s = f->getHeaderLine( idx );
	if( s != 0 )
	{
		RString res(s);
		return res.get();
	}

	//
	df1("No header line #%d to get!\n", idx );
	return R_NilValue;
}


	//
	//		Contig Identifiers / Chromosome Names (e.g. "1", "22", "chr2" )
	//
	//


/*!
**
**
**
*/
EXPORT SEXP	VCF_getContigNames( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getContigNames : argument 1 is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	int numseqs = f->getNumSequenceNames();
	//Rprintf("%d seqs\n",numseqs);
	
	//
	RString res;
	res.alloc( numseqs );
	for( int i=0; i < numseqs; i++ )
	{
		//Rprintf("#%d=%s\n",i,f->getSequenceName(i));
		res.set( f->getSequenceName(i), i );
	}
	
	//
	return res.get();
	
	//
}





/*!
**
**
**
*/
EXPORT SEXP	VCF_getNumContigs( SEXP vcfptr )
{
	//
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getNumContigs : argument is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	//
	RInteger	res( f->getNumSequenceNames() );
	return res.get();
}





		//
		//			VCF fields
		//
		//


/*!
**
**
**
*/
EXPORT  SEXP VCF_getFieldNames( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getFieldNames : argument 1 is not a VCF!\n");
		return R_NilValue;
	}
	
	//
	//
	RString st;
	st.alloc( f->getNumFields() );
	
	//
	for( unsigned int k=0; k < f->getNumFields(); k++ )
	{
		st.set( f->getFieldName( k ), k );
		//Rprintf("Token = '%s'\n",f->getFieldName( k ) );
	}
	return st.get();
	
	//
}

