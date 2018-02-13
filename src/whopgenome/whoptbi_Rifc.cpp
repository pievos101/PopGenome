/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		Tabix for R
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


SEXP		tabixhandle_attrname_filename = R_NilValue;

//*
//*			EXTERNS
//*



//*
//*			CODE
//*


/*!	Helper function to get a valid SEXP for the symbol "Tabix.filename"
**
*/
SEXP _internal_TabixGetAttrFilename( void )
{
	if( tabixhandle_attrname_filename == R_NilValue )
	{
		tabixhandle_attrname_filename = install("Tabix.filename");
	}
	return tabixhandle_attrname_filename;
}


/*!
**
**
**
*/
static void tabix_finalize(SEXP extPtr)
{
	df1("TABIX FINALIZE!\n");
	whop_tabix* f = (whop_tabix*)R_GetExtPtr( extPtr, "TabixHandle" );
	if( f )
	{
		df1("Finalizing TabixHandle!\n");
		tabix_close( extPtr );
		df1("Successfully finalized TabixHandle!\n");
	}
	else
	{
		df1("Could not finalize potential TabixHandle!\n");
	}

	//
}


/*!
**
**
*/
EXPORT	SEXP	tabix_open( SEXP filename )
{
	//
	//
	if(! ( isString(filename) && (length(filename) == 1) ) )
	{
		error("tabix_open : filename is not a single string!");
		return R_NilValue;
	}

	//open tabix with filename
	//	if failed: return NULL
	//
	whop_tabix * f = new whop_tabix( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df1("tabix_open : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}

	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df1("tabix_open : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return R_NilValue;
	}

	//
	df1("tabix_open : opened file '%s' is a Tabix-indexed file!\n",CHAR(STRING_ELT(filename, 0)));
	
	//
	//
	SEXP res;
	PROTECT(
		res = R_MakeExternalPtr( f, install("TabixHandle"), R_NilValue)
	);
	if( res == R_NilValue )
	{
		df1("tabix_open : could not create external pointer SEXP!\n");
		return res;
	}
	
	//
	R_RegisterCFinalizerEx(res, tabix_finalize, Rboolean_TRUE );
	
	//
	setAttrib( res, _internal_TabixGetAttrFilename(), filename );
	
	UNPROTECT( 1 );
	return res;
}



/*!
**
**
*/
EXPORT	SEXP	tabix_close( SEXP tabix )
{
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( 0 == f )
	{
		df1("tabix_close : parameter is not a TabixHandle or nil!\n");
		return RBool::False();
	}
	
	//
	R_ClearExternalPtr(tabix);
	
	//
	delete f;
	
	//
	return RBool::True();
}





/*!	Reopens a tabix file if the TabixHandle got stale (e.g. after crash and restart of R)
**
**
**
**
**
*/
EXPORT SEXP	tabix_reopen( SEXP tabix )
{
	//
	if( false == RType::IsExtPtr( tabix ) || ( strcasecmp( RExtPtr::getTag(tabix), "TabixHandle" ) != 0 ) )
	{
		df1("tabix_reopen : parameter is not an externalptr TabixHandle!\n");
		return RBool::False();
	}
	
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( f != 0 )
	{
		return RBool::True();
	}
	
	//
	SEXP filename = getAttrib( tabix, _internal_TabixGetAttrFilename() );
	
	//
	//open tabix with filename
	//	if failed: return NULL
	//
	f = new whop_tabix( CHAR(STRING_ELT(filename, 0)) );
	if( f == 0 )
	{
		df1("tabix_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}
	
	//	if loading failed, return with error
	//
	if( false == f->isValid() )
	{
		delete f;
		f=0;
		df1("tabix_reopen : Could not open file '%s' as tabix-indexed!\n",CHAR(STRING_ELT(filename, 0)));
		return RBool::False();
	}

	R_SetExternalPtrAddr( tabix, f );
	
	//
	return RBool::True();
}



/*!
**
**
*/
EXPORT	SEXP	tabix_setRegion( SEXP tabix, SEXP tid, SEXP begin, SEXP endpos )
{
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	//
	if( false == RString::isStr( tid ) )
	{
		df0("tabix_setregion : 'tid', is not a string!\n");
		return RBool::False();
	}
	
	//
	return f->setRegion( RString::get(tid), RNumeric::getInt(begin), RNumeric::getInt(endpos) ) ? 
				RBool::True() : RBool::False();

	//
}



/*!
**
**
*/
EXPORT	SEXP	tabix_getRegion( SEXP tabix )
{
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	return R_NilValue;
}






/*!
**
**
**
*/
EXPORT	SEXP tabix_restartRegion( SEXP tabix )
{
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( 0 == f )
	{
		Rprintf("tabix_restartRegion : parameter 1 needs to be valid TabixHandle!\n");
		return RBool::False();
	}
	
	if( f->restartRegion() )
		return RBool::True();
	
	return RBool::False();
	
}



/*!
**
**
*/
EXPORT	SEXP	tabix_readLine( SEXP tabix )
{
	//
	whop_tabix * f = (whop_tabix*)R_GetExtPtr( tabix , "TabixHandle" );
	if( 0 == f )
	{
		return R_NilValue;
	}
	
	//
	const char * s = f->readNextLine();
	if( s != 0 )
	{
		//
		RString resvec( s );
		return resvec.get();
	}
	return R_NilValue;
}


