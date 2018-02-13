/*
**
**		WhopGen
**
**		WHOle genome Population genetics with popGEN
**
**
**		controlling region of VCF
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



//*
//*			EXTERNS
//*



//*
//*			CODE
//*



/*!	Set the region from which the following calls to readNextLine/parseNextLine should
**		return lines.
**
**	- tid STRSXP
**
**
*/
EXPORT SEXP	VCF_setRegion( SEXP vcfptr, SEXP tid, SEXP beg, SEXP end )
{

	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_setRegion : argument is not a VCF!\n");
		return RBool::False();
	}

	//
	if( false == RString::isStr( tid ) )
	{
		Rprintf("VCF_setRegion : argument 1, 'tid', is not a string!\n");
		return RBool::False();
	}

	//
	//
	const char* tidcstr = RString::get(tid);
	int frompos = RNumeric::getInt(beg);
	int topos = RNumeric::getInt(end);
	if( frompos <= 0 || topos <= 0 )
	{
		Rprintf("VCF_setRegion : unexpected values for start(%d), end(%d)\n",frompos,topos);
		return RBool::False();
	}
	
	//
	bool setregsuccess = f->setRegion( tidcstr, frompos, topos );
	return setregsuccess ? RBool::True() : RBool::False();

	//
}






/*!
**
**
**
*/
EXPORT	SEXP VCF_getCurrentRegionTid( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getCurrentRegionTid : parameter 1 needs to be valid VCFhandle!\n");
		return R_NilValue;
	}
	
	//
	const char * regtid = f->getRegionTid();
	if( regtid == 0 )
		return R_NilValue;
	
	RString res( regtid );
	return res.get();
}






/*!
**
**
**
*/
EXPORT	SEXP VCF_getCurrentRegionBegin( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getCurrentRegionBegin : parameter 1 needs to be valid VCFhandle!\n");
		return R_NilValue;
	}
	
	RInteger res( f->getRegionBegin() );
	return res.get();
}






/*!
**
**
**
*/
EXPORT	SEXP VCF_getCurrentRegionEnd( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getCurrentRegionEnd : parameter 1 needs to be valid VCFhandle!\n");
		return R_NilValue;
	}
	
	RInteger res( f->getRegionEnd() );
	return res.get();	
}





/*!	Return a "contigid:beginpos-endpos" formatted string describing the currently set region
**
**
*/
EXPORT SEXP VCF_getRegion( SEXP vcfptr )
{
	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getRegion :: Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}
	
	const char * regtid = f->getRegionTid();
	unsigned int regbeg = f->getRegionBegin();
	unsigned int regend = f->getRegionEnd();

	{
		char	buffer[128];
		snprintf( buffer, sizeof(buffer)-2, "%s:%u-%u",regtid,regbeg,regend);
		buffer[ sizeof(buffer)-1 ] = 0;
		RString res( buffer );
		return res.get();
	}
	
	
	//
	return R_NilValue;
}






/*!	Restart the region so that the next read-call will return the first entry in the region again
**
**
**
*/
EXPORT	SEXP VCF_restartRegion( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_restartRegion : parameter 1 needs to be a valid VCFhandle!\n");
		return R_NilValue;//RBool::False();
	}
	
	if( f->restartRegion() )
		return RBool::True();
	
	return RBool::False();
	
}






/*!	Returns TRUE if all VCF-lines in the region have been read
**
**
**
*/
EXPORT	SEXP VCF_eor( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_eor : parameter 1 needs to be a valid VCFhandle!\n");
		return R_NilValue;//RBool::False();	//FIXME : maybe TRUE is a valid option here?
	}
	
	if( f->eor() )
		return RBool::True();
	
	return RBool::False();
	
}


