/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		Tabix wrapper class implementation
**
**
**

	TODO:




**
*/




//*
//*			INCLUDES
//*

#include	"whop_tabix_internal.h"


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



//*
//*			EXTERNS
//*



//*
//*			CODE
//*





/*!
**
*/
bool			whop_tabix::isValid( void )
{
	return(
		( tabix != 0 )
			&&
		( index != 0 )
			&&
		( num_seqnames!=0 )
	);
}





/*!
**
*/
const char*	whop_tabix::getSequenceName( unsigned int idx )
{
	if( idx >= num_seqnames )
		return 0;
	if( sequence_names == 0 )
		return 0;
	return sequence_names[idx];
}



/*!	Parse header lines of opened tabix file and call the specified functor for each line.
**		Returns TRUE as soon as functor returns true, does not continue parsing
**		Returns false if functor did not return true for any header line
**
*/
bool			whop_tabix::parseHeader( ParseFunctor &f )
{
	
	//	validate arguments
	//
	// comment out due to CRAN warnings #FIXME	start
	/*
	if( this == 0 )
	{
		Rprintf("whop_tabix::parseheader : NULL vcf*!\n");
		return false;
	}

	if( &f == 0 ){
		Rprintf("whop_tabix::parseheader : NULL functor!\n");
		return false;
	}
	*/
	// end fixme 
	// end of quick fix CRAN  

	//
	for( unsigned int i=0 ; i < header_lines.size() ; i++ )
	{
		const char * s = header_lines[i].c_str();
		int len = strlen(s);
		if( true == f(s,len) )
		{
			return true;
		}
	}
	
	//
	//Rprintf("VCF_parseheader false[%s]\n",s);
	return false;
}







/*!
**
*/
const char*	whop_tabix::getHeaderLine( unsigned int index )
{
	if( index >= header_lines.size() || index < 0 )
		return 0;
	return header_lines[ index ].c_str();
}
