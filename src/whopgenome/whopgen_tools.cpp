/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		Miscelleanous utility functions module
**
**
**

	TODO:

- parse fields of TabSeparatedValues-string
	- embed class methods of TSVParser



**
*/

//*
//*			INCLUDES
//*


#include	"whopgen_common.h"

	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>



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


TSVParser::TSVParser() : strbegin(0),num_fields(0),field_offsets(0)
{
}


TSVParser::TSVParser( const char * str ) : strbegin(str),num_fields(0),field_offsets(0)
{
	//	count the number of fields
	//
	int res=0, p=0;
	while( strbegin[p] != 0 )
	{
		if( strbegin[p] == '\t' )
			res++;
		p++;
	}
	num_fields = res+1;
	
	//
	//
//	Rprintf("%d fields\n",num_fields);
	field_offsets = (unsigned int*)malloc( sizeof(int) * num_fields );
	if( field_offsets == 0 )
	{
		Rprintf("Failed to malloc %d bytes!\n",sizeof(int) * num_fields);
		num_fields = 0;
		strbegin = 0;
		return;
	}
	
	//@FIXME handle malloc failure here!
	
	//	copy the offsets of the fields inside the string
	//
	p=0;
	field_offsets[0] = 0;
	res=1;
	while( strbegin[p] != 0 )
	{
		if( strbegin[p] == '\t' )
			field_offsets[res++]=p+1;
		p++;
	}

	//
}




/*!
**
**
*/
bool	TSVParser::getField( char* buf, unsigned int buflen, unsigned int idx )
{
	
	//
	if( buf == 0 )
	{
		Rprintf("(!!) TSVParser::getField : buf == 0!\n");
		return false;
	}
	//
	if( buflen < 1 )
	{
		Rprintf("(!!) TSVParser::getField : buflen == 0!\n");
		return false;
	}
	//
	if( idx >= num_fields )
	{
		Rprintf("(!!) TSVParser::getField : idx > num.Fields !\n");
		return false;
	}
	//
	if( field_offsets == 0 )
	{
		Rprintf("(!!) TSVParser::getField : field_offsets == 0 !\n");
		return false;
	}
	
	//
	//
	int offs = field_offsets[idx];
	if( offs < 0 )
	{
		Rprintf("(!!) TSVParser::getField : field_offsets == 0 !\n");
		return false;
	}
	
	//
	//
	const char * tokenaddr = &strbegin[ offs ];
	
	//
	//
	unsigned int i=0;
	for( ; i < buflen-1; i++ )
	{
		if( tokenaddr[i] == '\t' || tokenaddr[i] == 0 )
			break;
		buf[i] = tokenaddr[i];
	}
	buf[i]=0;
	
	//
	return true;
	
	//
}




/*!	Count the number of (tab-separated) fields in the given line
**
**
*/
unsigned int	TSVParser::numFields( void )
{
	return num_fields;
}
