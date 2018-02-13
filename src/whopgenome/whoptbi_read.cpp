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



//*
//*			EXTERNS
//*



//*
//*			CODE
//*


/*!	Just read the next line, 
**
*/
const char*	whop_tabix::readNextLine( void )
{
	int len=0;
	if( iter != 0 )
	{
		const char * s = ti_read( tabix, iter, &len);
		bEOR = (s==0);
		return s;
	}
	bEOR=true;
	return 0;
}


/*!	Read next line from VCF and build an index-list into the string
**		with the offsets of each field.
**
**	- minimal preprocessing to quickly access a specific field in the line
**	- First field is always at offset 0
**
*/
bool			whop_tabix::parseNextLine( void )
{
	//
	if( iter == 0 )
		return false;
	
	//@TODO remove
	if( field_offsets == 0 )
	{
		Rprintf("No field_offsets !\n");
		return false;
	}
	
	//
	//
	current_line = ti_read( tabix, iter, &current_line_len);
	if( current_line == 0 )
	{
		bEOR = true;
		return false;
	}

	//	count number of fields
	//
	unsigned int numfields = 0 ;
	unsigned int strpos = 0;
	for( int tmplen = current_line_len ; current_line[strpos] != 0 && tmplen > 0; strpos++, tmplen-- )
	{
		if( current_line[strpos] == '\t' )
		{
			numfields++;
		}
	}

		//@TODO eliminate these checks by catching lines with differing numbers of fields

	//
	if( numfields > field_offsets_size )
	{
		Rprintf("(!!) whop_tabix::parseNextLine : ERROR : %d > %d : new numfields > field_offsets_size!\n",numfields,field_offsets_size);
		return false;
	}

	//	find the starting indices of each field
	//
	field_offsets[0]=0;
	strpos = 0;
	unsigned int	idxpos = 1;
	for( int tmplen = current_line_len ; current_line[strpos] != 0 && tmplen > 0; strpos++, tmplen-- )
	{
		//
		if( current_line[strpos] == '\t' )
		{
			//
			if( idxpos > field_offsets_size )
			{
				Rprintf("(!!) whop_tabix::parseNextLine : ERROR : More fields in this line than expected! (%d>%d)\n",idxpos,last_num_fields);
				return false;
			}
			
			//
			field_offsets[idxpos]=strpos+1;
			idxpos++;
		}
	}
	
	//
	last_num_fields = numfields;

	//
	return true;
}


/*!	Copies a field of the last parseNextLine()'d line into a buffer and ends it with a null-byte
**
**
**
*/
bool			whop_tabix::copyField( unsigned int fieldidx, char*buffer, unsigned int maxbuflen )
{
	if( current_line == 0 )
		return false;
	if( field_offsets == 0 )
		return false;
	if( fieldidx >= last_num_fields )
		return false;
	if( buffer == 0 )
		return false;
	if( maxbuflen < 1 )
		return false;

	//
	int fldoffs = field_offsets[ fieldidx ];
	const char * str = &current_line[ fldoffs ];
	int i=0;
	for( ; i < (maxbuflen-1); i++ )
	{
		if( str[i] == '\t' || str[i] == 0 )
			break;
		buffer[i] = str[i];
	}
	buffer[i]=0;
		
	return true;
}


/*!	Returns a pointer to the text of field nr. <fieldidx>
**
**	- field-text ends with a tab-character, not a null-byte !
**	- fieldidx is 0-based : first element = 0
**	- requires a previous parseNextLine() call and therefore, a region having been set
**
*/
const char*	whop_tabix::getFieldPtr( unsigned int fieldidx )
{
	if( current_line == 0 )
	{
		Rprintf("(!!) whop_tabix::getFieldPtr : did not read a line of data from the Tabix-file yet!\n");
		return 0;
	}
	if( field_offsets == 0 )
	{
		Rprintf("(!!) whop_tabix::getFieldPtr : did not determine field-offsets for this Tabix-File!\n");
		return 0;
	}
	if( fieldidx > last_num_fields )
	{
		//fieldidx >= last_num_fields <--> %d >= %d!\n",fieldidx , last_num_fields);
		Rprintf("(!!) whop_tabix::getFieldPtr : requested field %d but only %d present!\n", fieldidx , last_num_fields);
		return 0;
	}
	//
	int fldoffs = field_offsets[ fieldidx ];
	const char * str = &current_line[ fldoffs ];
	return str;
}

