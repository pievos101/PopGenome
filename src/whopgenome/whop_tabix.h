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


	//
#include	<string>
#include	<vector>

	//
#include	"../tabix/tabix.h"

//*
//*			DEFINES
//*


	//
	//
#define		CHROM				0
#define		POS					1
#define		ID					2
#define		REF					3
#define		ALT					4
#define		QUAL				5
		//optional
#define		FILTER				6
#define		INFO				7
#define		FORMAT				8


//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*


/*!	Class wrapping working with VCF files
**
**
**
*/
class whop_tabix
{
	//

public:
	//TODO : embed methods : open, close, parseheader, parseline, ...
	
	//
	//
	whop_tabix();
	whop_tabix(const char * filename );
	~whop_tabix();
	
	//-
	//-		METHODS
	//-
	
	
	//	Open + Close
	//
	bool			open( const char * filename );
	bool			isValid( void );
	bool			eor( void ){ return bEOR;	}

	//	Header
	//
	bool			parseHeader( class ParseFunctor &f );
	const char*	getHeaderLine( unsigned int index );

	//	Regions
	//
	bool			setRegion( const char* tid, int begin, int end );
	bool			restartRegion( void );	//restart region from beginning
	int				lastReadPos( void );
	const char*	getRegionTid( void );
	int				getRegionBegin( void );
	int				getRegionEnd( void );
	
	//	Parsing
	//
	const char*	readNextLine( void );	//just read
	const char*	getCurrentLine( void ){ return current_line; }
	bool			parseNextLine( void );	//minimal parsing to access fields
	bool			copyField( unsigned int fieldidx, char*buffer, unsigned int maxbuflen );
	const char*	getFieldPtr( unsigned int fieldidx );
	unsigned int	numParsedFields( void ){	return last_num_fields; }
	
	//	Sequence Names
	//
	unsigned int	getNumSequenceNames( void ){	return num_seqnames;	}
	const char*	getSequenceName( unsigned int idx );
	
	//-
	//-		DATA
	//-
	
protected:

	//	Tabix file format data
	//
	tabix_t*				tabix;
	ti_index_t*				index;
	ti_iter_t				iter;
	

	bool						bEOR;

	//
	//
	std::string					currentTid;
	int							currentBegin;
	int							currentEnd;

	//	Header data
	//
	std::vector<std::string>	header_lines;
	
	//	Sequencename data
	//
	unsigned int				num_seqnames;
	char const					**sequence_names;
	
	//	transient data while parsing a TSV-line
	//
	unsigned int				last_num_fields;	//
	unsigned int				*field_offsets;		//	offsets to the substrings of each TSV-field
	unsigned int				field_offsets_size;	//	how many entries are ready in field_offsets ?
	const char					*current_line;		//	start-address of the string containing the line
	int							current_line_len;	//	length of string at current_line

private:

	//
};

