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



//*
//*			DEFINES
//*



//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*


/*!	Helper class to parse tab-serparated value (TSV) lines of text
**
**
*/
class TSVParser {
public:
	TSVParser();
	TSVParser( const char * str );
	
	//
	//
	unsigned int	numFields( void );
	bool			getField( char* buf, unsigned int buflen, unsigned int idx );

	//
	//
protected:
	const char*	strbegin;
	unsigned int	num_fields;
	unsigned int	*field_offsets;
};




/*!	Virtual base class to parse lines extracted from text files
**
*/
class ParseFunctor {
public:
	virtual	bool	operator()( const char * s, int len ){return false;};
};







inline int parseDecInt( char * &s )
{
	int res=0;
	while( s[0] >= '0' && s[0] <= '9' )
	{
		res = (res*10)+(s[0]-'0');
		s++;
	}
	return res;
}

inline void skipDecInt( char * & s )
{
	while( s[0] >= '0' && s[0] <= '9' )
	{
		s++;
	}
}
