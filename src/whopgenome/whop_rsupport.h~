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
//#include <stdio.h>
//#include <vector>
//#include <string>

#include <vector>
#include <string>
#include <strings.h>

//*
//*			DEFINES
//*




//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*

/*!
**
**
*/
class RType
{
public:
	static const char*	TypeName( SEXP v ){	return Rf_type2char( TYPEOF(v) ); }
	static bool			IsType( SEXP v, const char*name ){ return( strcasecmp( TypeName(v), name ) == 0 ); }
	
	static	bool			IsSymbol( SEXP v ){	return Rf_isSymbol(v);	}
	static	bool			IsExtPtr( SEXP v ){ return (TYPEOF(v) == EXTPTRSXP) ; }
};


/*!	Wraps around R EXTPTRSEXP, simplifies working with external pointers in R
**
**
*/
class RExtPtr
{
public:
	RExtPtr();
	RExtPtr(SEXP v, const char* tagname );
/*
	 void *R_ExternalPtrAddr(SEXP s);
     SEXP R_ExternalPtrTag(SEXP s);
     SEXP R_ExternalPtrProtected(SEXP s);
     void R_ClearExternalPtr(SEXP s);
     void R_SetExternalPtrAddr(SEXP s, void *p);
     void R_SetExternalPtrTag(SEXP s, SEXP tag);
     void R_SetExternalPtrProtected(SEXP s, SEXP p);
*/

	//
	//

	SEXP		getTag( void );
	SEXP		getProtected( void );
	void*		getPointer( void );
	
	static	bool			IsExtPtr( SEXP v ){ return (TYPEOF(v) == EXTPTRSXP) ; }
	//
	//
	static	const char*	getTag( SEXP v )
	{
		if( false == RType::IsExtPtr(v) )
			return (const char*)0;
		SEXP r = ::R_ExternalPtrTag( v );
		PROTECT( r = coerceVector( r, STRSXP ) );
		const char * res = CHAR( STRING_ELT( r, 0 ) );
		UNPROTECT(1);
		return res;
	}

protected:
	SEXP		extptr;
};


/*!	Wraps around R STRSEXP, simplifies working with strings
**
*/
class RString
{
public:

	//
	RString();
	RString(unsigned int numelems){ alloc(numelems);	}
	RString(const char * str);
//	RString(SEXP s);
	~RString();
	
	//!
	static	bool			isStr( SEXP var );
	static	const char*	get( SEXP var , unsigned int idx=0 );
	static	unsigned int	length( SEXP s );
	
	//!
			bool			alloc( int len=1 );
			SEXP			get( void );
			SEXP			getElem( int elemidx=0 );
			bool			set( const char * str , int idx=0 );
			unsigned int	length( void );
	
	//
protected:
	int		_length;
	SEXP	_value;
};



/*!
**
**
**
*/
class RNumeric
{
public:
	static	bool	isInt( SEXP v );
	static	bool	isFloat( SEXP v );
	static	int		getInt( SEXP v );
	static	float	getFloat( SEXP v );
};



/*!
**
**
**
*/
class RInteger
{
public:
	RInteger();
	RInteger(int value);
	~RInteger();
	SEXP	get( void );
	int		getAsInt( void );
	void	set( int v );
protected:
	SEXP	_value;
};



/*!
**
**
**
*/
extern SEXP	_const_false,_const_true;
class RBool
{
public:
	RBool();
	RBool( bool value );
	RBool( SEXP sexpvalue );
	~RBool();
	
	//
	void		set( bool v );
	SEXP		get( void ){	return _value;}
	bool		getValue( void );
	
	//
	
	//
	static	bool		getValue( SEXP v );
	
	static	SEXP		True( void )
	{
		//if( _const_true == R_NilValue )
		{
			PROTECT( _const_true = allocVector( LGLSXP , 1 ) );
			INTEGER( _const_true )[0] = 1;
			UNPROTECT(1);
		}
		return _const_true;
	}
	static	SEXP		False( void )
	{
		//if( _const_false == R_NilValue )
		{
			PROTECT( _const_false = allocVector( LGLSXP , 1 ) );
			INTEGER( _const_false )[0] = 0;
			UNPROTECT(1);
		}
		return _const_false;
	}

	//
protected:
			SEXP		_value;
};


/*!
**
**
**
*/
class RMatrix
{
public:
	RMatrix();
	RMatrix( SEXP mat );
	RMatrix( int SEXPTYPE, int numcols, int numrows );
	~RMatrix();
	
	//
	//
	
	//
	bool			alloc( int SEXPTYPE, int numcols, int numrows );
	void			dealloc( void );
	bool			isValid( void );
	
	//
	SEXP			get( void );
	void			set( SEXP );
	
	//
	unsigned int	numRows( void );
	unsigned int	numCols( void );
	int				getType( void );

	SEXP			getRowNames( void ){ return getRowNames(_matvar);	}
	SEXP			getColNames( void ){ return getColNames(_matvar);	}
	
	bool			setRowNames( SEXP rownames ){ return setRowNames(_matvar,rownames);	}
	bool			setColNames( SEXP colnames ){ return setColNames(_matvar,colnames);	}
	
	bool			setRowNames( std::vector<std::string>& rownames ){ return setRowNames(_matvar,rownames);	}
	bool			setColNames( std::vector<std::string>& colnames ){ return setColNames(_matvar,colnames);	}
	
	//
	int*			getIntPtr( void );
	double*			getDoublePtr( void );
	SEXP*			getStrPtr( void );
	
	//
	void			clearD( double clearvalue=0.0 );
	void			clearI( int clearvalue=0 );
	
	//
	//		--	static methods
	//
	static	bool			isMatrix( SEXP m );
	static	unsigned int	numRows( SEXP m );
	static	unsigned int	numCols( SEXP m );
	
	static	SEXP	getRowNames( SEXP m );
	static	SEXP	getColNames( SEXP m );
	
	static	bool	setRowNames( SEXP m , SEXP rownames );
	static	bool	setColNames( SEXP m , SEXP colnames );
	
	static	bool	setRowNames( SEXP m , std::vector<std::string>& rownames );
	static	bool	setColNames( SEXP m , std::vector<std::string>& colnames );
	
protected:
	SEXP		_matvar;
	bool		self_alloced;
};



/*!
**
**
**
*/
class RList
{
public:
	RList();
	RList( int length );
	~RList();
	
	//
	//
	
	//
	//
	static	bool	isList( SEXP m );
	static int		length( SEXP m );
	
protected:
	SEXP		_listvar;
	bool		self_alloced;
};




//*
//*			DATA
//*



//*
//*			EXTERNS
//*



void	setListElement(SEXP list, char *str, SEXP value);
SEXP	getListElement(SEXP list, const char *str);
SEXP	getvar(SEXP name, SEXP rho);
SEXP	createList( SEXP r );

	//
	//
void*			R_GetExtPtr( SEXP var , const char* expectedname );
const char*	R_GetExtPtrTag( SEXP var );

	//	data.frame
	//
SEXP	createDataFrame( unsigned int numrows, unsigned int numcolumns,
							const char** rownames, const char**columnnames );


