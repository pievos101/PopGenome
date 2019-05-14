/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		R support module
**
**
**

	TODO:

- easy create STRSXP vector
- easy set string in STRSXP-vector




**
*/

//*
//*			INCLUDES
//*


#include	"whopgen_common.h"
/*
	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>
*/


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




//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//**		RExtPtr
//***************************************************************



//!
void*	R_GetExtPtr( SEXP var , const char* expectedname )
{
	try
	{
		//
		//
		if( false == ( (TYPEOF(var) == EXTPTRSXP) ) )
		{
			df1("R_GetExtPtr : Not an extptrsxp!\n");
			return 0;
		}

		//
		//
		void * addr = R_ExternalPtrAddr( var );
		if( addr == 0 )
		{
			ONDBG df1("R_GetExtPtr : ExtPtr is Null!\n");
			return 0;
		}
		
		//
		//
		const char * tag = R_GetExtPtrTag( var );
		if( (0 == tag) || (strcmp( tag , expectedname ) != 0) )
		{
			ONDBG df1("R_GetExtPtr : ExtPtrTag is not '%s'!\n",expectedname);
			return 0;
		}
		
		//
		//
		ONDBG Rprintf("GetExtPTR : f=%x\n",addr);
		return addr;
		
		//
	}
	catch( const char* s )
	{
		df1("Exception trying to get extptr:\n\t%s\n",s);
	}
	catch( ... )
	{
		df1("Generic exception trying to get extptr address\n");
	}
	
	return 0;
}

//!
const char* R_GetExtPtrTag( SEXP var )
{
	try
	{
		//
		//
		if( false == ( (TYPEOF(var) == EXTPTRSXP) ) )
		{
			df1("not a extptr!\n");
			return 0;
		}

		//
		//
		SEXP tag = R_ExternalPtrTag(var);
		if( tag == R_NilValue )
		{
			return 0;
		}
		
		//
		//
		SEXP nam;
		PROTECT( nam = coerceVector( tag , STRSXP ) );
		if( nam == R_NilValue )
		{
			UNPROTECT( 1 );
			return 0;
		}

		//
		//
		const char * res = CHAR(STRING_ELT(nam,0));
		UNPROTECT( 1 );
		return res;
	}
	catch( const char* s )
	{
		df1("Exception trying to get extptr tag:\n\t%s\n",s);
	}
	catch( ... )
	{
		df1("Generic exception trying to get extptr tag\n");
	}
	//UNPROTECT( 1 );
	return 0;
}




//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//**		RString
//***************************************************************



/*!
**
*/
RString::RString() : _length(0),_value(R_NilValue)
{
}



/*!
**
*/
RString::RString(const char * str) : _length(0),_value(R_NilValue)
{
	alloc(1);
	set(str,0);
}



/*!
**

RString::RString(SEXP s) : _length(0),_value(R_NilValue)
{
	if( true == isStr( s ) )
	{
		_value = s;
		_length = ::Rf_length( s );
	}
}
*/


/*!
**
*/
RString::~RString()
{
	if( _value != R_NilValue )		//required to prevent garbage-collection when
		UNPROTECT(1);				//	time-consuming tasks are done until return of the using function
}



/*!
**
*/
unsigned int	RString::length( void )
{
	if( false == ::isString( _value ) )
		return 0;
	return( _length );
}




/*!
**
*/
unsigned int	RString::length( SEXP s )
{
	if( false == ::isString( s ) )
		return 0;
	return ::Rf_length( s );
}



/*!
**
*/
bool	RString::isStr( SEXP var )
{
	return ::isString( var );
}



/*!
**
*/
const char*	RString::get( SEXP var , unsigned int idx )
{
	if( false == ::isString( var ) )
		return 0;
	if( idx >= RString::length(var) )
		return 0;
	return CHAR(STRING_ELT( var, idx ));
}



/*!
**
*/
bool	RString::alloc( int len )
{
	if( len <= 0 )
	{
		df1("RString::alloc: Tried 0-length alloc on R-String!\n");
		return false;
	}
	if( _value != R_NilValue )	//make sure we unprotect the old value if we're going to alloc a new one
		UNPROTECT(1);
	PROTECT(_value = allocVector(STRSXP, len));
	//NOTE: cannot unprotect the value right now, otherwise it may become gc'd by R before this DLL returns
	_length=len;
	return true;
}



/*!
**
*/
SEXP	RString::get( void )
{
	return _value;
}



/*!
**
*/
SEXP	RString::getElem( int elemidx )
{
	if( (elemidx >= _length) || (elemidx < 0) || (_value == R_NilValue) )
		return R_NilValue;
	RString res;
	res.alloc( 1 );
	res.set( CHAR(STRING_ELT( _value, elemidx )), 0 );
	return res.get();
}



/*!
**
*/
bool	RString::set( const char * str , int idx ){
	if( _value == R_NilValue )
		return false;
	if( idx >= _length || idx < 0 )
		return false;
	SET_STRING_ELT( _value, idx, mkChar(str) );
	return true;
}


//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//**		RNumeric
//***************************************************************


/*! 
**
*/
bool	RNumeric::isInt( SEXP v )
{
	return isInteger( v );
}

/*! 
**
*/
bool	RNumeric::isFloat( SEXP v )
{
	return isReal( v );
}

/*! 
**
*/
int		RNumeric::getInt( SEXP v )
{
	if( isInt( v ) )
	{
		return INTEGER(v)[0];
	}
	else if( isFloat( v ) )
	{
		return (int)( REAL(v)[0] );
	}
	else if( RString::isStr( v ) )
	{
		const char * vstr = RString::get( v, 0 );
		if( vstr )
		{
			double res = atof(vstr);
			return (int)( res );
		}
		else
			Rprintf("(!!) RNumeric::getInt : cannot read an integral number from an empty string!\n");
	}
	else
		Rprintf("(!!) RNumeric::getInt : Trying to get an integer from a non-numeric datatype!\n");

	return 0;
}

/*! 
**
*/
float	RNumeric::getFloat( SEXP v )
{
	if( isInt( v ) )
	{
		return (float)( INTEGER(v)[0] );
	}
	else if( isFloat( v ) )
	{
		return REAL(v)[0];
	}
	else if( RString::isStr( v ) )
	{
		const char * vstr = RString::get( v, 0 );
		if( vstr )
		{
			double res = atof(vstr);
			return float( res );
		}
		else
			Rprintf("(!!) RNumeric::getInt : cannot read a floating-point number from an empty string!\n");
	}
	else
		Rprintf("(!!) RNumeric::getInt : Trying to get a floating-point number from a non-numeric datatype!\n");

	return 0.f;
}



//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//**		RBool
//***************************************************************

SEXP		_const_true = R_NilValue;
SEXP		_const_false = R_NilValue;

/*!
**
**
*/
RBool::RBool( void ) : _value(R_NilValue)
{
}



/*!
**
**
*/
RBool::RBool( bool value ) : _value(R_NilValue)
{
	RBool();
	set( value );
}

/*!
**
**
*/
RBool::RBool( SEXP sexpvalue ) : _value(sexpvalue)
{
	RBool();
}


/*!
**
**
*/
RBool::~RBool()
{
	if( _value != R_NilValue ){
		UNPROTECT(1);
	}
}


/*!
**
**
*/
void	RBool::set( bool v ){
	if( _value == R_NilValue ){
		PROTECT( _value = allocVector(LGLSXP,1 ) );
	}
	INTEGER(_value)[0] = v?1:0;
        //UNPROTECT(1);
}


/*!
**
**
*/
bool	RBool::getValue( void ){
	int intvalue = INTEGER(_value)[0];
	bool res = ( intvalue != 0 );
	return res;
}


/*!
**
**
*/
bool	RBool::getValue( SEXP v ){
	int intvalue = INTEGER(v)[0];
	bool res = ( intvalue != 0 );
	return res;
}



//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//**		RInteger
//***************************************************************


/*!
**
**
*/
RInteger::RInteger( int value ) : _value(R_NilValue){
	set( value );
}


/*!
**
**
*/
RInteger::~RInteger(){
	UNPROTECT(1);
}


/*!
**
**
*/
void	RInteger::set( int v ){
	if( _value == R_NilValue ){
		PROTECT( _value = allocVector(INTSXP,1 ) );
	}
	INTEGER(_value)[0] = v;
        //UNPROTECT(1);
}


/*!
**
**
*/
SEXP	RInteger::get( void ){
	return _value;
}


/*!
**
**
*/
int		RInteger::getAsInt( void ){
	if( _value == R_NilValue )
		return -1;
	return INTEGER(_value)[0];
}
	



//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************



		
RMatrix::RMatrix() : _matvar(R_NilValue),self_alloced(false)
{
}
		
RMatrix::RMatrix(SEXP mat) : _matvar(R_NilValue),self_alloced(false)
{
	if( isMatrix(mat) )
		set( mat );
	else
		ONDBG df1("RMatrix::RMatrix(SEXP) : param not a matrix!\n");
}

RMatrix::RMatrix(int type, int w, int h ) : _matvar(R_NilValue),self_alloced(false)
{
	RMatrix::alloc( type, w, h );
}

RMatrix::~RMatrix()
{
	dealloc();
}
	
	
//
bool		RMatrix::alloc( int sexptype, int numcols, int numrows )
{
	dealloc();
	PROTECT( _matvar = allocMatrix( sexptype, numcols , numrows ) );
	self_alloced = (_matvar!=R_NilValue);
	return self_alloced;
}

//!
void		RMatrix::dealloc( void )
{
	if( (_matvar !=  R_NilValue) && self_alloced )
	{
		self_alloced = false;
		_matvar = R_NilValue;
		UNPROTECT(1);
	}
}

//!
bool		RMatrix::isValid( void )
{
	if( false == isMatrix(_matvar) )
		return false;
	return (_matvar != R_NilValue );
}

//!
SEXP		RMatrix::get( void )
{
	return _matvar;
}

//!
void		RMatrix::set( SEXP m )
{
	//@@ FIXME is <m> really a matrix ??
	dealloc();
	_matvar=m;
}

//!
unsigned int			RMatrix::numRows( void )
{
	SEXP matdims = getAttrib( _matvar, R_DimSymbol);
	return INTEGER( matdims )[0];
}


//!
unsigned int			RMatrix::numCols( void )
{
	SEXP matdims = getAttrib( _matvar, R_DimSymbol);
	return INTEGER( matdims )[1];
}


//!
int			RMatrix::getType( void )
{
	return TYPEOF(_matvar);
}


//!
int*		RMatrix::getIntPtr( void )
{
	if( _matvar == R_NilValue )
		return 0;
	if( getType() != INTSXP )
		return 0;
	return( INTEGER(_matvar ) );
}


//!
SEXP*		RMatrix::getStrPtr( void )
{
	if( _matvar == R_NilValue )
		return 0;
	if( getType() != STRSXP )
		return 0;
	return 0;//( STRING_ELT(_matvar,0) );
}



//!
double*	RMatrix::getDoublePtr( void )
{
	if( _matvar == R_NilValue )
		return 0;
	if( getType() != REALSXP )
		return 0;
	return( REAL(_matvar ) );
}


//!
void		RMatrix::clearD( double clearvalue )
{
	unsigned int numelems = numRows() * numCols();
	double * ptr = getDoublePtr();
	if( ptr == 0 )
		return;	//@@ TODO warn
	while( numelems-- )
	{
		ptr[0] = clearvalue;
		ptr++;
	}
}



//!
void		RMatrix::clearI( int clearvalue )
{
	unsigned int numelems = numRows() * numCols();
	int * ptr = getIntPtr();
	if( ptr == 0 )
		return;	//@@ TODO warn
	while( numelems-- )
	{
		ptr[0] = clearvalue;
		ptr++;
	}
}



//!
bool		RMatrix::isMatrix( SEXP m )
{
	//@@ FIXME enough to declare this a matrix ? class=="matrix" not necessary ?
	SEXP matdims = getAttrib( m, R_DimSymbol);
	return ( LENGTH(matdims) == 2 );
}



//!
unsigned int	RMatrix::numRows( SEXP m )
{
	SEXP matdims = getAttrib( m, R_DimSymbol);
	return INTEGER( matdims )[0];
}



//!
unsigned int	RMatrix::numCols( SEXP m )
{
	SEXP matdims = getAttrib( m, R_DimSymbol);
	return INTEGER( matdims )[1];
}

		//	get row names

//!	Returns a string vector of row names or NULL
SEXP		RMatrix::getRowNames( SEXP m )
{
	SEXP matdimnames = getAttrib( m, R_DimNamesSymbol);
	if( matdimnames == R_NilValue )
	{
		//TODO error message
		return R_NilValue;
	}
	return VECTOR_ELT( matdimnames, 0 );
}
	

//!	Returns a string vector of column names or NULL
SEXP		RMatrix::getColNames( SEXP m )
{
	SEXP matdimnames = getAttrib( m, R_DimNamesSymbol);
	if( matdimnames == R_NilValue )
	{
		//TODO error message
		return R_NilValue;
	}
	return VECTOR_ELT( matdimnames, 1 );
}

//!
bool		RMatrix::setRowNames( SEXP m , SEXP rownames )
{
	
	//
	if( m == R_NilValue || rownames == R_NilValue )
	{
		//TODO error message
		return false;
	}
	
	//
	if( false == RString::isStr( rownames ) )
	{
		//TODO error message
		return false;
	}
	
	//
	if( numRows(m) != (unsigned)Rf_length(rownames) )
	{
		df1("RMatrix::setColNames : Vector length mismatch: %d matrix rows != %d names!\n",numRows(m) ,Rf_length(rownames) );
		return false;
	}
	
	//
	SEXP matdimnames = getAttrib( m, R_DimNamesSymbol);
	if( matdimnames == R_NilValue )
	{
		PROTECT( matdimnames = allocVector( VECSXP, 2 ) );
		SET_VECTOR_ELT( matdimnames, 0, rownames );
		UNPROTECT( 1 );
		return true;
	}
	
	//
	SET_VECTOR_ELT( matdimnames, 0, rownames );
	return true;
}

//!
bool		RMatrix::setColNames( SEXP m , SEXP colnames )
{
	
	//
	if( m == R_NilValue || colnames == R_NilValue )
	{
		//TODO error message
		return false;
	}
	
	//
	if( false == RString::isStr( colnames ) )
	{
		//TODO error message
		return false;
	}
	
	//
	if( numCols(m) != (unsigned)Rf_length(colnames) )
	{
		df1("RMatrix::setColNames : Vector length mismatch: %d matrix rows != %d names!\n",numCols(m) ,Rf_length(colnames) );
		return false;
	}
	
	//
	//
	SEXP matdimnames = getAttrib( m, R_DimNamesSymbol);
	if( matdimnames == R_NilValue )
	{
		PROTECT( matdimnames = allocVector( VECSXP, 2 ) );
		SET_VECTOR_ELT( matdimnames, 1, colnames );
		UNPROTECT( 1 );
		return true;
	}
	
	//
	SET_VECTOR_ELT( matdimnames, 1, colnames );
	return true;
}

//!
bool		RMatrix::setRowNames( SEXP m , std::vector<std::string>& rownames )
{
	if( m == R_NilValue )
	{
		//TODO error message
		return false;
	}
	if( numRows(m) != rownames.size() )
	{
		df1("RMatrix::setColNames : Vector length mismatch: %d matrix rows != %d names!\n",numCols(m),rownames.size() );
		return false;
	}
	
	//
	//
	SEXP rownames_sexp;
	PROTECT( rownames_sexp = allocVector( STRSXP, rownames.size() ) );
	for( unsigned int i=0; i < rownames.size(); i++ )
	{
		SET_STRING_ELT( rownames_sexp, i, mkChar( rownames[i].c_str() ) );
	}
	
	//
	//
	bool res = setRowNames( m , rownames_sexp );
	UNPROTECT( 1 );
	return res;
}

bool		RMatrix::setColNames( SEXP m , std::vector<std::string>& colnames )
{
	if( m == R_NilValue )
	{
		//TODO error message
		return false;
	}
	if( numCols(m) != colnames.size() )
	{
		df1("RMatrix::setColNames : Vector length mismatch: %d matrix rows != %d names!\n",numCols(m),colnames.size() );
		return false;
	}
	
	//
	//
	SEXP colnames_sexp;
	PROTECT( colnames_sexp = allocVector( STRSXP, colnames.size() ) );
	for( unsigned int i=0; i < colnames.size(); i++ )
	{
		SET_STRING_ELT( colnames_sexp, i, mkChar( colnames[i].c_str() ) );
	}
	
	//
	//
	bool res = setColNames( m , colnames_sexp );
	UNPROTECT( 1 );
	return res;
}


//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************




//!
RList::RList() : _listvar(R_NilValue),self_alloced(false)
{
}


//!
RList::RList( int len ) : _listvar(R_NilValue),self_alloced(false)
{
	if( len > 0 )
	{
		PROTECT( _listvar = allocList( len ) );
		self_alloced = true;
	}
}


//!
RList::~RList()
{
	if( self_alloced )
		UNPROTECT(1);
}



bool	RList::isList( SEXP m )
{
	SEXP cls;
	PROTECT( cls = getAttrib( m , R_ClassSymbol ) );
	ONDBG Rprintf("RList::isList : class is '%s'\n", RString::get(cls,0) );
	int r = TYPEOF( m );
	ONDBG Rprintf("RList::isList : typeof = %d (LIST=%d,INT=%d,STR=%d,VEC=%d,REAL=%d)\n", r,LISTSXP,INTSXP,STRSXP,VECSXP,REALSXP);
	return ::isList( m );
}


int		RList::length( SEXP m )
{
	return ::Rf_length( m );
}



//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************

/*
**
**
*/
SEXP	createList( SEXP r )
{
	SEXP res;
	PROTECT( res = allocList( 5 ) );
        UNPROTECT(1);
	return res;
}

		
/*! get the list element named str, or return NULL
**
*/
SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	for (R_len_t i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
	/*enables us to say

       double g;
       g = REAL(getListElement(a, "g"))[0];
	*/
}

void setListElement(SEXP list, char *str, SEXP value)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      SET_VECTOR_ELT(list, i, value);
      break;
    }
}

/*! 
**
*/
SEXP getvar(SEXP name, SEXP rho)
 {
	 SEXP ans;
 
	 if(!isString(name) || length(name) != 1)
		 error("name is not a single string");
	 if(!isEnvironment(rho))
		 error("rho should be an environment");
	 ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);
	 Rprintf("first value is %f\n", REAL(ans)[0]);
	 return(R_NilValue);
 }
  




