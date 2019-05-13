/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		R support module - data.frame
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


/*
**
**
*/
SEXP	createDataFrame( unsigned int numrows, unsigned int numcolumns, const char** rownames, const char**columnnames )
{
//	SEXP dataFrameSym = ::Rf_install( "data.frame"); // cannot be gc()ed once in symbol table
//    SEXP emptydataframe = ::Rf_eval( ::Rf_lang1( dataFrameSym ), R_GlobalEnv );
	int unprotcount = 0;

    //
    SEXP newlist = PROTECT( allocVector( VECSXP, numcolumns ) );
		unprotcount++;
    
    //	allocate vectors for the columns
    //
    for( unsigned int i=0; i < numcolumns; i++ )
    {
		SEXP newlist1 = PROTECT( allocVector( STRSXP, numrows ) );
		unprotcount++;
		//SEXP newlist2 = PROTECT( allocVector( STRSXP, 2 ) );
		SET_VECTOR_ELT( newlist, i, newlist1 ); 
	}
    
    //SET_VECTOR_ELT( newlist, 1, newlist2 ); 
    
    //	set column names
    //
    SEXP newnames = PROTECT( allocVector( STRSXP, numcolumns ) );
		unprotcount++;
    if( columnnames )
    {
		for( unsigned int i=0; i < numrows; i++ )
		{
			SET_STRING_ELT( newnames, i, Rf_mkChar( columnnames[i] ) ) ;
		}
	}
	else
	{
		char buf[64];
		for( unsigned int i=0; i < numrows; i++ )
		{
			sprintf(buf,"%d",i);
			SET_STRING_ELT( newnames, 0, Rf_mkChar( buf ) ) ;
		}
	}
    Rf_setAttrib( newlist , R_NamesSymbol, newnames );
    
    //	set row names
    //
	SEXP newrownames = PROTECT( allocVector( STRSXP, numrows ) );
		unprotcount++;
    if( rownames )
    {
		for( unsigned int i=0; i < numrows; i++ )
		{
			SET_STRING_ELT( newrownames, i, Rf_mkChar( rownames[i] ) ) ;
		}
		//SET_STRING_ELT( newrownames, 1, Rf_mkChar("2") ) ;
	}
	else
	{
		char buf[64];
		for( unsigned int i=0; i < numrows; i++ )
		{
			sprintf(buf,"%d",i);
			SET_STRING_ELT( newrownames, 0, Rf_mkChar( buf ) ) ;
		}
	}
    Rf_setAttrib( newlist , R_RowNamesSymbol, newrownames );

    //
    //
    Rf_setAttrib( newlist , R_ClassSymbol, Rf_mkString("data.frame") );
    
    //
    UNPROTECT(unprotcount);
    //UNPROTECT(3);
    //
    return newlist;
}

