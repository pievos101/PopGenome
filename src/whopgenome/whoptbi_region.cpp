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


/*!
**
*/
bool			whop_tabix::setRegion( const char* tid, int begin, int end )
{
	//
	//
	// quick fix for CRAN #FIXME (commented out) start
	//
	/*
        if( this == 0 )
	{
		Rprintf("(!!) whop_tabix::setregion called on this==0!\n");
		return 0;
	}
	*/
	// end fixme 

	//
	//
	if( tabix == 0 )
	{
		Rprintf("(!!) whop_tabix::setregion called on this==0!\n");
		return 0;
	}

	//
	//
	ti_iter_t newiter = ti_query( tabix , tid, begin, end );
	bEOR = (newiter != 0 );
	
	
	//
	if( newiter != 0 )
	{
		ti_iter_destroy( iter );
		iter = newiter;
		
		//
		currentTid = tid;
		currentBegin = begin;
		currentEnd = end;

		//
		return true;
	}
	
	//
	Rprintf("whop_tabix::setRegion : '%s' %d - %d NOT SET! (tabix=%x)\n",tid,begin,end,tabix);
	return false;
}


/*!	
**
*/
bool			whop_tabix::restartRegion( void )	//rest region to initial settings
{
	if( tabix )
	{
		iter = ti_query( tabix , currentTid.c_str(), currentBegin, currentEnd );
		bEOR = (iter != 0 );
		return( iter != 0 );
	}
	return false;
}



/*!
**
*/
const char*	whop_tabix::getRegionTid( void )
{
	return currentTid.c_str();
}



/*!
**
*/
int				whop_tabix::getRegionBegin( void )
{
	return currentBegin;
}



/*!
**
*/
int				whop_tabix::getRegionEnd( void )
{
	return currentEnd;
}

