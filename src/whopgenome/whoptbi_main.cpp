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




/*!
**
*/
whop_tabix::whop_tabix(){

	//
	tabix = 0;
	index = 0;
	iter = 0;

	//
	num_seqnames = 0;
	sequence_names = 0;

	//
	last_num_fields = 0;
	field_offsets_size = 0;
	field_offsets = 0;
	current_line = 0;

	//
	//currentTid = ;
	currentBegin = 0;
	currentEnd = 0;

	//
}



/*!
**
*/
whop_tabix::whop_tabix(const char * filename ){

	//
	tabix = 0;
	index = 0;
	iter = 0;

	//
	num_seqnames = 0;
	sequence_names = 0;

	//
	last_num_fields = 0;
	field_offsets = 0;
	field_offsets_size = 0;
	current_line = 0;

	//
	//currentTid = 0;
	currentBegin = 0;
	currentEnd = 0;

	//
	this->open( filename );

	//
}



/*!
**
*/
whop_tabix::~whop_tabix()
{
	ONDBG Rprintf("whop_tabix DTOR: %lx\n",(uint64_t)tabix);

	//
	if( sequence_names != 0 )
	{
		free( sequence_names );	//were allocated with calloc() in tabix
		sequence_names=0;
		num_seqnames = 0;
	}

	//
	if( tabix != 0 )
	{
		ti_close( tabix );
		tabix = 0;
	}

	//
	if( iter != 0 )
	{
		ti_iter_destroy(iter);
		iter = 0;
	}

	//
	if( field_offsets != 0 )
	{
		last_num_fields = 0;
		field_offsets_size = 0;
		free( field_offsets );
		field_offsets = 0;
	}

	//
}



/*!
**
*/
bool			whop_tabix::open( const char * filename )
{
	//
	ti_iter_t hdriter=0;
	
	ONDBG Rprintf("whoptabix_open\n");

	//
	//
	try
	{
		//	load tabix-index of file
		//
		ONDBG Rprintf("ti_open(%s,0)...\n",filename,0);
		tabix = ti_open( filename , 0 );
		ONDBG Rprintf("ti_open=%d,OK\n",tabix);
		if( tabix == 0 )
			throw "whop_tabix::open : Failed to open tabix index file";

		//		load index
		//
		ONDBG Rprintf("ti_index_load(%s)",filename);
		index = ti_index_load( filename );
		ONDBG Rprintf("OK\n");
		if( index == 0 )
			throw "whop_tabix::open : index load failed!";

		//
		//
		ONDBG Rprintf("ti_query(%x,0,1,99999999)",tabix);
		iter = ti_query( tabix, 0, 1, 99999999 );//FIXME best way to have a catch-all iterator?
		ONDBG Rprintf("OK\n");
		if( 0 == iter )
			throw "whop_tabix::open : Failed to create whole-file iterator!";
		bEOR = (iter != 0 );

		//	iterator that includes header lines
		//
		ONDBG Rprintf("ti_query(%p,0,0,0)",tabix);
		hdriter = ti_query(tabix,0,0,0);
		ONDBG Rprintf("OK\n");
		if( 0 == hdriter )
			throw "whop_tabix::open : Failed to create header-iterator!";

		//	read header lines until the ##CHROM
		//
		int len=0;
		const char	*s;
		ONDBG Rprintf("ti_read(%p,%p,%p)*x",tabix,hdriter,&len);
		while(  (s=ti_read( tabix, hdriter, &len ))!=0 && (s[0]=='#') )
		{
			//Rprintf("[%s]\n",s);
			header_lines.push_back( s );
			
			if( s[1] == 'C' && s[2] == 'H' )	//stop parsing lines once the '#CHROM	POS	ID	REF	ALT	QUAL...' line is found
				break;							//		..helps not skipping the first data-line if no region was set
		}
		ONDBG Rprintf("OK\n");

		//	get sequence names (first column, usually chromosome name)
		//
		num_seqnames=0;
		sequence_names = ti_seqname( index, (signed*)&num_seqnames);

		ONDBG Rprintf("return TRUE from whoptabix_open\n");
#if 0
		//
		//
		for( unsigned int i = 0; i < num_seqnames; i++ )
		{
			printf("Range of nucleotides for '%s':\n",sequence_names[i] );
			
			//
			int min = 0;
			int max = 100*1000*1000, okmax=0, failmax=cmax;
			const char * tid = &sequence_names[i];
			ti_iter_t newiter = 0;
			int testiters=20;
			
			//	find maximum
			//
			while( testiters-- > 0 )
			{
				newiter = ti_query( tabix , tid, max, max+1 );
				
				// memorize which value works as maximum
				//
				if( newiter == 0 )	//half of max
				{
					failmax = max;
				}
				else
					okmax = max;
				
				// if there is too much difference between known-to-fail and known-to-work maximum
				//
				if( (failmax - okmax) > 1 )
				{
					if( newiter == 0 )
						
				}
				
				
			}
		}//...for( all contig-ids )
#endif

		return true;

		//
	}
	catch( const char * errmsg )
	{
		Rprintf("Caught exception inside whop_tabix::open('%s'):\n\t'%s'\n",filename,errmsg);
	}
	catch( ... )
	{
		Rprintf("Caught generic exception inside whop_tabix::open('%s')\n",filename);
	}

	//
	//
	if( sequence_names != 0 )
	{
		free( sequence_names );	//were allocated with calloc() in tabix
		sequence_names=0;
		num_seqnames = 0;
	}
	if( tabix != 0 )
		ti_close( tabix );
		tabix=0;
	if( hdriter )
		ti_iter_destroy( hdriter );
		hdriter=0;
	if( iter )
		ti_iter_destroy( iter );
		iter=0;
	//FIXME : free index ?
	
	Rprintf("return FALSE from whoptabix_open\n");

	return false;
}




