/*
**
**		dynstorage.h
**
**
**			class to store strings of not a priori known lengths and number with as little memory overhead as possible
**
**
**
*/


#ifndef _DYNSTORAGE_H_INCLUDED_
#define _DYNSTORAGE_H_INCLUDED_

/*!
**
**
*/
class dynstorage			//
{
public:

	//
	dynstorage(unsigned int inblocksize ){
		blocksize=inblocksize;
		currentblock = 0;
		currentblockpointer=0;
		currentblockindex=0;
		top=0;
		enlarge();
	}

	//! ONLY CALL IF CURRENT BLOCK IS EXHAUSTED - OTHERWISE REMAINING BYTES ARE INVALID DATA AND WILL BE CONSIDERED VALID!
	void	enlarge( void )
	{
		chunkheader * newptr;
		if( currentblock == 0 || currentblock->next == 0 )
		{

			newptr = (chunkheader*)malloc( blocksize+sizeof(chunkheader)-1 );
											//-1 because chunkheader contains 1 databyte
			if( newptr == 0 )
				throw "dynstorage : failed malloc to enlarge !\n";
			if( currentblock )
				currentblock->next = newptr;
			newptr->next = 0;
		}
		else
			newptr = currentblock->next;

		//switch to next block
		//
		currentblock = newptr;
		currentblockpointer = &currentblock->databegin[0];
		currentblockindex = 0;

		//
		//
		if( top == 0 )
			top = newptr;

		//
	}//enlarge

	//
	struct chunkheader {
		chunkheader		*next;
		char			databegin[1];
	};

	//! Returns whether the given string at <seqname> is stored here
	int		hasname( const char * seqname )
	{
		int		numm=0;
		char	cc=1;
		int		i=0;
		int		si=0;
		bool	res=false;
		int		ni=0;

		while( get(i,cc) )
		{
			//
			if( cc )
			{
//				Rprintf("[%c]",cc);
				if( seqname[si++] == cc )
				{
					numm++;
				}
			}
			else
			{
//				Rprintf("(%d <=> %d)\n",si,numm);
				if( seqname[si] == 0 && numm==si )
					return ni;
				si=0;
				numm=0;
				ni++;
			}

			i++;
		}

		return -1;
	}


	//! Get character number <index> in the datablock
	inline	bool	get( int index, char&c )
	{
		int numchunk = index/blocksize;
		int	numbyte = index%blocksize;
		chunkheader * p = top;
		//printf("get[%d:chunk %d/byte %d]",index,numchunk,numbyte);
		while( p && numchunk-- )
		{
			p = p->next;
		}

		//
		if( p == 0 || (p==currentblock && numbyte >= currentblockindex) )
		{
			c=0;
			return false;
		}

		//
		c = p->databegin[numbyte];
		//printf("-->%c (%d)\n",c,c);

		//
		return true;

		//
	}

	//! Store <numbytes> characters stored at <p>
	inline void		serialise( char * p, int numbytes=1 )
	{

		//
		//
		while( numbytes > 0 )
		{

			// determine the maximum number of bytes we can copy
			int numbytescopyable = numbytes>(blocksize-currentblockindex) ? (blocksize-currentblockindex) : numbytes;

			//
			memcpy(currentblockpointer,p, numbytescopyable);

			//
			currentblockpointer += numbytescopyable;
			currentblockindex += numbytescopyable;
			p += numbytescopyable;
			numbytes -= numbytescopyable;

			if( currentblockindex >= blocksize )
				enlarge();

		}


		//
	}//serialise( buffer )

	//
	inline	void	reset(void)
	{
		chunkheader * p = top;
		while( p )
		{
			memset( &p->databegin[0], 0, blocksize );
			p = p->next;
		}
		currentblock=top;
		currentblockpointer=&top->databegin[0];
		currentblockindex=0;
	}


	//!
	inline	void	serialise( char c )
	{
//		printf("ser[%02d:%d,%d] -> '%s'\n",c,currentblockindex,blocksize,&top->databegin[0]);
		if( currentblockindex >= blocksize )
			enlarge();
		*currentblockpointer = c;
		currentblockpointer++;
		*currentblockpointer=0;
		currentblockindex++;
	}//serialise( char )

	//
	chunkheader			*top;
	chunkheader			*currentblock;
	char				*currentblockpointer;
	unsigned int		currentblockindex;

	unsigned int		blocksize;
};



#endif

