/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		Whopgen - Reading data per line
**
**
**

	TODO:


	METHODS:
		SEXP	VCF_getNextSNP( SEXP vcfptr )
		SEXP	VCF_getChrom( SEXP vcfptr )
		SEXP	VCF_getPos( SEXP vcfptr )
		SEXP	VCF_getID( SEXP vcfptr )
		SEXP	VCF_getRef( SEXP vcfptr )
		SEXP	VCF_getAlt( SEXP vcfptr )
		SEXP	VCF_getQual( SEXP vcfptr )
		SEXP	VCF_getFilter( SEXP vcfptr )
		SEXP	VCF_getInfo( SEXP vcfptr )
		SEXP	VCF_getInfoField( SEXP vcfptr, SEXP fieldnam )
		SEXP	VCF_getFormat( SEXP vcfptr )
		SEXP	VCF_getSample( SEXP vcfptr, SEXP stridx )
		SEXP	VCF_readLineRaw( SEXP vcfptr )
		SEXP	VCF_readLineTSV( SEXP vcfptr )
		SEXP	VCF_readLineDF( SEXP vcfptr )
		inline bool isBiallelic( const char * _REF, const char* _ALT )
		SEXP VCF_getBial( SEXP vcfptr, SEXP mat )
	-	----------------------------------------------------------
		SEXP	VCF_countSNPs( SEXP vcfptr )
		SEXP	VCF_countBiallelicSNPs( SEXP vcfptr )
	-	----------------------------------------------------------
	
		TODO : rename for more consistency
	-	----------------------------------------------------------
		SEXP	VCF_readLineRaw( SEXP vcfptr )
		SEXP	VCF_readLineTSV( SEXP vcfptr )
		SEXP	VCF_readLineDF( SEXP vcfptr )
		SEXP	VCF_readLineBial( SEXP vcfptr, SEXP mat )
		SEXP	VCF_readLineNextSNP( SEXP vcfptr )

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

//#define		ALLOW_NONSTANDARD_VCF_GT


#define		__INTERNAL_RETURN_FIELD_OF_LINE( fld )\
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );\
	if( 0 == f )	{		return R_NilValue;		}\
	RString resvec;\
	char	buffer[512];\
	if( f->copyField( fld, &buffer[0], sizeof(buffer)-2 ) )\
	{\
		resvec.alloc( 1 );\
		resvec.set( &buffer[0] );\
		return resvec.get();\
	}\
	return R_NilValue;\




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
**
**
*/
EXPORT	SEXP	VCF_countBiallelicSNPs( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	//
	int numsnps=0;
	while( f->parseNextLine() )
	{

		//if( _internal_isBiallelicSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ) )
		
		//
		const char * rp = f->getFieldPtr( REF );
		const char * ap = f->getFieldPtr( ALT );
		if( rp ==0 || ap == 0 )
			continue;
		else if( rp[1] == '\t' && ap[1] == '\t' )
			numsnps++;
	}

	RInteger res( numsnps );
	return res.get();

}//..VCF_countBiallelicSNPs




/*!	
**
**
**
*/
EXPORT	SEXP	VCF_countSNPs( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	//
	int numsnps=0;
	while( f->parseNextLine() )
	{
		if( _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ) )
			numsnps++;
	}

	RInteger res( numsnps );
	return res.get();

}//..VCF_countSNPs




/*!	
**
**
**
*/
EXPORT	SEXP	VCF_parseNextSNP( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	//
	while( f->parseNextLine() )
	{
		if( _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ) )
			return RBool::True();
	}

	return RBool::False();
}




/*!	
**
**
**
*/
EXPORT	SEXP	VCF_parseNextLine( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}
	
	//
	//
	if( f->parseNextLine() )
		return RBool::True();

	return RBool::False();
}






/*!	Return the contents of the CHROM field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getChrom( SEXP vcfptr )
{
	__INTERNAL_RETURN_FIELD_OF_LINE( CHROM );
}




/*!	Return the contents of the POS field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getPos( SEXP vcfptr )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );\
	if( 0 == f )	{		return R_NilValue;		}\
	const char * posptr = f->getFieldPtr( POS );
	if( posptr == 0 )
		return R_NilValue;
	int posasint = atoi( posptr );
	RInteger res( posasint );
	return res.get();
}




/*!	Return the contents of the ID field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getID( SEXP vcfptr )
{
	__INTERNAL_RETURN_FIELD_OF_LINE( ID );
}




/*!	Return the contents of the REF field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getRef( SEXP vcfptr )
{
	__INTERNAL_RETURN_FIELD_OF_LINE( REF );
}




/*!	Return the contents of the ALT field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getAlt( SEXP vcfptr )
{
	__INTERNAL_RETURN_FIELD_OF_LINE( ALT );
}




/*!	Return the contents of the QUAL field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getQual( SEXP vcfptr )
{
	//TODO : make sure this field is set ! (some fields are optional in VCF)
	__INTERNAL_RETURN_FIELD_OF_LINE( QUAL );
}




/*!	Return the contents of the FILTER field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getFilter( SEXP vcfptr )
{
	//TODO : make sure this field is set ! (some fields are optional in VCF)
	__INTERNAL_RETURN_FIELD_OF_LINE( FILTER );
}




/*!	Return the contents of the INFO field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getInfo( SEXP vcfptr )
{
	//TODO : make sure this field is set ! (some fields are optional in VCF)
	__INTERNAL_RETURN_FIELD_OF_LINE( INFO );
}




/*!
**
**
**
*/
EXPORT	SEXP	VCF_getInfoField( SEXP vcfptr, SEXP fieldnam )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )	{		return R_NilValue;		}
	
	//
	const char* infofieldptr = f->getFieldPtr( INFO );
	const char* fieldnamstr = RString::get( fieldnam );
	
	//
	if( infofieldptr && fieldnamstr )
	while( true )
	{
		
		if( *infofieldptr == 0 || *infofieldptr == '\t' )
			break;
		
		//
		//
		int j=0;
		for( ; infofieldptr[j] != 0 && infofieldptr[j] != '\t'; j++ )
		{
		//	Rprintf("j=%d : nam %c <=> %c info\n", j, fieldnamstr[j] , infofieldptr[j] );
			if( fieldnamstr[j] != infofieldptr[j] )
				break;
		}
		infofieldptr+=j;

		//	success if current char is =  OR ; as well as fieldnamstr is 0
		//
		if( fieldnamstr[j] == 0 && ('='  == *infofieldptr || ';' == *infofieldptr) )
		{
			infofieldptr++;
		
			//	determine length of result string
			//
			int infolength=0;
			for( ; infofieldptr[infolength] != ';' && infofieldptr[infolength]!=0 && infofieldptr[infolength] != '\t'; infolength++ )
				//Rprintf("%c",infofieldptr[infolength])
				;
			
			//	create an R string variable and return it
			//
			{
				char*	buf = new char[infolength+2];
				int i;
				for( i=0; i < infolength; i++ )
					buf[i] = infofieldptr[i];
				buf[i]=0;
				
				//
				RString res;
				res.alloc(1);
				res.set( buf );
				
				//
				delete [] buf;
				
				//
				return res.get();
			}
			
			//
		}
		else
		{
			//	skip until the next subfield of the INFO field
			//
			for( infofieldptr+=j; *infofieldptr != ';'; infofieldptr++ )
			{
				if( *infofieldptr == 0 || *infofieldptr == '\t' )
					return R_NilValue;
			}
			infofieldptr++;
		}

	}//..while forever

	return R_NilValue;
}




/*!	Return the contents of the FORMAT field of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getFormat( SEXP vcfptr )
{
	//TODO : make sure this field is set ! (some fields are optional in VCF)
	__INTERNAL_RETURN_FIELD_OF_LINE( FORMAT );
}




/*!	Return the contents of one of the sample fields of the last getNextSNP'd/getNextLine'd/getNextInDel'd line
**
**
**
*/
EXPORT	SEXP	VCF_getSample( SEXP vcfptr, SEXP stridx )
{	
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )	{		return R_NilValue;		}
	
	//	make sure sample index is an integer
	//
	if( RNumeric::isInt( stridx ) == false )
		return R_NilValue;
	
	//	get as int and make sure its inside the valid range
	//
	int stridxint = RNumeric::getInt( stridx );
	if( stridxint < 0 )
		return R_NilValue;
	stridxint += f->getFirstSampleFieldIndex();
	
	//
	if( f->getNumFields() < (unsigned)stridxint )
		return R_NilValue;
	
	//	get a copy of the field, turn it into an R string and return it
	//
	RString resvec;
	char	buffer[512];
	if( f->copyField( stridxint, &buffer[0], sizeof(buffer)-2 ) )
	{
		resvec.alloc( 1 );
		resvec.set( &buffer[0] );
		return resvec.get();
	}
	
	//
	return R_NilValue;
}





/*!	
**
**
**
*/
EXPORT SEXP VCF_readLineRaw( SEXP vcfptr, SEXP str )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}

	const char * s = f->readNextLine();
	if( s != 0 )
	{
		SET_STRING_ELT( str, 0, mkChar(s) );
		return RBool::True();
	}
	
	//
	return RBool::False();
}

/*!	
**
**
**
*/
EXPORT SEXP VCF_readLineRawFiltered( SEXP vcfptr, SEXP str )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return RBool::False();
	}

	//
	//
	RString resvec;
	int runs=1;
	while( f->parseNextLine() )
	{
		//
		if( runs-- <= 0 )
			break;
		
		//
		if( false == filterLine( f ) )
		{
			Rprintf("Line has been filtered away\n");
			continue;
		}
		
		//
		//
		const char * s = f->getCurrentLine();
		if( s != 0 )
		{
			SET_STRING_ELT( str, 0, mkChar(s) );
			return RBool::True();
		}
		else
			return RBool::False();
			
	
	}//...while

	
	//
	return RBool::False();
}












		//-----------------------------------------------------------------------------		TSV
		
		
		
		
		
		
		
		
		
		
		

/*!
**
**
**
**
*/
EXPORT SEXP VCF_readLineTSV( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_readLineTSV : need VCFhandle as parameter!\n");
		return R_NilValue;
	}

	//
	//
	const char * s = f->readNextLine();
	if( s != 0 )
	{
		//
		RString resvec;
		TSVParser pars(s);
		
		//
		int numfields = pars.numFields();
		
		//
		if( resvec.alloc(numfields) )
		{
			char	tokenbuffer[512];
			int idx=0;
			for( unsigned int i=0; i < pars.numFields() ; i++ )
			{
				if( pars.getField(tokenbuffer,500,i) )
				{
					resvec.set( tokenbuffer, idx );
					idx++;
				}
			}
			
			//
			//
			return resvec.get();
		}
		else
			df0("Could not allocate a R string vector with %d elements!\n",numfields);

	}

	//
	return R_NilValue;
    //
}




/*!
**
**
**
**
*/
EXPORT SEXP VCF_readLineTSVFiltered( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		df0("VCF_readLineTSV : need VCFhandle as parameter!\n");
		return R_NilValue;
	}

	//
	//
	RString resvec;
	int runs=1;
	while( f->parseNextLine() )
	{
		//
		if( runs-- <= 0 )
			break;
		
		//
		if( false == filterLine( f ) )
		{
			//Rprintf("Line has been filtered away\n");
			continue;
		}
		
		unsigned int numfields = f->numParsedFields();
		//Rprintf("%d parsed fields\n", numfields);
		
		//
		if( resvec.alloc(numfields) == false )
		{
			df0("Could not allocate a R string vector with %d elements!\n",numfields);
			return R_NilValue;
		}

		char	tokenbuffer[512];
		int		idx=0;
		for( unsigned int i=0; i < numfields ; i++ )
		{
			if( f->copyField(i,tokenbuffer,sizeof(tokenbuffer)-2) )
			{
				resvec.set( tokenbuffer, idx );
				idx++;
			}
		}
		return resvec.get();
	
	}//...while

	//
	return R_NilValue;
    //
}



#if 0

	//	Due to the complexities of constructing data.frames in C, this code is
	//		excluded from the build.
	//	For speed reasons it is recommended to use the combination of 1x VCF_getFieldNames() to
	//		get the names for each data field and <n> times VCF_readLineTSV to get a vector containing a line's data on all fields or
	//		or the combination of VCF_parseNextLine()/SNP() and the VCF_getXXX() functions
	//


		//-----------------------------------------------------------------------------		Data.Frame


/*!	Reads a line 
**
**
**
**
*/
EXPORT SEXP VCF_readLineDF( SEXP vcfptr )
{
	Rprintf("TODO : due to complexities of data.frames, this function is not yet properly implemented!\n");
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return R_NilValue;
	}

	//
	unsigned int numfields = f->getNumFields();

	//
	//
	if( f->parseNextLine() )
	{
		char buf[1024];

		if( numfields == f->numParsedFields() )
		{
			//
			for( unsigned int i=0; i < numfields; i++ )
			{
				const char * fieldnam = f->getFieldName(i);
				if( f->copyField(i,buf,sizeof(buf)-2) )
					df0("'%s'='%s'\n", fieldnam,buf);
				else
					df0("Problem getting field %d (%s!\n",i,fieldnam);
			}
			
			//
		}
		else
			df0("Parsed %d fields in line but expected %d!\n",f->numParsedFields(),numfields);
	}
	else
		df0("Could not parse line!\n");
	
	return R_NilValue;
    //
}











/*!	Reads a line 
**
**
**
**
*/
EXPORT SEXP VCF_readLineDFFiltered( SEXP vcfptr )
{
	Rprintf("TODO : due to complexities of data.frames, this function is not yet properly implemented!\n");
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		return R_NilValue;
	}

	//
	//
	RString resvec;
	int runs=1;
	while( f->parseNextLine() )
	{
		//
		if( runs-- <= 0 )
			break;
		
		//
		if( false == filterLine( f ) )
		{
			//Rprintf("Line has been filtered away\n");
			continue;
		}
		
		//	construct data.frame ...
		//
		unsigned int numfields = f->numParsedFields();
		char buf[1024];

		//
		for( unsigned int i=0; i < numfields; i++ )
		{
			const char * fieldnam = f->getFieldName(i);
			if( f->copyField(i,buf,sizeof(buf)-2) )
				df0("'%s'='%s'\n", fieldnam,buf);
			else
				df0("Problem getting field %d (%s!\n",i,fieldnam);
		}

		//
		return R_NilValue;
	
	}//...while

	return R_NilValue;
    //
}





#endif
