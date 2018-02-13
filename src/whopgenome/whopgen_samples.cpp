/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		main module - dll exports
**
**
**

	TODO:

- get list of samples
	getSampleNames( ..filename.. )
- 



**
*/

//*
//*			INCLUDES
//*

#include	"whopgen_common.h"


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
**
**
*/
EXPORT  SEXP VCF_selectSamples( SEXP vcfptr, SEXP samplesvec )
{
	vcff * f=0;

	//
	guard;

	//
	f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_selectSamples : parameter 1 needs to be valid VCFhandle!\n");
		return RBool::False();
	}
	
	//
	if( RString::isStr( samplesvec ) == false || Rf_length(samplesvec) < 1 )
	{
		Rprintf("VCF_selectSamples : parameter 2 needs to be a vector of strings!\n");
		return RBool::False();
	}

	//	Match the vector of sample names against the samples present
	//		in the VCF file
	//	and create a list of TSV-field-indices to quickly look up the
	//		wanted samples in the TSV-lines
	//
	{
		//unsigned int numsamples = f->getNumSamples();
		unsigned int samplesbeginindex = f->getFirstSampleFieldIndex();
		unsigned int numfields = f->getNumFields();
		unsigned int veclen = RString::length(samplesvec);
		
		//
		//
		unsigned int curstridx = 0;
		f->resetSampleSelection();
		for( unsigned int i=0; i < veclen; i++ )
		{
			const char* vecentry = RString::get( samplesvec, i );
		
			//
			//
			for( unsigned int j = samplesbeginindex; j < numfields; j++ )
			{
				const char * fieldname = f->getFieldName(j);
				if( strcmp( fieldname , vecentry ) == 0 )
				{
					ONDBG df1("Sample %d matches '%s' from input(%s)!\n",j,fieldname,vecentry);
					if( false == f->selectSample( j ) )
					{
						Rprintf("FAILED : %d wanted samples, adding field-index %d/%d\n",f->num_wanted_samples,j,f->getNumFields() );
						//throw "TOO MANY SAMPLES TO SELECT!";
						i=veclen;//break outer for-loop too
					}
					curstridx++;
					break;
				}

			}//..for all VCF sample names
			
			//	make sure we don't select more samples than are stored in the file
			//
			if( curstridx > f->getNumSamples() )
			{
				ONDBG Rprintf("at %d/%d elems, broke setting selected samples (allowed=%d,absolute max=%d)!\n" ,i,veclen,f->getNumSamples(),f->getNumFields());
				break;
			}

		}//..for all elements of the sample-names-input-vector
		
		//
		f->num_wanted_samples=curstridx;
		f->wanted_samples[ f->getNumSamples() ] = ((unsigned)-1);
		f->wanted_samples[ f->getNumFields()-1 ] = ((unsigned)-1);
		ONDBG df1("%d , %d!\n",f->num_wanted_samples,curstridx);
		
		//
	}

	//
	return RBool::True();
	
	//
	unguards;
	if( f )
		f->resetSampleSelection();
	return RBool::False();
	
}//..VCF_selectSamples





/*!
**
**
**
*/
EXPORT SEXP	VCF_getSelectedSamples( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getSelectedSamples : parameter 1 needs to be valid VCFhandle!\n");
		return RBool::False();
	}
	
	if( f->num_wanted_samples < 1 )
	{
		ONDBG Rprintf("VCF_getSelectedSamples : no samples selected!\n");
		return R_NilValue;
	}

	RString res;
	res.alloc( f->num_wanted_samples );
	
	//
	unsigned int strnbegidx = f->getFirstSampleFieldIndex();
	for( unsigned int k=0; k < f->num_wanted_samples; k++ )
	{
		const char * samplefieldname = f->getFieldName( k + strnbegidx );
		res.set( samplefieldname, k );
	}

	return res.get();
}//..VCF_getSelectedSamples





/*!
**
**
**
*/
EXPORT  SEXP VCF_getSampleNames( SEXP vcfptr )
{
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("VCF_getSampleNames : parameter 1 needs to be valid VCFhandle!\n");
		return RBool::False();
	}
	
	//
	//
	RString st;

	//
	int		samplenamebeginindex=7,
			numfields = f->getNumFields();
	for( int k=samplenamebeginindex; k < numfields; k++ )
	{
#define		CMP_FIELD_NAME( n )			(strcasecmp( fieldnam , #n ) == 0)
		const char * fieldnam = f->getFieldName( k );
		if( (k==0) && CMP_FIELD_NAME( CHROM ) ){	samplenamebeginindex++; continue;}
		if( (k==1) && CMP_FIELD_NAME( POS ) ){	samplenamebeginindex++; continue;}
		if( (k==2) && CMP_FIELD_NAME( ID ) ){	samplenamebeginindex++; continue;}
		if( (k==3) && CMP_FIELD_NAME( REF ) ){	samplenamebeginindex++; continue;}
		if( (k==4) && CMP_FIELD_NAME( ALT ) ){	samplenamebeginindex++; continue;}
		if( (k==5) && CMP_FIELD_NAME( QUAL ) ){	samplenamebeginindex++; continue;}
		if( (k==6) && CMP_FIELD_NAME( FILTER ) ){	samplenamebeginindex++; continue;}
		if( (k==7) && CMP_FIELD_NAME( INFO ) ){	samplenamebeginindex++; continue;}
		if( (k==8) && CMP_FIELD_NAME( FORMAT ) ){	samplenamebeginindex++; continue;}
	}
	
	//
	ONDBG df1("sample names begin at index %d\n",samplenamebeginindex );
	ONDBG df1("%d sample names\n",numfields - samplenamebeginindex );
	
	//
	unsigned int unprot=0;
	SEXP _value;
	PROTECT(_value = allocVector(STRSXP, numfields - samplenamebeginindex));
	unprot++;
	for( int k=samplenamebeginindex; k < numfields; k++ )
	{
		const char * str = f->getFieldName(k);
		if(0==str){
			Rprintf("(!!) Unexpected error condition: getFieldName(k=%d) == 0x%08x!\n",k,str);
			break;
		}
		SET_STRING_ELT( _value, k - samplenamebeginindex, mkChar( str ) );
	}
	
	//
	UNPROTECT( unprot );
	return _value;
	//
}//..VCF_getSampleNames


