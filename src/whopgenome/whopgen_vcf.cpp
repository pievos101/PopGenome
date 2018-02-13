/*
**
**		WhopGen
**
**		WHOle genome population genetics with PopGENome
**
**
**		VCF module
**
**
**

	TODO:

- parse header for ##fileformat=VCFv4.1 or ##fileformat=VCFv4.0
- extract sample names
- match a SEXP-STRSXP-vector of sample names with existing sample names
	- print errors if names are not found or double
	- order: like in given list or ordered like in the vcf file ?
	- produce a vector with table fields to extract (by sample name)
- 




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

struct tbi_file {
	tabix_t			*t;
	ti_conf_t		*idxconf;
	ti_iter_t		hdriter;
	ti_iter_t		it;
};

//*
//*			CLASSES
//*



//*
//*			DATA
//*

	//
tabix_t				*t=0;
struct timeval		tv1,
					tv2;

	//
const ti_conf_t		*idxconf;

	//
ti_iter_t			hdriter;		//iterator to get the header lines
ti_iter_t			it;				//iterator to parse the chromosome-regions of interest

	//
vector<string>		indnames;



//*
//*			EXTERNS
//*



//*
//*			CODE
//*




/*!
**
*/
vcff::vcff()
{
	ONDBG Rprintf("VCFF CTOR()@%08x\n",this);
	
	//
	num_fields=0;
	current_line=0;
	sample_begin_index=8;
	wanted_samples=0;
	num_wanted_samples=0;
	
	//TODO clear filters
	num_rules_used=0;
	num_fieldnames_used=0;
	
	//
}




/*!
**
*/
vcff::vcff(const char * filename )
{
	ONDBG Rprintf("VCFF CTOR('%s')@%08x\n",filename,this);
	
	//
	num_fields=0;
	current_line=0;
	sample_begin_index=8;
	wanted_samples=0;
	num_wanted_samples=0;
	
	//
	num_rules_used=0;
	num_fieldnames_used=0;
	
	//
	open( filename );
}




/*!
**
*/
vcff::~vcff()
{
	ONDBG Rprintf("VCFF DTOR (%x)[ws=%x,fo=%x]\n",this,wanted_samples,field_offsets);
	num_fields=0;
	current_line=0;
	num_wanted_samples=0;
	if( wanted_samples )
	{
		free( wanted_samples );
		wanted_samples = 0;
	}
	if( field_offsets )
	{
		free( field_offsets );
		field_offsets = 0;
	}
}




/*
**
*/
bool			vcff::testfunc( void )
{
	//
	int numnames=getNumSequenceNames();
	Rprintf("got %d seqnames!\n",numnames);
	for( int i=0; i < numnames; i++ )
	{
		Rprintf("SeqName#%d='%s'\n",i,getSequenceName(i));
	}
	return true;
}



/*!	Opens the file 'filename' as VCF file.
**		Returns a vcff object if file exists and has the VCF fileformat header line
**		Returns NULL if file could not be found, is not tabix-indexed or does not have the ##fileformat=VCF line
**
**	Reads all header-lines into a header_lines, a vector<string>
**
*/
bool			vcff::open( const char * filename )
{
	guard;

	//
	//
	if( whop_tabix::open( filename ) == false )
	{
		Rprintf("vcff::open : could not open tabix-index!\n");
		return false;
	}
	
	//	check for ##fileformat=VCFv4. (ignore whether followed by 0\n or 1\n for now)
	//
	const char * hdrline0 = getHeaderLine(0);
	if( (hdrline0 == 0) || (strncmp( hdrline0 ,"##fileformat=VCFv4.",19)!= 0) )
	{
		Rprintf("vcff::open : Not a VCF [%s]!\n",getHeaderLine(0));
		return false;
	}
	
	//	find and parse header-line defining the variant-line format (#CHROM	POS	...)
	//
	{
		class VcfParseSampleNames : public ParseFunctor {
		public:
			VcfParseSampleNames() : fmtline(0){}
			bool operator()(const char*s,int l)
			{
				//Rprintf("Formatline[%d]='%s'!!\n",l,s);
				if( (s==0) || (l<10) )
					return false;
				if( strncmp( s, "#CHROM\tPOS", 10 ) == 0 )
				{
					fmtline=(char*)s;
					return true;
				}
				return false;
			}
			char	*fmtline;
		} vpsn;

		//
		//	if the header line defining the line format of the variants was found...
		//
		if( parseHeader( vpsn ) )
		{
		
			//
			//
			TSVParser		tsv( vpsn.fmtline+1	);	//skip leading # of header line
			char			tokenbuf[256];
			unsigned int	i=0;
			sample_begin_index = 6;
			
			//
			//
			for( ; i < tsv.numFields(); i++ )
			{
				if( tsv.getField( &tokenbuf[0] , sizeof(tokenbuf)-1, i ) )
				{
					//Rprintf("Token = '%s'\n",tokenbuf);
					field_names.push_back( std::string( tokenbuf) );
					
					//	make sure the prescribed fields are there
					//
					if( ( i <= 8 ) && ( i >= 0 ) )
					{
						if( (i == 0) && ( strcmp(tokenbuf,"CHROM")!= 0 ) )	break;
						if( (i == 1) && ( strcmp(tokenbuf,"POS")!= 0 ) )	break;
						if( (i == 2) && ( strcmp(tokenbuf,"ID")!= 0 ) )		break;
						if( (i == 3) && ( strcmp(tokenbuf,"REF")!= 0 ) )	break;
						if( (i == 4) && ( strcmp(tokenbuf,"ALT")!= 0 ) )	break;
						if( (i == 5) && ( strcmp(tokenbuf,"QUAL")!= 0 ) )	break;
						if( i == 6 ){
							if( strcmp(tokenbuf,"FILTER")!= 0 ) break;
							else	sample_begin_index = 7;//i+1;
						}
						if( i == 7 ){
							if( strcmp(tokenbuf,"INFO")!= 0 ) break;
							else	sample_begin_index = 8;//i+1;
						}
						if( i == 8 ){
							if( strcmp(tokenbuf,"FORMAT")!= 0 ) break;
							else	sample_begin_index = 9;//i+1;
						}
					}//...if( one of the first 9 fields )

				}//...if( copied field from tsv-line )

			}//...for( each field in the tsv-line )

			//	determine the number of fields and the number of samples in this VCF
			//
			num_fields = field_names.size();
			num_samples = getNumFields() - getFirstSampleFieldIndex();

			//	Allocate buffer to memorize in-string start offsets of each tab-separated field
			//
			if( num_fields > 0 )
			{
				field_offsets = (unsigned int*)malloc( sizeof(unsigned int) * num_fields );
				field_offsets[0] = 0;
				field_offsets_size = num_fields;
				if( field_offsets == 0 )
					throw "vcff::open : failed to allocate buffer to memorize field offsets!";
			}
			else
			{
				Rprintf("whopgen::vcff::open : unexpected # of fields in TSV (%d<=0)!\n",num_fields);
				field_offsets=0;
			}
			
			//	Allocate buffer to memorize which samples were selected
			//
			//		use num_fields instead of num_samples to allocate a large enough buffer
			wanted_samples = (unsigned int*)malloc( sizeof( unsigned int ) * num_fields );
			if( wanted_samples == 0 )
				throw "vcff::open : failed to allocate buffer for sample-selection!";
			
			//
			Rprintf("vcff::open : file opened, contains %d samples\n",num_samples);
			
			//
		}//...if parseheader( vspn )
		else
		{
			Rprintf("vcff::open : could not find format-defining header line!\n");
			return false;
		}

		//
	}//...{

	return true;
	
	//
	unguards;
	if( field_offsets )	free( field_offsets );	field_offsets=0;
	if( wanted_samples ) free( wanted_samples ); wanted_samples=0;
	return false;
}




/*!
**
**
*/
unsigned int	vcff::getNumFields( void ){
	return field_names.size();
}




/*!
**
**
*/
const char*	vcff::getFieldName( unsigned int idx ){
	if( idx >= field_names.size() )
		return 0;
	return field_names[ idx ].c_str();
}




/*!
**
**
*/
unsigned int	vcff::getNumSamples( void )
{
	return num_samples;
}




/*!
**
**
*/
unsigned int	vcff::getFirstSampleFieldIndex( void )
{
	return sample_begin_index;
}





#if 0

/*!	Returns one of the header lines and iterates by using 
**		the "headerindex" member variable
**
*/
const char*	vcff::parseHeader( void )
{
	
	//	validate arguments
	//
	if( this == 0 )
	{
		Rprintf("VCF_parseheader : NULL vcf*!\n");
		return 0;
	}
	
	//
	if( vcf_tabixed == 0 )
	{
		Rprintf("VCF_parseheader : NULL vcf->tabix-index!\n");
		return 0;
	}

	//
	//
	if( headerindex < header_lines.size() )
	{
		
		//
		const char	*s=header_lines[headerindex].c_str();
		//Rprintf("VCF_parseheader line='%s'!!\n",s);
		headerindex++;
		
		//
		if( (s!= 0) && (s[0]=='#') )
		{
			return s;
		}
		
		//
	}
	
	//
	//Rprintf("VCF_parseheader false[%d]\n",headerindex);
	resetheader();
	return 0;
}



/*!
**
*/
const char*	vcff::parseregion( void )
{
	//
	if( this == 0 )
	{
		Rprintf("vcff::parseregion called on this==0!\n");
		return 0;
	}

	//
	if( vcftabix == 0 )
	{
		Rprintf("vcff::parseregion called on vcftabix==0!\n");
		return 0;
	}

	//
	if( variantiter == 0 )
	{
		Rprintf("vcff::parseregion called on variantiter==0!\n");
		return 0;
	}
	
	
	//
	int len=0;
	const char * s = ti_read(vcftabix, variantiter, &len);
	bEOR = (s==0);
	return s;

	//
}
#endif

/*!
**
*/

//-------------------------------------------------------


