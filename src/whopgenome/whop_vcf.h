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

//*
//*			DEFINES
//*


enum compare_how_enum {
	DOES_EXIST	=	0,	//subfield exists?
	INT_CMP,			//read int and compare
	INT_CMP_OO,
	INT_CMP_OC,
	INT_CMP_CO,
	INT_CMP_CC,
	FLT_CMP,			//
	FLT_CMP_OO,
	FLT_CMP_OC,
	FLT_CMP_CO,
	FLT_CMP_CC,
	NUM_CMP_OPS
};

enum do_what_enum {
	DW_NOP		=	0x00,
	
	DW_SKIP		=	0x01,			//skip if false, dont test further
	DW_KEEP		=	0x02,			//keep if true, dont test further
	
	DW_SKIP_NOT	=	0x81,			//skip if true, dont test further
	DW_KEEP_NOT	=	0x82,			//keep if false, dont test further

	NUM_DO_WHATS
};

#define		MAX_NUM_RULES			5
#define		MAX_NUM_SUBFIELDS		3

//*
//*			STRUCTS
//*

//
//	entry lookup table
//

struct filtered_fields {
	char		fieldnam[64];
	bool		valid;
	union {
		char	*fldptr;
		int		intval;
		float	fltval;
	} data;
};

//	filter rules
//

struct filter_entry_t {
	int		get_what;				//index into info_ffl array : name and string-dataptr of a subfield
	int		compare_how;			//one of some constants describing comparison operation
	int		do_what;				//
	int		i1,i2;
	float	f1,f2;
};

//*
//*			CLASSES
//*



/*!	Class wrapping working with VCF files
**
**
**
*/
class vcff : public whop_tabix
{
public:
	
	//-
	//-		CTORs + DTOR
	//-

	vcff();
	vcff(const char * filename );
	~vcff();
	
	//-
	//-		METHODS
	//-
	
	//
	//
	bool			open( const char * filename );
	
	//
	//
	unsigned int	getNumFields( void );
	const char*	getFieldName( unsigned int idx );
	unsigned int	getNumSamples( void );
	unsigned int	getFirstSampleFieldIndex( void );
	
	//
	unsigned int	lastReadPos( void ){
		return _last_readline_pos;
	}
	
	//
	void			resetSampleSelection( void )
	{
		num_wanted_samples=0;
	}

	bool			selectSample( unsigned int samplefieldidx )
	{
		if( samplefieldidx >= num_fields ){	Rprintf("sampleidx=%d >= %d=num_fields",samplefieldidx,num_fields); return false;}
		if( num_wanted_samples >= num_samples ){ Rprintf("num_wanted_samples=%d >= %d=num_samples",num_wanted_samples,num_samples); return false;}
		wanted_samples[ num_wanted_samples ] = samplefieldidx;
		wanted_samples[ num_wanted_samples+1 ] = -1;
		num_wanted_samples++;
		return true;
	}

	//
	//
	bool			testfunc( void );

	
	//-
	//-		DATA
	//-
	
	//	Filter rules
	//
	int							num_rules_used;
	int							num_fieldnames_used;
	struct filter_entry_t		ruleset[ MAX_NUM_RULES ];
	struct filtered_fields	fieldset[ MAX_NUM_SUBFIELDS ];
	
	//	Parse-related (transient) data
	//

		//	samples to extract
	
	unsigned int				*wanted_samples;	//
	unsigned int				num_wanted_samples;


	//-
	//-		DATA
	//-

protected:
	
	//	VCF specific data
	//
	unsigned int				num_fields;			//number of fields in VCF/TSV file
	unsigned int				num_samples;		//number of samples in VCF
	std::vector<std::string>	field_names;	//CHROM,POS,ID,REF,ALT,... from the "#CHROM	POS	ID	REF	ALT..." line
	unsigned int				sample_begin_index;
	unsigned int				_last_readline_pos;

private:

	//
};

