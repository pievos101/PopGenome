/*
**
**		WhopGen
**
**		WHOle genome population genetics with popGEN
**
**
**		VCF Filtering functions
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

	//	R includes
	//
#include	<R.h>
#include	<Rinternals.h>
#include	<R_ext/Rdynload.h>



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
bool	get_subfield( filtered_fields& a, vcff* f )
{
	if( a.valid )
		return true;

	//
	const char * info = f->getFieldPtr( INFO );
	if( info == 0 )
	{
		Rprintf("field %d is %x\n",INFO,info);
		return false;
	}

//	Rprintf("INFO=[[");
//	for( int i=0; i < 80 && info[i] != '\t'; i++ )
//		Rprintf("%c",info[i] );
//	Rprintf("]]\n");

	//
	const char *fnam = a.fieldnam;

	//
//	Rprintf("fieldnam='%s'\n",fnam);
//	return false;

	//
	//
	int cnt=500;
	while( cnt-- )
	{
		int i=0;
		for( ; fnam[i] == *info ; i++ )
		{
			info++;
		}
//		Rprintf("last cmp was : %c <-> %c\n",fnam[i],*info);

		//
		//
		if( fnam[i] == 0 )
		{
			//<fnam[i] == 0 : has fieldname ended too ?>
			if( *info == '=' )
			{
				//
				a.data.fldptr = (char*)info+1;
				a.valid=true;

				//success, get ptr behind =
				//
			//	Rprintf("success:[[");
			//	for( int i=0; i < 40 && a.data.fldptr[i] != '\t'; i++ )
			//		Rprintf("%c",a.data.fldptr[i] );
			//	Rprintf("]]\n");

				//
				return true;
			}
			else
			if( *info == ';' || *info == '\t' || *info == 0 )	//fieldname is last entry and has no params ?( e.g. AC=200;FLAG_XYZ\t )
			{
				//success, set ptr to 0
			//	Rprintf("success: no data associated[[");
			//	for( int i=0; i < 40 && info[i] != '\t'; i++ )
			//		Rprintf("%c",info[i] );
			//	Rprintf("]]\n");
				a.data.fldptr = 0;
				a.valid=true;
				return true;
			}
			//else
			//{
			//	// fnam at end, but fieldname not yet
			//}
		}

		//fail, find next ';'
		//
//		Rprintf("	subfield didnt match, find next:(");
		for( ; *info != ';' ; info++ )
		{
//			Rprintf("%c",*info);
			if( *info == '\t' || *info == 0 )
				return false;
		}
		info++;
//		Rprintf("\n");

		//
	}

	//
	return false;
}

//---------------------------------------------

/*!	Applies the ruleset stored in the given VCFhandle and returns TRUE if this line passes
**
**
**
*/
bool	filterLine( vcff* f )
{
	int numfes = f->num_rules_used;

	//	new line to filter, invalidate all data for subfields of INFO
	//
	for( int i=0; i < numfes; i++ )
		f->fieldset[i].valid = false;

	//
	//
	for( int i=0; i < numfes; i++ )
	{

		// get ptr to data of field
		//
		Rprintf("FE#%d: get %d (%s)\n",i,f->ruleset[i].get_what, f->fieldset[ f->ruleset[i].get_what ].fieldnam );
		if( get_subfield( f->fieldset[ f->ruleset[i].get_what ] , f ) == false )
		{
			Rprintf("	could not get field!\n");
			//TODO : what to do if a subfield of INFO could not be found ?
			//	depends on f->fieldset : test for MUSTEXIST bit
			//
			continue;
		}

		//
		bool cmpresult = false;
		int iv;
		float fv;

		//	got ptr to subfield, not compare it
		//
		switch( f->ruleset[i].compare_how )
		{
			case DOES_EXIST:
				cmpresult = true;
				break;

			//	integer comparisons
			//
			case INT_CMP:
				iv = atoi( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				Rprintf("intval=%d\n",iv);
				cmpresult = (f->ruleset[i].i1 == iv );
				break;
			case INT_CMP_OO:
				iv = atoi( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].i1 < iv ) && (iv < f->ruleset[i].i2));
				break;
			case INT_CMP_OC:
				iv = atoi( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].i1 < iv ) && (iv <= f->ruleset[i].i2));
				break;
			case INT_CMP_CO:
				iv = atoi( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].i1 <= iv ) && (iv < f->ruleset[i].i2));
				break;
			case INT_CMP_CC:
				iv = atoi( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].i1 <= iv ) && (iv <= f->ruleset[i].i2));
				break;			//

			//	float comparisons
			//
			case FLT_CMP:
				fv = atof( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				Rprintf("fltval=%f\n",fv);
				cmpresult = (f->ruleset[i].f1 == fv );
				break;
			case FLT_CMP_OO:
				fv = atof( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				Rprintf("fltval=%f, range (%f,%f)\n",fv,f->ruleset[i].f1,f->ruleset[i].f2);
				cmpresult = ((f->ruleset[i].f1 < fv ) && (fv < f->ruleset[i].f2));
				break;
			case FLT_CMP_OC:
				fv = atof( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].f1 < fv ) && (fv <= f->ruleset[i].f2));
				break;
			case FLT_CMP_CO:
				fv = atof( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].f1 <= fv ) && (fv < f->ruleset[i].f2));
				break;
			case FLT_CMP_CC:
				fv = atof( f->fieldset[ f->ruleset[i].get_what ].data.fldptr );
				cmpresult = ((f->ruleset[i].f1 <= fv ) && (fv <= f->ruleset[i].f2));
				break;
			default:
				Rprintf("compare-how : unknown op!");
				break;
		}

		Rprintf("cmpresult=%s\n",cmpresult?"TRUE":"FALSE");

		//
		//
		switch( f->ruleset[i].do_what )
		{

			//	skip line if this filter matches
			//
			case DW_SKIP_NOT://	=	0x81,			//skip if true, dont test further
				Rprintf("action=skipnot..");
				cmpresult=!cmpresult;
			case DW_SKIP://		=	0x01,			//skip if false, dont test further
				Rprintf("action=skip\n");
				if( cmpresult )
					return false;
				break;

			//	keep line if this filter matches
			//
			case DW_KEEP_NOT://	=	0x82,			//keep if false, dont test further
				Rprintf("action=keepnot..");
				cmpresult=!cmpresult;
			case DW_KEEP://		=	0x02,			//keep if true, dont test further
				Rprintf("action=keep\n");
				if( cmpresult )
					return true;
				break;

			//
			case DW_NOP:
				Rprintf("NOP\n");
			default:
				Rprintf("UNKNOWN!\n");
				break;
		}

		//
	}//...for all rules

	return true;
}








/*!	Remove all filters
**
**
**
*/
EXPORT SEXP VCF_clearFilters( SEXP vcfptr )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}

	//
	f->num_rules_used = 0;
	f->num_fieldnames_used = 0;

	//
	for( int i=0; i < MAX_NUM_RULES; i++ )
	{
		f->ruleset[i].get_what = 0;
		f->ruleset[i].compare_how = 0;
		f->ruleset[i].do_what = 0;
		f->ruleset[i].i1 =
		f->ruleset[i].i2 = 0;
		f->ruleset[i].f1 =
		f->ruleset[i].f2 = 0.0f;
	}

	//
	for( int i=0; i < MAX_NUM_SUBFIELDS; i++ )
	{
		for( unsigned int j=0; j < sizeof( f->fieldset[0].fieldnam ); j++ )
			f->fieldset[i].fieldnam[j] = 0;
		f->fieldset[i].valid = false;
		f->fieldset[i].data.fldptr = 0;
		f->fieldset[i].data.intval = 0;
		f->fieldset[i].data.fltval = 0;
	}

	//
	return R_NilValue;
}


/*!	Add a filter
**
**
*/
EXPORT SEXP VCF_addFilter( SEXP vcfptr, SEXP fieldnam, SEXP cmptype, SEXP action, SEXP arg1, SEXP arg2 )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}

	//
	//
	Rprintf("used %d rules!\n",f->num_rules_used);
	if( f->num_rules_used >= MAX_NUM_RULES )
	{
		Rprintf("Too many rules already set!\n");
		return R_NilValue;
	}

	//
	//
	const char * fieldname = RString::get( fieldnam , 0 );
	if( fieldname == 0 )
	{
		Rprintf("Fieldname empty!\n");
		return R_NilValue;
	}

	//	check comparison type
	//
	int cmp = INTEGER( cmptype )[0];
	if( cmp < 0 || cmp > 10 )
	{
		Rprintf("cmptype %d not within [0,10]!\n",cmp);
		return R_NilValue;
	}
	Rprintf("cmptype=%d\n",cmp);

	//	check action type
	//
	int act = INTEGER( action )[0];
	if( ! (( act >= 0 && act <= 2) || (act == 0x81 || act == 0x82 )) )
	{
		Rprintf("acttype %d not valid!\n",act);
		return R_NilValue;
	}
	Rprintf("acttype=%d\n",act);

	//
	//
	filter_entry_t * fe = &f->ruleset[f->num_rules_used];

	fe->compare_how = cmp;			//how to compare

	//
	//
	int i1=0,i2=0;
	float f1=0,f2=0;
	if( cmp >= 1 && cmp <= 5 )
	{
		i1 = INTEGER( arg1 )[0];
		i2 = INTEGER( arg2 )[0];
	}
	else if( cmp >= 7 && cmp <= 10 )
	{
		f1 = REAL( arg1 )[0];
		f2 = REAL( arg2 )[0];
	}

	Rprintf("i %d,%d   f %f,%f\n",i1,i2,f1,f2);


	fe->i1 = i1;		//integer/float values to compare against
	fe->i2 = i2;
	fe->f1 = f1;
	fe->f2 = f2;


	Rprintf("%d\n",f->num_fieldnames_used);

	//-----
	//	add fieldname
	//
	int get_what_index=0;
	for( ; get_what_index < f->num_fieldnames_used; get_what_index++ )
	{
		Rprintf("%x\n",&(f->fieldset[get_what_index].fieldnam[0]));
		if( strcmp( fieldname, &(f->fieldset[get_what_index].fieldnam[0]) ) == 0 )
		{
			Rprintf("match at %d\n",get_what_index);
			break;
		}
	}

	//	if none of the fieldnames matches the one we want to filter against, add a new one
	//
	if( get_what_index >= f->num_fieldnames_used )
	{
		//
		//
		Rprintf("fieldname '%s' not yet found!\n",fieldname);

		//	check if there is room for one more subfield name
		//
		if( f->num_fieldnames_used >= MAX_NUM_SUBFIELDS )
		{
			Rprintf("Cannot use more fieldnames!\n");
			return R_NilValue;
		}

		//	add a new subfield name
		//
		filtered_fields*ffl = &(f->fieldset[get_what_index]);
		strcpy( ffl->fieldnam , fieldname );


		//
		f->num_fieldnames_used++;

		//
	}

	//

	//	add rules
	//
	fe->get_what = get_what_index;		//what subfieldname to look up

	fe->do_what = act;

	f->num_rules_used++;

	//
	return R_NilValue;
}


/*!
**
*/
EXPORT SEXP VCF_describeFilterConfig( SEXP vcfptr )
{

	//	get vcf
	//
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f )
	{
		Rprintf("Parameter not a VCFhandle EXTPTR!\n");
		return R_NilValue;
	}


	//-----
	const char * cmptype_to_string[] = {
		"EXISTS",
		"INT==",
		"INT(,)",
		"INT(,]",
		"INT[,)",
		"INT[,]",
		"FLT==",
		"FLT()",
		"FLT(]",
		"FLT[)",
		"FLT[]",
		0
	};

	//
	//
	Rprintf("Filtering rules:\n---------\n");
	for( int i=0; i < f->num_rules_used; i ++ )
	{
		filter_entry_t * t = &f->ruleset[i];
		Rprintf("#%d\t",i);
		Rprintf("%s(%d)\t", f->fieldset[ t->get_what ].fieldnam, t->get_what );
		Rprintf("%s(%d)\t", cmptype_to_string[ t->compare_how ], t->compare_how );

		//
		Rprintf("(%d)",t->do_what);
		switch( t->do_what )
		{
			case DW_SKIP_NOT:
				Rprintf("skipnot\t");
				break;
			case DW_SKIP:
				Rprintf("skip\t");
				break;
			case DW_KEEP_NOT:
				Rprintf("keepnot\t");
				break;
			case DW_KEEP:
				Rprintf("keep\t");
				break;
			case DW_NOP:
				Rprintf("NOP\t");
				break;
			default:
				Rprintf("UNKNOWN:%d!\t",t->do_what);
				break;
		}

		//
		Rprintf("\n");
	}
	Rprintf("---------\n");

	//
	return R_NilValue;
}



/*!
**
**
*/
EXPORT	SEXP	VCF_isSNP( SEXP vcfptr )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f ){	Rprintf("Parameter not a VCFhandle EXTPTR!\n");	return R_NilValue;	}
	if( _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ))
		return RBool::True();
	return RBool::False();
}


/*!
**
**
*/
EXPORT	SEXP	VCF_isInDel( SEXP vcfptr )
{
	vcff * f = (vcff*)R_GetExtPtr( vcfptr , "VCFhandle" );
	if( 0 == f ){	Rprintf("Parameter not a VCFhandle EXTPTR!\n");	return R_NilValue;	}
	if( ! _internal_isSNP( f->getFieldPtr( REF ), f->getFieldPtr( ALT ) ))
		return RBool::True();
	return RBool::False();
}


