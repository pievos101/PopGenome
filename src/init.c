#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ap_pop_ancestral_C(SEXP);
extern SEXP ap_pop_C(SEXP);
extern SEXP Ccompare(SEXP, SEXP);
extern SEXP C_get_sfreqh_C(SEXP, SEXP);
extern SEXP code_nucs(SEXP);
extern SEXP combnsum2_C(SEXP, SEXP);
extern SEXP combnsum_C(SEXP);
extern SEXP compute_FREQ_C(SEXP);
extern SEXP compute_FREQOUT_C(SEXP);
extern SEXP count_congruent(SEXP);
extern SEXP FASTA_getNextIndividual(SEXP, SEXP);
extern SEXP FASTA_open(SEXP, SEXP, SEXP);
extern SEXP find_lines_GFF(SEXP, SEXP);
extern SEXP find_lines_GFF_Human2(SEXP, SEXP);
extern SEXP find_lines_SNP(SEXP, SEXP);
extern SEXP find_windowC(SEXP, SEXP, SEXP, SEXP);
extern SEXP fittingGFFC(SEXP, SEXP);
extern SEXP get_dim_fasta(SEXP);
extern SEXP get_gff_info_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_ind_fasta(SEXP, SEXP, SEXP);
extern SEXP makeBialMatrix(SEXP);
extern SEXP makeBialMatrixinclude(SEXP);
extern SEXP myReadVCFC(SEXP);
extern SEXP my_unique_C(SEXP);
extern SEXP pimpMatrix(SEXP, SEXP);
extern SEXP polyC(SEXP);
extern SEXP polyCinclude(SEXP);
extern SEXP R2_between_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP R2_C(SEXP, SEXP, SEXP);
extern SEXP R2_C_plus(SEXP, SEXP, SEXP, SEXP);
extern SEXP readdna(SEXP);
extern SEXP split_VCF_scaffolds(SEXP, SEXP);
extern SEXP VCF_getContigNames(SEXP);
extern SEXP VCF_getSampleNames(SEXP);
extern SEXP VCF_open(SEXP);
extern SEXP VCF_readIntoCodeMatrix(SEXP, SEXP);
extern SEXP VCF_readIntoCodeMatrixdiploid2(SEXP, SEXP);
extern SEXP VCF_selectSamples(SEXP, SEXP);
extern SEXP VCF_setRegion(SEXP, SEXP, SEXP, SEXP);
extern SEXP verify_ancestral_C(SEXP);
extern SEXP whichbigger_C(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ap_pop_ancestral_C",             (DL_FUNC) &ap_pop_ancestral_C,             1},
    {"ap_pop_C",                       (DL_FUNC) &ap_pop_C,                       1},
    {"Ccompare",                       (DL_FUNC) &Ccompare,                       2},
    {"C_get_sfreqh_C",                 (DL_FUNC) &C_get_sfreqh_C,                 2},
    {"code_nucs",                      (DL_FUNC) &code_nucs,                      1},
    {"combnsum2_C",                    (DL_FUNC) &combnsum2_C,                    2},
    {"combnsum_C",                     (DL_FUNC) &combnsum_C,                     1},
    {"compute_FREQ_C",                 (DL_FUNC) &compute_FREQ_C,                 1},
    {"compute_FREQOUT_C",              (DL_FUNC) &compute_FREQOUT_C,              1},
    {"count_congruent",                (DL_FUNC) &count_congruent,                1},
    {"FASTA_getNextIndividual",        (DL_FUNC) &FASTA_getNextIndividual,        2},
    {"FASTA_open",                     (DL_FUNC) &FASTA_open,                     2},
    {"find_lines_GFF",                 (DL_FUNC) &find_lines_GFF,                 2},
    {"find_lines_GFF_Human2",          (DL_FUNC) &find_lines_GFF_Human2,          2},
    {"find_lines_SNP",                 (DL_FUNC) &find_lines_SNP,                 2},
    {"find_windowC",                   (DL_FUNC) &find_windowC,                   4},
    {"fittingGFFC",                    (DL_FUNC) &fittingGFFC,                    2},
    {"get_dim_fasta",                  (DL_FUNC) &get_dim_fasta,                  1},
    {"get_gff_info_C",                 (DL_FUNC) &get_gff_info_C,                 4},
    {"get_ind_fasta",                  (DL_FUNC) &get_ind_fasta,                  3},
    {"makeBialMatrix",                 (DL_FUNC) &makeBialMatrix,                 1},
    {"makeBialMatrixinclude",          (DL_FUNC) &makeBialMatrixinclude,          1},
    {"myReadVCFC",                     (DL_FUNC) &myReadVCFC,                     1},
    {"my_unique_C",                    (DL_FUNC) &my_unique_C,                    1},
    {"pimpMatrix",                     (DL_FUNC) &pimpMatrix,                     2},
    {"polyC",                          (DL_FUNC) &polyC,                          1},
    {"polyCinclude",                   (DL_FUNC) &polyCinclude,                   1},
    {"R2_between_C",                   (DL_FUNC) &R2_between_C,                   4},
    {"R2_C",                           (DL_FUNC) &R2_C,                           3},
    {"R2_C_plus",                      (DL_FUNC) &R2_C_plus,                      4},
    {"readdna",                        (DL_FUNC) &readdna,                        1},
    {"split_VCF_scaffolds",            (DL_FUNC) &split_VCF_scaffolds,            2},
    {"VCF_getContigNames",             (DL_FUNC) &VCF_getContigNames,             1},
    {"VCF_getSampleNames",             (DL_FUNC) &VCF_getSampleNames,             1},
    {"VCF_open",                       (DL_FUNC) &VCF_open,                       1},
    {"VCF_readIntoCodeMatrix",         (DL_FUNC) &VCF_readIntoCodeMatrix,         2},
    {"VCF_readIntoCodeMatrixdiploid2", (DL_FUNC) &VCF_readIntoCodeMatrixdiploid2, 2},
    {"VCF_selectSamples",              (DL_FUNC) &VCF_selectSamples,              2},
    {"VCF_setRegion",                  (DL_FUNC) &VCF_setRegion,                  4},
    {"verify_ancestral_C",             (DL_FUNC) &verify_ancestral_C,             1},
    {"whichbigger_C",                  (DL_FUNC) &whichbigger_C,                  2},
    {NULL, NULL, 0}
};

void R_init_PopGenome(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

