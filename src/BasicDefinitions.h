// BasicDefinitions.h: constants and enumerations for GERMLINE

#ifndef BASICDEFINITIONS_H
#define BASICDEFINITIONS_H
// #define VERBOSE
#include <vector>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;

// where we are in the sequence (markerset & physical)
extern unsigned int position_ms;
extern unsigned int num_sets;
extern size_t num_samples;

extern unsigned long num_matches;

// parameters
extern double MIN_MATCH_LEN;
extern int MARKER_SET_SIZE;
extern int MAX_ERR_HOM;
extern int MAX_ERR_HET;
extern int ERR_W;
extern bool PRINT_MATCH_HAPS;
extern bool ROI;
extern bool HAPLOID;
extern bool HAP_EXT;
extern bool WIN_EXT;
extern bool ALLOW_HOM;
extern bool HOM_ONLY;
extern bool SILENT;
extern bool DEBUG;
extern bool BINARY_OUT;
extern bool REDUCE;
extern bool NO_SUFFIX;

class SNPs;
class Individuals;
extern SNPs ALL_SNPS;
extern Individuals ALL_SAMPLES;
extern ofstream MATCH_FILE; //bmatch filr
//extern ofstream MATCH_FILE2; //match file
extern string PEDFILE;
extern string MAPFILE;
extern string OUTFILE;
extern double IBD_THRESHOLD;
enum ErrorType{RECOMB=0,MI=1};
const int HET=2;
const int MIS=9;
const short PAR_M = 0;
const short PAR_F = 1;
extern bool IBD;
// type for file format
enum FileFormat{HAPS,PED,HM};
// type for families
enum FamilyType{NOPARENT,SINGLEPARENT,TWOPARENT};
// type for sex
enum Sex{MALE,FEMALE,UNKNOWN};
// type for nucleotides
enum Nucleotide{A=1,C=2,G=3,T=4};
// type for chromosomes types
enum ChromosomeType{TRANS=0,UNTRANS=1,MISSING=-1};
// type for chromosome ids
enum ChromosomeID{ONE, TWO, THREE, FOUR, FIVE, SIX, SEVEN, EIGHT, NINE, TEN, ELEVEN, TWELVE,
                THIRTEEN, FOURTEEN, FIFTEEN, SIXTEEN, SEVENTEEN, EIGHTEEN, NINETEEN, TWENTY,
				TWENTYONE, TWENTYTWO, X, Y};
#endif

// end of BasicDefinitions.h
