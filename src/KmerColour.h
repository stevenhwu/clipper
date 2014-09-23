#ifndef KmerColour_h
#define KmerColour_h

#include <string>
#include <stddef.h>

const int kErrorCode = 100;//TODO: with errorCode, fix or flexible? define only == or >=

typedef unsigned char KmerColour;
//typedef int kmer_colour;

const size_t kSizeOfKmerColour = sizeof(KmerColour);

char* kmer_colour_pattern_string(KmerColour *colour, char *seq, int length, int error_code=kErrorCode);

void rev_colour(KmerColour *colour, int length);
int count_colour(KmerColour *colour, int length);
int number_of_colour(KmerColour colour);
int number_of_colour2(KmerColour *colour);
int number_of_colour3(KmerColour &colour);

class KmerColourUtil {
public:

//	static int number_of_colour_c(KmerColour colour);
	static int number_of_colour_c(KmerColour colour) {
		int count = 0;
		while (colour) {
			count += colour & 1;
			colour = colour >> 1;
		}

		return count;
	}

	static void get_all_colour_count2(KmerColour *colour_seq, int colour_len,
				int *all_colour, int error_code = kErrorCode);

	static int number_of_colour_s(KmerColour colour);

	static void get_all_colour_count(KmerColour *colour_seq, int colour_len,
			int *all_colour, int error_code = kErrorCode);

	static int append_colour(KmerColour* left_colour_traversal,
			long long len_left, KmerColour* contig_colour, int &colour_len);

	static int colour_table(std::string &report, KmerColour *colour,
			int colour_len, int max_colour_count);

	static std::string summary(std::string &report, KmerColour *kmer_colour,
			int colour_len);

	static void summary_colour_code(std::string &report, KmerColour *kmer_colour,
			int colour_len);

	static void summary_colour_count(std::string &report, KmerColour *kmer_colour,
			int colour_len);

};

class KmerColourSummary {

private:
	int length;
	KmerColour* kmer_colour;
	int* colour_count;
	int* delta_colour_count;

public:
	KmerColourSummary(KmerColour *colour, int colour_len): kmer_colour(colour), length(colour_len){}
	void full_summary();

};



namespace KmerColourN {

int number_of_colour_n(KmerColour colour);
//static int number_of_colour_ns(KmerColour colour);
}
//
//
//int number_of_colour_s(KmerColour colour)
//void get_all_colour(KmerColour *colour_seq, int colour_len, int *all_colour);
//int append_colour(KmerColour* left_colour_traversal, long long len_left,
//		KmerColour* contig_colour, int &colour_len);
//}

//
//extern int sizeKmer;
//extern kmer_type kmerMask;
//extern uint64_t nsolids;
//
//
//int NT2int(char nt);
//int revcomp_int(int nt_int);
//kmer_type  codeSeed(char *seq, int sizeKmer, kmer_type kmerMask);
//kmer_type  codeSeed(char *seq);
//kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read);
//kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask);
//kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read);
//kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask);
//unsigned char  code_n_NT(char *seq, int nb);
//unsigned char  code4NT(char *seq);
//
//uint64_t revcomp(uint64_t x);
//uint64_t revcomp(uint64_t x, int size);
//
//#ifdef _largeint
//LargeInt<KMER_PRECISION> revcomp(LargeInt<KMER_PRECISION> x);
//LargeInt<KMER_PRECISION> revcomp(LargeInt<KMER_PRECISION> x, int size);
//#endif
//#ifdef _ttmath
//ttmath::UInt<KMER_PRECISION> revcomp(ttmath::UInt<KMER_PRECISION> x);
//ttmath::UInt<KMER_PRECISION> revcomp(ttmath::UInt<KMER_PRECISION> x, int size);
//#endif
//#ifdef _LP64
//__uint128_t revcomp(__uint128_t x);
//__uint128_t revcomp(__uint128_t x, int size);
//#endif
//
//int code2seq ( kmer_type code,char *seq);
//int code2seq ( kmer_type code,char *seq, int sizeKmer, kmer_type kmerMask);
//int code2nucleotide( kmer_type code, int which_nucleotide);
//
//kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp);
//kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential);
//kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential, int sizeKmer, kmer_type kmerMask);
//
//// compute the next kmer w.r.t forward or reverse strand, e.g. for ACTG (revcomp = CAGT)
//// it makes sure the result is the min(kmer,revcomp_kmer)
//// indicates if the result is the revcomp_kmer by setting *strand
//// examples:
//// next_kmer(ACTG,A,&0)=CTGA with strand = 0 (because revcomp=TCAG);
//// next_kmer(ACTG,A,&1)= (revcomp of ACTG + A = CAGT+A = ) AGTA with strand = 0 (because revcomp = TACT)
//kmer_type next_kmer(kmer_type graine, int added_nt, int *strand);
//
//void revcomp_sequence(char s[], int len);
//
//
//kmer_type  codeSeed_bin(char *seq);
//
//kmer_type  codeSeedRight_bin(char *seq, kmer_type  val_seed, bool new_read);
//
//kmer_type  codeSeedRight_revcomp_bin(char *seq, kmer_type  val_seed, bool new_read);
//
//kmer_type extractKmerFromRead_bin(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool use_compressed);
//
//char* print_kmer(kmer_type kmer); // debugging
//char* print_kmer(kmer_type kmer, int sizeKmer, kmer_type kmerMask); // debugging
//

#endif
