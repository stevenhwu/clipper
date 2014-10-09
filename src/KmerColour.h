#ifndef KmerColour_h
#define KmerColour_h

#include <string>
#include <stddef.h>
#include <unordered_map>
#include <map>
//#include <utility>
#include <vector>

const int kErrorCode = 100;//TODO: with errorCode, fix or flexible? define only == or >=
//TODO: KmerColour unsign char 0-255, 8 species. dynamicly change to unsign short/int/long? for more species?

typedef unsigned char KmerColour; //TODO char -> assume no more than 8 for now
//typedef int kmer_colour;

const size_t kSizeOfKmerColour = sizeof(KmerColour);

int kmer_colour_pattern_string(KmerColour *colour, int length, char *seq,
		int error_code = kErrorCode);



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
	static int kmer_colour_pattern_string(KmerColour *colour, int length,
			char *seq, int error_code = kErrorCode);

	static void rev_colour(KmerColour *colour, int length);
	static int count_colour(KmerColour *colour, int length);

	static int number_of_colour(KmerColour colour);

	static void get_all_colour_count(KmerColour *colour_seq, int colour_len,
			int *all_colour, int error_code = kErrorCode);

	static int append_colour(KmerColour* left_colour_traversal,
			long long len_left, KmerColour* contig_colour, int &colour_len);

	static int colour_table(std::string &report, KmerColour *colour,
			int colour_len, int max_colour_count) __attribute__ ((deprecated));

	static std::string summary(std::string &report, KmerColour *kmer_colour,
			int colour_len) __attribute__ ((deprecated));

	static void summary_colour_code(std::string &report,
			KmerColour *kmer_colour, int colour_len)
					__attribute__ ((deprecated));

	static void summary_colour_count(std::string &report,
			KmerColour *kmer_colour, int colour_len)
					__attribute__ ((deprecated));

};

class KmerColourSummary {


struct colourCounter{



};
private:
	KmerColour* kmer_colour;
	int length;
	int max_colour;

	int valid_length;
	double coverage;
	std::vector<int> colour_count;
	std::vector<int> delta_colour_count;
	std::vector<char> colour_code;//	std::string colour_code;

	std::unordered_map<int, int> colour_counter;
	std::unordered_map<int, int> colour_counter_percentage;
	std::map<int, int> delta_colour_pattern;
	std::map<KmerColour, int> delta_kmer_colour_pattern;

	std::map<std::pair<KmerColour, int>, int> delta_pattern;

	bool keep = true;
	std::string reason;

	char** matrix;

	void preprocess_kmer_colour();

//	int* colour_count;
//	char* colour_code;

public:
	KmerColourSummary(KmerColour *colour, int colour_len, int max_colour);
//	:
//			kmer_colour(colour), length(colour_len), max_colour(max_colour) {
//		printf("Constructor in header\n");
//	}
	~KmerColourSummary();
	void full_summary();


	void summary_colour_code(std::string &report);
	void summary_colour_count(std::string &report);
	void summary_stat(std::string &report);
	int colour_table(std::string &report);

//	void KmerColourUtil::summary_colour_count(std::string &report,
//			KmerColour *kmer_colour, int colour_len) {
//
//		int *all_colour = (int *) malloc(sizeof(int)*colour_len);
//	//	int all_colour[colour_len];
//		get_all_colour_count(kmer_colour, colour_len, all_colour);
//		report.append("ColourCount:");
//		for (int i = 0; i < colour_len; ++i) {
//			report+= std::to_string(all_colour[i]) ;
//		}
//		report += "\n";
//	}


//	all_colour_count(int *all_colour, int error_code = kErrorCode);

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
