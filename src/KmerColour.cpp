#include <cstdio>
#include <cstring>

#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif 

#include <stdio.h>
#include <assert.h>
#include <algorithm> // for min
#include "KmerColour.h"

#include <iostream>
//using namespace std;

//extern int sizeKmer;
//uint64_t nsolids = 0;

char* print_kmer_colour_pattern(KmerColour *colour, char *seq, int length) {
//	printf("Length:%d\n", length);
//	int p = 0;
	for (int i = 0; i < length; ++i) {
		seq[i] = colour[i] + 48;
	}
	seq[length] = '\0';

//		int j = (int) ((colour[i]));

//		p += snprintf(seq+p,length-p,"%x",colour[i]+3);

//old testing
//		seq[i] = (char) +(colour[i]);
//		int iaoeu = +(colour[i]);
//		cout <<  static_cast<unsigned>(colour[i]);
//		cout << to_string(colour[i]);

//		cout <<  colour[i];
//		cout <<  +colour[i];
//		printf("%u\t%u\n", colour[i], seq[i]);
//		printf("SEQ:%d\t%s\n", p, seq);

//	cout << endl;
//	printf("SEQ:%s\n", seq);
//	printf("SEQ:%s\n", seq);
	return seq;

}

void rev_colour(KmerColour *colour, int length) {

	KmerColour temp;
	for (int i = 0; i < length / 2; ++i) {
		temp = colour[length - i - 1];
		colour[length - i - 1] = colour[i];
		colour[i] = temp;
	}
}

int count_colour(KmerColour *colour, int length) {
	int count = 0;
	for (int i = 0; i < length; ++i) {
		count += number_of_colour(colour[i]);
//		printf("%u\n",colour);
	}
	return count;
}

int number_of_colour(KmerColour colour) {
	int count = 0;
	colour = 2;
	printf("%u\n", colour);
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
		printf("%u\n", colour);
	}

	return count;
}

int number_of_colour2(KmerColour *colour) {
	int count = 0;
//	*colour = 2;
	printf("%u\n", *colour);
	while (*colour) {
		count += *colour & 1;
		*colour = *colour >> 1;
		printf("%u\n", *colour);
	}

	return count;
}

int number_of_colour3(KmerColour &colour) {
	int count = 0;
//	*colour = 2;
	printf("%u\n", colour);
	printf("%x\n", &colour);
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
//		printf("%u\n",colour);
	}

	return count;
}

//using namespace KmerColourS;

int KmerColourC::number_of_colour_s(KmerColour colour) {
	int count = 0;
//		colour = 20;
	//	printf("%u\n",colour);
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
		//		printf("%u\n",colour);
	}

	return count;
}

int KmerColourN::number_of_colour_n(KmerColour colour) {
	int count = 0;
//		colour = 20;
	//	printf("%u\n",colour);
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
		//		printf("%u\n",colour);
	}

	return count;
}

void KmerColourC::get_all_colour(KmerColour *colour, int colour_len, int *all_colour){

	for (int i = 0; i < colour_len; ++i) {
		if(colour[i]>8){
			all_colour[i] = 0;
		}
		else{
			all_colour[i] = KmerColourC::number_of_colour_s(colour[i]);
		}
	}
}

int KmerColourC::append_colour(KmerColour* left_colour_traversal, long long len_left,
		KmerColour* contig_colour, int &colour_len) {
	memcpy(contig_colour + colour_len, left_colour_traversal, len_left); // contig = revcomp(left_traversal) //TODO: need to reverse as well??
	colour_len += len_left;
	return colour_len;
}





// debug only: convert a kmer_type to char*
//char debug_kmer_buffer[1024];
//char* print_kmer(kmer_type kmer)
//{
//    return print_kmer(kmer,sizeKmer,kmerMask);
//}
//
//char* print_kmer(kmer_type kmer, int sizeKmer, kmer_type kmerMask)
//{
//    code2seq(kmer,debug_kmer_buffer, sizeKmer, kmerMask);
//
//    return debug_kmer_buffer;
//}
