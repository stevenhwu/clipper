#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif 

#include <stdio.h>
#include <assert.h>
#include <algorithm> // for min
#include "KmerColour.h"

#include <iostream>
#include <unordered_map>
#include <map>
#include <memory>
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
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
	}

	return count;
}

//using namespace KmerColourS;

int KmerColourC::number_of_colour_s(KmerColour colour) {
	int count = 0;
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
	}

	return count;
}

void KmerColourC::get_all_colour(KmerColour *colour, int colour_len, int *all_colour){

	for (int i = 0; i < colour_len; ++i) {
		all_colour[i] = KmerColourC::number_of_colour_s(colour[i]);
	}
}

int KmerColourC::append_colour(KmerColour* left_colour_traversal, long long len_left,
		KmerColour* contig_colour, int &colour_len) {
	memcpy(contig_colour + colour_len, left_colour_traversal, len_left); // contig = revcomp(left_traversal) //TODO: need to reverse as well??
	colour_len += len_left;
	return colour_len;
}

int KmerColourC::colour_table(std::string &report, KmerColour *colour, int colour_len, int max_colour_count){

	report.clear();
//	report.
	KmerColour temp = colour[0];
	int max = 0;
	while (temp) {
		temp = temp >> 1;
		max++;
	}
	char matrix[max][colour_len];

	for (int i = 0; i < colour_len; ++i) {

		KmerColour temp = colour[i];
		for (int j = 0; j < max; ++j) {
			if(temp & 1){
				matrix[j][i] = '*';
			}
			else{
				matrix[j][i] = '.';
			}
			temp = temp >> 1;
		}
	}
	for (int j = 0; j < max; ++j) {
		for (int i = 0; i < colour_len; ++i) {
			printf("%c",matrix[j][i]);
			matrix[j][i]='0';
		}
		printf("\n");
	}

	return 0;


}

std::string KmerColourC::summary (std::string &report, int *all_colour, int colour_len){

	char temp_char[10000];
	int update_len = colour_len;
//	std::unordered_map<int, int> counter;
	std::map<int, int> counter;

	int no_of_swap = 0;
	int no_of_swap_exclude_n = 0;
	int temp_colour = all_colour[0];

	const int exclude = 0;
	int temp_colour_exclude = exclude;

	for (int i = 0; i < colour_len; ++i) {
		counter[all_colour[i]]++;
		if( temp_colour != all_colour[i] ){
			no_of_swap++;
		}
		temp_colour = all_colour[i];

		if (all_colour[i] != exclude) {
			if( temp_colour_exclude != all_colour[i] && temp_colour_exclude!=exclude){
				no_of_swap_exclude_n ++;
			}
			temp_colour_exclude = all_colour[i];
		}



	}

//	printf("SWAP:%d %d\n",no_of_swap, no_of_swap_exclude_zero);

//	typedef std::unordered_map<int, int> MapIntInt;
//	for (MapIntInt::iterator it = counter.begin(); it != counter.end(); ++it ) {
//		printf("%d:%d\n",it->first, it->second);
//	}

	if (counter.count(0)){
		update_len -= counter[0];
	}
	std::vector<int> keys;
	double mean = 0;
	double mean_length = 0;
	for (auto &it : counter){
		if (it.first != 0) {
			keys.push_back(it.first);
			mean += it.second*it.first;
			mean_length += it.second;
		}
		sprintf(temp_char, "No. of Colour:%d  Count:%d/%d Percentage:%.3f%% \n", it.first, it.second,
				update_len, 100 * it.second / (double) update_len);
		printf(temp_char);
		report.append(temp_char);

	}
	mean /= update_len;
	mean_length /= keys.size(); //TODO: actually, it's not quite correct, 1,3 == 2,3 here!!
	report += "Total length exclude (" + std::to_string(exclude) + "):" + std::to_string(update_len) +"\n";
	report += "Number of swap:" + std::to_string(no_of_swap) +"\n";
	report += "Number of swap exclude (" + std::to_string(exclude) + "):"
			+ std::to_string(no_of_swap_exclude_n) + "\n";
	report += "Mean colour count per site:" + std::to_string(mean) + "\n";
	report += "Mean length per colour count:" +std::to_string(mean_length) + "\n";
//	printf("\nSIZE: %zu\t%zu\n", report.size(), report.max_size());
//	report.append("Mean colour count per site:");


//	printf("%s\n%x\n", report.data(), pt);
	report += "==========\n";
	return report;

}


int KmerColourN::number_of_colour_n(KmerColour colour) {
	int count = 0;
		while (colour) {
			count += colour & 1;
			colour = colour >> 1;
		}
	return count;
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
