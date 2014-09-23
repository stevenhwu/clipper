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
#include <cstdlib>
#include <stdlib.h>
//using namespace std;

//extern int sizeKmer;
//uint64_t nsolids = 0;

char* kmer_colour_pattern_string(KmerColour *colour, char *seq, int length, int error_code) {
//	printf("Length:%d\n", length);
//	int p = 0;
	for (int i = 0; i < length; ++i) {
		if(colour[i]< error_code){
			seq[i] = colour[i] + 48;
		}
		else{
			seq[i] = '#';//change to something later
		}
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

int KmerColourUtil::number_of_colour_s(KmerColour colour) {
	int count = 0;
	while (colour) {
		count += colour & 1;
		colour = colour >> 1;
	}

	return count;
}

void KmerColourUtil::get_all_colour_count(KmerColour *colour, int colour_len, int *all_colour, int error_code){

	for (int i = 0; i < colour_len; ++i) {
		if(colour[i] < error_code){
			all_colour[i] = KmerColourUtil::number_of_colour_s(colour[i]);
		}
		else{
			all_colour[i] = 0;
		}
	}
}

int KmerColourUtil::append_colour(KmerColour* left_colour_traversal, long long len_left,
		KmerColour* contig_colour, int &colour_len) {
	memcpy(contig_colour + colour_len, left_colour_traversal, len_left); // contig = revcomp(left_traversal) //TODO: need to reverse as well??
	colour_len += len_left;
	return colour_len;
}

int KmerColourUtil::colour_table(std::string &report, KmerColour *colour, int colour_len, int max_colour_count){

//	report.clear();
//	report.
	KmerColour temp = colour[0];
	int max = 0;
	while (temp) {
		temp = temp >> 1;
		max++;
	}
	max=max_colour_count;
	char matrix[max][colour_len];
	printf("MAX_COLOUR:%d\n",max);
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
void KmerColourUtil::summary_colour_code (std::string &report, KmerColour *kmer_colour, int colour_len){
	char colour_seq[colour_len];
	kmer_colour_pattern_string(kmer_colour, colour_seq, colour_len);
	report.append("ColourCode:").append(colour_seq).append("\n");

}

void KmerColourUtil::summary_colour_count(std::string &report, KmerColour *kmer_colour, int colour_len){

	int *all_colour = (int *) malloc(sizeof(int)*colour_len);
//	int all_colour[colour_len];
	get_all_colour_count(kmer_colour, colour_len, all_colour);
	report.append("ColourCount:");
	for (int i = 0; i < colour_len; ++i) {
		report+= std::to_string(all_colour[i]) ;
	}
	report += "\n";
}

std::string KmerColourUtil::summary (std::string &report, KmerColour *kmer_colour_code, int colour_len){

	summary_colour_code(report, kmer_colour_code, colour_len);
	summary_colour_count(report, kmer_colour_code, colour_len);
	int *all_colour_count = (int *) malloc(sizeof(int)*colour_len);
	get_all_colour_count(kmer_colour_code, colour_len, all_colour_count);


	char temp_char[10000];
	int update_len = colour_len;
//	std::unordered_map<int, int> counter;
	std::map<int, int> counter;
	std::map<int, int> delta_colour_pattern;
	std::map<KmerColour, int> kmer_colour_pattern;

	int no_of_swap = 0;
	int no_of_swap_exclude_n = 0;
	const int exclude = 0;

	int last_colour = all_colour_count[0];
	int last_colour_exclude = exclude;
	KmerColour last_kmer_colour = kmer_colour_code[0];

//	int delta_count[colour_len];
	std::string delta_string;
	int colour_delta = 0;
	for (int i = 0; i < colour_len; ++i) {
		int temp_count = all_colour_count[i];
		KmerColour temp_kmer = kmer_colour_code[i];
		counter[temp_count]++;
		if(temp_count > kErrorCode){
			update_len--;
		}
		else{
			if (last_colour != temp_count) {
				no_of_swap++;
//				printf("%d_%s:%s\n",delta_count[i], std::to_string(delta_count[i]).data(), delta_string.data());
//				itoa(t, delta_count, 10);
			}


			if (temp_count != exclude) {
				if (last_colour_exclude != temp_count
						&& last_colour_exclude != exclude) {
					no_of_swap_exclude_n++;

					colour_delta = temp_count-last_colour_exclude;
				}
				else{
					colour_delta = 0;
				}

//				delta_count[i] = KmerColourUtil::number_of_colour_s(delta_kmer);

//				delta_colour_pattern[colour_delta]++;
//				delta_string += std::to_string(colour_delta);

				last_colour_exclude = temp_count;
			}
			colour_delta = temp_count-last_colour;
			delta_colour_pattern[colour_delta]++;
			delta_string += std::to_string(colour_delta);
			last_colour = temp_count;

			KmerColour delta_kmer = last_kmer_colour ^ temp_kmer;
			kmer_colour_pattern[delta_kmer]++;
			last_kmer_colour = temp_kmer;
		}

	}
//	char* delta_string = reinterpret_cast<char*>(delta_count);

//	printf("SWAP:%d %d\n",no_of_swap, no_of_swap_exclude_n);

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
		if (it.first != 0 && it.first < kErrorCode) {
			keys.push_back(it.first);
			mean += it.second*it.first;
			mean_length += it.second;
//		}
		sprintf(temp_char, "ColourCode:%d  Count:%d/%d Percentage:%.3f%% \n", it.first, it.second,
				update_len, 100 * it.second / (double) update_len);
//		printf(temp_char);
		report.append(temp_char);
		}

	}
	mean /= update_len;
	mean_length /= (keys.size() * update_len); //TODO: actually, it's not quite correct, 1,3 == 2,3 here!!
	report += "Total length exclude (" + std::to_string(exclude) + "):"
			+ std::to_string(update_len) + "(" + std::to_string(colour_len)
			+ ")\n";
	report += "Number of swap:" + std::to_string(no_of_swap) +"\n";
	report += "Number of swap exclude (" + std::to_string(exclude) + "):"
			+ std::to_string(no_of_swap_exclude_n) + "\n";

	report += "Swap Pattern: " + delta_string + "\n";
	report += "Mean colour count per site :" + std::to_string(mean) + "\n";
	report += "Mean length per colour count (delete later):" +std::to_string(mean_length) + "\n";

	for(auto d : delta_colour_pattern){
		sprintf(temp_char, "Change:%d  Count:%d\n", d.first, d.second);
		report.append(temp_char);
	}

	for(auto k : kmer_colour_pattern){
		sprintf(temp_char, "pattern:%d  Count:%d\n", k.first, k.second);
		report.append(temp_char);
	}


//	printf("\nSIZE: %zu\t%zu\n", report.size(), report.max_size());
//	report.append("Mean colour count per site:");


//	printf("%s\n", report.data());
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
