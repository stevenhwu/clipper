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

int kmer_colour_pattern_string(KmerColour *colour, int length, char *seq, int error_code) {
//	printf("Length:%d\n", length);
	int valid_length = 0;
	for (int i = 0; i < length; ++i) {
		if(colour[i]< error_code){
			seq[i] = colour[i] + 48;
			valid_length++;
		}
		else{
			seq[i] = '#';//change to something later
		}
	}
//	printf("%d %d\n",length, valid_length);
//	if(valid_length < length){
	seq[length] = '\0';
//	}
	return valid_length;
//	seq[length] = '\0';

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
//	return seq;

}

void KmerColourUtil::rev_colour(KmerColour *colour, int length) {

	KmerColour temp;
	for (int i = 0; i < length / 2; ++i) {
		temp = colour[length - i - 1];
		colour[length - i - 1] = colour[i];
		colour[i] = temp;
	}
}

int KmerColourUtil::count_colour(KmerColour *colour, int length) {
	int count = 0;
	for (int i = 0; i < length; ++i) {
		count += KmerColourUtil::number_of_colour(colour[i]);
//		printf("%u\n",colour);
	}
	return count;
}



int KmerColourUtil::number_of_colour(KmerColour colour) {
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
			all_colour[i] = KmerColourUtil::number_of_colour(colour[i]);
		}
		else{
			all_colour[i] = 0;
		}
	}
}

int KmerColourUtil::append_colour(KmerColour* colour_traversal, long long len,
		KmerColour* contig_colour, int &colour_len) {
	memcpy(contig_colour + colour_len, colour_traversal, len);
	colour_len += len;
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
	char colour_seq[colour_len+1];
	kmer_colour_pattern_string(kmer_colour, colour_len, colour_seq);
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

KmerColourSummary::KmerColourSummary(KmerColour *colour, int colour_len, int max_colour):
		kmer_colour(colour), length(colour_len), max_colour(max_colour) {

	colour_count.reserve(length);
	delta_colour_count.reserve(length);
	colour_code.reserve(length);

	matrix = new char*[max_colour];
	for (int i = 0; i < max_colour; ++i) {
		matrix[i] = new char[length];
	}
//	auto matrix = new char[max_colour][colour_len];
	preprocess_kmer_colour();

//	colour_code = new char[length+1];
//	colour_count = new int[length];
}
void KmerColourSummary::preprocess_kmer_colour(){


	valid_length = 0;
	for (int i = 0; i < length; ++i) {
		if(kmer_colour[i]< kErrorCode ){
			colour_code.push_back(kmer_colour[i] + 48);
			colour_count.push_back(KmerColourUtil::number_of_colour(kmer_colour[i]) );
			valid_length++;
		}
		else{
			colour_code.push_back('#');//change to something later
			colour_count.push_back(0);
		}
	}
	colour_code.push_back('\0');
}


KmerColourSummary::~KmerColourSummary(){
//	delete[] colour_code;
//	delete[] colour_count;
	for (int i = 0; i < max_colour; ++i) {
		delete[] matrix[i];
	}
}



void KmerColourSummary::summary_stat(std::string &report){

		char temp_char[10000];
		int no_of_colour_swap = 0;
		int no_of_kmer_swap = 0;

		int i = 0;
		int last_colour = colour_count[i];
		KmerColour last_kmer_colour = kmer_colour[i];
		do{
			last_kmer_colour = kmer_colour[i];
			last_colour = colour_count[i];
			i++;
		}while(last_kmer_colour >= kErrorCode);

	//	int delta_count[colour_len];
		std::string delta_string;
		int delta_colour = 0;
		for (int i = 0; i < length; ++i) {
			int temp_count = colour_count[i];
			KmerColour temp_kmer = kmer_colour[i];
			colour_counter[temp_count]++;

			if(temp_kmer < kErrorCode){
				if (last_colour != temp_count) {
					no_of_colour_swap++;
	//				printf("%d_%s:%s\n",delta_count[i], std::to_string(delta_count[i]).data(), delta_string.data());
	//				itoa(t, delta_count, 10);
				}
				if(last_kmer_colour != temp_kmer){

					no_of_kmer_swap++;
				}

//				if (temp_count != exclude) {
//					if (last_colour_exclude != temp_count
//							&& last_colour_exclude != exclude) {
//						no_of_swap_exclude_n++;
//
//						colour_delta = temp_count-last_colour_exclude;
//					}
//					else{
//						colour_delta = 0;
//					}
//
//	//				delta_count[i] = KmerColourUtil::number_of_colour_s(delta_kmer);
//
//	//				delta_colour_pattern[colour_delta]++;
//	//				delta_string += std::to_string(colour_delta);
//
//					last_colour_exclude = temp_count;
//				}
				delta_colour = temp_count-last_colour;
				delta_colour_pattern[delta_colour]++;
				delta_string += std::to_string(delta_colour);
				last_colour = temp_count;

				KmerColour delta_kmer = last_kmer_colour ^ temp_kmer;
				delta_kmer_colour_pattern[delta_kmer]++;
				last_kmer_colour = temp_kmer;

				auto pattern = std::make_pair (delta_kmer, delta_colour);
				delta_pattern[pattern]++;
			}

		}
	//	char* delta_string = reinterpret_cast<char*>(delta_count);

	//	printf("SWAP:%d %d\n",no_of_swap, no_of_swap_exclude_n);

	//	typedef std::unordered_map<int, int> MapIntInt;
	//	for (MapIntInt::iterator it = counter.begin(); it != counter.end(); ++it ) {
	//		printf("%d:%d\n",it->first, it->second);
	//	}

		std::vector<int> keys;
		double mean = 0;
		double mean_length = 0;

		double max_percentage = 0;
		for (int i = 1; i <= max_colour; ++i) {
			if(colour_counter.count(i)){
				keys.push_back(i);
				auto value = colour_counter.at(i);
				mean += i*value;
				mean_length += value;
				double percentage = 100 * value/ (double) valid_length;
				sprintf(temp_char, "ColourCount:%d  Count:%d/%d Percentage:%.3f%% \n", i, value,
									valid_length, percentage);
				if(i==1 && percentage > 50){
					keep = false;
					reason+= "Single colour > 50%: " + std::to_string(percentage) +" \n";
				}
				else if(percentage > max_percentage){
					max_percentage = percentage;
				}
			}
			else{
				sprintf(temp_char, "ColourCount:%d  Count:%d/%d Percentage:0%% \n", i, 0,
										valid_length);
			}
			report.append(temp_char);

		}
		if(max_percentage<0.5){
			keep = false;
			reason.append("Don't have ?good pattern?, max_percentage < 50%\n");
		}




//		for (auto &it : counter){
//			if (it.first != 0 && it.first < kErrorCode) {
//				keys.push_back(it.first);
//				mean += it.second*it.first;
//				mean_length += it.second;
//
//			sprintf(temp_char, "ColourCount:%d  Count:%d/%d Percentage:%.3f%% \n", it.first, it.second,
//					valid_length, 100 * it.second / (double) valid_length);
//	//		printf(temp_char);
//			report.append(temp_char);
//			}
//
//		}
		mean /= valid_length;
		coverage = (mean / (double) max_colour);
		mean_length /= (keys.size() * valid_length); //TODO: actually, it's not quite correct, 1,3 == 2,3 here!!
		report += "Total length: "+ std::to_string(valid_length) + " (" + std::to_string(length)
				+ ")\n";
		report += "Number of colour swap:" + std::to_string(no_of_colour_swap) +"\n";
		report += "Number of  kmer  swap:" + std::to_string(no_of_kmer_swap) +"\n";

		report += "Swap Pattern: " + delta_string + "\n";
		report += "Mean colour count per site :" + std::to_string(mean) + "/"
				+ std::to_string(max_colour) + "\n";
		report += "Mean length per colour count (delete later):" +std::to_string(mean_length) + "\n";

		if(coverage< 0.5 ){
			keep = false;
			reason += "coverage < 0.5: "+ std::to_string(coverage) +"\n";
		}
		for(auto d : delta_colour_pattern){
			sprintf(temp_char, "delta colour:%d  Count:%d\n", d.first, d.second);
			report.append(temp_char);
		}

		for(auto k : delta_kmer_colour_pattern){
			sprintf(temp_char, "delta kmer:%d  Count:%d\n", k.first, k.second);
			report.append(temp_char);
		}

		for(auto p : delta_pattern){
			sprintf(temp_char, "delta pattern:%d %d Count:%d\n", p.first.first, p.first.second, p.second);
						report.append(temp_char);

		}

		for(auto it : colour_counter){
			double p = it.second/ (double) valid_length;
			if(!keep && p>0.6 && it.first>1){
				keep = true;
			reason += "Overwrite previous decision, colour_count >0.6: Colour:"
					+ std::to_string(it.first) +" = " +std::to_string(p) + "%\n";

			}
		}
	//	printf("\nSIZE: %zu\t%zu\n", report.size(), report.max_size());
	//	report.append("Mean colour count per site:");


	//	printf("%s\n", report.data());
		report += "==========\n";
		if(keep){
			report +="Keep!!!\n";
		}
		else{
			report += "Remove!!!\n";

		}
		report +="Reason:" + reason;
//		char* keep0 = keep ? "Keep!!!!\n" : "Remove!!!!\nReason:"+reason.data();
//
//		report.append(keep0);
//		report += " Reason: +"\n";

}


void KmerColourSummary::summary_colour_code(std::string &report) {
	report.append("ColourCode: ").append(colour_code.data()).append("\n");
}

void KmerColourSummary::summary_colour_count(std::string &report) {
	report.append("ColourCount:");
	for (auto t : colour_count) {
		report += std::to_string(t);
	}
	report.append("\n");
}


int KmerColourSummary::colour_table(std::string &report){

//	report.clear();
//	report.


	for (int i = 0; i < length; ++i) {

		KmerColour temp = kmer_colour[i];
//		if(temp> kErrorCode){
//			for (int j = 0; j < max_colour; ++j) {
//				matrix[j][i] = '#';
//			}
//		}
//		else{
			for (int j = 0; j < max_colour; ++j) {
				if(temp & 1){
					matrix[j][i] = '*';
				}
				else{
					matrix[j][i] = '.';
				}
				temp = temp >> 1;
			}
//		}
	}
	for (int j = 0; j < max_colour; ++j) {
		for (int i = 0; i < length; ++i) {
//			printf("%c",matrix[j][i]);
			report += matrix[j][i];
		}
		report += "\n";
//		printf("\n");
	}

	return 0;


}
/* TODO:??
 * keep or remove??
ColourCount:1  Count:0/58 Percentage:0%
ColourCount:2  Count:58/58 Percentage:100.000%
ColourCount:3  Count:0/58 Percentage:0%
ColourCount:4  Count:0/58 Percentage:0%
ColourCount:5  Count:0/58 Percentage:0%
ColourCount:6  Count:0/58 Percentage:0%

ColourCount:1  Count:11/106 Percentage:10.377%
ColourCount:2  Count:95/106 Percentage:89.623%
ColourCount:3  Count:0/106 Percentage:0%
ColourCount:4  Count:0/106 Percentage:0%
ColourCount:5  Count:0/106 Percentage:0%
ColourCount:6  Count:0/106 Percentage:0%

ColourCount:1  Count:7/806 Percentage:0.868%
ColourCount:2  Count:297/806 Percentage:36.849%
ColourCount:3  Count:372/806 Percentage:46.154%
ColourCount:4  Count:130/806 Percentage:16.129%
ColourCount:5  Count:0/806 Percentage:0%
ColourCount:6  Count:0/806 Percentage:0%
Total length: 806 (808)

ColourCount:1  Count:21/517 Percentage:4.062%
ColourCount:2  Count:37/517 Percentage:7.157%
ColourCount:3  Count:421/517 Percentage:81.431%
ColourCount:4  Count:38/517 Percentage:7.350%
ColourCount:5  Count:0/517 Percentage:0%
ColourCount:6  Count:0/517 Percentage:0%


*/
