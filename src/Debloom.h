#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <cmath> // for log2f
#include <algorithm> // for max
#include <unistd.h> // for truncate

#ifndef DEBLOOM_H
#define DEBLOOM_H

#include "Bank.h"
#include "Bloom.h"
#include "Kmer.h"
#include "Hash16.h"
#include "Utils.h"
#include "Traversal.h"
#include "Memory.h"

using namespace std;


typedef ListSet FPSet; // list-based 
//typedef HashSet FPSet; // hash-based

// GUS: see comment in Debloom.cpp, where false_positive is been declared
extern Set *false_positives;

int debloom(int order, int max_memory);
void end_debloom_partition(bool last_partition);

Set *dummy_false_positives();
Set *load_false_positives();
Set *load_false_positives_cascading4();

void print_size_summary(FPSet *fp);
void print_size_summary(FPSetCascading4 *fp);  

uint64_t get_FPSetCascading4_size (FPSetCascading4 *fp);

class DebloomUtils{
public:
	static int debloom_partition(char *solid_kmer_partition_file, int max_memory);
	static void end_debloom_partition_multi_files(bool last_partition, char* solid_kmer_partition_file);
	static Set *load_false_positives_cascading4_partition(char* solid_kmer_partition_file);
};

#endif
