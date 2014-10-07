#ifndef SORTINGCOUNTPARTITIONS_H
#define SORTINGCOUNTPARTITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <sys/stat.h> // for S_IRWXU etc
#include <vector>
#include <sys/statvfs.h> // to determine available disk space
#include <dirent.h> // to clear the temp directory
#include <cstdint>


#include "Bank.h"
//#include "Kmer.h"
#include "Utils.h"
//#include "OAHash.h"




using namespace std;
typedef uint32_t uint_abundance_t;


class SortingCountPartitions {

	static bool clear_cache;// = false; // clear file cache from memory (for timing only)
	static float load_factor;// = 0.7;
	static int optimism;// = 1; // optimism == 1 mean that we garantee worst case the memory usage, any value above assumes that, on average, a k-mer will be seen 'optimism' times


public:

	static void sorting_count_partitions(Bank *Sequences,
			char solid_kmer_partition_file[][Utils::MaxFileNameLength],
			int max_memory, int max_disk_space, int nb_splits,int verbose=0);

	static void splitBinaryFile(Bank* infile,
			char split_fasta_file[][Utils::MaxFileNameLength],
			int nb_splits);


private:
	static void sorting_count_partitions_core(Bank *Sequences, char* solid_kmer_partition_file,
				int max_memory, int max_disk_space, int split_index, int verbose=0);
	static uint32_t estimate_nb_partitions(uint64_t volume, uint32_t &nb_passes, int max_memory, int verbose=0);


};
//extern int optimism;

#endif
