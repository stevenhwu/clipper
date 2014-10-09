/*
 * Assembler.h
 *
 *  Created on: Oct 8, 2014
 *      Author: sw167
 */

#ifndef ASSEMBLER_H_
#define ASSEMBLER_H_

#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "Bank.h"
#include "Bloom.h"
#include "Debloom.h"
#include "Kmer.h"
#include "KmerColour.h"
#include "Memory.h"
#include "Set.h"
#include "Terminator.h"
#include "Traversal.h"
#include "Utils.h"

class Assembler{

private:

	static const long long maxlen=10000000;
	const uint kMinContigSize  = 2*sizeKmer+1;

	char *left_traversal;
	char *right_traversal;
	char *side_traversal;

	KmerColour *left_colour_traversal;
	KmerColour *right_colour_traversal;
	KmerColour *side_colour_traversal;

	char *contig;
	KmerColour *contig_colour;

	FILE * file_assembly;
    FILE * file_colour_assembly;


	BinaryBank *solid_kmers_colour;
    Bloom* bloom;
	Set* false_positives;
	BranchingTerminator *terminator;
	MonumentTraversal *traversal;

	int debug = 0;
	//FIXME
//	uint max_colour_count = Reads->nb_files;
	uint max_colour_count = 3;

	//deal with these later
    char *assemble_only_one_region = NULL; // debugging, set to a ASCII kmer to activate, NULL to desactivate
    bool LOAD_BRANCHING_KMERS=false; // debugging
    bool DUMP_BRANCHING_KMERS=false;


public:

	void init();
	void run();

	Assembler(int debug);
	Assembler(char *solid_kmer_partition_file, int max_memory, int debug=0);
	~Assembler();

	void runOLD();
};




#endif /* ASSEMBLER_H_ */
