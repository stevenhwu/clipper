/*
 * Memory.h
 *
 *  Created on: Sep 26, 2014
 *      Author: sw167
 */

#ifndef MEMORY_H_
#define MEMORY_H_

#include <stddef.h>
#include "Bloom.h"
#include "Terminator.h"
#include "Traversal.h"
#include "Debloom.h"
#include <cstdint>
#include <string>


double bits_to_MB(double bits);
double bits_to_KB(double bits);


class MemoryMonitor{

private:
	static int parseLine(char* line);

public:
	static int getValue();
	static void printValue(const std::string message = "");

	static size_t getPeakRSS( );
	static size_t getCurrentRSS( );

};

class MemoryUtils{

public:
	static uint64_t estimate_memory_number_only(char* solid_kmers_colour_file,
			int max_memory);

	static uint64_t estimate_memory(char* solid_kmers_colour_file, int max_memory);
	static uint64_t estimate_memory(BinaryBank* solid_kmers_colour,
			Bloom* bloom, Set* false_positives, BranchingTerminator *terminator);

//	static long int estimate_memory(BinaryBank* solid_kmers_colour,
//				Bloom* bloom, Set* false_positives, Terminator terminator);




};


#endif /* MEMORY_H_ */
