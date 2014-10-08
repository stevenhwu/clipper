/*
 * Memory
 *
 *  Created on: Sep 26, 2014
 *      Author: sw167
 */
#include "Memory.h"
//#include <stddef.h>
//#include <cstdint>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cinttypes>

#include "OAHash.h"
#include "Kmer.h"
#include "KmerColour.h"
#include "Utils.h"
#include "Bank.h"
#include "Set.h"

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif


int MemoryMonitor::parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}


int MemoryMonitor::getValue(){ //Note: this value is in KB!
	sleep(3);
	FILE* file = fopen("/proc/self/status", "r"); //proc/curproc in FreeBSD, and different format
	int result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

void MemoryMonitor::printValue(const std::string message){

#if defined(__FreeBSD__)
printf("MemoryMonitor: in FreeBSD, TODO\n");
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
printf("MemoryMonitor: %d KB (%s) \n",getValue(), message.c_str());
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
printf("MemoryMonitor: in UNIX/MAC, TODO\n");
#else
printf("Nothing??:\n");
#endif

}
/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t MemoryMonitor::getPeakRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
	/* AIX and Solaris ------------------------------------------ */
	struct psinfo psinfo;
	int fd = -1;
	if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
		return (size_t)0L;		/* Can't open? */
	if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
	{
		close( fd );
		return (size_t)0L;		/* Can't read? */
	}
	close( fd );
	return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
	return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
	/* Unknown OS ----------------------------------------------- */
	return (size_t)0L;			/* Unsupported. */
#endif
}





/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t MemoryMonitor::getCurrentRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
	/* OSX ------------------------------------------------------ */
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;		/* Can't access? */
	return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
	return (size_t)0L;			/* Unsupported. */
#endif
}

uint64_t MemoryUtils::estimate_memory_number_only(char* solid_kmers_colour_file,
		int max_memory) {

	BinaryBank *solid_kmers_colour = new BinaryBank(
				return_file_name(solid_kmers_colour_file),
				kSizeOfKmerType + kSizeOfKmerColour, 0);
	off_t nbElements = solid_kmers_colour->nb_elements();


	Bloom* bloom = bloom_create_bloo1_partition((BloomCpt *) NULL,
			solid_kmers_colour_file, false);
	printf("Estimate bloom %lu\n", (uint64_t) (nbElements*NBITS_PER_KMER));
	printf("Actual   bloom %lu\n", bloom->tai);

	DebloomUtils::debloom_partition(solid_kmers_colour_file, max_memory);
	Set* false_positives = DebloomUtils::load_false_positives_cascading4_partition(
			solid_kmers_colour_file);
	double fpRate = powf(0.62, NBITS_PER_KMER);
//	(double)powf((double)0.62, (double)NBITS_PER_KMER))
	int nbFP = (int) ceilf(nbElements*7*fpRate);
	printf("Estimate FP: rate%f\nB2:%lu %lu\n"
				"B3:%lu %lu\n"
				"B4:%lu %lu\n"
				"T4:%lu %lu\n",fpRate,
				(uint64_t) nbFP,
				(uint64_t) ceilf(nbFP * NBITS_PER_KMER),

				(uint64_t) ceilf(nbElements*fpRate),
				(uint64_t) ceilf(nbElements*fpRate * NBITS_PER_KMER),

				(uint64_t) ceilf(nbFP*fpRate),
				(uint64_t) ceilf(nbFP*fpRate * NBITS_PER_KMER),

				(uint64_t) ceilf(nbElements*fpRate*fpRate),
				(uint64_t) ceilf(nbElements*fpRate*fpRate) * FPSet::bits_per_element );



	//Can do better on this
	uint64_t genome_size_zero = 0;
//	BranchingTerminator *terminator = new BranchingTerminatorColour(
//			solid_kmers_colour, genome_size_zero, bloom, false_positives);

	BranchingTerminatorColour terminator  (
				solid_kmers_colour, genome_size_zero, bloom, false_positives);
	uint64_t memory = estimate_memory(solid_kmers_colour, bloom, false_positives, &terminator); //, terminator);

	delete bloom;
	delete false_positives;
	solid_kmers_colour->close();
	delete solid_kmers_colour;
//	delete terminator;
//	printf("Max:%d\nMax FP:%d", max_memory*8*1024*1024, * 4 * powf(0.6,11)))

	printf("End estimating memory\n============================\n");
	return memory;
}

uint64_t MemoryUtils::estimate_memory(char* solid_kmers_colour_file,
		int max_memory) {
	MemoryMonitor::printValue("Estimating1:");
	Bloom* bloom = bloom_create_bloo1_partition((BloomCpt *) NULL,
			solid_kmers_colour_file, false);

	DebloomUtils::debloom_partition(solid_kmers_colour_file, max_memory);

	Set* false_positives = DebloomUtils::load_false_positives_cascading4_partition(
			solid_kmers_colour_file);

//	exit(-1);
	//Can do better on this
//	BinaryBank *solid_kmers_colour = new BinaryBank(
//			return_file_name(solid_kmers_colour_file),
//			kSizeOfKmerType + kSizeOfKmerColour, 0);
	 BinaryBank solid_kmers_colour(
				return_file_name(solid_kmers_colour_file),
				kSizeOfKmerType + kSizeOfKmerColour, 0);
	uint64_t genome_size_zero = 0;


	BranchingTerminatorColour terminator  (
				&solid_kmers_colour, genome_size_zero, bloom, false_positives);

	uint64_t memory = estimate_memory(&solid_kmers_colour, bloom, false_positives, &terminator); //, terminator);
	MemoryMonitor::printValue("Estimating2:");

	delete bloom;
	delete false_positives;
	solid_kmers_colour.close();
//	delete solid_kmers_colour;

//	delete &terminator;
	MemoryMonitor::printValue("Estimating after:");
//	printf("Max:%d\nMax FP:%d", max_memory*8*1024*1024, * 4 * powf(0.6,11)))

	printf("End estimating memory\n============================\n");
	return memory;
}

//estimated_BL1_freesize

uint64_t MemoryUtils::estimate_memory(BinaryBank* solid_kmers_colour,
		Bloom* bloom, Set* false_positives, BranchingTerminator *terminator) {

//		    BinaryBank *solid_kmers_colour = new BinaryBank(return_file_name(solid_kmers_colour_file), kSizeOfKmerType+kSizeOfKmerColour, 0);
//	BranchingTerminator *terminator;
//	uint64_t genome_size_zero = 0;
//	terminator = new BranchingTerminatorColour(solid_kmers_colour,
//					genome_size_zero, bloom, false_positives);
//	MonumentTraversal *traversal = new MonumentTraversal(bloom,false_positives,terminator);

//	    traversal->set_maxlen(maxlen);
//	    traversal->set_max_depth(500);
//	    traversal->set_max_breadth(20);
	//    traversal->SetSolidKmersColour(solid_kmers_colour, max_memory);
	off_t nbElements = solid_kmers_colour->nb_elements();

	printf("Estimate bloom %lu\n", (uint64_t) (nbElements*NBITS_PER_KMER));
	printf("Actual   bloom %lu\n", bloom->tai);

	double fpRate = powf(0.62, NBITS_PER_KMER);
//		printf("nbElement %llu. should be (3354114)\n", nbElements);
//		int nbFP = (int) ceilf(nbElements*0.5);// 1-log(2) = 0.3
	int nbFP = (int) ceilf(nbElements*7*fpRate);
//	nbFP = 1109063;
	printf("Estimate FP: rate%f\nB2:%lu %lu\n"
			"B3:%lu %lu\n"
			"B4:%lu %lu\n"
			"T4:%lu %lu\n",fpRate,
			(uint64_t) nbFP,
			(uint64_t) ceilf(nbFP * NBITS_PER_KMER),

			(uint64_t) ceilf(nbElements*fpRate),
			(uint64_t) ceilf(nbElements*fpRate * NBITS_PER_KMER),

			(uint64_t) ceilf(nbFP*fpRate),
			(uint64_t) ceilf(nbFP*fpRate * NBITS_PER_KMER),

			(uint64_t) ceilf(nbElements*fpRate*fpRate),
			(uint64_t) ceilf(nbElements*fpRate*fpRate) * FPSet::bits_per_element );

	uint64_t fp_total = (uint64_t) ceilf(nbFP * NBITS_PER_KMER + nbElements * fpRate * NBITS_PER_KMER
		+ nbFP * fpRate * NBITS_PER_KMER
		+ nbElements * fpRate * fpRate * FPSet::bits_per_element);

	printf("Estimate FP %lu\n", fp_total);
	printf("Actual   FP %lu\n", false_positives->get_total_memory());

	printf("Estimate terminator??? %lu %lu %lu (how many inster??\n", (uint64_t) ceilf(nbElements*fpRate), (uint64_t) ceilf(nbElements*fpRate*4), (uint64_t) ceilf(nbElements*fpRate * 10) );
	printf("Actual   terminator %lu\n", terminator->get_total_memory());
	uint64_t est_terminator =  (uint64_t) ceilf(nbElements*fpRate*4);
	uint64_t total_estimate = (nbElements*NBITS_PER_KMER) + fp_total + est_terminator;
	printf("\n====================\nMemory usage:\n");
	uint64_t min_solid = (1 * solid_kmers_colour->nb_elements()
			* OAHashColour::size_entry());
//	off_t nbElements = solid_kmers_colour->nb_elements();

	printf("Bloom structure: %" PRId64" (Divided up bloom, recalculate FP&Terminator\n",
			bloom->tai);
	printf("Terminator: %" PRId64"\n", terminator->get_total_memory());
	printf("false_positives: %" PRId64"\n",
			false_positives->get_total_memory());
	double scale =  (double) powf( (double)0.62, (double)NBITS_PER_KMER) ;
	printf("Testing estimaitng FP size:%f %f %f %f\n", (nbElements * scale), (nbElements * 4 * scale), scale, NBITS_PER_KMER);
	printf("Soldi_kmer_colour: %" PRId64" to %" PRId64"\n", min_solid,
			min_solid * 2);
	printf("MonumentTraversal: %d bloom+debloom+terminator + (colour)\n",0);

	uint64_t total_bits = bloom->tai + terminator->get_total_memory()
			+ false_positives->get_total_memory();
	printf("Initial max: %" PRId64" %.3f\n", total_bits, bits_to_KB(total_bits));
	printf("Estimate   : %" PRId64" %.2f\n", total_estimate, bits_to_KB(total_estimate));
	return total_bits;
}






double bits_to_MB(double bits)
{
  return bits / 8LL/1024LL/1024LL;
}
double bits_to_KB(double bits)
{
  return bits / 8LL/1024LL;
}
