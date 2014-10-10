#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>

#include <iterator>
#include <iostream>
using namespace std;




#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)

#include "Bank.h"
#include "Hash16.h"
#include "Set.h"
#include "Pool.h"
#include "Bloom.h"
#include "Debloom.h"
#include "Utils.h"
#include "SortingCount.h"
#include "Terminator.h"
#include "Kmer.h"
#include "Traversal.h"
#include "rvalues.h" // for 4bloom
#include "Memory.h"
#include "SortingCountPartitions.h"
#include "Assembler.h"


int max_memory; // the most memory one should alloc at any time, in MB
int order = 0; // deblooming order; 0 = debloom everything; 1 = don't debloom 1-node tips (experimental, untested, shouldn't work);// (made extern int in Traversal.h)
//ALWAYS 0 NOW!!

bool FOUR_BLOOM_VERSION = true;

 // shortcuts to go directly to assembly using serialized bloom and serialized hash
int START_FROM_SOLID_KMERS=0; // if = 0, construct the fasta file of solid kmers, if = 1, start directly from that file
int LOAD_FALSE_POSITIVE_KMERS=0; // if = 0, construct the fasta file of false positive kmers (debloom), if = 1, load that file into the hashtable
int NO_FALSE_POSITIVES_AT_ALL=0; // if = 0, normal behavior, if = 1, don't load false positives (will be a probabilistic de bruijn graph)
int max_disk_space = 0;// let dsk decide

int max_colour_count;

inline void assemble(Bloom* bloo1, Set* false_positives)
{

    //////-------------------------------------------------------------------------------------------
    fprintf (stderr,"______________________________________________________ \n");
    fprintf (stderr,"___________ Assemble from bloom filter _______________ \n");
    fprintf (stderr,"______________________________________________________ \n\n");

    //////-------------------------------------------------------------------------------------------


    long long len_left = 0;
    long long len_right = 0;
    long long contig_len =0;
    long long maxlen=10000000;

    char *left_traversal  = (char *) malloc(maxlen*sizeof(char));
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    KmerColour *left_colour_traversal  = (KmerColour *) malloc(maxlen*sizeof(KmerColour));
    KmerColour *right_colour_traversal = (KmerColour *) malloc(maxlen*sizeof(KmerColour));

    char *contig          = (char *) malloc(2*(maxlen+sizeKmer)*sizeof(char));
    KmerColour *contig_colour   = (KmerColour *) malloc(2*(maxlen+sizeKmer)*sizeof(KmerColour));

    kmer_type kmer;

    long long nbContig =0;
    long long totalnt=0;
    long long max_contig_len=0;
    long long mlenleft=0,mlenright=0;
    int64_t NbBranchingKmer=0;
    char kmer_seq[sizeKmer+1];
    FILE * file_assembly = fopen(return_file_name(assembly_file),"w+");
    FILE * file_colour_assembly = fopen(return_file_name(assembly_colour_file),"w+");
    BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);
    BinaryBank *solid_kmers_colour = new BinaryBank(return_file_name(solid_kmers_colour_file), kSizeOfKmerType+kSizeOfKmerColour, 0);

    char colour_seq[1000];

    STARTWALL(assembly);

    char *assemble_only_one_region = NULL; // debugging, set to a ASCII kmer to activate, NULL to desactivate
    bool LOAD_BRANCHING_KMERS=false; // debugging
    bool DUMP_BRANCHING_KMERS=false;

    BranchingTerminator *terminator;

    if (LOAD_BRANCHING_KMERS)
    {printf("LOA:%d\n",LOAD_BRANCHING_KMERS);
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),false);
        terminator = new BranchingTerminator(BranchingKmers,SolidKmers, bloo1,false_positives);
        BranchingKmers->close();
    }
    else
        terminator = new BranchingTerminator(SolidKmers,genome_size, bloo1,false_positives);

    if (DUMP_BRANCHING_KMERS)
    {printf("DUMP:%d\n",DUMP_BRANCHING_KMERS);
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),true);
        terminator->dump_branching_kmers(BranchingKmers);
        BranchingKmers->close();
    }
    printf("Check boolean:%i\t%i\n", LOAD_BRANCHING_KMERS, DUMP_BRANCHING_KMERS);

#ifdef UNITIG
    SimplePathsTraversal *traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    fprintf (stderr,"_________________Assembling in Unitig mode ..._____________________ \n\n");
#else
    MonumentTraversal *traversal = new MonumentTraversal(bloo1,false_positives,terminator);
#endif
    //RandomBranchingTraversal *traversal = new RandomBranchingTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(maxlen);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    traversal->SetSolidKmersColour(solid_kmers_colour, max_memory);

    while (terminator->next(&kmer))
    {
        // keep looping while a starting kmer is available from this kmer
		// everything will be marked during the traversal()'s
		kmer_type starting_kmer;
		code2seq(kmer,kmer_seq); // convert
//		printf("StartWhile, init Kmer:%li\t%s\n",kmer, kmer_seq);// Varified! kmer's matched seq from the original creation
		while (traversal->find_starting_kmer(kmer,starting_kmer))
//		while (traversal->find_starting_kmer_inside_simple_path(kmer,starting_kmer))
		{
		    code2seq(starting_kmer,kmer_seq); // convert starting kmer to nucleotide seq
		    KmerColour kmer_colour = traversal->GetColour(starting_kmer);

//		    printf("Starting_kmer:%lu %s",starting_kmer, kmer_seq);
            if (assemble_only_one_region != NULL)
            {
                kmer_type dummy;
                starting_kmer = extractKmerFromRead(assemble_only_one_region,0,&kmer,&dummy,false);
            }

            // right extension
//            len_right = traversal->traverse(starting_kmer, right_traversal, 0);
			len_right = traversal->traverse_colour(starting_kmer, right_traversal, right_colour_traversal, 0);
            mlenright= max(len_right,mlenright);
            int debug=1;
            if(debug>1){
            	printf("RightSeq:%lld\t%s\n", len_right, right_traversal);
//            	printf("RightColour:");
//            	for (int i = 0; i < len_right; ++i) {
//            		printf("%u ",right_colour_traversal[i]);
//				}
            	KmerColourUtil::kmer_colour_pattern_string(right_colour_traversal, len_right, colour_seq);
				printf("RightColour:%s\n",colour_seq);
            }

            // left extension, is equivalent to right extension of the revcomp
//            len_left = traversal->traverse(starting_kmer, left_traversal, 1);
            len_left = traversal->traverse_colour(starting_kmer, left_traversal,
            										left_colour_traversal, 1);
            mlenleft= max(len_left,mlenleft);

            // form the contig

//            printf("before Rev:%s\n",left_traversal);
            revcomp_sequence(left_traversal,len_left);
            KmerColourUtil::rev_colour(left_colour_traversal, len_left);

//            printf("after Rev:%s\n",left_traversal);
            strcpy(contig,left_traversal); // contig = revcomp(left_traversal)
	        strcat(contig,kmer_seq);//               + starting_kmer
            strcat(contig,right_traversal);//           + right_traversal
			contig_len=len_left+len_right+sizeKmer;


            int colour_len = 0;
            KmerColour sep_colour = kErrorCode+1;// output with %x, so anything greater than 100;
			colour_len = KmerColourUtil::append_colour(left_colour_traversal, len_left,
					contig_colour, colour_len);
			if(debug){
				KmerColourUtil::append_colour(&sep_colour, 1, contig_colour,
						colour_len);
			}
//            memset(contig_colour+pt_len, (int) kmer_colour, kSizeOfKmerColour*sizeKmer);
//            pt_len += sizeKmer;

			KmerColourUtil::append_colour(&kmer_colour, 1, contig_colour,
					colour_len);

			if(debug){
				KmerColourUtil::append_colour(&sep_colour, 1, contig_colour,
										colour_len);
			}

//            memcpy(contig_colour+colour_len, right_colour_traversal, len_right);
//			colour_len += len_right;
			KmerColourUtil::append_colour(right_colour_traversal, len_right,
					contig_colour, colour_len);


            if(debug>1){
            	printf("LeftSeq:%lld\t%s\n", len_left, left_traversal);
//            	printf("LeftColour:");
//				for (int i = 0; i < len_left; ++i) {
//					printf("%u ",left_colour_traversal[i]);
//				}
//				printf("\n");
            	KmerColourUtil::kmer_colour_pattern_string(left_colour_traversal, len_left, colour_seq);
				printf("LeftColour:%s\n",colour_seq);
				printf("Kmer:%s\n",kmer_seq);
				printf("KmerColour:%u\n",kmer_colour);
				printf("Contig:%lld\t%s\n",contig_len ,contig);
//				printf("Colour:");
//				for (int i = 0; i < pt_len; ++i) {
//					printf("%x", contig_colour[i]);
//				}

//				printf("Colour:%d\t%s\n\n",pt_len+len_right ,contig_colour);


            }

            std::string report("==========Summary==========\n");
//			KmerColourUtil::summary(report, contig_colour, colour_len);
//			KmerColourUtil::colour_table(report, contig_colour, colour_len, max_colour_count);
//			printf("%s", report.data());

			KmerColourSummary kcs(contig_colour, colour_len, max_colour_count);
			kcs.summary_colour_code(report);
			kcs.summary_colour_count(report);
			kcs.summary_stat(report);
			kcs.colour_table(report);
//			printf("%s", report.data());
//			delete &kcs;
//			delete &kcs;
//			printf("================END======================\n\n\n");
			// save the contig
            if(contig_len >= MIN_CONTIG_SIZE)//TODO: add colour info here
            {
                max_contig_len = max(max_contig_len,contig_len);
                fprintf(file_assembly,">%lli__len__%lli \n",nbContig,contig_len);
                fprintf(file_assembly,"%s\n",contig);

                fprintf(file_colour_assembly,">%lli__len__%lli \n",nbContig,contig_len);
				fprintf(file_colour_assembly,"%s\n",contig);
////				fprintf(file_colour_assembly,"%s\n",contig_colour);
//				for (int i = 0; i < colour_len; ++i) {
//					fprintf(file_colour_assembly, "%d", all_colour[i]);
//				}
				fprintf(file_colour_assembly,"%s\n",report.data());
                nbContig++;
                totalnt+=contig_len;
            }
            if (assemble_only_one_region != NULL)
                break;
//            printf("Memory: %zu %zu\n", getPeakRSS(), getCurrentRSS() );
//exit(-1);
        }
//		printf("Done while look is assemble()\n");
//fclose(file_assembly);
//fclose(file_colour_assembly);
//exit(-2);
        NbBranchingKmer++;
        if ((NbBranchingKmer%300)==0) fprintf (stderr,"%cLooping through branching kmer n° %" PRId64 "/ %" PRId64 " total nt   %lli" ,13,NbBranchingKmer,terminator->nb_branching_kmers,totalnt );

        if (nbContig > 0 && assemble_only_one_region != NULL)
            break;

    }


    fprintf (stderr,"\n Total nt assembled  %lli  nbContig %lli\n",totalnt,nbContig);
    fprintf (stderr,"\n Max contig len  %lli (debug: max len left %lli, max len right %lli)\n",max_contig_len,mlenleft,mlenright);

    STOPWALL(assembly,"Assembly");
    MemoryMonitor::printValue("big here??");
    free(left_traversal);
    free(right_traversal);
    free(left_colour_traversal);
    free(right_colour_traversal);
    free(contig);
    free(contig_colour);

    fclose(file_assembly);
    fclose(file_colour_assembly);
    SolidKmers->close();
    solid_kmers_colour->close();
    delete SolidKmers;
    delete solid_kmers_colour;
    delete terminator;
    delete traversal;

    MemoryMonitor::printValue("why not decrease??");
    printf("===========DONE=========EXIT========\n");

//exit(-9);
}


void divided_solid_kmer_colour_into_partitions(int nb_partitions, char solid_kmer_partition_file[][Utils::MaxFileNameLength]){

    printf("Split into loop:%d\n", nb_partitions);
//    int nb_partitions = 3;

//	char redundant_colour_filename[nb_partitions][256];
    BinaryBank *solid_kmers_colour_partition[nb_partitions];
    long distinct_kmers_per_partition[nb_partitions];
	for (int p=0;p<nb_partitions;p++)
	{
//		sprintf(solid_kmer_partition_file[p] ,"%s_partition_%d", solid_kmers_colour_file, p);
		printf("Divide %s into %s\n", solid_kmers_colour_file, solid_kmer_partition_file[p]);
		solid_kmers_colour_partition[p] = new BinaryBank(
				return_file_name(solid_kmer_partition_file[p]),
				kSizeOfKmerType + kSizeOfKmerColour, 1);
		distinct_kmers_per_partition[p]=0;
	}

	BinaryBank *solid_kmers_colour = new BinaryBank(return_file_name(solid_kmers_colour_file), kSizeOfKmerType+kSizeOfKmerColour, 0);

	kmer_type kmer;
	KmerColour kmer_colour;
	solid_kmers_colour->rewind_all();
	int total = 0;
	int total_sub = 0;

	for (int p=0;p<nb_partitions;p++){

    	while( solid_kmers_colour-> read_kmer_colour(&kmer, &kmer_colour)){
    		kmer_type kmer_hash  = AbstractOAHash::static_hashcode(kmer);
    		int index = kmer_hash % nb_partitions;
    		distinct_kmers_per_partition[index]++;
			solid_kmers_colour_partition[index]-> write(&kmer, kSizeOfKmerType);
			solid_kmers_colour_partition[index]-> write(&kmer_colour, kSizeOfKmerColour);

//			Testing mulit subparts method -- FAIL
//    		int temp = kmer_hash % 10;
//    		for (int i = 0; i < nb_partitions; ++i) {
//    			if(i != temp){
//    				int index = i;
//    				distinct_kmers_per_partition[index]++;
//    				solid_kmers_colour_partition[index]-> write(&kmer, kSizeOfKmerType);
//    				solid_kmers_colour_partition[index]-> write(&kmer_colour, kSizeOfKmerColour);
//    				total_sub++;
//    			}
//			}
    		total ++;

    	}

    }
	solid_kmers_colour->close();
	delete solid_kmers_colour;
	int total2 = 0;
	for (int p=0;p<nb_partitions;p++) {
		printf("Size for p%d:%ld\n",p, distinct_kmers_per_partition[p]);
		total2 += distinct_kmers_per_partition[p];
		solid_kmers_colour_partition[p]->close();
		delete solid_kmers_colour_partition[p];

	}
	printf("Total:%d\t%d\t%d\t%d\n",total, total2, total_sub, total/10*9*nb_partitions);
}

int preprocess_arg(int argc, char *argv[]){

    if(argc <  6)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s fasta_file kmer_size min_abundance estimated_genome_size prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_genome_size is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    bool FOUR_BLOOM_VERSION = true;

     // shortcuts to go directly to assembly using serialized bloom and serialized hash
    int START_FROM_SOLID_KMERS=0; // if = 0, construct the fasta file of solid kmers, if = 1, start directly from that file
    int LOAD_FALSE_POSITIVE_KMERS=0; // if = 0, construct the fasta file of false positive kmers (debloom), if = 1, load that file into the hashtable
    int NO_FALSE_POSITIVES_AT_ALL=0; // if = 0, normal behavior, if = 1, don't load false positives (will be a probabilistic de bruijn graph)
    int max_disk_space = 0;// let dsk decide
    for (int n_a = 6; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"--original") == 0)
    	    FOUR_BLOOM_VERSION = false;

        if (strcmp(argv[n_a],"--dont-count")==0)
            START_FROM_SOLID_KMERS = 1;

        if (strcmp(argv[n_a],"--dont-debloom")==0)
            LOAD_FALSE_POSITIVE_KMERS = 1;

        if (strcmp(argv[n_a],"--just-assemble")==0)
        {
            START_FROM_SOLID_KMERS = 1;
            LOAD_FALSE_POSITIVE_KMERS = 1;
        }

        if (strcmp(argv[n_a],"--titus-mode")==0)
            NO_FALSE_POSITIVES_AT_ALL = 1;


        if (strcmp(argv[n_a],"-d")==0)
            max_disk_space = atoi(argv[n_a+1]);


        if (strcmp(argv[n_a],"-maxc")==0)
        	max_couv = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"--le-changement")==0)
            {printf("c'est maintenant!\n");exit(0);}
    }


    // kmer size
    sizeKmer=27; // let's make it even for now, because i havnt thought of how to handle palindromes (dont want to stop on them)
    if(argc >=  3)
    {
        sizeKmer = atoi(argv[2]);
        if (sizeKmer%2==0)
        {
            sizeKmer-=1;
            printf("Need odd kmer size to avoid palindromes. I've set kmer size to %d.\n",sizeKmer);
        }
        if (sizeKmer>((int)sizeof(kmer_type)*4))
        {
            printf("Max kmer size on this compiled version is %zu\n", sizeof(kmer_type)*4);
            exit(1);
        }
    }

    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    double lg2 = log(2);

    if (sizeKmer > 128)
    {
        FOUR_BLOOM_VERSION = false;
        printf("Reverted to single Bloom filter implementation for k>128\n");
    }

    /*
     *
     * m = -nln(p) / (ln2)^2
     */
    if (!FOUR_BLOOM_VERSION)
      NBITS_PER_KMER = log(16*sizeKmer*(lg2*lg2))/(lg2*lg2); // needed to process argv[5]
    else
      NBITS_PER_KMER = rvalues[sizeKmer][1];

    // solidity
    min_abundance =NNKS;
    if(argc >=  4)
    {
        min_abundance = atoi(argv[3]);
        if (min_abundance==0) min_abundance=1; // min abundance can't be 0
    }


   if(argc >=  5)
    {
       genome_size  = atoll(argv[4]);
//       int estimated_bloom_size = max( (int)ceilf(log2f(genome_size * NBITS_PER_KMER )), 1);
        uint64_t estimated_bloom_size = (uint64_t) (genome_size * NBITS_PER_KMER);

       uint64_t estimated_nb_FP =  (uint64_t)(genome_size * 4 * powf(0.6,11)); // just indicative

//       max_memory = max( (1LL << estimated_bloom_size)/8LL /1024LL/1024LL, 1LL );
        max_memory =  max((int64_t) estimated_bloom_size/8LL /1024LL/1024LL,1LL);

      printf("estimated values: nbits Bloom %" PRId64 ", nb FP %" PRId64 ", max memory %i MB\n",estimated_bloom_size,estimated_nb_FP,max_memory);

    }

    // output prefix
    if(argc >=  6)
    {
        strcpy(Utils::outfile_prefix,argv[5]);
    }
    return 0;

}

int old_minia(int argc, char *argv[]){

    START_FROM_SOLID_KMERS = 0; //TODO change back to 0 later
    LOAD_FALSE_POSITIVE_KMERS =0; //TODO: change back to 0 later
    NO_FALSE_POSITIVES_AT_ALL = 0;//TODO: change back to 0 later
    STARTWALL(0);
printf("argv[1]:%s\n", argv[1]);
    Bank *Reads = new Bank(argv[1]);

    // counter kmers, write solid kmers to disk
    if (!START_FROM_SOLID_KMERS)
    {
        int verbose = 1; //default 0
        bool write_count = false; //default false
        bool skip_binary_conversion = false; //default false
printf("==========START_FROM_SOLID_KMERS\n");
		sorting_count(Reads, Utils::outfile_prefix, max_memory, max_disk_space, write_count,
				verbose, skip_binary_conversion);
	}

    max_colour_count = Reads->nb_files;

    // debloom, write false positives to disk, insert them into false_positives

    if (! LOAD_FALSE_POSITIVE_KMERS)
    {
    	printf("==========LOAD_FALSE_+ve_KMERS\n");//Don't think we need to add colour to this section, double check LATER
        debloom(order, max_memory);
    }


	int LOAD_BLOOM_FROM_DUMP = 0;
    if(!LOAD_BLOOM_FROM_DUMP){} //TODO later
    Bloom* bloo1 = bloom_create_bloo1((BloomCpt *)NULL, false);

    Set* false_positives;
    FOUR_BLOOM_VERSION = 1;//Save memory
	if (!NO_FALSE_POSITIVES_AT_ALL) { //TODO: deal with this later, use dummy_false_positivies()
		printf("===============!NO_FALSE_+ve_KMERS\n");
		// load false positives from disk into false_positives
		if (!FOUR_BLOOM_VERSION){
			printf("===============NoFourBloom\n");
			false_positives = load_false_positives();//Memory: nb_fp * 64 bits
		}
		else{
			printf("===============FourBloomVersion\n");
			false_positives = load_false_positives_cascading4();
		}
	} else {
		printf("=============else NO_FALSE_+ve_KMERS: titus mode: no FP's\n");
		// titus mode: no FP's
		false_positives = dummy_false_positives();
	}
	assemble(bloo1, false_positives);

    STOPWALL(0,"Total");

    delete Reads;
    return 0;
}




void test_memory_partitions(int nb_partitions, char solid_kmer_partition_file[][Utils::MaxFileNameLength]){
//	long long maxlen=10000000;
//	BinaryBank *solid_kmers_colour_partition[nb_partitions];

	uint64_t memory[nb_partitions];
	for (int p=0;p<nb_partitions;p++) {
		memory[p] = MemoryUtils::estimate_memory(solid_kmer_partition_file[p], max_memory);//, bloo1, false_positives);
	}
	printf("Summary memory usage for all partitions:\n");
	for (int p=0;p<nb_partitions;p++) {
		printf("P%d: %lu %.3fKB %.2fMB\n", p, memory[p], bits_to_KB(memory[p]), bits_to_MB(memory[p]));
	}

}
void test_partitions(int nb_partitions, char solid_kmer_partition_file[][Utils::MaxFileNameLength]){
//	MemoryMonitor::printValue("Partition");
	STARTWALL(assembly);
	for (int p=0;p<nb_partitions;p++) {
//		assemble_partition(solid_kmer_partition_file[p]);
		Assembler assemble(solid_kmer_partition_file[p], max_memory,0);
//		Assembler assemble(4);
			assemble.run();
		MemoryMonitor::printValue("End each partition");
		break;
	}
	STOPWALL(assembly, "Assembly");
	MemoryMonitor::printValue("End Partition");
}


int main(int argc, char *argv[])
{
#ifdef OSX
printf("OSX:\n");
#elif defined(__FreeBSD__)
printf("FreeBSD\n");
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
printf("Linux\n");
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
printf("Some unix/apple");
#else
printf("Nothing??:\n");
#endif

	MemoryMonitor::printValue("Init");

	preprocess_arg(argc, argv);

    max_memory = 500;
    max_disk_space = 1000;

	int nb_splits=1;
	char solid_kmer_partition_file[nb_splits][1024];
	Utils::initilise_partition_names(solid_kmer_partition_file, nb_splits);

//	SortingCountPartitions::sorting_count_partitions(argv[1],
//			solid_kmer_partition_file, max_memory, max_disk_space, nb_splits);

//	test_memory_partitions(nb_splits, solid_kmer_partition_file);
//	test_memory_partitions_using_number_only(nb_splits, solid_kmer_partition_file);

	test_partitions(nb_splits, solid_kmer_partition_file);

//	bool divide_kmers = false;
//		if(divide_kmers){
//			divided_solid_kmer_colour_into_partitions(nb_splits, solid_kmer_partition_file);
//	}
exit(-1);
	old_minia(argc, argv);
	return 0;
}
