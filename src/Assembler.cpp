/*
 * Assembler.cpp
 *
 *  Created on: Oct 8, 2014
 *      Author: sw167
 */


#include "Assembler.h"

Assembler::Assembler(int debug) :
		debug(debug) {}

Assembler::Assembler(char *solid_kmer_partition_file, int max_memory, int debug) :
		debug(debug) {
//
//	init();
//}
//
//void Assembler::init(){
	if(debug){
		printf("debug mode:%d\n", debug);
	}
//	left_traversal  = (char *) malloc(maxlen*sizeof(char));
//	right_traversal = (char *) malloc(maxlen*sizeof(char));
//	left_colour_traversal  = (KmerColour *) malloc(maxlen*sizeof(KmerColour));
//	right_colour_traversal = (KmerColour *) malloc(maxlen*sizeof(KmerColour));
//	if(debug>3){
//		MemoryMonitor::printValue("after both side");
//	}
	side_traversal  = (char *) malloc(maxlen*sizeof(char));
	side_colour_traversal  = (KmerColour *) malloc(maxlen*sizeof(KmerColour));
	if(debug>3){
		MemoryMonitor::printValue("after side only");
	}

	contig          = (char *) malloc(2*(maxlen+sizeKmer)*sizeof(char));
	contig_colour   = (KmerColour *) malloc(2*(maxlen+sizeKmer)*sizeof(KmerColour));
	if(debug>3){
		MemoryMonitor::printValue("after contigs");
	}

    char temp_filename[2014];
    sprintf(temp_filename, "%s.%s", solid_kmer_partition_file, assembly_file);
    file_assembly = fopen(return_file_name(temp_filename),"w+");
    sprintf(temp_filename, "%s.%s", solid_kmer_partition_file, assembly_colour_file);
    file_colour_assembly = fopen(return_file_name(temp_filename),"w+");

	solid_kmers_colour = new BinaryBank(
			return_file_name(solid_kmer_partition_file),
			kSizeOfKmerType + kSizeOfKmerColour, 0);
	if(debug>2){
		MemoryMonitor::printValue("after solidkmer");
	}
    bloom = bloom_create_bloo1_partition((BloomCpt *) NULL,
    		solid_kmer_partition_file, false);
    if(debug>2){
    	MemoryMonitor::printValue("after bloom");
    }
    false_positives = DebloomUtils::create_false_positives_cascading4_partition(
    		solid_kmer_partition_file, max_memory);
    if(debug>2){
    	MemoryMonitor::printValue("after FP");
    }

    uint64_t genome_size  =10000000;
	terminator = new BranchingTerminatorColour(solid_kmers_colour, genome_size,
			bloom, false_positives);
	if(debug>2){
		MemoryMonitor::printValue("after terminator");
	}


    traversal = new MonumentTraversal(bloom,false_positives,terminator);
//    SimplePathsTraversal *traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
//    RandomBranchingTraversal *traversal = new RandomBranchingTraversal(bloo1,false_positives,terminator);

    traversal->set_maxlen(maxlen);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    if(debug>2){
    	MemoryMonitor::printValue("before loop, before solid_kmer_colour in hash");
    }
    traversal->SetSolidKmersColour(solid_kmers_colour, max_memory);
    if(debug>2){
    	MemoryMonitor::printValue("before loop");
    }

}
Assembler::~Assembler(){

//    MemoryMonitor::printValue("Start Destructor!");

//    free(left_traversal);
//    free(right_traversal);
//    free(left_colour_traversal);
//    free(right_colour_traversal);
    free(contig);
    free(contig_colour);
//    MemoryMonitor::printValue("Free some");

    fclose(file_assembly);
    fclose(file_colour_assembly);
//    MemoryMonitor::printValue("Free some more");

    solid_kmers_colour->close();
    delete solid_kmers_colour;
//    MemoryMonitor::printValue("maybe the min here");

    delete bloom;
	delete false_positives;
    delete terminator;
    delete traversal;

//    MemoryMonitor::printValue("Done!! why not decrease??");

}

void Assembler::run(){
//printf("Start run\n");

	//////-------------------------------------------------------------------------------------------
	fprintf (stderr,"______________________________________________________ \n");
	fprintf (stderr,"_______ Assemble partition from bloom filter _________ \n");
	fprintf (stderr,"______________________________________________________ \n\n");

	//////-------------------------------------------------------------------------------------------
//	MemoryMonitor::printValue("init assemble");


	STARTWALL(assembly);
//printf("a\n");
	long long len_left = 0;
	long long len_right = 0;
	long long contig_len = 0;
	long long nbContig = 0;
	long long totalnt = 0;
	long long max_contig_len = 0;
	long long max_lenleft = 0, max_lenright = 0;
//	printf("b\n");

	int64_t NbBranchingKmer = 0;
	char kmer_seq[sizeKmer + 1];
//	char colour_seq[1000000];//crash at 10000000//;
//	http://www.cs.nyu.edu/exact/core/doc/stackOverflow.txt

	kmer_type kmer;
	while (terminator->next(&kmer)) {

		// keep looping while a starting kmer is available from this kmer
		// everything will be marked during the traversal()'s
		kmer_type starting_kmer;
		while (traversal->find_starting_kmer(kmer, starting_kmer)) {
			int colour_len = 0;
			code2seq(starting_kmer,kmer_seq);
			KmerColour sep_colour = kErrorCode+1;// output with %x, so anything greater than 100;

			KmerColour kmer_colour = traversal->GetColour(starting_kmer);

			len_left = traversal->traverse_colour(starting_kmer, side_traversal,
													side_colour_traversal, 1);
			max_lenleft= max(len_left,max_lenleft);

			revcomp_sequence(side_traversal,len_left);
			KmerColourUtil::rev_colour(side_colour_traversal, len_left);
			strcpy(contig,side_traversal); // contig = revcomp(left_traversal)
			colour_len = KmerColourUtil::append_colour(side_colour_traversal, len_left,
					contig_colour, colour_len);
			if (debug > 1) {
				printf("LeftSeq:%lld\t%s\n", len_left, side_traversal);
				char colour_seq[len_left+1];
				KmerColourUtil::kmer_colour_pattern_string(
						side_colour_traversal, len_left, colour_seq);
				printf("LeftColour:%s\n", colour_seq);
				printf("Kmer:%s\n", kmer_seq);

				printf("KmerColour:%u\n", kmer_colour);
			}

			if(debug){
				KmerColourUtil::append_colour(&sep_colour, 1, contig_colour,
						colour_len);
			}


			strcat(contig,kmer_seq);//               + starting_kmer

			len_right = traversal->traverse_colour(starting_kmer, side_traversal, side_colour_traversal, 0);
			max_lenright= max(len_right,max_lenright);

			strcat(contig,side_traversal);//           + right_traversal

//			KmerColour kmer_colour = traversal->GetColour(starting_kmer);
			KmerColourUtil::append_colour(&kmer_colour, 1, contig_colour,
					colour_len);

			if(debug){
				KmerColourUtil::append_colour(&sep_colour, 1, contig_colour,
										colour_len);
			}

			KmerColourUtil::append_colour(side_colour_traversal, len_right,
					contig_colour, colour_len);

			contig_len=len_left+len_right+sizeKmer;

			if (debug > 1) {
				printf("RightSeq:%lld\t%s\n", len_right, side_traversal);
				char colour_seq[len_right+1];
				KmerColourUtil::kmer_colour_pattern_string(
						side_colour_traversal, len_right, colour_seq);
				printf("RightColour:%s\n", colour_seq);

				printf("Contig:%lld\t%s\n", contig_len, contig);
			}

			KmerColourSummary kcs(contig_colour, colour_len, 3);
			kcs.summary(1);
			if(debug>0){
				printf("%s", kcs.get_report());
			}

			// save the contig
			if(contig_len >= kMinContigSize)//TODO: add colour info here
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
				fprintf(file_colour_assembly,"%s\n",kcs.get_report());
				nbContig++;
				totalnt+=contig_len;
			}
			if (assemble_only_one_region != NULL)
				break;

//exit(-1);
		}
	//		printf("Done while look is assemble()\n");
		NbBranchingKmer++;
		if ((NbBranchingKmer% 5000)==0) fprintf (stderr,"%cLooping through branching kmer n° %" PRId64 "/ %" PRId64 " total nt   %lli" ,13,NbBranchingKmer,terminator->nb_branching_kmers,totalnt );

		if (nbContig > 0 && assemble_only_one_region != NULL)
			break;

	}

	fprintf(stderr, "\n Total nt assembled  %lli  nbContig %lli\n", totalnt,
			nbContig);
	fprintf(stderr,
			"\n Max contig len  %lli (debug: max len left %lli, max len right %lli)\n",
			max_contig_len, max_lenleft, max_lenright);

	STOPWALL(assembly, "Assembly");
	printf("===========DONE=========EXIT========\n");
//	exit(-9);
}

/*
void Assembler::runOLD(){
//printf("Start run\n");

	//////-------------------------------------------------------------------------------------------
	fprintf (stderr,"______________________________________________________ \n");
	fprintf (stderr,"_______ Assemble partition from bloom filter _________ \n");
	fprintf (stderr,"______________________________________________________ \n\n");

	//////-------------------------------------------------------------------------------------------
//	MemoryMonitor::printValue("init assemble");


	STARTWALL(assembly);
//printf("a\n");
	long long len_left = 0;
	long long len_right = 0;
	long long contig_len = 0;
	long long nbContig = 0;
	long long totalnt = 0;
	long long max_contig_len = 0;
	long long max_lenleft = 0, max_lenright = 0;
//	printf("b\n");

	int64_t NbBranchingKmer = 0;
	char kmer_seq[sizeKmer + 1];
//	char colour_seq[1000000];//crash at 10000000//;
//	http://www.cs.nyu.edu/exact/core/doc/stackOverflow.txt

	kmer_type kmer;
	while (terminator->next(&kmer)) {

		// keep looping while a starting kmer is available from this kmer
		// everything will be marked during the traversal()'s
		kmer_type starting_kmer;
		while (traversal->find_starting_kmer(kmer, starting_kmer)) {

			if (assemble_only_one_region != NULL)
			{
				kmer_type dummy;
				starting_kmer = extractKmerFromRead(assemble_only_one_region,0,&kmer,&dummy,false);
			}

			// right extension
//            len_right = traversal->traverse(starting_kmer, right_traversal, 0);
			len_right = traversal->traverse_colour(starting_kmer, right_traversal, right_colour_traversal, 0);
			max_lenright= max(len_right,max_lenright);

			if(debug>1){
				printf("RightSeq:%lld\t%s\n", len_right, right_traversal);
            	printf("RightColour:");
            	char colour_seq[len_right];
				KmerColourUtil::kmer_colour_pattern_string(right_colour_traversal, len_right, colour_seq);
				printf("RightColour:%s\n",colour_seq);
			}
//			MemoryMonitor::printValue("debug");

			// left extension, is equivalent to right extension of the revcomp
//            len_left = traversal->traverse(starting_kmer, left_traversal, 1);
			len_left = traversal->traverse_colour(starting_kmer, left_traversal,
													left_colour_traversal, 1);
			max_lenleft= max(len_left,max_lenleft);

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

			KmerColour kmer_colour = traversal->GetColour(starting_kmer);
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
				char colour_seq[len_left];
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
			if(debug>0){
				printf("%s", report.data());
			}
//			delete &kcs;
//			delete &kcs;
//			printf("================END======================\n\n\n");
			// save the contig
			if(contig_len >= kMinContigSize)//TODO: add colour info here
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

//exit(-1);
		}
	//		printf("Done while look is assemble()\n");
		NbBranchingKmer++;
		if ((NbBranchingKmer% 5000)==0) fprintf (stderr,"%cLooping through branching kmer n° %" PRId64 "/ %" PRId64 " total nt   %lli" ,13,NbBranchingKmer,terminator->nb_branching_kmers,totalnt );

		if (nbContig > 0 && assemble_only_one_region != NULL)
			break;

	}

	fprintf(stderr, "\n Total nt assembled  %lli  nbContig %lli\n", totalnt,
			nbContig);
	fprintf(stderr,
			"\n Max contig len  %lli (debug: max len left %lli, max len right %lli)\n",
			max_contig_len, max_lenleft, max_lenright);

	STOPWALL(assembly, "Assembly");
	printf("===========DONE=========EXIT========\n");
//	exit(-9);
}


*/
