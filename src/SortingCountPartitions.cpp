#include "SortingCountPartitions.h"
//#include "inttypes.h"
#include <sys/resource.h> // for getrlimit()
#include <unistd.h>
#include <dirent.h>
#include <sys/statvfs.h>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>

#if OMP
#include "omp.h"
#endif

#include "KmerColour.h"
#include "Memory.h"


#define SINGLE_BAR 1

bool SortingCountPartitions::clear_cache = false; // clear file cache from memory (for timing only)
float SortingCountPartitions::load_factor = 0.7;
int SortingCountPartitions::optimism = 1; // optimism == 1 mean that we garantee worst case the memory usage, any value above assumes that, on average, a k-mer will be seen 'optimism' times


//FORCE it to be the hybird = false, using_hashing=true, compressed_read=false
//bool hybrid_mode = false;
//bool use_hashing = true; // use hashing instead of sorting (better control of memory)
//bool use_compressed_reads = false;//true;default TRUE // true; // write compressed read file //TODO: change default

void SortingCountPartitions::sorting_count_partitions(Bank *Sequences,
		char solid_kmer_partition_file[][Utils::MaxFileNameLength], int max_memory,
		int max_disk_space, int nb_splits,  int verbose){

	char split_fasta_file[nb_splits][Utils::MaxFileNameLength];

	splitBinaryFile(Sequences, split_fasta_file, nb_splits);

	for (int p = 0; p < nb_splits; ++p) {
		printf("Process:%s\n",return_file_name(split_fasta_file[p]));
		Bank *read_split = new Bank(return_file_name(split_fasta_file[p]));
		sorting_count_partitions_core(read_split, solid_kmer_partition_file[p], max_memory, max_disk_space, p);
		read_split->close();
		delete read_split;
	}
}

// main k-mer counting function, shared between minia and dsk
// verbose == 0 : stderr progress bar
// verbose >= 1 : print basic status
// verbose >= 2 : print extra partition information
// write_count == True: include kmer count in results file, in that form:
//           - save kmer count for each kmer in the resulting binary file
//           - the very first four bytes of the result file are the kmer length
void SortingCountPartitions::sorting_count_partitions_core(Bank *Sequences,
		char* solid_kmer_partition_file,
		int max_memory, int max_disk_space, int split_index, int verbose ) {

	//verbose=1;
	//verbose=2;

    // create a temp dir from the prefix
    char temp_dir[1024];
    sprintf(temp_dir,"%s_temp", Utils::outfile_prefix);

    // clear the temp folder (needs to be done before estimating disk space)
    DIR*            dp;
    struct dirent*  ep;
    char            p_buf[512] = {0};
    dp = opendir(temp_dir);
    while ( (dp != NULL) && ((ep = readdir(dp)) != NULL)) {
        sprintf(p_buf, "%s/%s", temp_dir, ep->d_name);
        remove(p_buf);
    }
    if(dp != NULL)
        closedir(dp);

    if (max_disk_space == 0)
    {
        // default max disk space
        struct statvfs buffer ;
        char current_path[1000];
        getcwd(current_path,sizeof(current_path));
        // int ret =
        statvfs(current_path, &buffer);
        uint32_t available = (uint32_t)(((double)buffer.f_bavail * (double)buffer.f_bsize) / 1024.0 / 1024.0);
        printf("Available disk space in %s: %d MB\n",current_path,available); // not working in osx (is that a TODO then?)
        uint32_t input_size = max(1, (int)(( (double)(Sequences->filesizes) ) / 1024.0 / 1024.0));
        max_disk_space = min(available/2, input_size);
        // half of total disk space or input file size, why min?? only make sense if ava<input

        if (max_disk_space == 0) // still 0? for osx
			max_disk_space = 10000; // = default for osx

    } 



    int nb_threads=1;
    
#if OMP
    use_compressed_reads = true;
    nb_threads = 8;
    max_memory /= nb_threads;
    max_memory = max (max_memory,1);
#endif

    // estimate number of iterations
    uint64_t volume = Sequences->estimate_kmers_volume(sizeKmer);
    uint32_t nb_passes = ( volume / max_disk_space ) + 1;
    uint32_t nb_partitions = estimate_nb_partitions(volume, nb_passes, max_disk_space);
 // volume / (sizeof(kmer_type)*4)   is approx size of read file stored in binary, read nb_passes -1 times
//    uint64_t total_IO =   volume * 2LL * 1024LL*1024LL   ;// in bytes  +   nb_passes * ( volume / (sizeof(kmer_type)*4) )    ; // in bytes


    BinaryBankConcurrent * redundant_partitions_file[nb_partitions];
    BinaryBankConcurrent * redundant_partitions_colour_file[nb_partitions];
    char redundant_filename[nb_partitions][256];
    char redundant_colour_filename[nb_partitions][256];
	for (uint32_t p = 0; p < nb_partitions; p++) {
		sprintf(redundant_filename[p], "%s/partition%d.redundant_kmers",temp_dir, p);
		sprintf(redundant_colour_filename[p], "%s/partition%d.redundant_kmers_colour", temp_dir, p);
	}

    int max_read_length = KMERSBUFFER_MAX_READLEN;

    kmer_type * kmer_table_seq = (kmer_type * ) malloc(sizeof(kmer_type)*max_read_length); ;
    KmerColour * kmer_table_col = (KmerColour * ) malloc(sizeof(kmer_type)*max_read_length); ;
    fprintf(stderr,"Sequentially counting ~%llu MB of kmers with %d partition(s) and %d"
    		" passes using %d thread(s), ~%d MB of memory and ~%d MB of disk space\n",
    		(unsigned long long)volume, nb_partitions,nb_passes, nb_threads,
    		max_memory * nb_threads, max_disk_space);

    STARTWALL(count);
    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


//    sprintf(solid_kmer_partition_file, "%s_%d", solid_kmer_partition_file, split_index);
	BinaryBankConcurrent * solid_kmers_colour = new BinaryBankConcurrent(
			return_file_name(solid_kmer_partition_file),
			kSizeOfKmerType + kSizeOfKmerColour , true, nb_threads);
//			kSizeOfKmerType , true, nb_threads);

    int64_t estimated_NbReads = Sequences->estimate_nb_reads();
    char * rseq;
    int readlen;
    KmerColour read_colour;
    int64_t * NbSolid_omp = (int64_t  *) calloc(nb_threads,sizeof(int64_t));
    //long total_kmers_per_partition[nb_partitions]; //guillaume probably commented it because updating this variable would require synchronization
    long distinct_kmers_per_partition[nb_partitions];



#if OMP
    uint64_t  **  histo_count_omp = (uint64_t  **) calloc(nb_threads,sizeof(uint64_t *));
    for(int ii=0;ii<nb_threads;ii++)
    {
        histo_count_omp[ii]= (uint64_t  *) calloc(10001,sizeof(uint64_t));
    }
#endif

    if (clear_cache)
    {
#ifdef OSX
        system("purge");
#else
        system("echo 3 > /proc/sys/vm/drop_caches");
#endif
    }
    
#if !SINGLE_BAR//Default on
    Progress progress;
    char message[1000];
    sprintf(message,"Counting kmers");
    progress.timer_mode=1;
    if (verbose == 0 ){
    	progress.timer_mode=0;
        progress.init(total_IO,message);//temp fix for this counter
    }
#endif
    
    // nb_passes = how many times we will traverse the whole reads file (has an influence on temp disk space)
    for (uint32_t current_pass = 0; current_pass < nb_passes; current_pass ++)
    {

        STARTWALL(debpass);
        STARTWALL(debw);
        Sequences->rewind_all();
        // partitioning redundant kmers
        for (uint32_t p=0;p<nb_partitions;p++)
        {
            redundant_partitions_file[p] =  new BinaryBankConcurrent (redundant_filename[p],sizeof(kmer_type),true, nb_threads);
            redundant_partitions_colour_file[p] =  new BinaryBankConcurrent (redundant_colour_filename[p], kSizeOfKmerColour, true, nb_threads);
            distinct_kmers_per_partition[p]=0;
        }

#if !SINGLE_BAR
        Progress progress;
        progress.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        char message[1000];
        sprintf(message,"Pass %d/%d, Step 1: partitioning",current_pass+1,nb_passes);
        if (verbose == 0 )
            progress.init(estimated_NbReads,message);
#endif
     
#if OMP
#pragma omp parallel if(use_compressed_reads)  num_threads(nb_threads)
#endif
        {//start OMP

#if OMP
            int tid = omp_get_thread_num();
#else
            int tid =0;
#endif

            int64_t  nbkmers_written =0;
            int64_t NbRead = 0;
            kmer_type * kmer_table ;

            while(1)
            {
				if(! Sequences->get_next_seq_colour_in_seq_name(&rseq, &readlen, &read_colour)) break;

				if(readlen > max_read_length) // realloc kmer_table_seq if needed
				{
					max_read_length = 2*readlen;
					kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
					kmer_table_col = (KmerColour * ) realloc(kmer_table_col,kSizeOfKmerColour*max_read_length);
				}


                int i;
                int nbkmers =readlen-sizeKmer+1;


				compute_kmer_table_from_one_seq_colour(readlen,rseq,kmer_table_seq, read_colour, kmer_table_col);
				nbkmers =readlen-sizeKmer+1;
				kmer_table = kmer_table_seq;
				NbRead++;

                nbkmers_written= 0;
                //compute the kmers stored in the buffer kmer_table

                for (i=0; i<nbkmers; i++)
                {
                    kmer_type lkmer;
                    KmerColour lkmer_colour;

                    // kmer = extractKmerFromRead(rseq,i,&graine,&graine_revcomp);
//                    printf("A1:%d_A2:%d\n", &kmer_table[i], &kmer_table_seq[i]); TRUE, why declare kmer_table??
//                    lkmer = kmer_table_seq[i];//kmer_table[i];
                    lkmer = kmer_table[i];
                    lkmer_colour = kmer_table_col[i];

                    kmer_type kmer_hash = AbstractOAHash::static_hashcode(lkmer);

                    // check if this kmer should be included in the current pass
                    if ((kmer_hash % nb_passes  ) != current_pass) 
                        continue;
                    kmer_type reduced_kmer = kmer_hash / nb_passes;

                    int p;// compute in which partition this kmer falls into

#ifdef _ttmath
                    (reduced_kmer % nb_partitions).ToInt(p);
#else
                    p = reduced_kmer % nb_partitions;
#endif

                    nbkmers_written++;
//                    printf("%u\n",lkmer_colour);
                    redundant_partitions_file[p]->write_element_buffered(&lkmer,tid); // save this kmer to the right partition file
                    redundant_partitions_colour_file[p]->write_element_buffered(&lkmer_colour,tid); // save this kmer to the right partition file
					// total_kmers_per_partition[p]++; // guillaume probably commented it because updating this variable would require synchronization

				}

#if !SINGLE_BAR //Default on
				if(verbose==0)
				{
					if (nb_threads == 1)
						progress.inc(nbkmers_written * sizeof(kmer_type));
					else
						progress.inc(nbkmers_written * sizeof(kmer_type),tid);
				}
#endif

#if !SINGLE_BAR

				if (nb_threads == 1)
					progress.set(NbRead);
				else
					progress.inc(10000,tid);
#else
				if (verbose>=2 && NbRead%print_table_frequency ==0)
					fprintf(stderr,
							"%cPass %d/%d, loop through reads to separate (redundant) kmers into partitions, processed %lluM reads out of %lluM",
							13, current_pass + 1, nb_passes,
							(unsigned long long) (NbRead / 1000 / 1000),
							(unsigned long long) (estimated_NbReads / 1000/ 1000));
#endif

            } //end while
        } // end OMP 

#if !SINGLE_BAR
        if (verbose == 0)
        {
            if (nb_threads == 1)
             progress.finish();
            else
              progress.finish_threaded();  // here only one thread
            
            sprintf(message,"Pass %d/%d, Step 2: computing kmer count per partition",current_pass+1,nb_passes);
            progress.init(nb_partitions+1,message);
        }
#endif
        
        if (verbose)fprintf(stderr,"Writing redundant kmers\n");

        if (verbose >= 2)
        {
            STOPWALL(debw,"Writing redundant kmers");
        }
        STARTWALL(debtri);

        // close partitions and open them for reading
		for (uint32_t p=0;p<nb_partitions;p++)
		{
			redundant_partitions_file[p]->close();
			redundant_partitions_file[p]->open(false);
			redundant_partitions_colour_file[p]->close();
			redundant_partitions_colour_file[p]->open(false);
		}


        // for better timing: clear the file cache, since the partitions may still be in memory, that's unfair to low mem machines
        if (clear_cache)
        {
#ifdef OSX
            system("purge");
#else
            system("echo 3 > /proc/sys/vm/drop_caches");
#endif
        }

        //quick and dirty parall with omp, testing
        //todo if we want omp and histo : separate histo_count tab per thread that needs to be merged at the end
        // TODO to guillaume: remove that todo above, because it is done, right?
#if OMP 
        //omp_set_numthreads(2);  //num_threads(2) //if(!output_histo) num_threads(nb_threads)
#pragma omp parallel for private (p)  num_threads(nb_threads)
#endif        

        // load, sort each partition to output solid kmers
        for (uint32_t p=0;p<nb_partitions;p++)
        {

            kmer_type lkmer;
            KmerColour lkmer_colour;

            int tid =0;
#if OMP
            tid = omp_get_thread_num();
#endif
#if !SINGLE_BAR//default #if SINGLE_BAR
	uint64_t nkmers_read=0;
#endif
			// hash partition and save to solid file
			OAHashColour hash(max_memory*1024LL*1024LL);


			while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
			{
				redundant_partitions_colour_file[p]-> read_element_buffered(&lkmer_colour) ;
				hash.increment_colour(lkmer, lkmer_colour);

#if !SINGLE_BAR//default #if SINGLE_BAR
				nkmers_read++;// Why not inside SINGLE_BAR?? MOVED
				if(verbose==0 && nkmers_read==print_table_frequency)
				{
					if (nb_threads == 1)
						progress.inc(nkmers_read*sizeof(kmer_type));
					else
						progress.inc(nkmers_read*sizeof(kmer_type),tid);
					nkmers_read=0;
				}
#endif
			}

			if (verbose >= 2)
				printf("Pass %d/%d partition %d/%d hash load factor: %0.3f\n",current_pass+1,nb_passes,p+1,nb_partitions,hash.load_factor());

			uint32_t local_counter = 0;
			hash.start_iterator();
			while (hash.next_iterator())
			{
				uint_abundance_t abundance = hash.iterator->value;
//				KmerColour colour = hash.iterator->colour;
				if (abundance >= nks && abundance <= max_couv) //&&colour > min_colour_count
				{

					solid_kmers_colour-> write_buffered(&(hash.iterator->key), kSizeOfKmerType, tid);
					solid_kmers_colour-> write_buffered(&(hash.iterator->colour), kSizeOfKmerColour, tid);
					local_counter++;
					NbSolid_omp[tid]++;
				}
				distinct_kmers_per_partition[p]++;
			}

            if (verbose >= 1)
				fprintf(stderr,
						"Pass %d/%d, loaded and sorted partition %d/%d, found %lld solid kmers so far, stored %u here.\n",
						current_pass + 1, nb_passes, p + 1, nb_partitions,
						(long long) (NbSolid_omp[tid]), local_counter );

			if (verbose >= 2)
				printf("Pass %d/%d partition %d/%d %ld distinct kmers\n",
						current_pass + 1, nb_passes, p + 1, nb_partitions,/*total_kmers_per_partition[p],*/
						distinct_kmers_per_partition[p]);
            
#if !SINGLE_BAR
            if (verbose == 0 && nb_threads==1)
                progress.inc(1);
            else if (verbose == 0 && nb_threads>1)
                progress.inc(1,tid);
#endif
            
            redundant_partitions_file[p]->close();
			redundant_partitions_colour_file[p]->close();

        } // end for partitions

#if OMP
        //merge histo
        if(output_histo)
        {
            for (int cc=1; cc<10001; cc++) {
                uint64_t sum_omp = 0;
                for(int ii=0;ii<nb_threads;ii++)
                {
                    sum_omp += histo_count_omp[ii][cc];
                }
                histo_count[cc] = sum_omp;
            }
        }
#endif
        
#if !SINGLE_BAR
        if (verbose == 0 && nb_threads == 1)
            progress.finish();
        else if (verbose == 0 && nb_threads > 1 )
            progress.finish_threaded();
#endif

        if (verbose >= 2)
        {
            STOPWALL(debtri,"Reading and sorting partitions");
            STOPWALL(debpass,"Pass total");
        }
       
		for (uint32_t p=0;p<nb_partitions;p++)
		{
			delete redundant_partitions_file[p];
			delete redundant_partitions_colour_file[p];
		}

    }//end pass

#if !SINGLE_BAR//Default on
    if (verbose == 0 && nb_threads == 1)
        progress.finish();//Needs better way to fix this counter
    else if (verbose == 0 && nb_threads > 1 )
        progress.finish_threaded();
#endif
    

#if OMP
    int64_t NbSolid=0;
    for(int ii=0;ii<nb_threads;ii++)
    {
        NbSolid += NbSolid_omp[ii];
    }
#else
    int64_t NbSolid = NbSolid_omp[0];
#endif


    solid_kmers_colour->close();
    delete solid_kmers_colour;

    free(kmer_table_seq);
    free(kmer_table_col);
    free(NbSolid_omp);

    printf("Saved %lld solid kmers in %s\n",(long long)NbSolid, return_file_name(solid_kmer_partition_file));

    rmdir(temp_dir);//TODO add it back

    STOPWALL(count,"Counted kmers");
    fprintf(stderr,"------------------ Counted kmers and kept those with abundance >=%i,     \n",nks);
} 

void SortingCountPartitions::splitBinaryFile(Bank* infile,
		char split_fasta_file[][Utils::MaxFileNameLength], int nb_splits){
//TODO, check file size before and after
	//TODO, need more names and setup proper naming convention here!!
printf("Split:%d\n", nb_splits);

	FILE *fasta_split[nb_splits];// = new Bank(return_file_name(solid_kmer_partition_file[p]) )
//FILE *solid_kmer_partition;// = NULL;
//FILE *file_false_positive_kmers =NULL;
	for (int p = 0; p < nb_splits; ++p) {
			sprintf(split_fasta_file[p] ,"split_%d.fasta", p);
	}
	for (int p = 0; p < nb_splits; ++p) {
		printf("prefix:%s\n", return_file_name(split_fasta_file[p]));

		fasta_split[p] = fopen(return_file_name(split_fasta_file[p]),"wb");
//		fasta_split = fopen(return_file_name(false_positive_kmers_file),"wb");
//		 = fopen(return_file_name(false_positive_kmers_file),"wb");
	}
	char * rseq;
	int readlen;
	KmerColour read_colour;
	infile->rewind_all();
	srand (time(NULL));
	while (1) {
		if (!infile->get_next_seq_colour(&rseq, &readlen, &read_colour))
			break; // read  original fasta file

		int p = rand() % nb_splits;
//		nbSplits;
//		int p = 0;
		fprintf(fasta_split[p],">%d\n", read_colour);
		fputs(rseq,fasta_split[p]);
		fprintf(fasta_split[p],"\n");


	}
	for (int p = 0; p < nb_splits; ++p) {
//		printf("prefix:%s\n", return_file_name(solid_kmer_partition_file[p]));
//		fasta_split[p].close();
		fclose(fasta_split[p]);
	}


//
//	char * pt_begin;
//	int idx = 0;
//	int64_t NbRead = 0;
//	Progress progress_conversion;
//	// progress_conversion.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
//	progress_conversion.init(estimated_NbReads,
//			"First step: Converting input file into Binary format");
//
//	binread = new BinaryReads(return_file_name(binary_read_file), true);
//
//	infile->rewind_all();
//
//	while (1) {
//		if (!Sequences->get_next_seq_colour(&rseq, &readlen, &read_colour))
//			break; // read  original fasta file
//
//		if (readlen > max_read_length) // realloc kmer_table_seq if needed
//				{
//			max_read_length = 2 * readlen;
//			kmer_table_seq = (kmer_type *) realloc(kmer_table_seq,
//					sizeof(kmer_type) * max_read_length);
//			kmer_table_col = (KmerColour *) realloc(kmer_table_col,
//					kSizeOfKmerColour * max_read_length);
//		}
//
//		pt_begin = rseq;
//
//		//should be ok
//		while (pt_begin < (rseq + readlen)) {
//			idx = 0; // start a new read
//
//			//skips NN
//			while (*pt_begin == 'N' && pt_begin < (rseq + readlen)) {
//				pt_begin++;
//			}
//			// goes to next N or end of seq
//			while ((pt_begin[idx] != 'N')
//					&& ((pt_begin + idx) < (rseq + readlen))) {
//				idx++;
//			}
//			//we have a seq beginning at  pt_begin of size idx  ,without any N, will be treated as a read:
//			binread->write_read(pt_begin, idx);
//			pt_begin += idx;
//			//HOW to get colour into this? This break into sequencs without 'file' information
//		}
//
//		// binread->write_read(rseq,readlen);
//
//		NbRead++;
//
//
//	}
//	progress_conversion.finish();
//	binread->close();
//
//}

}
uint32_t SortingCountPartitions::estimate_nb_partitions(uint64_t volume, uint32_t &nb_passes, int max_memory, int verbose){
    // loop to lower the number of partitions below the maximum number of simulatenously open files


	uint64_t volume_per_pass;
	uint32_t nb_partitions;

	// get max number of open files
	struct rlimit lim;
	uint32_t max_open_files = 1000;
	int err = getrlimit(RLIMIT_NOFILE, &lim);
	if (err == 0)
		max_open_files = lim.rlim_cur / 2;
    do
    {
        volume_per_pass = volume / nb_passes;
        nb_partitions = ( volume_per_pass / max_memory ) + 1;
        if(verbose){
        	printf("Initial estimation. volume:%lu nb_passes:%u Volume_per_pass:%lu nb_partitions:%u\n",
        			volume, nb_passes, volume_per_pass, nb_partitions);
        }
        // if partitions are hashed instead of sorted, adjust for load factor
        // (as in the worst case, all kmers in the partition are distinct and partition may be slightly bigger due to hash-repartition)
		nb_partitions = (uint32_t) ceil((float) nb_partitions / load_factor);
		nb_partitions = ((nb_partitions * OAHashColour::size_entry()) + sizeof(key_type)-1) / sizeof(key_type); // also adjust for hash overhead
		nb_partitions = max((int)(nb_partitions/(optimism+1)), 1);
		if (verbose){
			printf("Updated number of partitions for hash-based k-mer counting: %d\n",nb_partitions);
            printf("Estimate of number of partitions: %d, number of passes: %d\n",nb_partitions, nb_passes);
		}


        if (nb_partitions >= max_open_files)
        {
            if (verbose){
                printf("Number of partitions higher than max. number of open files (%d), need to increase the number of passes\n", max_open_files);
                printf("Change nb_passes:%d %d\n",err, max_open_files);
            }
            nb_passes++;
        }
        else
            break;
    }
    while (1);
    return nb_partitions;
}
