#include "SortingCount.h"
#include "inttypes.h"
#include <sys/resource.h> // for getrlimit()
#if OMP
#include "omp.h"
#endif

#define SINGLE_BAR 1

bool clear_cache = false; // clear file cache from memory (for timing only)

bool hybrid_mode = false;
bool use_hashing = true; // use hashing instead of sorting (better control of memory)
float load_factor = 0.7;
bool use_compressed_reads = false;//true ; // true; // write compressed read file
int optimism = 1; // optimism == 1 mean that we garantee worst case the memory usage, any value above assumes that, on average, a k-mer will be seen 'optimism' times

bool output_histo;


// main k-mer counting function, shared between minia and dsk
// verbose == 0 : stderr progress bar
// verbose >= 1 : print basic status
// verbose >= 2 : print extra partition information
// write_count == True: include kmer count in results file, in that form:
//           - save kmer count for each kmer in the resulting binary file
//           - the very first four bytes of the result file are the kmer length
void sorting_count(Bank *Sequences, char *prefix, int max_memory, int max_disk_space, bool write_count, int verbose, bool skip_binary_conversion)
{
int first =1;
    // create a temp dir from the prefix
    char temp_dir[1024];
    sprintf(temp_dir,"%s_temp",prefix);

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
    } 
    if (max_disk_space == 0) // still 0?
        max_disk_space = 10000; // = default for osx

    // estimate number of iterations
    uint64_t volume = Sequences->estimate_kmers_volume(sizeKmer);
    uint32_t nb_passes = ( volume / max_disk_space ) + 1;

    
    int nb_threads=1;
    
#if OMP
    use_compressed_reads = true;
    nb_threads = 8;
    max_memory /= nb_threads;
    max_memory = max (max_memory,1);
#endif
    
    // temp bugfix: don't use compressed reads for long reads
    if (Sequences->estimate_max_readlen() > 1000000)
        use_compressed_reads = false;
    
    
    uint64_t volume_per_pass;
    uint32_t nb_partitions;


    // loop to lower the number of partitions below the maximum number of simulatenously open files
    do
    {
        volume_per_pass = volume / nb_passes;
        nb_partitions = ( volume_per_pass / max_memory ) + 1;

        // if partitions are hashed instead of sorted, adjust for load factor
        // (as in the worst case, all kmers in the partition are distinct and partition may be slightly bigger due to hash-repartition)
        if (use_hashing)
        {
            nb_partitions = (uint32_t) ceil((float) nb_partitions / load_factor);
            nb_partitions = ((nb_partitions * OAHash::size_entry()) + sizeof(key_type)-1) / sizeof(key_type); // also adjust for hash overhead
            nb_partitions = max((int)(nb_partitions/(optimism+1)), 1);
            if (verbose)
                printf("Updated number of partitions for hash-based k-mer counting: %d\n",nb_partitions);
        }

        // round nb_partitions to mulitple of nthreads, for better perf
        //  nb_partitions = ((nb_partitions + nb_threads - 1) / nb_threads) * nb_threads;

        if (verbose)
            printf("Estimate of number of partitions: %d, number of passes: %d\n",nb_partitions, nb_passes);
        
        // get max number of open files
        struct rlimit lim;
        int max_open_files = 1000;
        int err = getrlimit(RLIMIT_NOFILE, &lim);
        if (err == 0)
            max_open_files = lim.rlim_cur / 2;
        if (nb_partitions >= max_open_files)
        {
            if (verbose)
                printf("Number of partitions higher than max. number of open files (%d), need to increase the number of passes\n", max_open_files);
            nb_passes++;
        }
        else
            break;
    }
    while (1);

 // volume / (sizeof(kmer_type)*4)   is approx size of read file stored in binary, read nb_passes -1 times
    uint64_t total_IO =   volume * 2LL * 1024LL*1024LL   ;// in bytes  +   nb_passes * ( volume / (sizeof(kmer_type)*4) )    ; // in bytes
    uint64_t temp_IO = 0;
    //if (nb_passes==1) use_compressed_reads=false;
    BinaryBankConcurrent * redundant_partitions_file[nb_partitions];
    BinaryBankConcurrent * redundant_partitions_colour_file[nb_partitions];
    char redundant_filename[nb_partitions][256];
    char redundant_colour_filename[nb_partitions][256];
    kmer_type kmer;
    int max_read_length = KMERSBUFFER_MAX_READLEN;
    kmer_type * kmer_table_seq = (kmer_type * ) malloc(sizeof(kmer_type)*max_read_length); ;
    KmerColour * kmer_table_col = (KmerColour * ) malloc(sizeof(kmer_type)*max_read_length); ;

    fprintf(stderr,"Sequentially counting ~%llu MB of kmers with %d partition(s) and %d passes using %d thread(s), ~%d MB of memory and ~%d MB of disk space\n", (unsigned long long)volume, nb_partitions,nb_passes, nb_threads, max_memory * nb_threads, max_disk_space);

    STARTWALL(count);

    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    BinaryBankConcurrent * SolidKmers = new BinaryBankConcurrent(return_file_name(solid_kmers_file),sizeof(kmer),true,nb_threads);
	BinaryBankConcurrent * solid_kmers_colour = new BinaryBankConcurrent(
			return_file_name(solid_kmers_colour_file),
			kSizeOfKmerType + kSizeOfKmerColour , true, nb_threads);
//			kSizeOfKmerType , true, nb_threads);

    if (write_count)
    {
        // write k-mer nbits as the first 4 bytes; and actual k-mer size as the next 4 bits
        uint32_t kmer_nbits = sizeof(kmer) * 8;
        SolidKmers->write_buffered(&kmer_nbits, 4,0);
        SolidKmers->write_buffered(&sizeKmer, 4,0);
        SolidKmers->flush(0);
    }
    int64_t estimated_NbReads = Sequences->estimate_nb_reads();
    char * rseq;
    int readlen;
    KmerColour read_colour;

    int64_t NbSolid = 0;
    int64_t * NbSolid_omp = (int64_t  *) calloc(nb_threads,sizeof(int64_t));
    //long total_kmers_per_partition[nb_partitions]; //guillaume probably commented it because updating this variable would require synchronization
    long distinct_kmers_per_partition[nb_partitions];
    uint64_t  * histo_count = (uint64_t  *) calloc(10001,sizeof(uint64_t));


#if OMP
    uint64_t  **  histo_count_omp = (uint64_t  **) calloc(nb_threads,sizeof(uint64_t *));
    for(int ii=0;ii<nb_threads;ii++)
    {
        histo_count_omp[ii]= (uint64_t  *) calloc(10001,sizeof(uint64_t));
    }
#endif
    

    
   
    //start by the conversion of the file to binary format
    BinaryReads *  binread = NULL;
    if (skip_binary_conversion)
    {
        binread = new BinaryReads(return_file_name(binary_read_file),false);
        binread->close();
    }

    if(use_compressed_reads && (!skip_binary_conversion))
    {
        char * pt_begin;
        int idx =0 ;
        int64_t NbRead = 0;
        Progress progress_conversion;
       // progress_conversion.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        progress_conversion.init(estimated_NbReads,"First step: Converting input file into Binary format");
        
        binread = new BinaryReads(return_file_name(binary_read_file),true);

        Sequences->rewind_all();

        while(1)
        {
            if(! Sequences->get_next_seq_colour(&rseq, &readlen, &read_colour)) break; // read  original fasta file


            if(readlen > max_read_length) // realloc kmer_table_seq if needed
            {
                max_read_length = 2*readlen;
                kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
                kmer_table_col = (KmerColour * ) realloc(kmer_table_col, kSizeOfKmerColour*max_read_length);
            }

            pt_begin = rseq;
            
            //should be ok
            while (pt_begin < (rseq+ readlen))
            {
                idx=0; // start a new read

                //skips NN
                while (*pt_begin =='N' && pt_begin < (rseq+ readlen))
                {
                    pt_begin ++;
                }
                // goes to next N or end of seq
                while ( (pt_begin[idx] !='N') &&  ((pt_begin +idx) < (rseq+ readlen))  )
                {
                    idx++;
                }
                //we have a seq beginning at  pt_begin of size idx  ,without any N, will be treated as a read:
                binread->write_read(pt_begin,idx);
                pt_begin += idx;

            }
            
            // binread->write_read(rseq,readlen);
            
            
            NbRead++;
            if ((NbRead%10000)==0)
            {
                progress_conversion.inc(10000);
            }
        }
        progress_conversion.finish();
        binread->close();

    }
    ///fin conversion

    if (clear_cache)
    {
#ifdef OSX
        system("purge");
#else
        system("echo 3 > /proc/sys/vm/drop_caches");
#endif
    }
    
    
    
#if SINGLE_BAR
    Progress progress;
    char message[1000];
    sprintf(message,"Counting kmers");
    progress.timer_mode=1;
    if (verbose == 0 )
        progress.init(total_IO,message);
#endif
    
    
    // nb_passes = how many times we will traverse the whole reads file (has an influence on temp disk space)
    for (uint32_t current_pass = 0; current_pass < nb_passes; current_pass ++)
    {
        
        if(use_compressed_reads ) //open binary reads for reading
            binread->open(false);
        
        STARTWALL(debpass);
        STARTWALL(debw);

        for (uint32_t p=0;p<nb_partitions;p++)
        {
            sprintf(redundant_filename[p],"%s/partition%d.redundant_kmers",temp_dir,p);
            sprintf(redundant_colour_filename[p],"%s/partition%d.redundant_kmers_colour",temp_dir,p);
            redundant_partitions_file[p] =  new BinaryBankConcurrent (redundant_filename[p],sizeof(kmer_type),true, nb_threads);
            redundant_partitions_colour_file[p] =  new BinaryBankConcurrent (redundant_colour_filename[p], kSizeOfKmerColour, true, nb_threads);
            distinct_kmers_per_partition[p]=0;
        }

        // partitioning redundant kmers
        
        Sequences->rewind_all();
#if !SINGLE_BAR
        Progress progress;
        progress.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        char message[1000];
        sprintf(message,"Pass %d/%d, Step 1: partitioning",current_pass+1,nb_passes);
        if (verbose == 0 )
            progress.init(estimated_NbReads,message);
#endif
     

        
        //current_pass> 0 &&
#if OMP
#pragma omp parallel if(use_compressed_reads)  num_threads(nb_threads)
#endif
        {
            int64_t  nbkmers_written =0;
            int tid =0;
            int64_t NbRead = 0;
            int64_t nread =0;
            int64_t tempread =0;
#if OMP

            tid = omp_get_thread_num();
#endif
            int nreads_in_buffer= 1000;
            KmersBuffer * kbuff =NULL;
            if(use_compressed_reads)
            {
                kbuff = new KmersBuffer (binread, 1000000,  nreads_in_buffer); //buffer size (in nb of kmers), seq per task // the buffer is per thread
                kbuff->binary_read_file = binread->binary_read_file;
            }

            kmer_type * kmer_table ;
            while(1)
            {

                //read the fasta file
                if(use_compressed_reads) // && current_pass>0
                {
                    nread = kbuff->readkmers();
                    if(nread < 0) break;
                    NbRead+= nread;
                    tempread+= nread;
                }
                else
                {
                    if(! Sequences->get_next_seq_colour(&rseq, &readlen, &read_colour)) break; // read  original fasta file
                    if(readlen > max_read_length) // realloc kmer_table_seq if needed
                    {
                        max_read_length = 2*readlen;
                        kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
                        kmer_table_col = (KmerColour * ) realloc(kmer_table_col,kSizeOfKmerColour*max_read_length);
                    }
                }

//                if(use_compressed_reads ) //write compressed read file at first pass //&& current_pass==0
//                    binread->write_read(rseq,readlen);

                int i;
                int nbkmers =readlen-sizeKmer+1;

                if( use_compressed_reads) //current_pass >0 &&
                {
                    nbkmers = kbuff->nkmers;
                    kmer_table = kbuff->kmers_buffer;
                   // printf("nb kmers read  %lli \n",nbkmers);
                 //   NbRead+= nreads_in_buffer;
                } 
                else //old fashion
                {
//                    compute_kmer_table_from_one_seq(readlen,rseq,kmer_table_seq, kmer_table_col);
                	compute_kmer_table_from_one_seq_colour(readlen,rseq,kmer_table_seq, read_colour, kmer_table_col);
                    nbkmers =readlen-sizeKmer+1;
                    kmer_table = kmer_table_seq;
                    NbRead++;
                }

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

                    // some hashing to uniformize repartition
                    kmer_type kmer_hash = lkmer ^ (lkmer >> 14);
                    kmer_hash = (~kmer_hash) + (kmer_hash << 18); 
                    kmer_hash = kmer_hash ^ (kmer_hash >> 31);
                    kmer_hash = kmer_hash * 21; 
                    kmer_hash = kmer_hash ^ (kmer_hash >> 11);
                    kmer_hash = kmer_hash + (kmer_hash << 6);
                    kmer_hash = kmer_hash ^ (kmer_hash >> 22);

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
                //NbRead++;
#if SINGLE_BAR
                if(verbose==0)
                {
                if (nb_threads == 1)
                    progress.inc(nbkmers_written * sizeof(kmer_type));
                else
                    progress.inc(nbkmers_written * sizeof(kmer_type),tid);
                }
#endif
             //   if ((NbRead%10000)==0)
                if(tempread> 10000)
                {
                    tempread -= 10000;
                    if (verbose)
                        fprintf (stderr,"%cPass %d/%d, loop through reads to separate (redundant) kmers into partitions, processed %lluM reads out of %lluM",13,current_pass+1,nb_passes,(unsigned long long)(NbRead/1000/1000),(unsigned long long)(estimated_NbReads/1000/1000));
#if !SINGLE_BAR
                    else
                        if (nb_threads == 1)
                            progress.set(NbRead);
                        else
                            progress.inc(10000,tid);
#endif
                }
            } //end while
            printf("kmer/type/colour_sizeof:%d\t%d\t%d\n", sizeof(kmer), sizeof(kmer_type), kSizeOfKmerColour);
            if(use_compressed_reads)
                delete kbuff;
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
        
        if (verbose)fprintf(stderr,"1\n");

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
        for (int p=0;p<nb_partitions;p++)
        {
            kmer_type lkmer;
            KmerColour lkmer_colour;
            
            bool use_hashing_for_this_partition = use_hashing;
//            printf("max mem %i MB  ,   parti size %i MB\n",max_memory,(redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type))/1024LL/1024LL);
            if(hybrid_mode)
            {
              //  printf("max mem %i MB  ,   parti size %i MB\n",max_memory,(redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type))/1024LL/1024LL);
                if(   (redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type)) <  (max_memory*1024LL*1024LL) )
                    use_hashing_for_this_partition = false;
                else
                    use_hashing_for_this_partition = true;
            }
            int tid =0;
#if OMP
            tid = omp_get_thread_num();
#endif
first = 1;
            if (use_hashing_for_this_partition)
            {
                // hash partition and save to solid file
                OAHash hash(max_memory*1024LL*1024LL);
                uint64_t nkmers_read=0;
                
                while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
                {
                	redundant_partitions_colour_file[p]-> read_element_buffered(&lkmer_colour) ;
                    hash.increment(lkmer, lkmer_colour);
if(first){
//	printf("%li\n", lkmer);
//	first = 0;
}
                    nkmers_read++;// Why not inside SINGLE_BAR??
#if SINGLE_BAR
                    if(verbose==0 && nkmers_read==10000)
                    {
                        if (nb_threads == 1)
                            progress.inc(nkmers_read*sizeof(kmer_type));
                        else
                            progress.inc(nkmers_read*sizeof(kmer_type),tid);
                        nkmers_read=0;
                    }
#endif
                }
                
                //single bar
                
first = 1;
                if (verbose >= 2)
                    printf("Pass %d/%d partition %d/%d hash load factor: %0.3f\n",current_pass+1,nb_passes,p+1,nb_partitions,hash.load_factor());
                int counter = 0;
                hash.start_iterator();
                while (hash.next_iterator())
                {
                    uint_abundance_t abundance = hash.iterator->value;
//                    printf("Hash:key:%ui Value:%li Colour:%li\n",hash.iterator->key, hash.iterator->value, hash.iterator->colour);
                    if(output_histo)
                    {
                        uint_abundance_t saturated_abundance;
                        saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
                        histo_count_omp[tid][saturated_abundance]++;
#else
                        printf("histo_count 0 1  2 %i %i %i \n",histo_count[0],histo_count[1],histo_count[2]);
                        
                        histo_count[saturated_abundance]++;
#endif
                    }

                    if (abundance >= nks && abundance <= max_couv)
                    {
                    	counter++;
//if(first){
//	first = 0;
//	printf("fkey: %li\n", hash.iterator->key);
//}
                		printf("K:%lu\t%hu\n",hash.iterator->key, hash.iterator->colour);

                        SolidKmers->write_element_buffered(&(hash.iterator->key),tid);
//                        Save key only, no need for value/colour
                        solid_kmers_colour-> write_buffered(&(hash.iterator->key), kSizeOfKmerType, tid);
//                        solid_kmers_colour->write_element_buffered(&(hash.iterator->key),tid);
                        solid_kmers_colour-> write_buffered(&(hash.iterator->colour), kSizeOfKmerColour, tid);
                        NbSolid_omp[tid]++;
                        if (write_count)
                            SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                        
                    }
                    distinct_kmers_per_partition[p]++;
                }printf("Total:%i\n\n", counter);
            }
            //In test dataset,
//            4729+4653+4642+4592 = 18616 (unique)
//            1252+1192+1165+1164 = 4773 (>min)
            else
            {//START non hash version
                // sort partition and save to solid file
                vector < kmer_type > kmers;
                uint64_t nkmers_read=0;
                
                
                
                while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
                {
                    kmers.push_back (lkmer);
                    nkmers_read++;
#if SINGLE_BAR
                    if(verbose==0 && nkmers_read==10000)
                    {
                        if (nb_threads == 1)
                            progress.inc(nkmers_read*sizeof(kmer_type));
                        else
                            progress.inc(nkmers_read*sizeof(kmer_type),tid);
                        nkmers_read=0;
                    }
#endif
                }
                
                
                sort (kmers.begin (), kmers.end ());
                
                kmer_type previous_kmer = *(kmers.begin ());
                uint_abundance_t abundance = 0;
                for (vector < kmer_type >::iterator it = kmers.begin (); it != kmers.end ();
                     it++)
                {
                    kmer_type current_kmer = *it;
                    
                    if (current_kmer == previous_kmer)
                        abundance++;
                    else
                    {
                        if(output_histo)
                        {
                            uint_abundance_t saturated_abundance;
                            saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
                            histo_count_omp[tid][saturated_abundance]++;
#else
                            histo_count[saturated_abundance]++;
#endif
                            
                        }

                        if (abundance >= nks  && abundance <= max_couv)
                        {
                            NbSolid_omp[tid]++;
                            SolidKmers->write_element_buffered(&previous_kmer,tid);
                            
                            if (write_count)
                                SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                        }
                        abundance = 1;
                        distinct_kmers_per_partition[p]++;
                    }
                    previous_kmer = current_kmer;
                }
                
                //last kmer
                distinct_kmers_per_partition[p]++;
                if(output_histo)
                {
                    uint_abundance_t saturated_abundance;
                    saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
#if OMP
                    histo_count_omp[tid][saturated_abundance]++;
#else
                    histo_count[saturated_abundance]++;
#endif
                    
                }
                if (abundance >= nks && abundance <= max_couv)
                {
                    NbSolid_omp[tid]++;
                    SolidKmers->write_element_buffered(&previous_kmer,tid);
                    
                    if (write_count)
                        SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                    
                }
                
            }//END non-hash verison, else{**}
            
            
            if (verbose >= 1)
                fprintf(stderr,"%cPass %d/%d, loaded and sorted partition %d/%d, found %lld solid kmers so far",13,current_pass+1,nb_passes,p+1,nb_partitions,(long long)(NbSolid_omp[tid]));
            
            if (verbose >= 2)
                printf("\nPass %d/%d partition %d/%d %ld distinct kmers\n",current_pass+1,nb_passes,p+1,nb_partitions,/*total_kmers_per_partition[p],*/distinct_kmers_per_partition[p]);
            
#if !SINGLE_BAR
            if (verbose == 0 && nb_threads==1)
                progress.inc(1);
            else if (verbose == 0 && nb_threads>1)
                progress.inc(1,tid);
#endif
            
            
            
            redundant_partitions_file[p]->close();
			redundant_partitions_colour_file[p]->close();
//            remove(redundant_filename[p]);
            
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

        if (verbose) fprintf(stderr,"\n");

        if (verbose >= 2)
        {
            STOPWALL(debtri,"Reading and sorting partitions");
            STOPWALL(debpass,"Pass total");

        }
       
        if(use_compressed_reads)
            binread->close();
        
        //delete
            for (uint32_t p=0;p<nb_partitions;p++)
            {
                delete redundant_partitions_file[p] ;
            }
	printf("End of pass %d END\n" ,current_pass+1);
    }
    //single bar
#if SINGLE_BAR
    if (verbose == 0 && nb_threads == 1)
        progress.finish();
    else if (verbose == 0 && nb_threads > 1 )
        progress.finish_threaded();
#endif
    
    if(output_histo)
    {
        FILE * histo_file = fopen(return_file_name(histo_file_name),"w");
        for (int cc=1; cc<10001; cc++) {
            fprintf(histo_file,"%i\t%llu\n",cc,(unsigned long long)(histo_count[cc]));
        }
        fclose(histo_file);
    }
    free(histo_count);

    NbSolid = NbSolid_omp[0];
#if OMP
    NbSolid=0;
    for(int ii=0;ii<nb_threads;ii++)
    {
        NbSolid += NbSolid_omp[ii];
    }
#endif

    SolidKmers->close();
    solid_kmers_colour->close();
    printf("\nSaved %lld solid kmers\n",(long long)NbSolid);
//    rmdir(temp_dir);//TODO add it back

    STOPWALL(count,"Counted kmers");
    fprintf(stderr,"\n------------------ Counted kmers and kept those with abundance >=%i,     \n",nks);
} 


