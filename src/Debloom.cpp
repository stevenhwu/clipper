#include "Debloom.h"

// GUS: the false positive set can be either FPSet or FPSetCascading4,
// both inheriting from Set. The choice is made by the function called
// to load the false positives: load_false_positives() or
// load_false_positives_cascading4().
Set *false_positives; 

uint64_t nbkmers_solid = 0, b1_size = 0; 

FILE * F_debloom_read;
FILE * F_debloom_write;
uint64_t n_false_positives=0;

Hash16 * hasht1;

//extern const int print_table_frequency;

void end_debloom_partition(bool last_partition)
{

    int value;
    char false_positive_kmer_char[sizeKmer+1];
    FILE *file_false_positive_kmers =NULL;
    kmer_type graine;

    /////////////////////////begin write files 
    rewind (F_debloom_read);
    rewind (F_debloom_write);

	#ifndef MINGW
	ftruncate(fileno(F_debloom_write), 0); //erase previous file 
	#else // tempfix? fileno is not accepted by mingw
	fclose(F_debloom_write);
	F_debloom_write = fopen(return_file_name("debloom2"),"wb+");
	#endif

	BinaryReads* file_false_positive_kmers_binary = NULL;
    if (last_partition)
    {   
        // write false positive kmers to fasta file
        file_false_positive_kmers = fopen(return_file_name(false_positive_kmers_file),"wb");
        char *false_positive_kmers_binary_file = (char *)"false_positive_kmers_binary";

        file_false_positive_kmers_binary = new BinaryReads(
				return_file_name(false_positive_kmers_binary_file), true);
		//TODO: change store FP into binary? reduce space?
    }

    n_false_positives = 0;
    while(fread(&graine, sizeof(graine),1, F_debloom_read)){

        if(hasht1->get(graine,&value)==0) //kmer not present == kmer not solid
        {
            n_false_positives ++;

            if (!fwrite(&graine, sizeof(graine), 1, F_debloom_write))
            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }


            if (last_partition)
            {
                code2seq(graine,false_positive_kmer_char);
                fprintf(file_false_positive_kmers,">fp\n");
                fputs(false_positive_kmer_char,file_false_positive_kmers);
                fprintf(file_false_positive_kmers,"\n");
                file_false_positive_kmers_binary->write_read(false_positive_kmer_char, sizeKmer);
            }
        }
        //else kmer is a true positive, do nothing

    }

    if (last_partition){
        fclose(file_false_positive_kmers);
    	file_false_positive_kmers_binary->close();
    }
} 

int debloom(int order, int max_memory)
{
    // read bloo1 from disk dump
    Bloom *bloo1 = bloom_create_bloo1((BloomCpt *)NULL);

    STARTWALL(pos);

    FILE * debloom_file = fopen(return_file_name("debloom"),"wb+");
    FILE * debloom_file_2 = fopen(return_file_name("debloom2"),"wb+");
    FILE * F_tmp;
    
    F_debloom_read = debloom_file;
    F_debloom_write = debloom_file_2;
	
    BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);
    
    uint64_t cc=0;
    kmer_type new_graine, kmer;
    int nt;
   
    uint64_t NbSolidKmer =0;
    // write all positive extensions in disk file
    while (SolidKmers->read_element(&kmer))
    {

        //8 right extensions   (4F and 4R); left extensions are redundant by revcomplementation
        for(nt=0; nt<4; nt++) 
        {
            int strand;
            for (strand = 0; strand < 2 ; strand++)
            {
                int current_strand = strand;
                new_graine = next_kmer(kmer,nt, &current_strand);

                if(bloo1->contains(new_graine)){   // extension is positive

                    // maybe do more lax deblooming; if it's a dead-end, it's no big deal, don't pass it to the false positive test
                    // what would have been needed if i decided to enable order>0 (but actually this won't happen): 
                    //  - better estimate of structure size in the presence of order>0 deblooming  
                    if (order == 1)  // this case just detects tips
                    {printf("ORDER==1");
                        bool is_linked = false;
                        for(int tip_nt=0; tip_nt<4; tip_nt++) 
                        {
                            int new_strand = current_strand;
                            kmer_type kmer_after_possible_tip = next_kmer(new_graine,tip_nt, &new_strand);
                            if(bloo1->contains(kmer_after_possible_tip))
                            {
                                is_linked = true;
                                break;
                            }
                        }
                        if (!is_linked)
                            continue; // it's a tip, because it's linked to nothing
                    }
    
                    if (order > 1) // general case. should work for order = 1, but i coded an optimized version above
                    { printf("ORDER>1");
                        Frontline frontline( new_graine, current_strand, bloo1, NULL, NULL, NULL);
                        while (frontline.depth < order)
                        {
                            frontline.go_next_depth();
                            if (frontline.size() == 0)
                                break;
                            // don't allow a breadth too large anywqy
                            if (frontline.size()> 10)
                                break;
                        }
                        if (frontline.size() == 0)
                            continue; // it's a deadend
                    }

                    if (!fwrite(&new_graine, sizeof(new_graine), 1, debloom_file))
                    {
                        printf("error: can't fwrite (disk full?)\n");
                        exit(1);
                    }
                    cc++;
                }

            }
        }
        NbSolidKmer++;
        if ((NbSolidKmer%print_table_frequency)==0) fprintf (stderr,"%c Writing positive Bloom Kmers %lld",13,NbSolidKmer);
    }
    nbkmers_solid =  NbSolidKmer; // GUS: it's global now

    fprintf(stderr,"\n%lli kmers written\n",cc);

    STOPWALL(pos,"Write all positive kmers");

    STARTWALL(deb);

    double bl1tai =  (double)bloo1->tai ;
//    printf("Missing var bl1tai:%f\n", bl1tai);
    delete bloo1;

    // now that bloo1 is deleted, initialize hasht1
    int NBITS_HT = max( (int)ceilf(log2f((0.1*max_memory*1024L*1024L)/sizeof(cell_ptr_t))), 1); // set hasht1 cells to occupy 0.1 * [as much mem as poss]
    hasht1 =new Hash16(NBITS_HT); 
    
    ////////////////////////////////////////////////////////////////   --find false positive, with hash table partitioning
    uint64_t max_kmer_per_part = (uint64_t) (0.8*max_memory*1024LL*1024LL /sizeof(cell<kmer_type>));
    //adapter taille ht en fonction
    

    printf("%d partitions will be needed\n",(int)(nbkmers_solid/max_kmer_per_part));

    NbSolidKmer =0;
    int numpart = 0;
    SolidKmers->rewind_all();

    // deblooming:
    // read the list of (non-redundant) solid kmers and load it, in chunks, into a hash table
    // at each pass, check all the positive extensions and keep those which are not indicated, by the current chunk, as solid kmers
    // at the end, only the positive extensions which are not solid are kept
    while (SolidKmers->read_element(&kmer))
    {
        hasht1->add(kmer);

        NbSolidKmer++;
        if ((NbSolidKmer%print_table_frequency)==0) fprintf (stderr,"%cBuild Hash table %lld",13,NbSolidKmer);

        if(hasht1->nb_elem >max_kmer_per_part) //end partition,  find false positives
        {
            fprintf(stderr,"End of debloom partition  %lli / %lld \n",hasht1->nb_elem,max_kmer_per_part);

            end_debloom_partition(false);

            //swap file pointers
            F_tmp = F_debloom_read;
            F_debloom_read = F_debloom_write;
            F_debloom_write = F_tmp;
            /////////end write files

            //reset hash table
            hasht1->empty_all();

            fprintf(stderr,"\n%lli false positives written , partition %i \n",n_false_positives,numpart);

            numpart++;
        } ///end partition


    }
    fprintf(stderr,"Nb kmers stored in the bloom table %lld\n",nbkmers_solid);


    ///////////////////////// last partition, will write all the FP's to the good file

    end_debloom_partition(true); 

    /////////end write files


    fprintf(stderr,"Total nb false positives stored in the Debloom hashtable %lli \n",n_false_positives);

    delete hasht1;


    STOPWALL(deb,"Debloom");
 
    // GUS: will use to output summary later
    b1_size = (uint64_t) bl1tai;
  
    fclose(debloom_file);
    fclose(debloom_file_2);
    SolidKmers->close();


    return 1;

}

uint64_t countFP(Bank *FalsePositives)
{
  char * rseq;
  int readlen;
  uint64_t nbFP = 0;

  while (FalsePositives->get_next_seq(&rseq,&readlen))
    nbFP++;
  
  FalsePositives->rewind_all();
  return nbFP;
}

Set *load_false_positives() 
{
    int64_t NbInsertedKmers = 0;
    char * rseq;
    int readlen;
    kmer_type kmer, graine, graine_revcomp;

    Bank *FalsePositives = new Bank(return_file_name(false_positive_kmers_file));


    // alloc false positives with the just the right estimated size

    uint64_t nbFP = countFP(FalsePositives);

    FPSet *fp = new FPSet(nbFP);
    
    while (FalsePositives->get_next_seq(&rseq,&readlen))
    {
        kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
                
        fp->insert(kmer);

        NbInsertedKmers++;

        if ((NbInsertedKmers%print_table_frequency)==0) fprintf (stderr,(char*)"%cInsert false positive Kmers in hash table %lld",13,NbInsertedKmers);
    }
    fp->finalize(); // always call this when finishing to create a FPSet

    fprintf (stderr,"\nInserted %lld false positive kmers in the hash structure.\n\n",NbInsertedKmers);

    print_size_summary(fp);

    return fp;
}


Set *dummy_false_positives() 
{
    FPSet *fp = new FPSet((uint64_t)1);
    return fp;
}

Set *load_false_positives_cascading4()
{

if (nbkmers_solid==0){
	BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);
	nbkmers_solid = SolidKmers->nb_elements();
//	b1_size = 20239964;
//	uint64_t estimated_BL1_freesize =  (uint64_t)(nbkmers_solid*NBITS_PER_KMER);

}
  int64_t NbInsertedKmers;
  char * rseq;
  int readlen;
  kmer_type kmer, graine, graine_revcomp;

  
  // **** Initialize B2, B3, B4 and T4 ****
  Bank *FalsePositives = new Bank(return_file_name(false_positive_kmers_file));
  uint64_t nbFP = countFP(FalsePositives);
  
  FPSetCascading4 *fp = new FPSetCascading4;
  
  fp->bloom2 = new Bloom((uint64_t)(nbFP * NBITS_PER_KMER));
  fp->bloom2->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  uint64_t estimated_T2_size = max((int)ceilf(nbkmers_solid * (double)powf((double)0.62, (double)NBITS_PER_KMER)), 1);
  uint64_t estimated_T3_size = max((int)ceilf(nbFP          * (double)powf((double)0.62, (double)NBITS_PER_KMER)) ,1);

  fp->bloom3 = new Bloom((uint64_t)(estimated_T2_size * NBITS_PER_KMER));
  fp->bloom3->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  fp->bloom4 = new Bloom((uint64_t)(estimated_T3_size * NBITS_PER_KMER));
  fp->bloom4->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));


  // **** Insert the false positives in B2 ****
  NbInsertedKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
    
    fp->bloom2->add(kmer);
    
    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive B2 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,"%cInsert false positive B2 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  printf("\nInserted %lld (estimated, %lld) kmers in B2.\n", NbInsertedKmers, nbFP);


  //  **** Insert false positives in B3 and write T2 
  int addKmers = 0;
  NbInsertedKmers = 0;
  FILE *T2_file = fopen(return_file_name("t2_kmers"), "w+"); // We will read this file later, when filling T4 
  BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer),0);
  while(SolidKmers->read_element(&kmer))
  {
    if (fp->bloom2->contains(kmer))
    {
      if (!fwrite(&kmer, sizeof(kmer), 1, T2_file))
      {
	printf("error: can't fwrite (disk full?)\n");
	exit(1);
      }

      fp->bloom3->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers% print_table_frequency)==0)
      fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  SolidKmers->close();

  printf("\nInserted %lld (estimated, %llu) kmers in B3.\n", addKmers, estimated_T2_size);

  
  // **** Insert false positives in B4 (we could write T3, but it's not necessary)
  FalsePositives = new Bank(return_file_name(false_positive_kmers_file));
  NbInsertedKmers = 0;
  addKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
    
    if (fp->bloom3->contains(kmer))
    {
      fp->bloom4->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive B4 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,"%cInsert false positive B4 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  printf("\nInserted %lld (estimated, %lld) kmers in B4.\n", addKmers, estimated_T3_size);
  

  // **** Count and insert false positives in T4
  rewind(T2_file);
  addKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
    if (fp->bloom4->contains(kmer))
      addKmers++;

  fp->false_positives = new FPSet(addKmers);
  rewind(T2_file);
  addKmers = 0;
  NbInsertedKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
  {  
    if (fp->bloom4->contains(kmer))
    {
      fp->false_positives->insert(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive T4 %lld",13,NbInsertedKmers);
  }
  fp->false_positives->finalize();
  fprintf (stderr,"%cInsert false positive T4 %lld", 13,NbInsertedKmers);
  fclose(T2_file);

  printf("\nInserted %lld (estimated, %lld) kmers in T4.\n\n", addKmers, (uint64_t)fp->false_positives->capacity());

  uint64_t total_size = get_FPSetCascading4_size(fp);
  fp->set_total_memory(total_size);
  print_size_summary(fp);

  return fp;
}

void print_size_summary(FPSet *fp)
{
  int bits_per_FP_element = FPSet::bits_per_element;

  uint64_t size_B1 = b1_size,
           size_T1 = fp->capacity() * FPSet::bits_per_element;
  double total_size = (double)(size_B1 + size_T1);

  fprintf(stderr,"Size of the Bloom table  : %.2lf MB\n", bits_to_MB(size_B1) );
  fprintf(stderr,"                           %.2lf bits / solid kmer\n", b1_size/(double)(nbkmers_solid) );  
  fprintf(stderr, "Size of the FP table     : %lli FP x %d bits =  %.2lf MB  \n", fp->capacity(), bits_per_FP_element, bits_to_MB((double)(size_T1)) );
  fprintf(stderr,"                                      actual implementation : %.2lf bits / solid kmer\n", size_T1/(double)nbkmers_solid);
  fprintf(stderr,"  assuming list of kmers, i.e. sizeof(kmer_type) bits / FP : %.2lf bits / solid kmer \n\n",(fp->capacity()*sizeof(kmer_type)*8LL)/(double)(nbkmers_solid));
  fprintf(stderr,"      Total %.2lf MB for %lld solid kmers  ==>  %.2lf bits / solid kmer\n\n", bits_to_MB(total_size), nbkmers_solid, total_size / nbkmers_solid);
}

void print_size_summary(FPSetCascading4 *fp)
{
  uint64_t size_B1 = b1_size,
    size_B2 = fp->bloom2->tai,
    size_B3 = fp->bloom3->tai,
    size_B4 = fp->bloom4->tai,
    size_T4 = fp->false_positives->capacity() * FPSet::bits_per_element;
  double total_size = (double)(size_B1 + size_B2 + size_B3 + size_B4 + size_T4);

  fprintf(stderr,"Size of the Bloom table (B1)  : %.2lf MB  %llu\n", bits_to_MB((double)size_B1), size_B1);
  fprintf(stderr,"Size of the Bloom table (B2)  : %.2lf MB  %llu\n", bits_to_MB((double)size_B2), size_B2);
  fprintf(stderr,"Size of the Bloom table (B3)  : %.2lf MB  %llu\n", bits_to_MB((double)size_B3), size_B3);
  fprintf(stderr,"Size of the Bloom table (B4)  : %.2lf MB  %llu\n", bits_to_MB((double)size_B4), size_B4);
  fprintf(stderr,"Size of the FP table (T4)     : %.2lf MB  %llu\n", bits_to_MB((double)size_T4), size_T4);
  fprintf(stderr,"      Total %.2lf MB for %lld solid kmers  ==>  %.2lf bits / solid kmer\n\n", bits_to_MB(total_size), nbkmers_solid, total_size / nbkmers_solid);
}

uint64_t get_FPSetCascading4_size (FPSetCascading4 *fp)
{
  uint64_t size_B1 = 0,//b1_size, size_B1 is NOT in FP
    size_B2 = fp->bloom2->tai,
    size_B3 = fp->bloom3->tai,
    size_B4 = fp->bloom4->tai,
    size_T4 = fp->false_positives->capacity() * FPSet::bits_per_element;
  uint64_t total_size = (size_B1 + size_B2 + size_B3 + size_B4 + size_T4);
  return total_size;
}


int debloom_partition(char *solid_kmer_partition_file, int max_memory)
{//ONLY work for order ==0 now.
    // read bloo1 from disk dump

    Bloom *bloo1 = bloom_create_bloo1_partition((BloomCpt *)NULL, solid_kmer_partition_file, true);
    //Load from dist

    STARTWALL(pos);

    FILE * debloom_file = fopen(return_file_name("debloom"),"wb+");
    FILE * debloom_file_2 = fopen(return_file_name("debloom2"),"wb+");
    FILE * F_tmp;

    F_debloom_read = debloom_file;
    F_debloom_write = debloom_file_2;

    BinaryBank *solid_kmer_colour = new BinaryBank(return_file_name(solid_kmer_partition_file),sizeof(kmer_type)+sizeof(KmerColour),0);

    uint64_t cc=0;
    kmer_type new_graine, kmer;
    KmerColour kmer_colour;
    int nt;

    uint64_t NbSolidKmer =0;
    // write all positive extensions in disk file
    uint64_t cc2=0;

    while (solid_kmer_colour->read_kmer_colour(&kmer, &kmer_colour))
    {
        //8 right extensions   (4F and 4R); left extensions are redundant by revcomplementation
        for(nt=0; nt<4; nt++)
        {
            int strand;
            for (strand = 0; strand < 2 ; strand++)
            {
                int current_strand = strand;
                new_graine = next_kmer(kmer,nt, &current_strand);
                cc2++;
                if(bloo1->contains(new_graine)){   // extension is positive
                    if (!fwrite(&new_graine, sizeof(new_graine), 1, debloom_file))
                    {
                        printf("error: can't fwrite (disk full?)\n");
                        exit(1);
                    }
                    cc++;
                }

            }
        }
        NbSolidKmer++;
        if ((NbSolidKmer%print_table_frequency)==0) fprintf (stderr,"%c Writing positive Bloom Kmers %lld",13,NbSolidKmer);
    }
    nbkmers_solid =  NbSolidKmer; // GUS: it's global now

    fprintf(stderr,"\n%lli kmers written. All possible %lli\n",cc, cc2);

    STOPWALL(pos,"Write all positive kmers");

    STARTWALL(deb);

    double bl1tai =  (double)bloo1->tai ;
//    printf("Missing var bl1tai:%f\n", bl1tai);
    delete bloo1;

    // now that bloo1 is deleted, initialize hasht1
    int NBITS_HT = max( (int)ceilf(log2f((0.1*max_memory*1024L*1024L)/sizeof(cell_ptr_t))), 1); // set hasht1 cells to occupy 0.1 * [as much mem as poss]
    hasht1 =new Hash16(NBITS_HT);

    ////////////////////////////////////////////////////////////////   --find false positive, with hash table partitioning
    uint64_t max_kmer_per_part = (uint64_t) (0.8*max_memory*1024LL*1024LL /sizeof(cell<kmer_type>));
    //adapter taille ht en fonction


    printf("%d partitions will be needed\n",(int)(nbkmers_solid/max_kmer_per_part));

    NbSolidKmer =0;
    int numpart = 0;
    solid_kmer_colour->rewind_all();

    // deblooming:
    // read the list of (non-redundant) solid kmers and load it, in chunks, into a hash table
    // at each pass, check all the positive extensions and keep those which are not indicated, by the current chunk, as solid kmers
    // at the end, only the positive extensions which are not solid are kept
    while (solid_kmer_colour->read_kmer_colour(&kmer, &kmer_colour))
    {

        hasht1->add(kmer);

        NbSolidKmer++;
        if ((NbSolidKmer%print_table_frequency)==0) fprintf (stderr,"%cBuild Hash table %lld",13,NbSolidKmer);

        if(hasht1->nb_elem >max_kmer_per_part) //end partition,  find false positives
        {
            fprintf(stderr,"End of debloom partition  %lli / %lld \n",hasht1->nb_elem,max_kmer_per_part);

            end_debloom_partition_multi_files(false, solid_kmer_partition_file);

            //swap file pointers
            F_tmp = F_debloom_read;
            F_debloom_read = F_debloom_write;
            F_debloom_write = F_tmp;
            /////////end write files

            //reset hash table
            hasht1->empty_all();

            fprintf(stderr,"\n%lli false positives written , partition %i \n",n_false_positives,numpart);

            numpart++;
        } ///end partition


    }
    fprintf(stderr,"Nb kmers stored in the bloom table %lld\n",nbkmers_solid);


    ///////////////////////// last partition, will write all the FP's to the good file

    end_debloom_partition_multi_files(true, solid_kmer_partition_file);


    /////////end write files


    fprintf(stderr,"Total nb false positives stored in the Debloom hashtable %lli \n",n_false_positives);

    delete hasht1;


    STOPWALL(deb,"Debloom");

    // GUS: will use to output summary later
    b1_size = (uint64_t) bl1tai;

    fclose(debloom_file);
    fclose(debloom_file_2);
    solid_kmer_colour->close();


    return 1;

}

void end_debloom_partition_multi_files(bool last_partition, char* solid_kmer_partition_file)
{

    int value;
    char false_positive_kmer_char[sizeKmer+1];
    FILE *file_false_positive_kmers =NULL;
    kmer_type graine;

    /////////////////////////begin write files
    rewind (F_debloom_read);
    rewind (F_debloom_write);

	#ifndef MINGW
	ftruncate(fileno(F_debloom_write), 0); //erase previous file
	#else // tempfix? fileno is not accepted by mingw
	fclose(F_debloom_write);
	F_debloom_write = fopen(return_file_name("debloom2"),"wb+");
	#endif

	BinaryReads* file_false_positive_kmers_binary = NULL;
    if (last_partition)
    {
        // write false positive kmers to fasta file
    	char temp_file[1024];
    	sprintf(temp_file,"%s_%s",solid_kmer_partition_file, false_positive_kmers_file);
        file_false_positive_kmers = fopen(return_file_name(temp_file),"wb");

        sprintf(temp_file,"%s_%s_binary",solid_kmer_partition_file, false_positive_kmers_file);
//        char *false_positive_kmers_binary_file = (char *)"false_positive_kmers_binary";

        file_false_positive_kmers_binary = new BinaryReads(
				return_file_name(temp_file), true);
		//TODO: change store FP into binary? reduce space?
    }

    n_false_positives = 0;
    while(fread(&graine, sizeof(graine),1, F_debloom_read)){

        if(hasht1->get(graine,&value)==0) //kmer not present == kmer not solid
        {
            n_false_positives ++;

            if (!fwrite(&graine, sizeof(graine), 1, F_debloom_write))
            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }


            if (last_partition)
            {
                code2seq(graine,false_positive_kmer_char);
                fprintf(file_false_positive_kmers,">fp\n");
                fputs(false_positive_kmer_char,file_false_positive_kmers);
                fprintf(file_false_positive_kmers,"\n");
                file_false_positive_kmers_binary->write_read(false_positive_kmer_char, sizeKmer);
            }
        }
        //else kmer is a true positive, do nothing

    }

    if (last_partition){
        fclose(file_false_positive_kmers);
    	file_false_positive_kmers_binary->close();
    }
}

Set *load_false_positives_cascading4_partition(char* solid_kmer_partition_file)
{

if (nbkmers_solid==0){
	BinaryBank *solid_kmer_colour = new BinaryBank(
			return_file_name(solid_kmer_partition_file),
			sizeof(kmer_type) + sizeof(KmerColour), 0);
	nbkmers_solid = solid_kmer_colour->nb_elements();
	solid_kmer_colour->close();
//	b1_size = 20239964;
//	uint64_t estimated_BL1_freesize =  (uint64_t)(nbkmers_solid*NBITS_PER_KMER);

}
  int64_t NbInsertedKmers;
  char * rseq;
  int readlen;
  kmer_type kmer, graine, graine_revcomp;
  KmerColour kmer_colour;

  char temp_file[1024];
  sprintf(temp_file,"%s_%s",solid_kmer_partition_file, false_positive_kmers_file);

  // **** Initialize B2, B3, B4 and T4 ****
//  Bank *FalsePositives = new Bank(return_file_name(false_positive_kmers_file));
  Bank *FalsePositives = new Bank(return_file_name(temp_file));
  uint64_t nbFP = countFP(FalsePositives);

  FPSetCascading4 *fp = new FPSetCascading4;

  fp->bloom2 = new Bloom((uint64_t)(nbFP * NBITS_PER_KMER));
  fp->bloom2->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  uint64_t estimated_T2_size = max((int)ceilf(nbkmers_solid * (double)powf((double)0.62, (double)NBITS_PER_KMER)), 1);
  uint64_t estimated_T3_size = max((int)ceilf(nbFP          * (double)powf((double)0.62, (double)NBITS_PER_KMER)) ,1);

  fp->bloom3 = new Bloom((uint64_t)(estimated_T2_size * NBITS_PER_KMER));
  fp->bloom3->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  fp->bloom4 = new Bloom((uint64_t)(estimated_T3_size * NBITS_PER_KMER));
  fp->bloom4->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  /*
   * NBITS_PER_KMER ~ m/n: predefinded
   * optimal k = m/n ln(2) ~~ 0.7*NBITS_PER_KMER
   * 2^(-k) ~~ 0.6185^(m/n)
   *
   * ln(2) ~ 0.6931472
   *
   *
   *  m = ceil((n * log(p)) / log(1.0 / (pow(2.0, log(2.0)))));
   *  k = round(log(2.0) * m / n);
   *
   *  ln p = - (m/n) ln(2)^2
   *
   *  m = (-n ln(p) ) / (ln(2)^2)
   *
   */



  // **** Insert the false positives in B2 ****
  NbInsertedKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);

    fp->bloom2->add(kmer);

    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive B2 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,"%cInsert false positive B2 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  printf("\nInserted %lld (estimated, %lld) kmers in B2.\n", NbInsertedKmers, nbFP);


  //  **** Insert false positives in B3 and write T2
  int addKmers = 0;
  NbInsertedKmers = 0;
  FILE *T2_file = fopen(return_file_name("t2_kmers"), "w+"); // We will read this file later, when filling T4
//  BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer),0);
//  while(SolidKmers->read_element(&kmer))
  BinaryBank *solid_kmer_colour = new BinaryBank(
  			return_file_name(solid_kmer_partition_file),
  			sizeof(kmer_type) + sizeof(KmerColour), 0);
  while(solid_kmer_colour->read_kmer_colour(&kmer, &kmer_colour) )
  {
    if (fp->bloom2->contains(kmer))
    {
      if (!fwrite(&kmer, sizeof(kmer), 1, T2_file))
      {
	printf("error: can't fwrite (disk full?)\n");
	exit(1);
      }

      fp->bloom3->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers% print_table_frequency)==0)
      fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  solid_kmer_colour->close();

  printf("\nInserted %lld (estimated, %llu) kmers in B3.\n", addKmers, estimated_T2_size);


  // **** Insert false positives in B4 (we could write T3, but it's not necessary)
  FalsePositives = new Bank(return_file_name(temp_file));
  NbInsertedKmers = 0;
  addKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);

    if (fp->bloom3->contains(kmer))
    {
      fp->bloom4->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive B4 %lld",13,NbInsertedKmers);
  }
  fprintf (stderr,"%cInsert false positive B4 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  printf("\nInserted %lld (estimated, %lld) kmers in B4.\n", addKmers, estimated_T3_size);


  // **** Count and insert false positives in T4
  rewind(T2_file);
  addKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
    if (fp->bloom4->contains(kmer))
      addKmers++;

  fp->false_positives = new FPSet(addKmers);
  rewind(T2_file);
  printf("====addMkersin T4:%d\n",addKmers);
  addKmers = 0;
  NbInsertedKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
  {
    if (fp->bloom4->contains(kmer))
    {
      fp->false_positives->insert(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    if ((NbInsertedKmers%print_table_frequency)==0)
      fprintf (stderr,"%cInsert false positive T4 %lld",13,NbInsertedKmers);
  }
  fp->false_positives->finalize();
  fprintf (stderr,"%cInsert false positive T4 %lld", 13,NbInsertedKmers);
  fclose(T2_file);

  printf("\nInserted %lld (estimated, %lld) kmers in T4.\n\n", addKmers, (uint64_t)fp->false_positives->capacity());

  uint64_t total_size = get_FPSetCascading4_size(fp);
  fp->set_total_memory(total_size);
  print_size_summary(fp);

  return fp;
}
