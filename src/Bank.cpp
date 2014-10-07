//
//  Bank.cpp
//
//  Created by Guillaume Rizk on 28/11/11.
//

//TEST

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <inttypes.h>
       #include <stdio.h>
       #include <stdlib.h>
       #include <string.h>
#include <cmath> // for log2f

#include "Bank.h"
#include "Kmer.h" // Bank (almost) doesn't need Kmer.h, but KmersBuffer certainly does
#include "lut.h"
#include <errno.h>
using namespace std;

off_t fsize(const char *filename) {
    struct stat st; 
    
    if (stat(filename, &st) == 0)
        return st.st_size;
    
    return -1; 
}

// just a macro to open file indexed by i
void Bank::open_stream(int i)
{
    buffered_file[i]->stream = gzopen(buffered_file[i]->fname,"r");
    if (buffered_file[i]->stream == NULL)
    {
        printf("error opening file: %s\n",buffered_file[i]->fname);
        exit(1);
    }
}
// and close it
void Bank::close_stream(int i)
{
    gzclose(buffered_file[i]->stream);
    buffered_file[i]->stream = NULL;
}

// the following functions are adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
inline bool rebuffer(buffered_file_t *bf)
{
    if (bf->eof)
        return false;
    bf->buffer_start = 0;
    bf->buffer_end = gzread(bf->stream, bf->buffer, BUFFER_SIZE);
    if (bf->buffer_end < BUFFER_SIZE)
        bf->eof = 1;
    if (bf->buffer_end == 0) 
        return false;
    return true;
}

inline signed char buffered_getc(buffered_file_t *bf)
{
    if (bf->buffer_start >= bf->buffer_end)
        if (! rebuffer(bf))
            return -1;
    return (signed char) ( bf->buffer[bf->buffer_start++] );
}

#define nearest_power_of_2(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

inline signed int Bank::buffered_gets(buffered_file_t *bf, variable_string_t *s, char *dret, bool append, bool allow_spaces)
{
    if (dret) *dret = 0;
    if (!append)
        s->length = 0;
    if (bf->buffer_start >= bf->buffer_end && bf->eof)
        return -1;
    while (1)
    {
        int i;
        if (bf->buffer_start >= bf->buffer_end)
            if (! rebuffer(bf))
                break;
        if (allow_spaces)
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                if (bf->buffer[i] == '\n') 
                    break;
        }
        else
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                // isspace() answers yes for ' ', \t, \n, \v, \f, \r
                if (isspace(bf->buffer[i]))
                    break;
        }
        if (s->max - s->length < (i - bf->buffer_start + 1))
        {
            s->max = s->length + (i - bf->buffer_start + 1);
            nearest_power_of_2(s->max);
            s->string = (char*)realloc(s->string,s->max);

        } 
        memcpy(s->string + s->length, bf->buffer + bf->buffer_start, i - bf->buffer_start);
        s->length += i - bf->buffer_start;
        bf->buffer_start = i + 1;
        if (i < bf->buffer_end)
        {
            if (dret)
                *dret = bf->buffer[i];
            break;
        }
    }
    if (s->string == NULL)
    {
        s->max = 256;
        s->string = (char*)calloc(256,1);
    }
    else if ( allow_spaces && s->length > 1 && s->string[s->length-1] == '\r')
        s->length--;
    s->string[s->length]= '\0';
    return s->length;
}

void Bank::rewind_all()
{
    for (int i=0; i<nb_files; i++)
    {
        if (buffered_file[i]->stream != NULL)
        {
            gzclose(buffered_file[i]->stream);
            buffered_file[i]->stream = NULL;
        }
        buffered_file[i]->last_char = buffered_file[i]->eof = buffered_file[i]->buffer_start = buffered_file[i]->buffer_end = 0;
    }
    index_file = 0;
    open_stream(index_file);
}

// THIS READS FASTQ or FASTA, compressed with gzip or not
// no limit on read length, allows multi-line reads
// returns true if a read was successfuly read
//         false if end of file
// adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
bool  Bank::get_next_seq_from_file(char **nseq, char **cheader, int *len, int *hlen, int file_id, KmerColour *col)
{
    signed char c;
    buffered_file_t *bf = buffered_file[file_id];
	if (col != NULL) {
		*col= (bf->file_colour);
	}

    if (bf->last_char == 0)
    {
        while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '@'); // go to next header
        if (c == -1)
            return false; // eof
        bf->last_char = c;
    }
    read->length = dummy->length = 0;

    if (buffered_gets(bf, header, (char *)&c, false, false) < 0) //ici
        return false; // eof
    if (c != '\n')
        buffered_gets(bf, dummy, NULL, true, true); // read header //dummy instead of header to stop before first space
    

    while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if (c == '\n')
            continue; // empty line
        read->string[read->length++] = c;
        buffered_gets(bf, read, NULL, true, true);
    }
    if (c == '>' || c == '@')
        bf->last_char = c;
    if (read->length + 1 >= read->max)
    {
        read->max = read->length + 2;
        nearest_power_of_2(read->max);
        read->string = (char*) realloc(read->string, read->max);
    }
    read->string[read->length] = '\0';
    if (c == '+') // fastq
    {
        if (dummy->max < read->max) // resize quality to match read length
        {
            dummy->max = read->max;
            dummy->string = (char*)realloc(dummy->string, dummy->max);
        }
        while ( (c = buffered_getc(bf)) != -1 && c != '\n'); // read rest of quality comment
        while (buffered_gets(bf, dummy, NULL, true, true) >= 0 && dummy->length < read->length); // read rest of quality
        bf->last_char = 0;
    }
    *len = read->length;
    *nseq = read->string;
    if (cheader && hlen)
    {
        *cheader = header->string;
        *hlen = header->length;
    }
    return true;
}

// wrapper
bool  Bank::get_next_seq_from_file(char **nseq, int *len, int file_id)
{
    return get_next_seq_from_file(nseq,NULL,len,NULL,file_id, NULL);
}

// wrapper
bool Bank::get_next_seq(char **nseq, char **cheader, int *len, int *hlen, KmerColour *col)
{
//	printf("NULL:%s\n", *col);
    bool success = get_next_seq_from_file(nseq,cheader,len,hlen,index_file, col);
    if (success)
        return true;
    
    // cycle to next file if possible
    if ( index_file < nb_files-1 )
    {
        close_stream(index_file);
        index_file++;
        open_stream(index_file);
        return get_next_seq(nseq,cheader, len,hlen, col);
    }
    return false;
}

// wrapper
bool Bank::get_next_seq(char **nseq, int *len)
{
  return get_next_seq(nseq,NULL,len,NULL, NULL); //&t | nullptr TODO: which is the best??
}


bool Bank::get_next_seq_colour(char **nseq, int *len, KmerColour *colour)
{
	bool nextSeq = get_next_seq(nseq, NULL, len, NULL, colour);
	return nextSeq;
}

bool Bank::get_next_seq_colour_in_seq_name(char **nseq, char **cheader, int *len, int *hlen, int file_id, KmerColour *kmer_colour)
{//Non-generic customised reader

	signed char c;
	buffered_file_t *bf = buffered_file[file_id];

	if (bf->last_char == 0)
	{
		while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '@'); // go to next header
		if (c == -1)
			return false; // eof
		bf->last_char = c;
	}
//	printf("last_char:%c\n", bf->last_char);
	if (bf->last_char == '>') {
		int code = buffered_gets(bf, colour, NULL, false, false);

		*kmer_colour = stoi(colour->string);
//		printf("this should be colour:%s %d %d =%c=\n", colour->string, *kmer_colour, aoeu);
		if(code==-1){
			return false;
		}
	}
	read->length = dummy->length = 0;

//	if (buffered_gets(bf, header, (char *)&c, false, false) < 0) //ici
//		return false; // eof
//	if (c != '\n'){
//		buffered_gets(bf, dummy, NULL, true, true); // read header //dummy instead of header to stop before first space
//	}
//	printf("AA header:%d dummy:%d\n", header->length, dummy->length);
	while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '+' && c != '@')
	{
		if (c == '\n')
			continue; // empty line
		read->string[read->length++] = c;
		buffered_gets(bf, read, NULL, true, true);
//		printf("in read;%s\n",read->string);
	}

	if (c == '>' || c == '@')
		bf->last_char = c;
	if (read->length + 1 >= read->max)
	{
		read->max = read->length + 2;
		nearest_power_of_2(read->max);
		read->string = (char*) realloc(read->string, read->max);
	}
	read->string[read->length] = '\0';
//	    if (c == '+') // fastq
//	    {
//	        if (dummy->max < read->max) // resize quality to match read length
//	        {
//	            dummy->max = read->max;
//	            dummy->string = (char*)realloc(dummy->string, dummy->max);
//	        }
//	        while ( (c = buffered_getc(bf)) != -1 && c != '\n'); // read rest of quality comment
//	        while (buffered_gets(bf, dummy, NULL, true, true) >= 0 && dummy->length < read->length); // read rest of quality
//	        bf->last_char = 0;
//	    }

	*len = read->length;
	*nseq = read->string;
	if (cheader && hlen)
	{
		*cheader = header->string;
		*hlen = header->length;
	}
	return true;

}

bool Bank::get_next_seq_colour_in_seq_name(char **nseq, int *len, KmerColour *col){
	bool nextSeq = get_next_seq_colour_in_seq_name(nseq, NULL, len, NULL, 0, col); //Assuming only 1 file
	return nextSeq;
}

// had to move the Bank(x,x) constructor to an init() to avoid calling a constructor inside the Bank(x) constructor
void Bank::init(char **fname, int nb_files_)
{
    int64_t i;
    nb_files = nb_files_;
    filesizes = 0;

    // open the reads file, don't know if it is a fasta/q file or a list of file names yet
    gzFile tempfile = gzopen(fname[0],"r");
    if (tempfile == NULL)
    {
        char *buffer = (char*)malloc(BUFSIZ);
        strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s   %s \n",buffer,fname[0]);
        free(buffer);
        exit(1);
    }
    char deb=(char)gzgetc(tempfile);

    char **nfname;// [MAX_NB_FILES][TAILLE_NOM];
    nfname = (char**) malloc(sizeof(char*)*MAX_NB_FILES);
    for(int jj=0; jj<MAX_NB_FILES; jj++ )
	nfname [jj] =  (char*) malloc(sizeof(char)*TAILLE_NOM);

    if(deb=='>' || deb=='@' || deb==EOF)
    { // file is a fasta/q file
        gzclose(tempfile);
    }
    else // file contains a list of file names
    {
        char* ret;
        gzungetc(deb,tempfile);
        printf("File %s starts with character \"%c\", hence is interpreted as a list of file names\n",fname[0],deb );
        int ii;
        // get the filenames
        for (ii=0; ii<MAX_NB_FILES ; ii++)
        {
            ret = gzgets(tempfile, nfname[ii], BUFFER_SIZE);
            if (ret != NULL) {
                // remove \r \n chars
                char *endline = strchr(nfname[ii], '\n');
                if (endline)
                    *endline='\0';
                endline = strchr(nfname[ii], '\r');
                if (endline)
                    *endline='\0';
            }
            else // no more filenames
                break;
        }
        printf("Reading %i read files (first 2 are..):\t%s\t%s\n",ii, nfname[0], nfname[1]);
        if(ii==MAX_NB_FILES)
            printf("Warning! using max number of read files (%i)\n",ii);

        nb_files = ii;
        fname = (char **) nfname;
        gzclose(tempfile);

    }

    // estimate total size of files
    for (i=0; i<nb_files; i++)
    {
        bool compressed = false;
        uint64_t estimated_filesize;

        if (strstr(fname[i],"gz") == (fname[i]+strlen(fname[i])-2) ) compressed=true;
        if (compressed)
            // crude hack, based on Quip paper reporting compression ratio (~0.3).
            // gzseek(SEEK_END) isn't supported. need to read whole file otherwise :/
            estimated_filesize = fsize(fname[i]) * 4;
        else
            estimated_filesize = fsize(fname[i]);

        filesizes += estimated_filesize;
    }

    // initialize the buffers
    buffered_file = (buffered_file_t**) malloc(sizeof(buffered_file_t *)*nb_files);
    for (i=0; i<nb_files; i++)
    {
        buffered_file[i] = (buffered_file_t *)calloc(1, sizeof(buffered_file_t));
        buffered_file[i]->buffer = (unsigned char*) malloc(BUFFER_SIZE);
        buffered_file[i]->fname = strdup(fname[i]);
        buffered_file[i]->file_colour = i;

    }
    rewind_all(); // initialize the get_next_seq iterator to the first file

    // init read and dummy (for readname and quality)
    read = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    dummy = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    header = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    colour = (variable_string_t*) calloc(1,sizeof(variable_string_t));

    //Moved malloc read-> string.
	read->max = 256;
	read->string = (char*) malloc(read->max);
	colour-> max = 8; //longer than this?? really?
	colour->string = (char*) malloc(colour->max);

    for(int jj=0; jj<MAX_NB_FILES; jj++ )
      free	(nfname [jj]);
    free(nfname);
}

Bank::Bank(char *fname0)
{
    char *fname[1] = { fname0 };
    init(fname, 1);
}


Bank::Bank(char **fname, int nb_files_)
{
    init(fname,nb_files_);
}

Bank::~Bank(){

    variable_string_t * to_free[4] = {read, dummy, header, colour};
    for (int i = 0; i < 4; i++)
    {
        if (to_free[i])
        {
            if (to_free[i]->string)
                free(to_free[i]->string);
            free(to_free[i]);
        }
    }
    close();
    for (int i=0; i<nb_files; i++)
    {
//    	printf("stirr %d %d\n", i, buffered_file[i]->stream);
//		if(buffered_file[i] != NULL){
//			int aoeu = gzclose(buffered_file[i]->stream);
//			printf("%d\n",aoeu);
//		}
////		aoeu;

        free(buffered_file[i]->buffer);
        free(buffered_file[i]->fname);
        free(buffered_file[i]);
    }
    free(buffered_file);
}

void Bank::close()
{
    for (int i=0; i<nb_files; i++){
    	if (buffered_file[i]->stream != NULL)
		{
    		gzclose(buffered_file[i]->stream);
		}
        buffered_file[i]->stream = NULL;
    }
}

// estimate the volume of all redundant kmers in the reads, if they were to be stored in 2bits
uint64_t Bank::estimate_kmers_volume(int k)
{

    char * rseq;
    int readlen;
    int NbRead = 0;
    //int kmer_nbits = std::max(64,(int)pow(2,ceilf(log2f(2*k)))); // Bank assumes that a kmer is stored in the smallest integer type (e.g. uint64_t or uint128_t) // not accurate anymore with _ttmath/_largeint
    int kmer_nbits = sizeof(kmer_type)*8;
    rewind_all();
    uint64_t volume = 0;

    while (get_next_seq(&rseq,&readlen))
    {
        if (readlen >= k)
            volume += (readlen-k+1) * (uint64_t) kmer_nbits;
        if (NbRead++ == 1000)
            break;
    }

    if ( gztell(buffered_file[index_file]->stream) == 0) // empty file
        return 1;
//printf("IN estimate_kmers_volume:%lu total_filesize:%lu gztell:%lu %d\n", volume, filesizes, gztell(buffered_file[index_file]->stream) ,index_file);
    volume = (uint64_t) ( ( (float) volume ) * ( ( (float)filesizes ) / ((float) gztell(buffered_file[index_file]->stream)) ) );

    volume = volume / 1024 /1024 /8; // put it in MB
    
    if (volume == 0)  // tiny files fix
        volume = 1;

    rewind_all();
    return volume;
}

// estimate the number of reads
uint64_t Bank::estimate_nb_reads()
{
    char * rseq;
    int readlen;
    int NbRead = 0;
    rewind_all();
    
    uint64_t volume = 0;
    while (get_next_seq(&rseq,&readlen))
    {
        volume += 1;
        if (NbRead++ == 1000)
            break;
    }

    if ( gztell(buffered_file[index_file]->stream) == 0) // empty file
        return 1;

    volume = (volume * filesizes) / gztell(buffered_file[index_file]->stream); // linear extrapolation from the first 1k reads 

    rewind_all();
    return volume;
}

// estimate maximum read length
// from the first 10000 reads of each file
int Bank::estimate_max_readlen()
{
    char * rseq;
    int readlen;
    rewind_all();
    int max_readlen = 0;
    uint64_t volume = 0;

    index_file = 0;
    
    while ( index_file < nb_files )
    {
        int NbRead = 0;
        while (get_next_seq_from_file(&rseq,NULL,&readlen,NULL,index_file, NULL))
        {
            max_readlen = max(readlen, max_readlen);
            if (NbRead++ == 10000)
                break;
        }
        index_file++;
    } 

    index_file = 0;
    return max_readlen;
}

void Bank::save_position()
{
    restore_index_file = index_file;
    restore_pos = gztell(buffered_file[index_file]->stream) - buffered_file[index_file]->buffer_end + buffered_file[index_file]->buffer_start;
}

void Bank::load_position()
{
    close_stream(index_file);
    index_file = restore_index_file;
    open_stream(index_file);
    gzseek(buffered_file[index_file]->stream, restore_pos, SEEK_SET);
    buffered_file[index_file]->eof = false;
    rebuffer(buffered_file[index_file]);
}


// BinaryBank: a binary file containing kmers

BinaryBankConcurrent::BinaryBankConcurrent(char *given_filename, int given_sizeElement, bool write, int given_nthreads) : BinaryBank(given_filename,given_sizeElement,write) 
{
    nthreads = given_nthreads;
    
    //free(buffer); buffer =NULL; //cannot do that
    bufferT = (void **) malloc(sizeof(void*) * nthreads);

    for (int i= 0; i< nthreads; i++)
    {
         ((void ** )bufferT)[i]= (void *) malloc( WRITE_BUFFER);

    }
    cpt_buffer_tid = (int  *)malloc(sizeof(int) * nthreads);
    memset (cpt_buffer_tid,0,sizeof(int) * nthreads);
}


void BinaryBankConcurrent::write_element_buffered( void *element, int tid)
{
    write_buffered(element, sizeElement, tid);    
}


void BinaryBankConcurrent::write_buffered( void *element, int size, int tid)
{
    write_buffered( element, size, tid, true);
}

void BinaryBankConcurrent::write_buffered( void *element, int size, int tid, bool can_flush)
{
    if(cpt_buffer_tid[tid]>= WRITE_BUFFER -100 && can_flush)
    {
        flush(tid);
    }
    
    char * buf_pt = ((char**) bufferT)[tid];    
    memcpy(buf_pt + cpt_buffer_tid[tid] , element, size);
    
    cpt_buffer_tid[tid]+=size;
    // cpt_buffer_tid[tid]++;
    

}



void BinaryBankConcurrent::flush(int tid)
{
    flockfile(binary_read_file);
    if (!fwrite( ((void **)bufferT)[tid], 1, cpt_buffer_tid[tid], binary_read_file))            
    {
        printf("error: can't fwrite (disk full?)\n");
        funlockfile(binary_read_file);
        exit(1);
    }
    cpt_buffer_tid[tid]=0;
    funlockfile(binary_read_file);

}


//should be called by only one of the threads
void BinaryBankConcurrent::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    for(int ii=0; ii< nthreads; ii++)
    {
        if(cpt_buffer_tid[ii])
        {
            if (!fwrite(((void **)bufferT)[ii], 1, cpt_buffer_tid[ii], binary_read_file))
          //      if (!fwrite(((void **)bufferT)[ii], sizeElement, cpt_buffer_tid[ii], binary_read_file))

            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }
        }
        cpt_buffer_tid[ii]=0;
    }
    
    fclose(binary_read_file);
}


BinaryBankConcurrent::~BinaryBankConcurrent()
{
    for (int i= 0; i< nthreads; i++)
    {
        free(((void ** )bufferT)[i]);
        ((void ** )bufferT)[i]=NULL;
    }
    free(bufferT);
    free(cpt_buffer_tid);

}



BinaryBank::BinaryBank(char *given_filename, int given_sizeElement, bool write) : sizeElement(given_sizeElement)
{
    strcpy(filename,given_filename);
    open(write);
    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
    buffer = (void *) malloc(given_sizeElement * buffer_size_nelem);
    cpt_buffer=0;
}


void BinaryBank::write_element( void *element)
{
  //  flockfile(binary_read_file);ca
   // fprintf(stderr,"write elem %lli \n",*(int64_t *)element);
    if (!fwrite(element, sizeElement, 1, binary_read_file))
    {
       // funlockfile(binary_read_file);
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
  //  funlockfile(binary_read_file);
}


void BinaryBank::write_element_buffered( void *element)
{
    
    if(cpt_buffer==buffer_size_nelem)
    {
        if (!fwrite(buffer, sizeElement, buffer_size_nelem, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    //((kmer_type *)buffer)[cpt_buffer]= *((kmer_type *)element);
    memcpy((unsigned char *)buffer + (cpt_buffer * sizeElement), element, sizeElement);
    cpt_buffer++;

 
    
}



size_t BinaryBank::read_element( void *element)
{
    return fread(element, sizeElement,1, binary_read_file);
}

size_t BinaryBank::read_element_buffered( void *element)
{
    if(cpt_buffer==0)
    {
        cpt_buffer=fread(buffer, sizeElement,buffer_size_nelem, binary_read_file);
        if (cpt_buffer==0) return 0;
    }
    memcpy(element, (unsigned char *)buffer + (cpt_buffer-1) * sizeElement, sizeElement);
    cpt_buffer --;
    return cpt_buffer+1; // nb remaining before read
}

// used to read/write raw information to the binary file (e.g. kmer count)

void BinaryBank::write( void *element, int size)
{
    if (!fwrite(element, size, 1, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
}

size_t BinaryBank::read( void *element, int size)
{

    return fread(element, size,1, binary_read_file);
}


void BinaryBank::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryBank::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    if(cpt_buffer)
    {
    if (!fwrite(buffer, sizeElement, cpt_buffer, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryBank::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",buffer,write,filename);
        free(buffer);
        exit(1);
    }

}

off_t BinaryBank::file_size()
{
  return fsize(filename);
}

off_t BinaryBank::nb_elements()
{
  return fsize(filename)/sizeElement;
}


size_t BinaryBank::read_kmer(void *element)
{
    return fread(element, kSizeOfKmerType, 1, binary_read_file);
}

size_t BinaryBank::read_kmer_colour(void *element_kmer, void* element_colour)
{
	size_t size = fread(element_kmer, kSizeOfKmerType, 1, binary_read_file);
	size += fread(element_colour, kSizeOfKmerColour, 1, binary_read_file);
    return size;
}

size_t BinaryBank::read_kmer_skip_colour( void *element)
{
	size_t size = fread(element, kSizeOfKmerType, 1, binary_read_file);
	fseek(binary_read_file, kSizeOfKmerColour, SEEK_CUR);
	return size;
}

size_t BinaryBank::read_colour( void *element)
{
    return fread(element, kSizeOfKmerColour, 1, binary_read_file);
}

BinaryBank::~BinaryBank()
{
    if(buffer!=NULL)
    {
        free (buffer); //buffer =NULL;
    }

}


/////////////class BinaryReads a file containing reads

BinaryReads::~BinaryReads()
{
    free (buffer); buffer = NULL;
}


BinaryReads::BinaryReads(char *given_filename,  bool write)
{
    read_write_buffer_size = BINREADS_BUFFER;
    strcpy(filename,given_filename);
    open(write);
    buffer = (unsigned char *) malloc(read_write_buffer_size*sizeof(unsigned char));
    cpt_buffer = 0;
}


void BinaryReads::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryReads::close()
{
    unsigned int block_size =0;
    //flush buffer
    if(cpt_buffer)
    {
//        printf("close :write block %i \n",cpt_buffer);
        block_size = cpt_buffer;
        fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryReads::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",buffer,write,filename);
        free(buffer);
        exit(1);
    }
    
}



//format is
// 32 bit integer = readlen,  then seq in binary
// then next read..
//32 bit len is overkill but simpler
//also makes buffer then write block with header : size of block to read, with n reads .... will allow large fread when reading this file ...
void BinaryReads::write_read(char * read, int readlen)
{
    int tai = readlen;
    unsigned char rbin;
    char * pt = read;
    unsigned int block_size = 0;
    
//    printf("write read %s %i / %i   readlen %i \n",pt , cpt_buffer,read_write_buffer_size,readlen);
    //todo : also flush to disk  sometimes (ie if very large buffer, to create smaller blocks..)
    if((cpt_buffer && cpt_buffer >= (read_write_buffer_size-readlen)) || cpt_buffer > 10000000 )  ////not enough space to store next read   true space is 4 + readlen/4 + rem
        //flush buffer to disk
    {
        
        block_size = cpt_buffer;
        
        //printf("write block %i\n",block_size);
        if(block_size) fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file)) // write a block, it ends at end of a read
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    //check if still not enough space in empty buffer : can happen if large read, then enlarge buffer
    if(read_write_buffer_size < readlen)
    {
        read_write_buffer_size = 2*readlen; // too large but ok
        buffer =  (unsigned char *) realloc(buffer,sizeof(unsigned char) * read_write_buffer_size);
    }

    memcpy(buffer+cpt_buffer,&readlen,sizeof(int));
    cpt_buffer+= sizeof(int);
    
    //fwrite( (void *) &readlen, sizeof(int), 1, binary_read_file);
//    printf("Buffer=%s=cpt=%d=\n", pt, cpt_buffer);
    
    for (tai=readlen; tai>=4  ; tai-=4)
    {
        rbin = code4NT(pt);
      //  fwrite((void *) &rbin, 1,1,binary_read_file );
        buffer[cpt_buffer]=rbin; cpt_buffer++;
//        printf("B:%d%c%c%c%c:%u ", tai,pt[0], pt[1], pt[2], pt[3], rbin);
        pt +=4;

    }
    
    //then remaining
    if(tai)
    {
        rbin = code_n_NT(pt,tai);
       // fwrite( (void *) &rbin,1,1,binary_read_file);
        buffer[cpt_buffer]=rbin; cpt_buffer++;
    }
}



void  compute_kmer_table_from_one_seq(int readlen, char * seq, kmer_type * kmer_table )  //,char * pkmer_table //pour remplissage table loc
{
    kmer_type graine = codeSeed(seq);
    kmer_type graine_revcomp = revcomp(graine);
    kmer_table[0] = min(graine,graine_revcomp);
    seq++;
    for (int i=1; i<readlen-sizeKmer+1; i++)
    {
        graine =   (graine * 4 + NT2int(seq[sizeKmer-1])) & kmerMask   ;
        graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask ;
        kmer_table[i] = min(graine,graine_revcomp);
        seq++;
    }
}

void  compute_kmer_table_from_one_seq_colour(int readlen, char * seq, kmer_type * kmer_table, KmerColour readColour, KmerColour *kmer_table_colour )  //,char * pkmer_table //pour remplissage table loc
{
    kmer_type graine = codeSeed(seq);
//    printf("SEQ:%s\n",seq);
    kmer_type graine_revcomp = revcomp(graine);
    kmer_table[0] = min(graine,graine_revcomp);
    kmer_table_colour[0] = readColour;
    seq++;
//    printf("kmer:%lu\n",kmer_table[0]);
    for (int i=1; i<readlen-sizeKmer+1; i++)
    {
        graine =   (graine * 4 + NT2int(seq[sizeKmer-1])) & kmerMask   ;
        graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask ;
        kmer_table[i] = min(graine,graine_revcomp);
        kmer_table_colour[i] = readColour;
//        printf("kmer:%lu\n",kmer_table[i]);
        seq++;
    }

}




////kmers buffer


KmersBuffer::KmersBuffer(BinaryReads *bfile, int  pbuffer_size, int nseq_task )
{
    
    read_write_buffer_size = BINREADS_BUFFER;
    buffer = ( char *) malloc(read_write_buffer_size*sizeof( char));
    cpt_buffer = 0;
    cpt_binSeq_read =0; binSeq_toread =0;
    max_read_length = KMERSBUFFER_MAX_READLEN;
    binfile = bfile;
    buffer_size = pbuffer_size;
    kmers_buffer =  (kmer_type *) malloc(sizeof(kmer_type) * buffer_size);
   // binSeq =  (char *) malloc(sizeof(char) * max_read_length); // no need to alloc ram for binse : will points to buffer
    binSeq_extended =  (char *) malloc(sizeof(char) * max_read_length);
    blocksize_toread =0;


    nseq_step = nseq_task;
    binary_read_file = bfile->binary_read_file;
    
}


void KmersBuffer::reset_max_readlen(int read_length)
{

    max_read_length = read_length;

  //  binSeq =  (char *) realloc(binSeq,sizeof(char) * max_read_length);
    binSeq_extended =  (char *) realloc(binSeq_extended,sizeof(char) * max_read_length);

}

 KmersBuffer::~KmersBuffer()
{
    free (kmers_buffer);
    free(buffer);
    //free(binSeq);
    free(binSeq_extended);
}


//now returns number of kmers read
int KmersBuffer::readkmers()
{
    
    int llen;
    int * len = & llen ;
    unsigned int block_size =0;
    
    //////reading new block from disk if needed
    // cpt_buffer == blocksize_toread    tells we finished reading previous buffer
    // (binSeq_toread  <= cpt_binSeq_read) tells  we finished reading the last sequence
    if(cpt_buffer == blocksize_toread  && (binSeq_toread  <= cpt_binSeq_read))
    {
        flockfile(binary_read_file);
        if( ! fread(&block_size,sizeof(unsigned int),1, binary_read_file)) //read block header
        {
            funlockfile(binary_read_file);
            return -1; // no more blocks to read
        }
        
        
        if(block_size >= read_write_buffer_size) // block buffer need to be enlarged
        {
            read_write_buffer_size = 2*block_size;
            buffer =  ( char *) realloc(buffer,sizeof( char) * read_write_buffer_size);
        }
        
        fread(buffer,sizeof( char),block_size, binary_read_file); // read a block of sequences into the buffer
        funlockfile(binary_read_file);
        cpt_buffer = 0;
        blocksize_toread = block_size;        
    }
    ///////////////////////
    
    
    
    
    //now parse the whole block in ram
    
    int i,j;
    int nchar;
    unsigned char fournt;
    
    nkmers = 0;
    int nseq_lues = 0;
    //cpt_buffer : how much we have already read in the buffer
    //blocksize_toread : how much there is to read in the buffer
    while(cpt_buffer < blocksize_toread || ( binSeq_toread  > cpt_binSeq_read)) //while work to do
    {
        
        if( binSeq_toread <= cpt_binSeq_read)// read new sequence if needed  //we  will put one sequence into binSeq_extended
        {
            memcpy(len,buffer+cpt_buffer,sizeof(int)); // the sequence length
            cpt_buffer += sizeof(int);
            nseq_lues ++;
            
            if( (*len) > max_read_length) reset_max_readlen((int)(1.2*(*len))); // resize memory for sequence if needed
            nchar = ((*len)+3)/4; // number of bytes used to encode the sequence in its binary format (4 nt per byte)

            binSeq = buffer + cpt_buffer; // point binseq to correct place //cpt_buffer ==  where we are now in the buffer
            cpt_buffer += nchar;
            
            // on disk data was encoded with 4 nucleotides per bytes,
            // here we expand one sequence  to one nucl per byte into binSeq_extended
            //nucleotides are still encoded in [0-3]
            j=0;
            for(i=0; i<nchar; i++)
            {
                fournt = binSeq[i];
                binSeq_extended[j+3]=fournt & 3; fournt = fournt >> 2; 
                binSeq_extended[j+2]=fournt & 3; fournt = fournt >> 2;
                binSeq_extended[j+1]=fournt & 3; fournt = fournt >> 2;
                binSeq_extended[j+0]=fournt & 3;
                j+=4;
            }
            binSeq_toread = *len-sizeKmer+1; // binSeq_toread tells how many kmers there are in this sequence
            cpt_binSeq_read = 0;  // tells how many kmers we have currently parsed in this sequence
                        
        }
        
        
        {
            // binSeq_extended = beginning of the sequence,
            //  cpt_binSeq_read = how much we have already read in this sequence (when kmers_buffer is full, we can halt parsing kmers (see below) in the middle of a sequence, so this value is not necessarily 0)
            char *seq = binSeq_extended+cpt_binSeq_read;  
            kmer_type graine;
            kmer_type graine_revcomp;
            if( binSeq_toread  > cpt_binSeq_read)
                // there are still unread kmers in this sequence, here we read the first one,
                // we put it in graine / graine_revcomp  and store it in the kmers_buffer
            {
                graine = codeSeed_bin(seq);
                graine_revcomp = revcomp(graine);
                
                if(nkmers>=buffer_size)
                {
                    return nkmers;
                }
                kmers_buffer[nkmers] = min(graine,graine_revcomp); nkmers++; cpt_binSeq_read ++;
                seq++;
            }
            
            while( binSeq_toread  > cpt_binSeq_read) //while there remains kmers to be read in this sequence
            {                
                graine =  (graine * 4 + (seq[sizeKmer-1])) & kmerMask ; //parse next nucleotide to construc the next kmer
                graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[(int)(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
                kmers_buffer[nkmers] = min(graine,graine_revcomp); nkmers ++; cpt_binSeq_read ++; //we store the kmer in the kmers_buffer
                seq++;
                if(nkmers>=buffer_size) //the kmers_buffer is full, we stop
                {                    
                    return nkmers;
                }
            }
            
        }
    }
    
    // we stop when we finished one block, or when kmers_buffer is full,
    // it can happen in the middle of a sequence : the next time we call readkmers we will have to continue
    // from where we stopped  in this sequence (counter cpt_binSeq_read tells us that) 

    //while buffer is non empty, we 'expand' a sequence into  binSeq_extended
    //then we parse binSeq_extended to store kmers in the kmers_buffer
    
    return nkmers;
    
    
}

 
