#ifndef OAHash_h
#define OAHash_h
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>

#include "Kmer.h"

#ifdef _largeint
#include "LargeInt.h"
typedef LargeInt<KMER_PRECISION> key_type;
#else
#ifdef _ttmath
#include "ttmath/ttmath.h"
typedef ttmath::UInt<KMER_PRECISION> key_type;
#else
#if (! defined kmer_type) || (! defined _LP64)
typedef uint64_t key_type;
#else
typedef kmer_type key_type;
#endif
#endif
#endif

class AbstractOAHash{

protected:

#ifdef _largeint
    inline uint64_t hashcode(LargeInt<KMER_PRECISION> elem);
#endif
#ifdef _ttmath
    inline uint64_t hashcode(ttmath::UInt<KMER_PRECISION> elem);
#endif
#ifdef _LP64
    uint64_t hashcode( __uint128_t elem);
#endif
    uint64_t hashcode( uint64_t elem);


    uint64_t hash_size;
    uint64_t nb_inserted_keys;
    size_t sizeof_element_pair;

//#pragma pack(1)

public:

//    struct element_pair
//	{
//		key_type key;
//	};
    virtual void start_iterator() = 0;
	virtual bool next_iterator() = 0;

	virtual void printstat();
	virtual uint64_t memory_usage();
    virtual float load_factor();


    virtual ~AbstractOAHash(){};

//    element_pair * find_slot(key_type key);

};

class OAHash : public AbstractOAHash{


private:
//#pragma pack(1)
	struct element_pair// : AbstractOAHash::element_pair
	{
		key_type key;
		uint32_t value;
	};
	element_pair* data;

protected:
    bool is_occupied(element_pair *element);

public:

    static int size_entry();
    
    // iterator functions:
    element_pair* iterator;
    void start_iterator();
    bool next_iterator();

    OAHash(uint64_t max_elements);
    ~OAHash();

    element_pair * find_slot(key_type key);
    bool has_key(key_type graine);

    void insert(key_type graine, int value);
    void increment(key_type graine);
    bool get( key_type graine, int * val);

};

class OAHashColour: public AbstractOAHash{

private:
//#pragma pack(16)
	struct element_pair
	{
		key_type key;
		uint32_t value;
		KmerColour colour;
	};
	element_pair* data;


protected:
    bool is_occupied(element_pair *element);
public:

    static int size_entry();

    element_pair *iterator;
	void start_iterator();
	bool next_iterator();

	OAHashColour(uint64_t max_elements);
	~OAHashColour();

    void insert_colour(key_type graine, KmerColour colour);
    void increment_colour(key_type graine, KmerColour colour);
    bool get_colour( key_type graine, KmerColour *colour);

    element_pair* find_slot(key_type key);
    bool has_key(key_type graine);
};

#endif

