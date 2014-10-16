UNAME := $(shell uname)

CC=g++
CFLAGS = -O3 -Wall -D_FILE_OFFSET_BITS=64 -std=c++11 # needed to handle files > 2 GB on 32 bits systems
LDFLAGS = -lz


CC=clang++
CFLAGS = -O3 -D_FILE_OFFSET_BITS=64 -std=c++11 -stdlib=libstdc++ #-I/usr/include/c++/4.8.3/ -I/usr/include/c++/4.8.3/x86_64-redhat-linux/

ifeq ($(UNAME), FreeBSD)
CC=g++48
CFLAGS = -O3 -Wall -D_FILE_OFFSET_BITS=64 -std=c++11 -D_GLIBCXX_USE_C99
#LDFLAGS= -L/usr/local/lib/gcc48/
LDFLAGS += -Wl,-rpath,/usr/local/lib/gcc48
endif

SOURCES=Pool.cpp Bank.cpp Bloom.cpp Hash16.cpp LargeInt.cpp \
	Kmer.cpp Terminator.cpp Traversal.cpp LinearCounter.cpp \
	Set.cpp Utils.cpp SortingCount.cpp Debloom.cpp OAHash.cpp \
	KmerColour.cpp Memory.cpp SortingCountPartitions.cpp \
	Assembler.cpp

TEST_SOURCES=KmerColour.cpp
#Pool.cpp Bank.cpp Bloom.cpp Hash16.cpp LargeInt.cpp \
#	Kmer.cpp Terminator.cpp Traversal.cpp LinearCounter.cpp \
#	Set.cpp Utils.cpp SortingCount.cpp Debloom.cpp OAHash.cpp \
#	KmerColour.cpp

EXEC=minia

SRC_DIR=src/
TEST_DIR=test/
OBJ_DIR=obj/

SRC=$(SOURCES:%=$(SRC_DIR)%)
OBJ=$(SOURCES:%.cpp=$(OBJ_DIR)%.o)

TEST_SRC=$(TEST_SOURCES:%.cpp=$(TEST_DIR)%Test.cpp)
TEST_OBJ=$(TEST_SOURCES:%.cpp=$(OBJ_DIR)%Test.o)
GTEST_DIR=/usr/include/gtest/

SRCEXT   = cpp
SRCDIR   = src/
OBJDIR   = obj/
BINDIR   = bin/
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)

all: $(EXEC)


ifeq ($(prof),1)
 CFLAGS+= -pg
endif

ifeq ($(unitig),1)
 CFLAGS+= -DUNITIG
endif


ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

ifeq ($(omp),1)
 CFLAGS=-O4  -fopenmp -DOMP=1
endif


k := 0$(k) # dummy k if not specified 
K_BELOW_32 := $(shell echo $(k)\<=32 | bc)
K_BELOW_64 := $(shell echo $(k)\<=64 | bc)
ARCH := $(shell getconf LONG_BIT) # detects sizeof(int)
USING_UINT128 := 0
largeintlib := 0

ifeq ($(K_BELOW_32),0)

    # use uint128 when k<=64 and 64-bit architecture
    ifeq ($(K_BELOW_64),1)
        ifeq ($(strip $(ARCH)),64)
            CFLAGS += -Dkmer_type=__uint128_t
            USING_UINT128 := 1
        endif
    endif
    
	# use a bigint library otherwise
    ifeq ($(USING_UINT128),0)
        largeintlib := largeint#ttmath
    endif
endif

# ttmath (now, largeint) is used when you type "make k=[kmer size]" with a kmer size longer than supported integer type,
ifeq ($(largeintlib),ttmath)
    KMER_PRECISION := $(shell echo \($(k)+15\)/16 | bc)
endif
ifeq ($(largeintlib),largeint)
    KMER_PRECISION := $(shell echo \($(k)+31\)/32 | bc)
endif
ifneq ($(largeintlib),0)
    CFLAGS += -D_$(largeintlib) -DKMER_PRECISION=$(KMER_PRECISION)
endif

.PHONY: minia

.DEFAULT: inc

all: $(EXEC)




minia: clean $(OBJDIR) $(OBJ) src/Minia.cpp
	$(CC) -o $@ $(OBJ) src/Minia.cpp $(CFLAGS) $(LDFLAGS)

inc: $(OBJ) src/Minia.cpp
	$(CC) -o minia $(OBJ) src/Minia.cpp $(CFLAGS) $(LDFLAGS)

clean:
	@rm -rf $(OBJ) $(TEST_OBJ)

$(OBJDIR):
	mkdir -p $(OBJDIR)
	 
obj/%.o: src/%.cpp src/%.h test/%.cpp
	echo "BOTH!!"
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CC) -o $@ -c $< $(CFLAGS) 

obj/%.o: src/%.cpp src/%.h
#	echo "SRC ONLY" 
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CC)  -o $@ -c $< $(CFLAGS) 

obj/%.o: test/%.cpp
#	echo "TEST ONLY"
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CC) -o $@ -c $< $(CFLAGS)   -I$(CURDIR)/


#%Test.o: clean %Test.cpp $(SRP)%.cpp $(SRP)%.h
#	echo  $(SRC_DIR)%.cpp $(SRC_DIR)%.h
#	$(CC) -o $@ -c $< $(CFLAGS) -I$(GTEST_DIR) -I$(CURDIR)
#	
#$(TEST_OBJ): $(TEST_SRC) $(OBJ)
#	echo $(TEST_SRC)
#	$(CC) -o $(TEST_OBJ) -c $(TEST_SRC) $(CFLAGS) -I$(GTEST_DIR) -I$(CURDIR)


test: $(OBJ) $(TEST_OBJ)
	echo "callTest" $(TEST_OBJ) $(OBJ) $(SRC_DIR)%.cpp $(SRC_DIR)%.h
#	$(CC) $(CFLAGS) $(TEST_OBJ) $(OBJ) -I$(CURDIR) $(CURDIR)/include/gtest_main.a -lpthread -lz -o your_test;
#	$(CC) $(CFLAGS) $(TEST_SRC) -I$(CURDIR) $(CURDIR)/include/gtest_main.a -lpthread -o your_test;
#	$(CC) $(CFLAGS) $(TEST_SRC) -I$(CURDIR) $(CURDIR)/include/libgtest_main.a $(CURDIR)/include/libgtest.a -lpthread -o your_test;
#	  g++ -I${GTEST_DIR}/include path/to/your_test.cc libgtest.a -o your_test
#	
debug2: clean $(TEST_OBJ)# %Test.cpp $(SRP)%.cpp $(SRP)%.h
	echo  $(SRC_DIR)%.cpp $(SRC_DIR)%.h
	echo $(CFLAGS)
	$(CC) -o $@ -c $< $(CFLAGS) -I$(GTEST_DIR) -I$(CURDIR)

debug: CFLAGS += -g -O0
debug: buildrepo $(OBJ) src/Minia.cpp
	$(CC) -o $@_minia $(OBJ) src/Minia.cpp $(CFLAGS)  $(LDFLAGS)


buildrepo:
	@$(call make-repo)


define make-repo
   for dir in .; \
   do \
	mkdir -p $(OBJ_DIR)/$$dir; \
   done
endef


# usage: $(call make-depend,source-file,object-file,depend-file)
define make-depend
  $(CC) -MM       \
        -MF $3    \
        -MP       \
        -MT $2    \
        $(CFLAGS) -I$(CURDIR)  \
        $1
endef




#g++ -o minia obj/KmerColourTest.o obj/Pool.o obj/Bank.o obj/Bloom.o obj/Hash16.o obj/LargeInt.o obj/Kmer.o obj/Terminator.o obj/Traversal.o obj/LinearCounter.o obj/Set.o obj/Utils.o obj/SortingCount.o obj/Debloom.o obj/OAHash.o obj/KmerColour.o /home/sw167/Postdoc2/Software/minia-1.6088/include/gtest_main.a -O4 -Wall -D_FILE_OFFSET_BITS=64 -std=c++11  -lz -lpthread
#g++ -g -std=c++11  obj/*.o /home/sw167/Postdoc2/Software/minia-1.6088/include/gtest_main.a -lpthread -o your_test -lz;
