CC=/fs/sz-user-supported/Linux-x86_64/packages/gcc-4.2.4/bin/g++
CFLAGS=-O3 -fopenmp
INCLUDEDIR=/fs/szannotation/src/boost_1_42_0
CPPUNIT_PATH = /fs/sz-user-supported/Linux-x86_64
EXE_FILES = correct count-kmers count-qmers reduce-kmers reduce-qmers trim build_bithash correct_stats
.PHONY: all clean

all: $(EXE_FILES)

clean:
	-rm $(EXE_FILES) *.o

correct: correct.cpp Read.o bithash.o
	$(CC) $(CFLAGS) correct.cpp Read.o bithash.o -o correct

count-kmers: count-kmers.cpp count.o
	$(CC) $(CFLAGS) count-kmers.cpp count.o -o count-kmers

count-qmers: count-qmers.cpp count.o
	$(CC) $(CFLAGS) count-qmers.cpp count.o -o count-qmers

reduce-kmers: reduce-kmers.cpp
	$(CC) $(CFLAGS) reduce-kmers.cpp -o reduce-kmers

reduce-qmers: reduce-qmers.cpp
	$(CC) $(CFLAGS) reduce-qmers.cpp -o reduce-qmers

trim: trim.cpp Read.o bithash.o
	$(CC) $(CFLAGS) trim.cpp Read.o bithash.o -o trim

build_bithash: build_bithash.cpp bithash.o
	$(CC) $(CFLAGS) build_bithash.cpp bithash.o -o build_bithash

correct_stats: stats.cpp
	$(CC) $(CFLAGS) stats.cpp -o correct_stats

Read.o: Read.cpp Read.h bithash.o
	$(CC) $(CFLAGS) -c Read.cpp

bithash.o: bithash.cpp bithash.h
	$(CC) $(CFLAGS) -c bithash.cpp

count.o: count.cpp
	$(CC) $(CFLAGS) -c count.cpp

#test: bithash.cpp bithash.h mainTest.cpp bithashTest.cpp bithashTest.h
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp ReadTest.cpp Read.cpp prefix_tree.cpp node.cpp -o test
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp bithashTest.cpp bithash.cpp -o test