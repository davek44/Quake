CC=/fs/sz-user-supported/Linux-x86_64/packages/gcc-4.2.4/bin/g++
CFLAGS=-g -fopenmp
CPPUNIT_PATH = /fs/sz-user-supported/Linux-x86_64

.PHONY: all clean

all: correct trim build_bithash

clean:
	-rm correct trim test *.o

correct: correct.cpp Read.o bithash.o
	$(CC) $(CFLAGS) correct.cpp Read.o bithash.o -o correct

trim: trim.cpp Read.o bithash.o
	$(CC) $(CFLAGS) trim.cpp Read.o bithash.o -o trim

build_bithash: build_bithash.cpp bithash.o
	$(CC) $(CFLAGS) build_bithash.cpp bithash.o -o build_bithash

Read.o: Read.cpp Read.h bithash.o
	$(CC) $(CFLAGS) -c Read.cpp

bithash.o: bithash.cpp bithash.h
	$(CC) $(CFLAGS) -c bithash.cpp

#test: bithash.cpp bithash.h mainTest.cpp bithashTest.cpp bithashTest.h
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp ReadTest.cpp Read.cpp prefix_tree.cpp node.cpp -o test
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp bithashTest.cpp bithash.cpp -o test