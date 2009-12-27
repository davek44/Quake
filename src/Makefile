CC=/fs/sz-user-supported/Linux-x86_64/packages/gcc-4.2.4/bin/g++
CFLAGS=-O3 -fopenmp
CPPUNIT_PATH = /fs/sz-user-supported/Linux-i686

correct: correct.cpp Read.o prefix_tree.o pt_node.o
	$(CC) $(CFLAGS) correct.cpp Read.o prefix_tree.o pt_node.o -o correct

Read.o: Read.cpp Read.h
	$(CC) $(CFLAGS) -c Read.cpp

prefix_tree.o: prefix_tree.cpp prefix_tree.h
	$(CC) $(CFLAGS) -c prefix_tree.cpp

pt_node.o: pt_node.cpp pt_node.h
	$(CC) $(CFLAGS) -c pt_node.cpp

all:
	$(CC) $(CFLAGS) correct.cpp Read.cpp prefix_tree.cpp pt_node.cpp -o correct

test:
	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp ReadTest.cpp Read.cpp prefix_tree.cpp pt_node.cpp -o test

clean:
	rm correct test *.o