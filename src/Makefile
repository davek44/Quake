CC=/fs/sz-user-supported/Linux-x86_64/packages/gcc-4.2.4/bin/g++
CFLAGS=-O3 -fopenmp
CPPUNIT_PATH = /fs/sz-user-supported/Linux-x86_64

#correct: correct.cpp Read.o prefix_tree.o node.o
#	$(CC) $(CFLAGS) correct.cpp Read.o prefix_tree.o node.o -o correct
correct: correct.cpp Read.o bithash.o
	$(CC) $(CFLAGS) correct.cpp Read.o bithash.o -o correct
#correct: correct.cpp Read.o dawg.o node.o
#	$(CC) $(CFLAGS) correct.cpp Read.o dawg.o node.o -o correct

trim: trim.cpp Read.o bithash.o
	$(CC) $(CFLAGS) trim.cpp Read.o bithash.o -o trim

Read.o: Read.cpp Read.h bithash.o
	$(CC) $(CFLAGS) -c Read.cpp

#prefix_tree.o: prefix_tree.cpp prefix_tree.h
#	$(CC) $(CFLAGS) -c prefix_tree.cpp

bithash.o: bithash.cpp bithash.h
	$(CC) $(CFLAGS) -c bithash.cpp

#dawg.o: dawg.cpp dawg.h
#	$(CC) $(CFLAGS) -c dawg.cpp

#node.o: node.cpp node.h
#	$(CC) $(CFLAGS) -c node.cpp

#test: bithash.cpp bithash.h mainTest.cpp bithashTest.cpp bithashTest.h
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp ReadTest.cpp Read.cpp prefix_tree.cpp node.cpp -o test
#	$(CC) $(CFLAGS) -I${CPPUNIT_PATH}/include -L${CPPUNIT_PATH}/lib -lcppunit -ldl mainTest.cpp bithashTest.cpp bithash.cpp -o test

clean:
	rm correct trim test *.o