
all: shuffleReads

CXX=g++
STATIC=
CCOPTS=-O2

htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

shuffleReads: ShuffleReads.o htslib/lib/libhts.a
	$(CXX) $(STATIC) $(CCOPTS) $^ -L $(PWD)/htslib/lib  -lhts -lz -lpthread -o $@ -Wl,-rpath,$(PWD)/htslib/lib

ShuffleReads.o: ShuffleReads.cpp $(HEADERS) htslib/lib/libhts.a
	$(CXX) $(CCOPTS) -c $< -I htslib/include 


