CXX=		g++
CXXFLAGS=	-std=c++11 -pipe -Wall -O3 -fopenmp
CC=			gcc
CFLAGS=		-O2 -Wall
CPPFLAGS=	-DNDEBUG
OBJS=		library.o jumpgate.o kthread.o bgzf.o hts.o sam.o

all:quartz

quartz:quartz.o $(OBJS)
		$(CXX) -D_GLIBCXX_PARALLEL -o $@ $(CXXFLAGS) $(OBJS) quartz.o -lz -lpthread

bgzf.o:bgzf.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT -o $@ bgzf.c

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(CXXFLAGS) -- *.cpp *.c)

clean:
		rm -f quartz *.o

# DO NOT DELETE

jumpgate.o: global.h jumpgate.h
library.o: global.h
quartz.o: global.h jumpgate.h
bgzf.o: bgzf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
