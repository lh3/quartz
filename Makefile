CXX=		g++
CXXFLAGS=	-std=c++11 -pipe -Wall -O3 -fopenmp
CPPFLAGS=	-DNDEBUG
OBJS=		library.o jumpgate.o

all:quartz

quartz:quartz.o $(OBJS)
		$(CXX) -D_GLIBCXX_PARALLEL -o $@ $(CXXFLAGS) $(OBJS) quartz.o

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(CXXFLAGS) -- *.cpp *.c)

clean:
		rm -f quartz *.o

# DO NOT DELETE

jumpgate.o: global.h jumpgate.h
library.o: global.h
quartz.o: global.h jumpgate.h
