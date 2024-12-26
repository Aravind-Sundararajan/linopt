# Author: Aravind "Sunny" Sundararajan
# Should be working on windows with mingw,
# other options are the linux subshell in Win10, and of course just use linux good luck!
# USAGE: make all, or make -f path/to/makefile all
# or individual components:
# ------------------------
default: all

# convertFile
# findBounds
# findMean
# mplrs

# compiler

CC := gcc

CXX := g++

mpicxx=mpicc

# build flags, debug options:
#-DLRS_QUIET <- comment that guy for debugging
CFLAGS := -g -O3  -pthread -std=c++11 -Wall -DLRS_QUIET
#

SRCDIR := ./src/

LIBDIR := ./lib/

BINDIR := ./bin/

TSTDIR := ./test/

LRSDIR := ./lrslib-071/

INCLUDES := -I$(LIBDIR) -I$(LRSDIR) -I$(SRCDIR)

FASTOBJ = base.o util.o vec.o mat.o kvp.o polytope.o optim.o redund.o mapper.o

LRSLIB = lrsgmp.o lrsdriver.o lrslib.o

SHRLIB = lrsgmp-shr.o lrsdriver-shr.o lrslib1-shr.o lrslib2-shr.o

TESTS = vectest mattest kvptest optimtest

DRIVERS = f2r

#################################################################################
#                                     OBJ                                       #
#################################################################################

r2f: $(LRSOBJ) $(SRCDIR)rat2float.c
	$(CC) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)rat2float.o ${SRCDIR}rat2float.c -Wno-write-strings -lgmp
	$(CC) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)rat2float.o  -o $(BINDIR)r2f $(LRSOBJ) -lgmp

base.o: $(SRCDIR)base.cpp $(SRCDIR)base.h
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrslib.o $(LRSDIR)lrslib.c -lgmp
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsdriver.c -lgmp
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrsgmp.o $(LRSDIR)lrsgmp.c -lgmp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)base.o ${SRCDIR}base.cpp -Wno-write-strings -lgmp 
	$(CXX) $(CFLAGS) $(INCLUDES) -DGMP -c $(SRCDIR)base.cpp -o $(LIBDIR)base.o  -lgmp 

util.o: $(SRCDIR)util.cpp $(SRCDIR)util.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)util.cpp -o $(LIBDIR)util.o

kvp.o: $(SRCDIR)kvp.cpp $(SRCDIR)kvp.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)kvp.cpp -o $(LIBDIR)kvp.o

vec.o: $(SRCDIR)vec.cpp $(SRCDIR)vec.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)vec.cpp -o $(LIBDIR)vec.o

mat.o: $(SRCDIR)mat.cpp $(SRCDIR)mat.h
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrslib.o $(LRSDIR)lrslib.c -lgmp
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsdriver.c -lgmp
	$(CC)  ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LRSDIR)lrsgmp.o $(LRSDIR)lrsgmp.c -lgmp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)mat.o ${SRCDIR}mat.cpp -Wno-write-strings -lgmp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)mat.cpp -o $(LIBDIR)mat.o -lgmp

polytope.o: $(SRCDIR)polytope.cpp $(SRCDIR)polytope.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)polytope.cpp -o $(LIBDIR)polytope.o

optim.o: $(SRCDIR)optim.cpp $(SRCDIR)optim.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)optim.cpp -o $(LIBDIR)optim.o

redund.o: $(SRCDIR)redund.cpp $(SRCDIR)redund.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)redund.cpp -o $(LIBDIR)redund.o

mapper.o: $(SRCDIR)mapper.cpp $(SRCDIR)mapper.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRCDIR)mapper.cpp -o $(LIBDIR)mapper.o


obj: $(FASTOBJ)

cleanobj:
	cd $(LIBDIR) && rm -f *.o

#################################################################################
#                                     TESTS                                     #
#################################################################################

vectest: $(LRSOBJ) $(TSTDIR)vec_unit_test.cpp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)vectest.o ${TSTDIR}vec_unit_test.cpp -Wno-write-strings
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)vectest.o -o $(BINDIR)vectest $(LRSOBJ) -lgmp

mattest: $(LRSOBJ) $(TSTDIR)mat_unit_test.cpp mat.o
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)mattest.o ${TSTDIR}mat_unit_test.cpp -Wno-write-strings
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)mattest.o -o $(BINDIR)mattest $(LRSOBJ) -lgmp

kvptest: $(LRSOBJ) $(TSTDIR)kvp_unit_test.cpp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)kvptest.o ${TSTDIR}kvp_unit_test.cpp -Wno-write-strings
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)kvptest.o -o $(BINDIR)kvptest $(LRSOBJ) -lgmp

optimtest: $(LRSOBJ) $(TSTDIR)optim_unit_test.cpp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)optimtest.o ${TSTDIR}optim_unit_test.cpp -Wno-write-strings
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)optimtest.o -o $(BINDIR)optimtest $(LRSOBJ) -lgmp

polytopetest: $(LRSOBJ) $(TSTDIR)polytope_unit_test.cpp
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP -c -o $(LIBDIR)polytopetest.o ${TSTDIR}polytope_unit_test.cpp -Wno-write-strings
	$(CXX) ${CFLAGS} $(INCLUDES) -DGMP $(LRSDIR)lrslib.o $(LRSDIR)lrsdriver.o $(LRSDIR)lrsgmp.o $(LIBDIR)polytopetest.o -o $(BINDIR)polytopetest $(LRSOBJ) -lgmp


tests: $(TESTS)

cleantest:
	cd $(BINDIR) && rm -f ${TESTS} *.exe

#################################################################################
#                                   DRIVERS                                     #
#################################################################################

f2r: $(SRCDIR)f2r.c
	$(CC) ${CFLAGS} -DGMP $(INCLUDES) -o $(BINDIR)f2r $(SRCDIR)f2r.c -lgmp

# #stats:
# #	$(CXX) ${CFLAGS} $(BINDIR)stats.cpp -I/user/include/boost -o $(BINDIR)stats

drivers: $(DRIVERS)

cleandrivers:
	cd $(BINDIR) && rm -f ${DRIVERS} *.exe

#################################################################################
#                                   MPLRS                                       #
#################################################################################

mplrs:
	cd $(LRSDIR) && $(MAKE) lrs mplrs all-shared

cleanmplrs:
	cd $(LRSDIR) && $(MAKE) clean

clean: cleanobj cleandrivers cleantest cleanmplrs

all: clean mplrs obj drivers tests
