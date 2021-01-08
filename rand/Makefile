GCC = g++       

WHERE=.
INCLUDE=${WHERE}/include
LD_PATH = ${WHERE}/bin
SRC = ${WHERE}/src

ROOTFLAGS = $(shell root-config --cflags)
CXXFLAGS = $(ROOTFLAGS) -std=c++11 -pthread -Wno-deprecated-declarations -m64
ROOTLIBS = $(shell root-config --libs)

RLIB = -I${INCLUDE} -I./ShUtil -I${ROOTSYS}/include

TARGET = LNATrand CRDMrand CRDMFLXrand

all:	${TARGET}

LNATrand: ${SRC}/LNATrand.cc
	g++ ${RLIB} $^ -o LNATRand ${CXXFLAGS} ${ROOTLIBS} -I./ShUtil
	mv ./LNATRand ${LD_PATH}

CRDMrand: ${SRC}/CRDMrand.cc
	g++ ${RLIB} $^ -o CRDMrand ${CXXFLAGS} ${ROOTLIBS} -I./ShUtil
	mv ./CRDMrand ${LD_PATH}

CRDMFLXrand: ${SRC}/CRDMFLXrand.cc
	g++ ${RLIB} $^ -o CRDMFLXrand ${CXXFLAGS} ${ROOTLIBS} -I./ShUtil
	mv ./CRDMFLXrand ${LD_PATH}

clean:
	rm -f *.o
	rm -f *~
