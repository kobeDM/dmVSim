WHERE=.
INCLUDE=${WHERE}/include
LD_PATH = ${WHERE}/bin
SRC = ${WHERE}/src

ROOTFLAGS = $(shell root-config --cflags)
CXXFLAGS = $(ROOTFLAGS) -Wno-deprecated-declarations -fPIC
ROOTLIBS = $(shell root-config --libs)

RLIB = -I${INCLUDE} -I./ShUtil

TARGET = LNATrand SimDMFlux SimNuclRecoil SimDMVelocity

all:	${TARGET}

LNATrand: ${SRC}/LNATrand.o ${SRC}/CRDMFunc.o
	g++ ${CXXFLAGS} -o LNATrand ${SRC}/LNATrand.o ${SRC}/CRDMFunc.o ${ROOTLIBS} -lgcc
	mv ./LNATrand ${LD_PATH}

SimDMVelocity: ${SRC}/SimDMVelocity.o ${SRC}/CRDMFunc.o
	g++ ${CXXFLAGS} -o SimDMVelocity ${SRC}/SimDMVelocity.o ${SRC}/CRDMFunc.o ${ROOTLIBS} -lgcc
	mv ./SimDMVelocity ${LD_PATH}

SimDMFlux: ${SRC}/SimDMFlux.o ${SRC}/CRDMFunc.o
	g++ ${CXXFLAGS} -o SimDMFlux ${SRC}/SimDMFlux.o ${SRC}/CRDMFunc.o ${ROOTLIBS} -lgcc
	mv ./SimDMFlux ${LD_PATH}

SimNuclRecoil: ${SRC}/SimNuclRecoil.o ${SRC}/CRDMFunc.o
	g++ ${CXXFLAGS} -o SimNuclRecoil ${SRC}/SimNuclRecoil.o ${SRC}/CRDMFunc.o ${ROOTLIBS} -lgcc
	mv ./SimNuclRecoil ${LD_PATH}

%.o: %.cc
	g++ ${CXXFLAGS} ${RLIB} -O4 -c -o $@ $<

clean:
	rm -f ${SRC}/*.o
	rm -f *~
