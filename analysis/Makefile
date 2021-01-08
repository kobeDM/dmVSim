TARGETS1=fillTree
TARGETS2=runAnalysis
TARGETS3=drawHist
TARGETS4=fillHist
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)
INCLUDE = -I ./include

.PHONY: all

all: ${TARGETS1} ${TARGETS2} ${TARGETS3} ${TARGETS4}

${TARGETS1}: ./src/${TARGETS1}.cxx
	g++ ./src/${TARGETS1}.cxx -o ./bin/${TARGETS1} ${CXXLIBS} $(CXXFLAGS) 

${TARGETS2}: ./src/${TARGETS2}.cxx
	g++ ./src/${TARGETS2}.cxx -o ./bin/${TARGETS2} ${CXXLIBS} $(CXXFLAGS) ${INCLUDE}

${TARGETS3}: ./src/${TARGETS3}.cxx
	g++ ./src/${TARGETS3}.cxx -o ./bin/${TARGETS3} ${CXXLIBS} $(CXXFLAGS)

${TARGETS4}: ./src/${TARGETS4}.cxx
	g++ ./src/${TARGETS4}.cxx -o ./bin/${TARGETS4} ${CXXLIBS} $(CXXFLAGS) ${INCLUDE}

.PHONY: clean

clean:
	@rm -f  ./bin/* \
        && echo "Make Clean --> all ok!"
