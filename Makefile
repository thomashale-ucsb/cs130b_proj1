CXX=g++

CXXFLAGS = -std=c++11 -Wall -Wextra -g

#-Wno-unused-parameter -Wno-unused-private-field

BINARIES=cvGen convexTest

all: ${BINARIES}

cvGen: main.o quickHull.o
	${CXX} $^ -o $@

convexTest: cvTest.o quickHull.o
	${CXX} $^ -o $@

clean:
	/bin/rm -f ${BINARIES} *.o