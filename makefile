SHELL = /bin/sh
CXX = mpicxx
CXXFLAGS = -Wall -O3 -g -std=c++11 -fopenmp -I../src -I/usr/local/include/ -I/home/becker/spwa4419/local/eigen/include/eigen3 -I/home/becker/spwa4419/Documents/spectra/include
binaries = Maxwell-TDSE.o Fields.o Coordinate_system.o Gas.o Maxwell.o Linear_algebra.o Complex.o TDSE.o

Maxwell-TDSE: $(binaries)
	$(CXX) $(CXXFLAGS) -o Maxwell-TDSE $(binaries)

Maxwell-TDSE.o: Maxwell-TDSE.cpp Fields.h Coordinate_system.h 
	$(CXX) $(CXXFLAGS) -c Maxwell-TDSE.cpp 

Gas.o: Coordinate_system.h Gas.h

Coordinate_system.o: Coordinate_system.h

Maxwell.o: Coordinate_system.h Fields.h Gas.h Maxwell.h

Fields.o: Coordinate_system.h Fields.h 

Linear_algebra.o: Complex.h Linear_algebra.h

Complex.o: Complex.h

TDSE.o: TDSE.h Complex.h Gas.h Coordinate_system.h Fields.h

.PHONY: clean
clean:
	rm -f $(binaries) *.o *.txt *.out 
