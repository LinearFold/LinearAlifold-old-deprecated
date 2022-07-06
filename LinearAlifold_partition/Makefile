
CC=g++
DEPS=src/bpp.cpp src/linearalifold_p.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
CFLAGS=-std=c++11 -O3
.PHONY : clean linearalifold_p
objects=bin/linearalifold_p

linearalifold_p: src/linearalifold_p.cpp $(DEPS) 
		mkdir -p bin
		$(CC) src/linearalifold_p.cpp $(CFLAGS) -Dlpv -o bin/linearalifold_p 
clean:
	-rm $(objects)
