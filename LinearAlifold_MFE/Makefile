
CC=g++
DEPS= src/Linearalifold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3
.PHONY : clean linearalifold
objects=bin/linearalifold

linearalifold: src/Linearalifold.cpp $(DEPS) 
		mkdir -p bin
		$(CC) src/Linearalifold.cpp $(CFLAGS) -Dlv -Dis_candidate_list -o bin/linearalifold 
clean:
	-rm $(objects)