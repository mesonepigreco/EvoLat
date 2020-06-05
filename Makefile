# Numpy location
GSLFLAGS = -lgsl -lgslcblas
CPPFLAGS = -std=c++11
LIBS=`pkg-config --libs libconfig++` -lm
GPP=g++

bin/wator.exe: src/main.cpp lib/parse_input.o lib/genome.o lib/wator.o
	$(GPP) -O3 -fPIC $(CPPFLAGS) src/main.cpp  -o bin/wator.exe lib/wator.o lib/parse_input.o lib/genome.o $(LIBS) 


lib/parse_input.o: src/parse_input.cpp src/parse_input.hpp
	$(GPP) -O3 -fPIC $(CPPFLAGS) -c src/parse_input.cpp -o lib/parse_input.o

lib/genome.o: src/genome.cpp src/genome.hpp 
	$(GPP) -O3 -fPIC $(CPPFLAGS)  -c src/genome.cpp -o lib/genome.o

lib/wator.o: src/waTor.cpp src/waTor.hpp
	$(GPP) -O3 -fPIC $(CPPFLAGS)  -c src/waTor.cpp -o lib/wator.o
