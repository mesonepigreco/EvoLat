# Numpy location
GSLFLAGS = -lgsl -lgslcblas


lib/parse_input.o: src/parse_input.cpp src/parse_input.hpp
	g++ -O3 -fPIC -c src/parse_input.cpp -o lib/parse_input.o

lib/genome.o: src/genome.cpp src/genome.hpp
	g++ -O3 -fPIC -c src/genome.cpp -o lib/genome.o

lib/wator.o: src/waTor.cpp src/waTor.hpp
	g++ -O3 -fPIC -c src/waTor.cpp -o lib/wator.o
