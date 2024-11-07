all: main

main: main.cpp BRDF.o Geometry.o WaveSim.o
	g++ -std=c++11 -I include main.cpp BRDF.o Geometry.o WaveSim.o -O3 -o main -lpthread

BRDF.o: BRDF.h BRDF.cpp
	g++ -std=c++11 -c -I include BRDF.cpp -O3 -o BRDF.o
	
Geometry.o: Geometry.h Geometry.cpp
	g++ -std=c++11 -c -I include Geometry.cpp -O3 -o Geometry.o

WaveSim.o: WaveSim.h WaveSim.cpp
	g++ -std=c++11 -c -I include WaveSim.cpp -O3 -o WaveSim.o

clean:
	rm -f main
	rm -f *.o