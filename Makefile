AFMSim: Main.o
	g++ Main.o -std=c++11 -lboost_thread -lboost_filesystem  -lboost_system -O2 -o AFMSim

Main.o: Main.cpp Main.h Class_Definitions.h
	g++ -std=c++11 -c -lboost_thread -lboost_filesystem  -lboost_system -O2 Main.cpp

clean:
	rm *.o 