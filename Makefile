All: AFMSim.o AFMSim clean 

AFMSim: AFMSim.o
	g++ AFMSim.o -std=c++11 -lboost_thread -lboost_filesystem -lboost_iostreams  -lboost_system -O2 -o AFMSim

AFMSim.o: AFMSim.cpp AFMSim.h Var.h gnuplot-iostream.h 
	g++ -std=c++11 -c -lboost_thread -lboost_filesystem -lboost_system -O2 AFMSim.cpp

clean:
	rm *.o
	