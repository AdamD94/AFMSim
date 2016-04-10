CFLAGS=-std=c++11 -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_system -Ofast

AFMSim: AFMSim.o
	g++ AFMSim.o $(CFLAGS) -o AFMSim

AFMSim.o: AFMSim.cpp AFMSim.h Var.h gnuplot-iostream.h 
	g++ $(CFLAGS) -c AFMSim.cpp

clean:
	rm *o