all:	
	g++ -w varhmm.cpp -O3 -std=c++0x -pthread -o varhmm
	g++ varhmmtrain.cpp -O3 -o varhmmtrain
