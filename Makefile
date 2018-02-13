CXX=g++ -std=c++14 -fopenmp -march=native -O2 -DNPRINT -DNDEBUG -O3 #-pg
HOME=/cluster/home/hmustafa
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lpthread -lgmp

all: metannot

wavelet_trie.o: wavelet_trie.hpp wavelet_trie.cpp Makefile
	$(CXX) $(CXXFLAGS) -c wavelet_trie.cpp

metannot.o: metannot_main.cpp Makefile
	$(CXX) $(CXXFLAGS) -c metannot_main.cpp

unix_tools.o: unix_tools.cpp unix_tools.hpp
	$(CXX) $(CXXFLAGS) -c unix_tools.cpp

metannot: Makefile metannot.o wavelet_trie.o unix_tools.o
	$(CXX) $(CXXFLAGS) -o metannot wavelet_trie.o metannot_main.o unix_tools.o $(LDFLAGS)

clean:
	rm metannot *.o
