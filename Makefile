CXX=g++ -std=c++14 -fopenmp -g -pg -march=native -O2 -DNPRINT -DNDEBUG #-O3
HOME=/cluster/home/hmustafa
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lboost_serialization -lpthread -lgmp

all: metannot

wavelet_trie.o: wavelet_trie.hpp wavelet_trie.cpp Makefile
	$(CXX) $(CXXFLAGS) -c wavelet_trie.cpp

metannot: Makefile metannot_main.cpp wavelet_trie.o
	$(CXX) $(CXXFLAGS) -o metannot wavelet_trie.o metannot_main.cpp $(LDFLAGS)

clean:
	rm metannot
