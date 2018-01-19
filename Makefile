CXX=g++ -std=c++14 -fopenmp -g -pg -march=native -O2
HOME=/cluster/home/hmustafa
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lboost_serialization -lpthread -lgmp

all: metannot

metannot: metannot_main.cpp wavelet_trie_pointer.hpp array_int.hpp
	$(CXX) $(CXXFLAGS) -o metannot metannot_main.cpp $(LDFLAGS)

clean:
	rm metannot
