CXX=g++ -std=c++11 -fopenmp -g
HOME=/cluster/home/hmustafa
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lboost_serialization -lpthread

all: metannot

metannot: metannot_main.cpp wavelet_trie_pointer.hpp array_int.hpp
	$(CXX) $(CXXFLAGS) -o metannot metannot_main.cpp $(LDFLAGS)

