CXX=g++ -std=c++11 -fopenmp -g -O3
HOME=~
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib -I$(HOME)/Applications/DYNAMIC/include -I$(HOME)/Applications/DYNAMIC/include/internal
LDFLAGS=-lsdsl -lboost_serialization -lpthread

all: rrr

rrr: rrr_wt.cpp wavelet_trie_pointer.hpp
	$(CXX) $(CXXFLAGS) -o metannot metannot_main.cpp $(LDFLAGS)

