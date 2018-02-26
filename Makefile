CXX=g++ -std=c++14 -fopenmp -march=native -DNPRINT -DNDEBUG -pthread -O3 #-pg -g -O2
HOME=/cluster/home/hmustafa
BLOOM=/cluster/home/hmustafa/metagraph/bloom/projects2014-metagenome/metagraph/external-libraries
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lpthread -lgmp -lboost_serialization

all: metannot pack_sd count_unique

wavelet_trie.o: wavelet_trie.hpp wavelet_trie.cpp Makefile
	$(CXX) $(CXXFLAGS) -c wavelet_trie.cpp

metannot.o: metannot_main.cpp Makefile
	$(CXX) $(CXXFLAGS) -c metannot_main.cpp

unix_tools.o: unix_tools.cpp unix_tools.hpp
	$(CXX) $(CXXFLAGS) -c unix_tools.cpp

metannot: Makefile wavelet_trie.o unix_tools.o metannot.o
	$(CXX) $(CXXFLAGS) -o metannot wavelet_trie.o metannot_main.o unix_tools.o $(LDFLAGS)

pack_sd: Makefile pack_sd.cpp
	$(CXX) $(CXXFLAGS) -o pack_sd pack_sd.cpp $(LDFLAGS)

count_unique: Makefile count_unique.cpp
	$(CXX) $(CXXFLAGS) -o count_unique count_unique.cpp $(LDFLAGS)

clean:
	rm metannot pack_sd count_unique *.o
