CXX=g++ -std=c++14 -fopenmp -march=native -DNPRINT -O2 -DNDEBUG -pthread -O3 #-g #-O2 -pg
HOME=/cluster/home/hmustafa
BLOOM=/cluster/home/hmustafa/metagraph/bloom/projects2014-metagenome/metagraph/external-libraries
CXXFLAGS=-Wall -I$(HOME)/local/include -L$(HOME)/local/lib
LDFLAGS=-lsdsl -lpthread -lgmp -lboost_serialization

all: metannot pack_sd count_unique unpack_commas metannot_merge

wavelet_trie.o: wavelet_trie.hpp wavelet_trie.cpp
	$(CXX) $(CXXFLAGS) -c wavelet_trie.cpp

metannot.o: metannot_main.cpp
	$(CXX) $(CXXFLAGS) -c metannot_main.cpp

unix_tools.o: unix_tools.cpp unix_tools.hpp
	$(CXX) $(CXXFLAGS) -c unix_tools.cpp

utils.o: utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

metannot: wavelet_trie.o unix_tools.o metannot.o utils.o
	$(CXX) $(CXXFLAGS) -o metannot wavelet_trie.o metannot_main.o unix_tools.o utils.o $(LDFLAGS)

metannot_merge: wavelet_trie.o metannot_merge.cpp utils.o
	$(CXX) $(CXXFLAGS) -o metannot_merge wavelet_trie.o utils.o metannot_merge.cpp $(LDFLAGS)

metannot_extract: wavelet_trie.o metannot_extract.cpp utils.o
	$(CXX) $(CXXFLAGS) -o metannot_extract wavelet_trie.o utils.o metannot_extract.cpp $(LDFLAGS)

pack_sd: pack_sd.cpp
	$(CXX) $(CXXFLAGS) -o pack_sd pack_sd.cpp $(LDFLAGS)

count_unique: count_unique.cpp
	$(CXX) $(CXXFLAGS) -o count_unique count_unique.cpp $(LDFLAGS)

unpack_commas: unpack_commas.cpp
	$(CXX) $(CXXFLAGS) -o unpack_commas unpack_commas.cpp $(LDFLAGS)

clean:
	rm metannot pack_sd count_unique *.o
