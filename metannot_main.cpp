#include<iostream>
#include<fstream>
#include<omp.h>
#include "wavelet_trie_pointer.hpp"
#include "array_int.hpp"

int main(int argc, char** argv) {

    array_int::array_int tester;
    bit_test(tester, 0);

    const char* njob = std::getenv("NJOBS");
    size_t n_jobs=1;
    if (njob != NULL) {
        n_jobs = atol(njob);
    }
    std::cout << "njobs: " << n_jobs << "\n";
    omp_set_num_threads(n_jobs);
    if (argc < 6) {
        std::cerr << "Too few arguments\n";
        return 1;
    }
    size_t batch = atol(argv[1]); //batch size for storing in RAM
    size_t cbatch = atol(argv[2]); //batch size for wavelet trie construction
    bool reconst = atoi(argv[3]);
    std::vector<WaveletTrie::WTR*> wtrs;
    std::vector<std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*>* > ivs;
    WaveletTrie::intvec_t riv;
    for (int i=4;i<argc-1;++i) {
        std::ifstream pfile(argv[i]);
        bool debug=(std::getenv("DEBUG") != NULL);
        std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*> iv;
        if (pfile.good()) {
            std::cout << "Reading\n";
            std::cout << "Constructing";
            if (batch > 0 || cbatch > 0)
                std::cout << " in parts";
            std::cout << "\n";
            while (!pfile.eof()) {
                for (size_t i=0;!pfile.eof() && (batch>0 ? i<batch : true);i+=cbatch) {
                    ivs.push_back(new std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*>(WaveletTrie::construct(pfile, cbatch, MAXNUM, debug)));
                    if (reconst) {
                        riv.insert(riv.end(), ivs.back()->first->begin(), ivs.back()->first->end());
                    }
                }
                size_t k=wtrs.size();
                wtrs.resize(wtrs.size()+ivs.size());
                #pragma omp parallel
                #pragma omp single
                {
                    #pragma omp taskloop
                    for (size_t i=0;i<ivs.size();++i) {
                        wtrs[k+i] = new WaveletTrie::WTR(ivs[i]->first, 0, ivs[i]->second->size(), ivs[i]->first->size());
                        delete ivs[i]->second;
                        delete ivs[i];
                        std::cout << "." << std::flush;
                    }
                }
                ivs.clear();
            }
            std::cout << "\n";
        } else {
            std::cerr << "Bad input file\n";
            return 1;
        }
        pfile.close();
        std::cout << "\n";
    }
    if (wtrs.size() > 1) {
        std::cout << "Merging " << wtrs.size() << "\n";
        for (size_t i=0;i<wtrs.size();++i) {
            if (wtrs[i] == NULL) {
                std::cerr << "Fail " << i << "\n";
                assert(false);
            }
        }
        WaveletTrie::WTR* wtr;
        #pragma omp parallel
        #pragma omp single nowait
        wtr = new WaveletTrie::WTR(WaveletTrie::merge(wtrs.begin(), wtrs.end()));
        delete wtrs[0];
        wtrs[0] = wtr;
        if (reconst) {
            assert(wtrs[0]->reconstruct(riv));
            std::cout << "Passed merged test\n";
        }
    }
    std::ofstream ofile(argv[argc-1]);
    if (ofile.good()) {
        std::cout << "Serializing\n";
        wtrs[0]->serialize(ofile);
    } else {
        std::cerr << "Bad output file\n";
        return 1;
    }
    ofile.close();
    std::cout << "\nCleaning up\n";
    for (size_t i=0;i<wtrs.size();++i) {
        delete wtrs[i];
    }

    return 0;
}

