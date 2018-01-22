#include<iostream>
#include<fstream>
#include<omp.h>
#include "wavelet_trie_pointer.hpp"

int main(int argc, char** argv) {

    const char* njob = std::getenv("NJOBS");
    size_t n_jobs=1;
    if (njob != NULL) {
        n_jobs = atol(njob);
    }
    std::cout << "njobs: " << n_jobs << "\n";
    std::cout << "size wtr: " << sizeof(WaveletTrie::WTR) << "\n";
    std::cout << "size Node: " << sizeof(WaveletTrie::Node) << "\n";
    omp_set_num_threads(n_jobs);
    if (argc < 6) {
        std::cerr << "Too few arguments\n";
        return 1;
    }
    size_t batch = atol(argv[1]); //batch size for storing in RAM
    size_t cbatch = atol(argv[2]); //batch size for wavelet trie construction
    bool reconst = atoi(argv[3]);
    //std::vector<WaveletTrie::WTR*> wtrs;
    //std::vector<std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*>* > ivs;
    std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*> civ;
    WaveletTrie::intvec_t riv;
    //std::vector<WaveletTrie::WTR*> wtr(2);
    //wtr[0] = new WaveletTrie::WTR();
    WaveletTrie::WTR wtr;
    for (int i=4;i<argc-1;++i) {
        std::ifstream pfile(argv[i]);
        std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*> iv;
        if (pfile.good()) {
            std::cout << "Reading\n";
            std::cout << "Constructing";
            if (batch > 0 || cbatch > 0)
                std::cout << " in parts";
            std::cout << "\n";
            while (!pfile.eof()) {
                for (size_t i = 0; !pfile.eof(); i += cbatch) {
                    civ = WaveletTrie::construct(pfile, cbatch, MAXNUM);
                    if (reconst) {
                        riv.insert(riv.end(), civ.first->begin(), civ.first->end());
                    }
                    //wtr[1] = new WaveletTrie::WTR(civ.first, civ.second->size(), civ.first->size());
                    //auto temp = new WaveletTrie::WTR(WaveletTrie::WTR::merge(wtr.begin(), wtr.end()));
                    //delete wtr[0];
                    //delete wtr[1];
                    #pragma omp parallel
                    #pragma omp single nowait
                    wtr.append(WaveletTrie::WTR(civ.first, civ.second->size(), civ.first->size()));
                    delete civ.first;
                    delete civ.second;
                    if (i % 10000 == 0)
                        std::cout << "." << std::flush;
                    if (i % 1000000 == 0)
                        std::cout << "\n" << std::flush;
                    //wtr[0] = temp;
                }
            /*
                for (size_t i=0;!pfile.eof() && (batch>0 ? i<batch : true);i+=cbatch) {
                    ivs.push_back(new std::pair<WaveletTrie::intvec_t*, WaveletTrie::intvec_t*>(WaveletTrie::construct(pfile, cbatch, MAXNUM)));
                    if (reconst) {
                        riv.insert(riv.end(), ivs.back()->first->begin(), ivs.back()->first->end());
                    }
                }
                std::vector<WaveletTrie::WTR*> wtrs_new(ivs.size());
                //size_t k=wtrs.size();
                //wtrs.resize(wtrs.size()+ivs.size());
                #pragma omp parallel
                #pragma omp single
                {
                    #pragma omp taskloop
                    for (size_t i=0;i<ivs.size();++i) {
                        wtrs_new[i] = new WaveletTrie::WTR(ivs[i]->first, ivs[i]->second->size(), ivs[i]->first->size());
                        if (ivs[i]->first) {
                            delete ivs[i]->first;
                        }
                        delete ivs[i]->second;
                        delete ivs[i];
                        std::cout << "." << std::flush;
                    }
                }
                wtrs.insert(wtrs.end(), std::make_move_iterator(wtrs_new.begin()), std::make_move_iterator(wtrs_new.end()));
                ivs.clear();
            */
            }
            std::cout << "\n";
        } else {
            std::cerr << "Bad input file\n";
            return 1;
        }
        pfile.close();
        std::cout << "\n";
    }
    /*
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
        wtr = new WaveletTrie::WTR(WaveletTrie::WTR::merge(wtrs.begin(), wtrs.end()));
        delete wtrs[0];
        wtrs[0] = wtr;
    }
    */
    if (reconst) {
        //assert(wtrs[0]->reconstruct(riv));
        assert(wtr.reconstruct(riv));
        std::cout << "Passed merged test\n";
    }
    std::ofstream ofile(argv[argc-1]);
    if (ofile.good()) {
        std::cout << "Serializing\n";
        wtr.serialize(ofile);
        //wtrs[0]->serialize(ofile);
    } else {
        std::cerr << "Bad output file\n";
        return 1;
    }
    ofile.close();
    std::cout << "\nCleaning up\n";
    //delete wtr[0];
    //for (size_t i=0;i<wtrs.size();++i) {
    //    delete wtrs[i];
    //}

    return 0;
}

