#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

#include "wavelet_trie.hpp"
#include "unix_tools.hpp"

typedef annotate::cpp_int cpp_int;

int main(int argc, char** argv) {
    assert(argc > 2);
    std::string line, digit;
    std::vector<cpp_int> nums_ref;

    //environment variable arguments
    const char *step_char = std::getenv("STEP");
    size_t step = step_char ? atoi(step_char) : -1llu;
    const char *test = std::getenv("TEST");
    const char *njobs = std::getenv("NJOBS");
    size_t n_jobs = 1;
    if (njobs) {
        n_jobs = atoi(njobs);
        omp_set_num_threads(n_jobs);
    }

    annotate::WaveletTrie *wtr = NULL;
    Timer timer;
    double runtime = 0, readtime = 0;
    for (int f = 1; f < argc - 1; ++f) {
        std::vector<cpp_int> nums;
        if (step < -1llu) {
            nums.reserve(step);
        }
        std::ifstream fin(argv[f]);

        std::cout << "Compressing " << f << std::endl;
        timer.reset();
        while (std::getline(fin, line)) {
            std::istringstream sin(line);
            nums.emplace_back(0);
            while (std::getline(sin, digit, ',')) {
                annotate::bit_set(nums.back(), std::stoi(digit));
            }
            readtime += timer.elapsed();
            if (nums.size() == step) {
                if (test != NULL) {
                    nums_ref.reserve(nums_ref.size() + step);
                    nums_ref.insert(nums_ref.end(), nums.begin(), nums.end());
                }
                timer.reset();
                if (!wtr) {
                    wtr = new annotate::WaveletTrie(nums.begin(), nums.end());
                } else {
                    wtr->insert(annotate::WaveletTrie(nums.begin(), nums.end()));
                }
                runtime += timer.elapsed();
                std::cout << "." << std::flush;
                nums.clear();
            }
            timer.reset();
        }
        if (nums.size()) {
            if (test != NULL) {
                nums_ref.reserve(nums_ref.size() + nums.size());
                nums_ref.insert(nums_ref.end(), nums.begin(), nums.end());
            }
            timer.reset();
            if (!wtr) {
                wtr = new annotate::WaveletTrie(nums.begin(), nums.end());
            } else {
                wtr->insert(annotate::WaveletTrie(nums.begin(), nums.end()));
            }
            runtime += timer.elapsed();
            std::cout << std::endl;
            get_RAM();
            nums.clear();
        }
        std::cout << std::endl;
        fin.close();
    }

    if (wtr)
        std::cout << "Num edges:\t" << wtr->size() << std::endl;
    std::cout << "Reading:\t" << readtime << std::endl;
    std::cout << "Compressing:\t" << runtime << std::endl;

    if (test != NULL) {
        timer.reset();
        std::cout << "Uncompressing:\t" << std::flush;
        assert(wtr->size() == nums_ref.size());
#ifndef NPRINT
        wtr->print();
#endif
        for (size_t i = 0; i < nums_ref.size(); ++i) {
            if (nums_ref.at(i) != wtr->at(i)) {
                std::cerr << "Fail at " << i << "\n";
                std::cerr << nums_ref.at(i) << "\n" << wtr->at(i) << "\n";
                exit(1);
            }
        }
        std::cout << timer.elapsed() << std::endl;
    }
    if (wtr) {
        timer.reset();
        std::cout << "Serializing:\t" << std::flush;
        std::ofstream fout(argv[argc - 1]);
        size_t num_nodes = wtr->serialize(fout);
        std::cout << timer.elapsed() << std::endl;
        std::cout << "Num nodes:\t" << num_nodes << std::endl;
        delete wtr;
    }
    std::cout << "Done\n";
    return 0;
}


