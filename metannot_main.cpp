#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>

#include "wavelet_trie.hpp"

typedef annotate::cpp_int cpp_int;

int main(int argc, char** argv) {
    assert(argc > 1);
    std::string line, digit;
    std::vector<cpp_int> nums_ref;

    //environment variable arguments
    const char *test = std::getenv("TEST");
    const char *njobs = std::getenv("NJOBS");
    size_t n_jobs = 1;
    if (njobs) {
        n_jobs = atoi(njobs);
        omp_set_num_threads(n_jobs);
    }

    annotate::WaveletTrie *wtr = NULL;
    for (int f = 1; f < argc - 1; ++f) {
        std::vector<cpp_int> nums;
        std::ifstream fin(argv[f]);

        std::cout << "Reading " << f << std::endl;
        while (std::getline(fin, line)) {
            std::istringstream sin(line);
            nums.emplace_back(0);
            while (std::getline(sin, digit, ',')) {
                annotate::bit_set(nums.back(), std::stoi(digit));
            }
        }
        fin.close();
        nums_ref.reserve(nums_ref.size() + nums.size());
        nums_ref.insert(nums_ref.end(), nums.begin(), nums.end());

        std::cout << "Compressing\n";
        const char *step_char = std::getenv("STEP");
        size_t step = step_char ? atoi(step_char) : nums.size();
        auto it = nums.begin();
        for (; it + step < nums.end(); it += step) {
            if (!wtr) {
                wtr = new annotate::WaveletTrie(it, it + step);
            } else {
                wtr->insert(annotate::WaveletTrie(it, it + step));
            }
            std::cout << "." << std::flush;
            if (it != nums.begin() && (static_cast<uint64_t>(it - nums.begin()) / step % 100 == 0)) {
                std::cout << std::endl;
            }
        }
        if (it != nums.end()) {
            if (!wtr) {
                wtr = new annotate::WaveletTrie(it, nums.end());
            } else {
                wtr->insert(annotate::WaveletTrie(it, nums.end()));
            }
            std::cout << "." << std::flush;
        }
        std::cout << std::endl;
    }

    if (test != NULL) {
        std::cout << "Uncompressing\n";
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
    }
    if (wtr) {
        std::cout << "Serializing\n";
        std::ofstream fout(argv[argc - 1]);
        wtr->serialize(fout);
        delete wtr;
    }
    std::cout << "Done\n";
    return 0;
}


