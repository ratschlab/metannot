#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>

#include "wavelet_trie.hpp"
#include "unix_tools.hpp"


//TODO: replace getenv with real argument parsing

typedef annotate::cpp_int cpp_int;

void serialize_vector_(std::ostream &out, const std::vector<cpp_int> &nums, size_t n_cols = 0) {
    n_cols += 64;
    for (auto it = nums.begin(); it != nums.end(); ++it) {
        size_t size = annotate::serialize(out, *it) + 8;
        if (size * 8 < n_cols) {
            char *zeros = (char*)malloc((n_cols + 7 - (size * 8)) >> 3);
            memset(zeros, 0, (n_cols + 7 - (size * 8)) >> 3);
            out.write(zeros, (n_cols + 7 - (size * 8)) >> 3);
            free(zeros);
        }
    }
}

uint64_t deserializeNumber(std::istream &in) {
    uint64_t n = 0;
    for (size_t i = 0; i < sizeof(n); ++i) {
        n = (n << 8) | in.get();
    }
    return n;
}

int main(int argc, char** argv) {
    if (argc <= 2) {
        std::cerr << "ERROR: please pass in at least two files: an input and an output." << std::endl;
        exit(1);
    }
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
    size_t seed = 0;

    //const char *dump_raw = std::getenv("DUMP");
    /*
    std::ofstream dout;
    if (dump_raw) {
        dump_cols = atoi(dump_raw);
        dout.open(std::string(argv[argc - 1]) + ".raw");
    }
    *///size_t dump_cols = 0;

    const char *read_comma = std::getenv("COMMA");

    const char *strmap = std::getenv("MAP");

    const char *shuf_seed = std::getenv("SHUF_SEED");
    if (shuf_seed)
        seed = atoi(shuf_seed);

    size_t set_bits = 0;
    size_t total_bits = 0;

    annotate::WaveletTrie *wtr = NULL;


    double runtime = 0;
    double readtime = 0;
    //double dumptime = 0;
    size_t num_rows = 0;
    for (int f = 1; f < argc - 1; ++f) {
        std::vector<cpp_int> nums;
        if (step < -1llu) {
            nums.reserve(step);
        }
        std::ifstream fin(argv[f]);
        if (!fin.good()) {
            std::cerr << "WARNING: file " << argv[f] << " bad." << std::endl;
            exit(1);
        }
        if (!read_comma && !strmap) {
            num_rows = deserializeNumber(fin);
        }
        std::unordered_map<std::string, std::set<size_t>> string_map;
        if (strmap) {
            boost::archive::binary_iarchive iarch(fin);
            size_t num_elements;
            std::cout << "Loading input" << std::endl;
            iarch & string_map;
            iarch & num_elements;
        }
        auto strmap_it = string_map.begin();

        std::cout << "Compressing " << argv[f] << std::endl;
        while (true) {
            if (read_comma) {
                //read from text index list
                if (!std::getline(fin, line)) {
                    break;
                }
                std::istringstream sin(line);
                nums.emplace_back(0);
                while (std::getline(sin, digit, ',')) {
                    annotate::bit_set(nums.back(), std::stoi(digit));
                    set_bits++;
                }
            } else if (strmap) {
                //read from serialized index set
                if (strmap_it == string_map.end())
                    break;
                nums.emplace_back(0);
                for (auto &index : strmap_it->second) {
                    annotate::bit_set(nums.back(), index);
                    set_bits++;
                }
                ++strmap_it;
            } else {
                if (!num_rows)
                    break;
                size_t size = deserializeNumber(fin);
                std::vector<uint64_t> row(size);
                for (auto it = row.begin(); it != row.end(); ++it) {
                    *it = deserializeNumber(fin);
                    set_bits += __builtin_popcountll(*it);
                }
                if (shuf_seed) {
                    //shuffle
                    std::srand(seed);
                    std::random_shuffle(
                            reinterpret_cast<uint8_t*>(row.data()),
                            reinterpret_cast<uint8_t*>(row.data() + row.size()),
                            [&](int i){ return std::rand()%i; });
                }
                total_bits += row.size() * 64;
                nums.emplace_back(0);
                mpz_import(nums.back().backend().data(), row.size(), -1, sizeof(row[0]), 0, 0, &row[0]);
#ifndef NDEBUG
                if (nums.back() != 0) {
                    for (size_t i = 0; i < row.size(); ++i) {
                        assert(((nums.back() >> (i * 64)) & -1llu) == row[i]);
                    }
                }
#endif
                num_rows--;
            }
            if (nums.size() == step) {
                if (test != NULL) {
                    nums_ref.reserve(nums_ref.size() + step);
                    nums_ref.insert(nums_ref.end(), nums.begin(), nums.end());
                }
                if (!wtr) {
                    wtr = new annotate::WaveletTrie(nums.begin(), nums.end());
                } else {
                    wtr->insert(annotate::WaveletTrie(nums.begin(), nums.end()));
                }
                std::cout << "." << std::flush;
                /*
                if (dump_raw) {
                    timer.reset();
                    serialize_vector_(dout, nums, dump_cols);
                    dumptime += timer.elapsed();
                }
                */
                nums.clear();
            }
        }
        if (nums.size()) {
            if (test != NULL) {
                nums_ref.reserve(nums_ref.size() + nums.size());
                nums_ref.insert(nums_ref.end(), nums.begin(), nums.end());
            }
            if (!wtr) {
                wtr = new annotate::WaveletTrie(nums.begin(), nums.end());
            } else {
                wtr->insert(annotate::WaveletTrie(nums.begin(), nums.end()));
            }
            std::cout << std::endl;
            get_RAM();
            /*
            if (dump_raw) {
                serialize_vector_(dout, nums, dump_cols);
            }
            */
            nums.clear();
        }
        std::cout << std::endl;
        fin.close();
    }

    std::cout << "Times:" << std::endl;
    std::cout << "Reading:\t" << readtime << std::endl;
    std::cout << "Compressing:\t" << runtime << std::endl;

    if (test != NULL) {
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
    }
    if (wtr) {
        std::cout << "Serializing:\t" << std::flush;
        std::ofstream fout(argv[argc - 1], std::ofstream::binary | std::ofstream::ate);
        if (!fout.good()) {
            std::cerr << "ERROR: bad file " << argv[argc - 1] << std::endl;
            exit(1);
        }
        auto stats = wtr->serialize(fout);
        std::cout << "Input:" << std::endl;
        std::cout << "Num edges:\t" << wtr->size() << std::endl;
        std::cout << "Total bits:\t" << total_bits << std::endl;
        std::cout << "Set bits:\t" << set_bits << std::endl;
        std::cout << "Wavelet trie:" << std::endl;
        std::cout << "Num leaves:\t" << stats << std::endl;
        std::cout << "Num bytes:\t" << fout.tellp() << std::endl;
        std::cout << "Bits/edge:\t" << (double)fout.tellp() * 8.0 / (double)wtr->size() << std::endl;
        if (shuf_seed) {
            std::cout << "Shuffle seed:\t" << seed << std::endl;
        }
        fout.close();
        delete wtr;
    }
    /*
    if (dump_raw) {
        dout.close();
        std::cout << "Raw dump:\t" << dumptime << std::endl;
    }
    */
    std::cout << "Done\n";
    return 0;
}


