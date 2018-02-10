#ifndef __WAVELET_TRIE___
#define __WAVELET_TRIE___

#include <iostream>
#include <fstream>
#include <sdsl/wavelet_trees.hpp>
#include <omp.h>
#include <algorithm>
#include <boost/multiprecision/gmp.hpp>

namespace annotate {
    typedef boost::multiprecision::mpz_int cpp_int;
    typedef cpp_int annot_t;
    typedef sdsl::bit_vector bv_t;
    typedef sdsl::rrr_vector<> rrr_t;
    typedef sdsl::sd_vector<> sd_t;

    typedef cpp_int alpha_t;
    typedef bv_t beta_t;
    typedef beta_t::rank_1_type rank1_t;
    typedef beta_t::rank_0_type rank0_t;

    //typedef boost::archive::binary_oarchive oarchive_t;
    //typedef boost::archive::binary_iarchive iarchive_t;

    bool bit_test(const cpp_int &a, const size_t &col);

    void bit_set(mpz_t &a_d, const size_t &col);

    void bit_set(cpp_int &a, const size_t &col);

    void bit_unset(cpp_int &a, const size_t &col);

    size_t next_bit(const cpp_int &a, const size_t &col);

    void clear_after(mpz_t &a_d, const size_t &col);

    void clear_after(cpp_int &a, const size_t &col);

    size_t msb(const cpp_int &a);

    size_t lsb(const cpp_int &a);

    template <typename Vector>
    bv_t insert_range(const Vector &target, const Vector &source, size_t i = 0);

    class Prefix {
        public:
            Prefix() { }
            Prefix(size_t col, bool allequal) : col(col), allequal(allequal) { }
            size_t col = -1llu;
            bool allequal = true;
    };

    class WaveletTrie {
        public:
            WaveletTrie();

            template <class Iterator>
            WaveletTrie(Iterator row_begin, Iterator row_end);

            //destructor
            ~WaveletTrie();

            cpp_int at(size_t i);

            size_t size();

            void insert(WaveletTrie&& wtr, size_t i = -1llu);

            void print();

            template <typename T>
            void insert(const T &a, size_t i = -1llu);

        private:
            class Node;
            class Leaf;
            Node* root;
    };

    class WaveletTrie::Node {
        friend class WaveletTrie;
        public:
            //empty constructor
            Node() { }

            Node(const size_t count);

            Node(const cpp_int &alpha, const size_t count);

            //destructor
            ~Node();

            //Copy constructor
            Node(const Node &that);

            //size_t count() { return 0; }

            //Constructor from array
            //col is number of bits to trim off left
            template <class Iterator>
            Node(const Iterator &row_begin, const Iterator &row_end,
                    const size_t &col = 0, Prefix prefix = Prefix());

            //swap
            void swap(Node&& that);

            size_t size() { return beta_.size(); }

            void fill_left(bool rightside = false);

            bool check(bool ind);

        protected:
            alpha_t alpha_ = 1;
            beta_t beta_;
            rank1_t rank1_;
            rank0_t rank0_;
            Node *child_[2] = {NULL, NULL};
            bool all_zero = false;
            size_t popcount = 0;

        private:

            template <class Iterator>
            static size_t next_different_bit_(const Iterator &a, const Iterator &b,
                    const size_t col = 0, size_t next_col = -1llu);

            static size_t next_different_bit_alpha(cpp_int a, cpp_int b);

            template <class Iterator>
            static Prefix longest_common_prefix(
                    const Iterator &row_begin, const Iterator &row_end, const size_t &col);

            template <class Vector>
            void set_beta_(const Vector &bv);

            int move_label_down_(size_t length);

            static bool overlap_prefix_(Node *curnode, Node *othnode);

            static void merge_beta_(Node *curnode, Node *othnode, size_t i = -1llu);

            bool is_leaf();

            template <class Container>
            static void push_child(Container &nodes, Node *curnode, Node *othnode, bool ind, const size_t i, std::string path);

    };

};

#endif // __WAVELET_TRIE___
