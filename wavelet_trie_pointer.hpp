#ifndef __WAVELET_TRIE__
#define __WAVELET_TRIE__

#include<iostream>
#include<fstream>
#include<sdsl/wavelet_trees.hpp>
#include<stack>
#include<vector>
#include<tuple>
#include<unordered_map>
#include<unordered_set>
#include<boost/functional/hash.hpp>
#include<map>
#include<omp.h>
#include<algorithm>
#include<parallel/algorithm>
#include "array_int.hpp"
#define MAXNUM (std::numeric_limits<std::size_t>::max())
#define TASKMIN 100

#include <boost/multiprecision/cpp_int.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include<functional>
#include<type_traits>


bool debug=false;


namespace WaveletTrie {
    typedef boost::multiprecision::cpp_int cpp_int;
    /*
    typedef bool (*bittest_t)(const cpp_int&, uint32_t i);
    bittest_t bit_test = boost::multiprecision::bit_test;
    typedef cpp_int& (*bitset_t)(cpp_int&, uint32_t i);
    bitset_t bit_set = boost::multiprecision::bit_set;
    bitset_t bit_unset = boost::multiprecision::bit_unset;
    typedef uint32_t (*sigbit_t)(const cpp_int&);
    sigbit_t msb = boost::multiprecision::msb;
    sigbit_t lsb = boost::multiprecision::lsb;
    */

    typedef array_int::array_int annot_t;
    //typedef cpp_int annot_t;
    typedef std::vector<annot_t> intvec_t;
    typedef sdsl::rrr_vector<> beta_t;
    typedef beta_t::rank_1_type rank1_t;
    typedef beta_t::rank_0_type rank0_t;
    typedef array_int::array_int alpha_t;
    //typedef cpp_int alpha_t;
    //typedef sdsl::sd_vector<> alpha_t;
    typedef sdsl::bit_vector bv_t;
    //common prefix, allequal, size of block
    typedef std::tuple<alpha_t, bool> prefix_t;
    typedef std::tuple<alpha_t, bool, intvec_t> prefix_trans_t;
    typedef boost::archive::binary_oarchive oarchive_t;
    typedef boost::archive::binary_iarchive iarchive_t;

    void clear_after(array_int::array_int &a, size_t i) {
        a.indices.erase(a.indices.lower_bound(i), a.indices.end());
    }

    void clear_after(cpp_int &a, size_t i) {
        //reference:
        //a &= (alpha_t(1) << i)-1;
        cpp_int temp=0;
        bit_set(temp, i);
        --temp;
        a &= temp;
        return;
//        if (i == 0 || a.backend().size() == 0) {
//            a = 0;
//            return;
//        }
//        size_t limb_size=sizeof(*a.backend().limbs())*8lu;
//        size_t start = ceil((double)i/(double)limb_size);
//        if (start*limb_size > i) {
//            if (start < a.backend().size() + 1)
//                *(a.backend().limbs()+start-1) &= (1lu << (limb_size + i + 1 - start*limb_size))-1;
//        }
//        if (start < a.backend().size()) {
//            std::fill(a.backend().limbs()+start, a.backend().limbs()+a.backend().size(), 0lu);
//        }
    }

    //end is the max iteration
    prefix_trans_t longest_prefix_transpose(size_t ivb, size_t ive, intvec_t::iterator begin, intvec_t::iterator end, bool sorted=false, bool debug=false) {
        bool allequal=true;
        intvec_t::iterator done=end;
        alpha_t prefix=1lu,temp;
        intvec_t newcols;
        if (ive-ivb) {
            auto it=begin;
            prefix=0lu;
            while (it != end) {
                //std::cout << it-begin << "\t" << end-begin << "\n";
                temp = (*it) >> ivb;
                clear_after(temp, ive-ivb);
                //temp &= (alpha_t(1) << (ive-ivb))-1;
                newcols.push_back(temp);
                if (debug)
                    std::cout << temp;
                if (done==end && temp != 0lu) {
                    if ((temp & (temp+1lu)) == 0lu && (ive-ivb == 1lu || temp > 1lu)) {
                        //all ones
                        bit_set(prefix, (size_t)(it-begin));
                    } else {
                        //end of common prefix
                        allequal = false;
                        if (debug)
                            std::cout << "<\n";
                        done = it;
                        //break;
                    }
                }
                if (debug)
                    std::cout << "\n";
                ++it;
            }
            //boost::multiprecision::bit_set(prefix, (size_t)(it-begin));
            bit_set(prefix, (size_t)(done-begin));
        }
        if (debug)
            std::cout << "\n";
        return std::make_tuple(prefix, allequal, newcols);
    }

    prefix_t longest_prefix(intvec_t::iterator ivb, intvec_t::iterator ive, size_t begin=0, size_t end=MAXNUM, bool sorted=false, bool debug=false) {
        std::pair<alpha_t, size_t> prefix(1,0);
        bool allequal=true;
        intvec_t::iterator i = ivb;
        if ((ive - ivb) > 0) {
            //if (*ivb != 0 && boost::multiprecision::msb(*ivb) >= end) {
            //    return std::make_tuple(1, allequal, 0);
            //}
            std::get<0>(prefix) = *ivb;
            std::get<1>(prefix) = MAXNUM;
            for (++i;i != ive;++i) {
            //for (size_t i=(sorted ? size-1 : 1);i<size;++i) {
                if (*i != std::get<0>(prefix)) {
                    //if (*i != 0 && boost::multiprecision::msb(*i) >= end)
                    //    break;
                    std::get<1>(prefix) = std::min((size_t)std::get<1>(prefix), (size_t)lsb(((*i) ^ (alpha_t)std::get<0>(prefix))));
                    //std::get<0>(prefix) &= (alpha_t(1) << std::get<1>(prefix))-1;
                    clear_after(std::get<0>(prefix), std::get<1>(prefix));
                    if (std::get<1>(prefix) == 0) {
                        break;
                    }
                }
            }
            if (std::get<1>(prefix) == MAXNUM) {
                std::get<1>(prefix) = (std::get<0>(prefix) != 0) ? msb((alpha_t)std::get<0>(prefix))+1 : 0;
            } else {
                allequal = false;
            }
            //set the last bit to indicate the length of the prefix
            bit_set(std::get<0>(prefix), std::get<1>(prefix));
        }
        return std::make_tuple(std::get<0>(prefix), allequal);

    }

    class Node {
        public:
            alpha_t alpha=1;
            beta_t beta;
            rank1_t rank1;
            rank0_t rank0;
            Node *child[2] = {NULL, NULL};
            Node(Node *that, bool copy=true);
            Node();
            void serialize(std::ofstream &os, std::ofstream &oas);
            size_t add(intvec_t::iterator ivb, intvec_t::iterator ive, size_t begin, size_t end, prefix_t prefix, bool sorted=false, bool debug=false);
            size_t add(size_t ivb, size_t ive, intvec_t::iterator begin, intvec_t::iterator end, bool sorted=false, bool debug=false);
            int move_label_down(size_t ol, size_t add=0);
            void append(Node *other, bool leftorright=false, size_t addme=0, size_t addother=0, bool debug=false);
            void print(std::ostream &os, bool recursive=false);
            void print(bool recursive=false);
            ~Node();
            void breakleaf(size_t addme);
    };

    Node::Node() {
    }

    Node::Node(Node *that, bool copy) {
        this->alpha = that->alpha;
        this->beta = that->beta;
        this->rank1 = that->rank1;
        this->rank0 = that->rank0;
        if (that->child[0]) {
            if (copy) {
                this->child[0] = new Node(that->child[0]);
            } else {
                this->child[0] = that->child[0];
                that->child[0] = NULL;
            }
        }
        if (that->child[1]) {
            if (copy) {
                this->child[1] = new Node(that->child[1]);
            } else {
                this->child[1] = that->child[1];
                that->child[1] = NULL;
            }
        }
    }

    bv_t copy_bits(array_int::array_int &abv, size_t length=0) {
        return bv_t();
    }

    //assume that limbs are 64 bits long
    bv_t copy_bits(cpp_int &cbv, size_t length=0) {
        if (cbv == 0) {
            return bv_t(1);
        }
        size_t limb_size = sizeof(*cbv.backend().limbs())*8lu;
        bv_t bv(cbv.backend().size()*limb_size);
        assert(cbv.backend().size() <= (bv.capacity()>>6));
        std::copy(cbv.backend().limbs(), cbv.backend().limbs()+cbv.backend().size(), bv.data());
        bv.resize(length == 0 ? msb(cbv)+1 : length);
        return bv;
    }

    array_int::array_int copy_bits(bv_t &sbv) {
        return array_int::array_int();
    }
    /*
    cpp_int copy_bits(bv_t &sbv) {
        cpp_int cbv=0;
        bit_set(cbv, sbv.size());
        assert((sbv.capacity()>>6) <= cbv.backend().size());
        std::copy(sbv.data(), sbv.data()+(sbv.capacity()>>6), cbv.backend().limbs());
        bit_unset(cbv, sbv.size());
        return cbv;
    }
    */

    void copy_bits(bv_t &target, beta_t &source, size_t start=0, size_t t_bs=15) {
        assert(target.size()-start >= source.size());
        //bv_t test = target;
        size_t i=0;
        for (;i+t_bs<source.size();i+=t_bs) {
            target.set_int(start+i, source.get_int(i, t_bs), t_bs);
        }
        target.set_int(start+i, source.get_int(i, source.size()-i), source.size()-i);
        /*
        //sanity check
        for (size_t i=0;i<source.size();++i)
            test[start+i]=source[i];
        if (test != target) {
            std::cerr << "start: " << start << "\tb_start: " << b_start << "\ttotalsize: " << totalsize << "\n";
            std::cerr << "other: " << source << "\n";
            std::cerr << "E: " << test << "\tO: " << target << "\n";
            assert(false);
        }
        */
    }

    void Node::serialize(std::ofstream &os, std::ofstream &oas) {
        //oarchive_t oa(os);
        //oarchive_t oaa(oas);
        //bv_t bv = copy_bits(this->alpha);
        bv_t bv = copy_bits(this->alpha);
        //bv_t bv(boost::multiprecision::msb(this->alpha));
        //for (size_t i=0;i<bv.size();++i) {
        //    bv[i] = boost::multiprecision::bit_test(this->alpha, i);
        //}
        //sdsl::sd_vector<> sd(bv);
        beta_t sd(bv);
        //oa & this->alpha;
        sd.serialize(os);
        this->beta.serialize(oas);
        this->rank1.serialize(oas);
        this->rank0.serialize(oas);
    }

    class WTR {
        public:
            bool sorted=false;
            size_t height=0;
            size_t size=0;
            size_t maxannot=0;
            Node *root=NULL;
            WTR();
            WTR(intvec_t &iv, bool transpose, size_t osize, size_t batch, bool sorted=false, bool debug=false);
            //WTR(intvec_t &iv, size_t begin, size_t end, size_t transpose, bool sorted=false, bool debug=false);
            //WTR(intvec_t &iv, intvec_t &tiv, size_t begin=0, size_t end=MAXNUM, bool sorted=false, bool debug=false);
            //WTR(intvec_t &cols, size_t row, size_t nrows, bool sorted=false, bool debug=false);
            ~WTR();
            annot_t at(size_t i);
            void serialize(std::ofstream &os);
            void serialize_concat(std::ofstream &os, std::ofstream &oas);
            void load(std::ifstream &is);
            void check_structure();
            void append(WTR &other, bool debug=false);
            bool reconstruct(intvec_t::iterator ivb, intvec_t::iterator ive);
            bool reconstruct(intvec_t &iv);
            void print(std::ostream &os);
            void print();
    };

    void WTR::append(WTR &other, bool debug) {
        //void Node::append(Node *other, bool leftorright, size_t addme, size_t addother, bool debug) {
        this->size += other.size;
        #pragma omp parallel
        #pragma omp single nowait
        this->root->append(other.root, 0, this->size, other.size, debug);
    }

    void WTR::serialize(std::ofstream &os) {
        std::stack<Node*> node_stack;
        if (this->root != NULL)
            node_stack.push(this->root);
        Node *cur;
        while (node_stack.size()) {
            cur = node_stack.top();
            node_stack.pop();
            cur->serialize(os, os);
            if (cur->child[0] != NULL)
                node_stack.push(cur->child[0]);
            if (cur->child[1] != NULL)
                node_stack.push(cur->child[1]);
        }
    }

    void WTR::serialize_concat(std::ofstream &os, std::ofstream &oas) {
        //oarchive_t oa(os);
        std::stack<Node*> node_stack;
        std::vector<alpha_t> alphas;
        std::vector<beta_t> betas;
        std::vector<size_t> alpha_l;
        std::vector<size_t> beta_l;
        //std::vector<rank1_t> rank1s;
        //std::vector<rank0_t> rank0s;
        if (this->root != NULL)
            node_stack.push(this->root);
        Node *cur;
        size_t asize=0;
        size_t bsize=0;
        while (node_stack.size()) {
            cur = node_stack.top();
            node_stack.pop();
            alphas.push_back(cur->alpha);
            betas.push_back(cur->beta);
            alpha_l.push_back(msb(cur->alpha));
            beta_l.push_back(cur->beta.size());
            asize += alpha_l.back();
            bsize += beta_l.back();
            //rank1s.push_back(cur->rank1);
            //rank0s.push_back(cur->rank0);
            if (cur->child[0] != NULL)
                node_stack.push(cur->child[0]);
            if (cur->child[1] != NULL)
                node_stack.push(cur->child[1]);
        }
        std::cout << asize << "\t" << bsize << "\n";
        bv_t bv_alpha(asize);
        bv_t bv_beta(bsize);
        size_t ai=0,bi=0;
        for (size_t i=0;i<alpha_l.size();++i) {
            assert(alpha_l[i] == msb(alphas[i]));
            for (size_t j=0;j<alpha_l[i];++j) {
                bv_alpha[ai++] = bit_test(alphas[i], j);
            }
            assert(beta_l[i] == betas[i].size());
            for (size_t j=0;j<beta_l[i];++j) {
                bv_beta[bi++] = betas[i][j];
            }
        }
        beta_t alpha = beta_t(bv_alpha);
        beta_t beta = beta_t(bv_beta);
        alpha.serialize(os);
        beta.serialize(oas);
    }

    Node::~Node() {
        if (this->child[0] != NULL) {
            delete this->child[0];
            this->child[0] = NULL;
        }
        if (this->child[1] != NULL) {
            delete this->child[1];
            this->child[1] = NULL;
        }
    }

    WTR::~WTR() {
        if (this->root != NULL) {
            delete this->root;
            this->root = NULL;
        }
    }

    void Node::breakleaf(size_t addme) {
        if (addme) {
            this->beta = beta_t(bv_t(addme));
            sdsl::util::init_support(this->rank1, &(this->beta));
            sdsl::util::init_support(this->rank0, &(this->beta));
            assert(this->child[0] == NULL);
            assert(this->child[1] == NULL);
            this->child[0] = new Node();
            this->child[0]->alpha = 1;
        }
    }

    size_t Node::add(intvec_t::iterator ivb, intvec_t::iterator ive, size_t begin, size_t end, prefix_t prefix, bool sorted, bool debug) {
        if (ive-ivb > 0 && end > begin) {
            if (debug)
                std::cout << ive - ivb << "\n";
            //prefix_t test = longest_prefix(ivb, ive, begin, end, sorted, debug);
            if (debug)
                std::cout << ive-ivb << "\t" << std::get<0>(prefix) << "\t" << std::get<1>(prefix) << "\n";
            this->alpha = std::get<0>(prefix);
            size_t length = msb(std::get<0>(prefix));
            //if (!length && prefix.second) {
            //    assert(prefix.first == 1);
            //} else {
            if (!std::get<1>(prefix)) {
                intvec_t* children[2] = {new intvec_t(), new intvec_t()};
                //std::cout << ive-ivb << "\t" << std::get<1>(prefix) << "\n";
                if (debug) {
                    assert(this->child[0] == NULL);
                    assert(this->child[1] == NULL);
                }
                this->child[0] = new Node();
                this->child[1] = new Node();
                //bv_t bv(civ.size());
                bv_t bv(ive-ivb);
                alpha_t temp = 0;
                bit_set(temp, end);
                --temp;
                prefix_t cprefix[2] = {std::make_pair(temp, true), std::make_pair(temp, true)};
                //prefix_t cprefix[2] = {std::make_pair(0, true), std::make_pair(0, true)};
                bool set[2] = {false, false};
                size_t lengths[2] = {MAXNUM, MAXNUM};
                //std::cout << "\n";
                for (size_t i=0;i<bv.size();++i) {
                    bv[i] = bit_test(*(ivb+i), length);
                    children[bv[i]]->push_back((*(ivb+i)) >> (length+1));
                    if (children[bv[i]]->back() != std::get<0>(cprefix[bv[i]])) {
                        if (!set[bv[i]]) {
                            std::get<0>(cprefix[bv[i]]) = children[bv[i]]->back();
                            set[bv[i]]=true;
                        } else {
                            lengths[bv[i]] = std::min(lengths[bv[i]], (size_t)lsb((children[bv[i]]->back()) ^ std::get<0>(cprefix[bv[i]])));
                            //std::get<0>(cprefix[bv[i]]) &= (alpha_t(1) << lengths[bv[i]])-1;
                            clear_after(std::get<0>(cprefix[bv[i]]), lengths[bv[i]]);
                        }
                    }
                    //std::cout << "p" << bv[i] << "\t" << children[bv[i]]->back() << "\t" << std::get<0>(cprefix[bv[i]]) << "," << std::get<1>(cprefix[bv[i]]) << "\t";
                }
                //std::cout << "\n";
                for (size_t i=0;i<2;++i) {
                    if (lengths[i] == MAXNUM) {
                        lengths[i] = (std::get<0>(cprefix[i]) != 0) ? msb((alpha_t)std::get<0>(cprefix[i]))+1 : 0;
                    } else {
                        std::get<1>(cprefix[i]) = false;
                    }
                    //set the last bit to indicate the length of the prefix
                    bit_set(std::get<0>(cprefix[i]), lengths[i]);
                    //prefix_t test = longest_prefix(children[i]->begin(), children[i]->end(), begin, end, sorted, debug);
                    //if (std::get<0>(test) != std::get<0>(cprefix[i]) || std::get<1>(test) != std::get<1>(cprefix[i])) {
                    //    std::cerr << "\nFail prefix" << i << "\niv\t";
                    //    for (auto it=children[i]->begin(); it!=children[i]->end(); ++it)
                    //        std::cerr << *it << "\t";
                    //    std::cerr << "\npre\t";
                    //    std::cerr << std::get<0>(test) << "\t" << std::get<0>(cprefix[i]) << "\nall\t";
                    //    std::cerr << std::get<1>(test) << "\t" << std::get<1>(cprefix[i]) << "\n";
                    //    assert(false);
                    //}
                }
                //if (debug)
                //    std::cout << children[bv[i]]->back() << ",";
                //if (debug)
                //    std::cout << "\t";
                beta_t beta(bv);
                this->beta = beta;
                sdsl::util::init_support(this->rank1, &(this->beta));
                sdsl::util::init_support(this->rank0, &(this->beta));
                //assert(this->rank1(this->beta.size()) == children[1]->size());
                //assert(this->rank0(this->beta.size()) == children[0]->size());
                assert(children[0]->size());
                assert(children[1]->size());
                #pragma omp task shared(cprefix) if(children[0]->size() > TASKMIN)
                this->child[0]->add(children[0]->begin(), children[0]->end(), begin, end, cprefix[0], sorted, debug);
                //#pragma omp taskwait
                //delete children[0];
                #pragma omp task shared(cprefix) if (children[1]->size() > TASKMIN)
                this->child[1]->add(children[1]->begin(), children[1]->end(), begin, end, cprefix[1], sorted, debug);
                //#pragma omp taskwait
                //delete children[1];
                #pragma omp taskwait
                {
                    delete children[0];
                    delete children[1];
                    for (size_t i=0;i<2;++i) {
                        assert(this->child[i] != NULL);
                        if (this->child[i]->alpha == 0) {
                            delete this->child[i];
                            this->child[i] = NULL;
                        }
                    }
                }
            } else {
                this->breakleaf(ive-ivb);
            }
            return ive-ivb;
        //} else {
            //std::cerr << "Input empty vector\n";
            //assert(false);
        } else {
            this->breakleaf(ive-ivb);
        }
        return 0lu;
    }

    size_t Node::add(size_t ivb, size_t ive, intvec_t::iterator begin, intvec_t::iterator end, bool sorted, bool debug) {
        if (ive > ivb && end > begin) {
            prefix_trans_t prefix = longest_prefix_transpose(ivb, ive, begin, end, sorted, debug);
            //ive = std::get<2>(prefix);
            this->alpha = std::get<0>(prefix);
            if (!std::get<1>(prefix)) {
                intvec_t* tchildren[2] = {new intvec_t(end-begin), new intvec_t(end-begin)};
                size_t prefix_length = msb(this->alpha);
                //assert(prefix_length < end-begin);
                if (debug) {
                    assert(this->child[0] == NULL);
                    assert(this->child[1] == NULL);
                }
                this->child[0] = new Node();
                this->child[1] = new Node();

                //copy over beta
                //alpha_t cbv = (*(begin+prefix_length) >> ivb) & ((alpha_t(1) << (ive-ivb))-1);
                alpha_t cbv = *(begin+prefix_length) >> ivb;
                clear_after(cbv, ive-ivb);
                //sanity check
                if (cbv == 0) {
                    std::cout << "Bad prefix\t" << ivb << "\t" << ive << "\n";
                    for (auto k=begin; k!= end;++k) {
                        std::cout << (((*k) >> ivb) & ((alpha_t(1) << (ive-ivb)-1)));
                        if (k == begin+prefix_length)
                            std::cout << "<";
                        std::cout << "\n";
                    }
                    assert(false);
                }
                bv_t bv = copy_bits(cbv, ive-ivb);
                //TODO: std::get<2>(prefix) has the new cols
                size_t tsizes[2]={0,0};
                size_t maxsize[2]={0,0};
                for (size_t i=0;i<bv.size();++i) {
                    for (auto j=begin+prefix_length+1;j!=end;++j) {
                        if (bit_test(*j, ivb+i)) {
                            maxsize[bv[i]] = std::max(maxsize[bv[i]], (size_t)(j-(begin+prefix_length+1)));
                            bit_set(tchildren[bv[i]]->operator[](j-(begin+prefix_length+1)), tsizes[bv[i]]);
                        }
                    }
                    tsizes[bv[i]]++;
                }
                assert(tsizes[0] + tsizes[1] == bv.size());
                assert(tsizes[0] < bv.size());
                tchildren[0]->resize(maxsize[0]+1);
                tchildren[1]->resize(maxsize[1]+1);
                beta_t beta(bv);
                this->beta = beta;
                sdsl::util::init_support(this->rank1, &(this->beta));
                sdsl::util::init_support(this->rank0, &(this->beta));
                #pragma omp task if (tsizes[0] > TASKMIN)
                this->child[0]->add(0, tsizes[0], tchildren[0]->begin(), tchildren[0]->end(), sorted, debug);
                #pragma omp taskwait
                delete tchildren[0];
                #pragma omp task if (tsizes[1] > TASKMIN)
                this->child[1]->add(0, tsizes[1], tchildren[1]->begin(), tchildren[1]->end(), sorted, debug);
                #pragma omp taskwait
                delete tchildren[1];
                #pragma omp taskwait
                {
                    for (size_t i=0;i<2;++i) {
                        assert(this->child[i] != NULL);
                        if (this->child[i]->alpha == 0) {
                            delete this->child[i];
                            this->child[i] = NULL;
                        }
                    }
                }
            } else {
                this->breakleaf(ive-ivb);
                //this->breakleaf(std::get<2>(prefix));
            }
            //return std::get<2>(prefix);
            return ive-ivb;
        }
        return 0;
    }


    int Node::move_label_down(size_t ol, size_t add) {
        size_t len = msb(this->alpha);
        if (ol > len) {
            this->alpha = 0;
            bit_set(this->alpha, ol);
            //this->alpha = alpha_t(1) << ol;
            return 3;
        }
        if (ol == len)
            return 1;
        Node *temp_node = new Node();
        temp_node->alpha = this->alpha >> (ol+1);
        temp_node->beta = this->beta;
        temp_node->child[0] = this->child[0];
        temp_node->child[1] = this->child[1];
        this->beta = beta_t(bv_t(this->beta.size() == 0 ? add : this->beta.size(), bit_test(this->alpha, ol)));
        sdsl::util::init_support(this->rank1, &(this->beta));
        sdsl::util::init_support(this->rank0, &(this->beta));
        //this->alpha = (this->alpha & ((alpha_t(1) << (ol))-1)) | (alpha_t(1) << ol);
        clear_after(this->alpha, ol);
        bit_set(this->alpha, ol);
        if (this->beta[0]) {
            this->child[1] = temp_node;
            this->child[0] = NULL;
        } else {
            this->child[0] = temp_node;
            this->child[1] = NULL;
        }
        assert(msb(this->alpha) == ol);
        return 0;
    }

    void Node::append(Node *other, bool leftorright, size_t addme, size_t addother, bool debug) {
        //debug=true;
        if (debug) {
            this->print(true);
            std::cout << "-\n";
        }
        if (other == NULL) {
            if (debug)
                std::cout << "-\n";
            //propagate on the right side
            if (this->beta.size() == 0) {
                if (this->alpha != 1) {
                    assert(false);
                } else {
                    assert(this->child[0] == NULL);
                    assert(this->child[1] == NULL);
                    if (debug)
                        std::cout << "--------\n";
                    return;
                }
            }
            this->move_label_down(lsb(this->alpha));
            bv_t bv(this->beta.size()+addother);
            size_t i=(leftorright ? 0 : addother);
            copy_bits(bv, this->beta, i);
            this->beta = beta_t(bv);
            sdsl::util::init_support(this->rank1, &(this->beta));
            sdsl::util::init_support(this->rank0, &(this->beta));
            if (debug) {
                this->print(true);
                std::cout << "--------\n";
            }
            if (this->child[0] == NULL) {
                if (debug)
                    std::cout << "Adding new leaf, since destroyed\n";
                this->child[0] = new Node();
                this->child[0]->alpha=1;
            }
            this->child[0]->append(other, leftorright, 0, addother);
            return;
        }
        if (debug) {
            other->print(true);
            std::cout << "\n";
        }

        if (this->beta.size() == 0) {
            assert(this->child[0] == NULL);
            assert(this->child[1] == NULL);
            if (this->alpha != 1) {
                this->breakleaf(addme);
                assert(this->beta.size() > 0);
                assert(this->child[0] != NULL);
                assert(this->child[0]->alpha == 1);
            }
        }
        if (other->beta.size() == 0) {
            assert(other->child[0] == NULL);
            assert(other->child[1] == NULL);
            if (other->alpha != 1) {
                other->breakleaf(addother);
                assert(other->beta.size() > 0);
                assert(other->child[0] != NULL);
                assert(other->child[0]->alpha == 1);
            }
        }
        if (other->beta.size() == 0) {
            if (other->alpha != 1) {
                std::cerr << "Bad leaf\n";
                std::cerr << other->alpha << "\n";
                assert(false);
            }
            if (this->beta.size() == 0) {
                if (debug)
                    std::cout << "Done with this branch";
                assert(this->alpha == 1 && other->alpha == 1);
                return;
            } else {
                if (debug)
                    std::cout << "Fixing left children, append to right\n";
                this->append(NULL, 1, 0, addother);
                return;
            }
        }
        if (this->beta.size() == 0) {
            assert(this->alpha == 1);
            assert(this->child[0] == NULL);
            assert(this->child[1] == NULL);
            this->alpha = other->alpha;
            this->beta = other->beta;
            this->rank1 = other->rank1;
            this->rank0 = other->rank0;
            if (other->child[0] != NULL)
                this->child[0] = new Node(other->child[0], false);
            if (other->child[1] != NULL)
                this->child[1] = new Node(other->child[1], false);
            if (debug)
                std::cout << "Fixing left children, append to left\n";
            this->append(NULL, 0, 0, addme);
            return;
        }
        if (debug)
            std::cout << "Fixing internal node\n";
        
        //std::cout << this->alpha << "\t" << other->alpha << "\n";
        intvec_t alphas;
        alphas.push_back(this->alpha);
        alphas.push_back(other->alpha);
        bit_unset(alphas[0], msb(alphas[0]));
        bit_unset(alphas[1], msb(alphas[1]));
        annot_t overlap = (alphas[0] ^ alphas[1]);
        size_t tlen=msb(this->alpha);
        size_t olen=msb(other->alpha);
        size_t ol=0;
        if (overlap == 0) {
            ol = std::min(tlen, olen);
        } else {
            ol = lsb(overlap);
            if (ol > tlen || ol > olen) {
                assert((ol > tlen) ^ (ol > olen));
                ol = std::min(tlen, olen);
            }
        }
        if (olen != tlen || overlap != 0) {
            if (debug)
                std::cout << "Different prefix, reshaping\n";
            assert(this->beta.size() > 0 || this->alpha == 1);
            assert(other->beta.size() > 0 || other->alpha == 1);
            this->move_label_down(ol, addme);
            other->move_label_down(ol, addother);
            assert(this->alpha == other->alpha);
        }
        if (debug)
            std::cout << "Merging\n";

        //fix alphas
        bv_t bv(this->beta.size()+other->beta.size());
        copy_bits(bv, this->beta);
        copy_bits(bv, other->beta, this->beta.size());
        //size_t i=0,j=0;
        //for (;i<this->beta.size();++i)
        //    bv[i] = this->beta[i];
        //for (;j<other->beta.size();++j)
        //    bv[i+j] = other->beta[j];
        if (debug) {
            std::cout << bv << "\n";
            this->print(true);
            std::cout << "-\n";
            other->print(true);
            std::cout << "--------\n";
        }
        this->rank0.set_vector(&(this->beta));
        other->rank0.set_vector(&(other->beta));
        this->rank1.set_vector(&(this->beta));
        other->rank1.set_vector(&(other->beta));
        if (this->child[0] != NULL) {
            #pragma omp task if(this->rank0(this->beta.size()) > TASKMIN || other->rank0(other->beta.size()) > TASKMIN)
            this->child[0]->append(other->child[0], 0, this->rank0(this->beta.size()), other->rank0(other->beta.size()));
        } else if (other->child[0] != NULL) {
            if (this->rank0(this->beta.size()) > 0) {
                std::cerr << this->beta << "\n";
                assert(false);
            }
            this->child[0] = new Node(other->child[0], false);
        }
        if (this->child[1] != NULL) {
            #pragma omp task if(this->rank1(this->beta.size()) > TASKMIN || other->rank1(other->beta.size()) > TASKMIN)
            this->child[1]->append(other->child[1], 1, this->rank1(this->beta.size()), other->rank1(other->beta.size()));
        } else if (other->child[1] != NULL) {
            if (this->rank1(this->beta.size()) > 0) {
                std::cerr << this->beta << "\n";
                assert(false);
            }
            this->child[1] = new Node(other->child[1], false);
        }
        #pragma omp taskwait
        this->beta = beta_t(bv);
        sdsl::util::init_support(this->rank1, &(this->beta));
        sdsl::util::init_support(this->rank0, &(this->beta));
    }

    Node* merge(std::vector<Node*>::iterator begin, std::vector<Node*>::iterator end) {
        if (end == begin)
            return NULL;
        if (end-begin == 1) {
            return *begin;
        }
        if (end-begin == 2) {
            (*begin)->append(*(begin+1), false);
            delete (*(begin+1));
            return *begin;
        }
        Node *a, *b;
        #pragma omp task shared(a)
        a = merge(begin, begin+((end-begin)/2));
        #pragma omp task shared(b)
        b = merge(begin+((end-begin)/2), end);
        #pragma omp taskwait
        a->append(b, false);
        delete b;
        return a;
    }

    WTR merge(std::vector<WTR*>::iterator begin, std::vector<WTR*>::iterator end) {
        std::vector<Node*> nodes(end-begin);
        auto jt = nodes.begin();
        size_t totalsize=0;
        for (auto it = begin; it != end; ++it) {
            assert(*it != NULL);
            *jt = (*it)->root;
            ++jt;
            totalsize += (*it)->size;
            (*it)->root = NULL;
        }
        WTR wtr;
        if (end > begin) {
            Node *a;
            #pragma omp parallel
            #pragma omp single nowait
            a = merge(nodes.begin(), nodes.end());
            wtr.size = totalsize;
            wtr.root = a;
        }
        return wtr;
    }


    //when transpose is true, osize is the number of rows, otherwise, it's the number of columns
    WTR::WTR(intvec_t &iv, bool transpose, size_t osize, size_t batch, bool sorted, bool debug) {
        //this->sorted = sorted;
        if (std::getenv("DEBUG") != NULL)
            debug=true;
        this->root = new Node();
        if (batch != 0 && iv.size() != 0) {
        size_t tsize=0;
        prefix_t prefix;
        //std::cout << "Computing alpha\n";
        if (transpose) {
            //prefix = longest_prefix_transpose(0, std::min(batch, osize), iv.begin(), iv.end(), sorted, debug);
            //std::cout << "Computing wavelet trie\n";
            #pragma omp parallel
            #pragma omp single nowait
            this->root->add(0, std::min(batch, osize), iv.begin(), iv.end(), sorted, debug);
            tsize = osize;
        } else {
            prefix = longest_prefix(iv.begin(), std::min(iv.begin()+batch, iv.end()), 0, osize, sorted, debug); 
            //std::cout << "Computing wavelet trie\n";
            #pragma omp parallel
            #pragma omp single nowait
            this->root->add(iv.begin(), std::min(iv.begin()+batch, iv.end()), 0, osize, prefix, sorted, debug); 
            tsize = iv.size();
        }
        if (debug) {
            this->print();
            std::cout << "-------------\n\n";
        }
        this->size = tsize;
        //Node *other;
        std::vector<Node*> others(tsize/batch,NULL);
        #pragma omp parallel for shared(others)
        for (size_t cursize = batch; cursize < tsize; cursize += batch) {
            //other = new Node();
            others[cursize/batch-1] = new Node();
            //std::cout << "Computing alpha\n";
            if (transpose) {
                //#pragma omp parallel
                //#pragma omp single nowait
                //other->add(cursize, std::min(cursize+batch, osize), iv.begin(), iv.end(), sorted, debug);
                //std::cout << "Computing wavelet trie\n";
                //prefix = longest_prefix_transpose(cursize, std::min(cursize+batch, osize), iv.begin(), iv.end(), sorted, debug);
                others[cursize/batch-1]->add(cursize, std::min(cursize+batch, osize), iv.begin(), iv.end(), sorted, debug);
            } else {
                //#pragma omp parallel
                //#pragma omp single nowait
                //other->add(iv.begin()+cursize, std::min(iv.begin()+cursize+batch, iv.end()), 0, osize, sorted, debug); 
                //std::cout << "Computing wavelet trie\n";
                prefix = longest_prefix(iv.begin()+cursize, std::min(iv.begin()+cursize+batch, iv.end()), 0, osize, sorted, debug); 
                others[cursize/batch-1]->add(iv.begin()+cursize, std::min(iv.begin()+cursize+batch, iv.end()), 0, osize, prefix, sorted, debug); 
            }
            if ((cursize/batch-1) % 100 == 0)
                std::cout << "\n";
            std::cout << "." << std::flush;
            
        }
        if (others.size() && others[0]) {
            std::cout << "\nMerging\n";
            Node *a;
            #pragma omp parallel
            #pragma omp single nowait
            a = merge(others.begin(), others.end());
            #pragma omp parallel
            #pragma omp single nowait
            this->root->append(a, false);
            delete a;
        }
        /*
        for (size_t i=0;i<others.size();++i) {
            #pragma omp parallel
            #pragma omp single nowait
            this->root->append(others[i], false);
            delete others[i];
        }
        */
        }
    }

    annot_t WTR::at(size_t i) {
        annot_t annot = 0;
        Node *curnode = this->root;
        size_t len=0;
        uint8_t curbit=0;
        if (curnode != NULL) {
            while (curnode->beta.size()) {
                annot |= curnode->alpha << len;
                //annot += curnode->alpha << len;
                len += msb(curnode->alpha);
                curbit = curnode->beta[i];
                if (curbit) {
                    curnode->rank1.set_vector(&(curnode->beta));
                    i = curnode->rank1(i);
                } else {
                    bit_unset(annot, len);
                    curnode->rank0.set_vector(&(curnode->beta));
                    i = curnode->rank0(i);
                }
                curnode = curnode->child[curbit];
                if (curnode == NULL) {
                    std::cerr << "Missing node\n";
                    assert(false);
                }
                len++;
            }
            assert(curnode != NULL);
            annot |= curnode->alpha << len;
            //annot += curnode->alpha << len;
            bit_unset(annot, msb(annot));
        }
        return annot;
    }

    bool WTR::reconstruct(intvec_t::iterator ivb, intvec_t::iterator ive) {
        annot_t annot=0;
        for (auto it=ivb; it!= ive; ++it) {
            annot = this->at(it-ivb);
            if (annot != *it) {
                std::cerr << "Fail " << it-ivb << " " << *it << " " << annot << "\n";
                return false;
            }
        }
        return true;
    }

    bool WTR::reconstruct(intvec_t &iv) {
        return this->reconstruct(iv.begin(), iv.end());
    }

    std::pair<intvec_t, intvec_t> construct(std::ifstream &pfile, size_t begin=0, size_t end=MAXNUM, bool debug=false) {
        intvec_t iv;
        std::string line, ed;
        annot_t cannot, ctemp;
        size_t edi = 0, maxedi=0;
        while (pfile >> line) {
            std::istringstream sfile(line);
            cannot = 0;
            while (std::getline(sfile, ed, ',')) {
                edi = atol(ed.c_str());
                if (debug)
                    std::cout << edi << ",";
                if (edi) {
                    maxedi = std::max(maxedi, edi);
                    ctemp = 0;
                    bit_set(ctemp, edi-1);
                    if (ctemp == 0) {
                        std::cerr << "Fail at " << edi << "\n";
                        assert(false);
                    }
                    cannot |= ctemp;
                    //ctemp = cannot | (annot_t(1) << (edi-1));
                    //if (ctemp == cannot) {
                    //    std::cerr << "Fail at " << edi << "\n";
                    //    assert(false);
                    //}
                    //cannot = ctemp;
                }
            }
            //if (edi >= begin && edi < end)
            iv.push_back(cannot);
            if (iv.size() == begin)
                break;
            if (debug)
                std::cout << "\t" << iv.back() << "\n";
        }
        if (debug)
            std::cout << "\n";
        return std::make_pair(iv, intvec_t(maxedi));
    }

    size_t popcount(const cpp_int &a) {
        size_t as = a.backend().size();
        size_t apc=0;
        auto al = a.backend().limbs();
        size_t i=0;
        for (;i<as;++i) {
            apc += __builtin_popcount(*al);
            ++al;
        }
        return apc;
    }

    size_t popcount(const array_int::array_int &a) {
        return a.indices.size();
    }

    bool midless(const annot_t &a, const annot_t &b, const size_t maxsize) {
        //all zeros
        if (a==0 || b==0)
            return a<b;
        //all ones
        if ((a & (a+1)) == 0)
            return 1;
        if ((b & (b+1)) == 0)
            return 0;
        //highest diversity
        /*
        size_t apc=0,al=MAXNUM;
        size_t bpc=0,bl=MAXNUM;
        bool lownum=0;
        annot_t ac=a,bc=b;
        if (ac > 0 || bc > 0) {
            while (true) {
                if (ac > 0) {
                    al = boost::multiprecision::lsb(ac);
                    boost::multiprecision::bit_unset(ac, al);
                    apc++;
                    if (bc == 0)
                        break;
                }
                if (bc > 0) {
                    bl = boost::multiprecision::lsb(bc);
                    boost::multiprecision::bit_unset(bc, bl);
                    bpc++;
                    if (ac == 0) {
                        if (bc>0)
                            bpc++;
                        break;
                    }
                }
                //check which one has the least set bit
                if (!lownum && al < MAXNUM && al > bl) {
                    lownum = 1;
                }
            }
        }
        */
        /*
        size_t as = a.backend().size(), bs = b.backend().size();
        size_t apc=0, bpc=0;
        auto al = a.backend().limbs();
        auto bl = b.backend().limbs();
        size_t i=0;
        for (;i<std::min(as,bs);++i) {
            apc += __builtin_popcount(*al);
            bpc += __builtin_popcount(*bl);
            ++al;
            ++bl;
        }
        if (as > bs) {
            for (;i<as;++i) {
                apc += __builtin_popcount(*al);
                ++al;
            }
        } else if (as > bs) {
            for (;i<bs;++i) {
                bpc += __builtin_popcount(*bl);
                ++bl;
            }
        }
        */
        size_t apc = popcount(a);
        size_t bpc = popcount(b);
        //if a tie, but the one with the largest indices first
        if (abs(maxsize/2-apc) == abs(maxsize/2-bpc)) {
            return lsb(a) > lsb(b);
            //return lownum;
        }
        return abs(maxsize/2-apc) < abs(maxsize/2-bpc);
    }

    bool intlex(const annot_t &a, const annot_t &b) {
        if (a==b)
            return 0;
        size_t clsb = lsb(a ^ b);
        return bit_test(a, clsb) == 0;
    }

    std::pair<intvec_t, intvec_t> construct_transpose(std::ifstream &pfile, size_t begin=0, size_t end=MAXNUM, bool sort=false, bool debug=false) {
        //std::map<size_t, annot_t> iv;
        intvec_t iv;
        std::string line, ed;
        //std::vector<size_t> edi;
        size_t edi=0;
        size_t maxsize=0;
        for (maxsize=0;pfile >> line;++maxsize) {
            std::istringstream sfile(line);
            //edi.clear();
            while (std::getline(sfile, ed, ',')) {
                //edi.push_back(atol(ed.c_str()));
                edi = atol(ed.c_str());
                if (edi) {
                    if (edi >= iv.size()) {
                        iv.resize(edi+1);
                    }
                    bit_set(iv.at(edi), maxsize);
                    assert(bit_test(iv[edi], maxsize));
                }
                //if (edi.back() < begin || edi.back() >= end)
                //    break;
            }
            /*
            for (auto it=edi.begin(); it!=edi.end();it++) {
                if (*it) {
                    if (*it >= iv.size())
                        iv.resize(*it+1);
                    boost::multiprecision::bit_set(iv.at(*it), maxsize);
                    assert(boost::multiprecision::bit_test(iv[*it], maxsize));
                }
            }
            */
        }
        //for (size_t i=0;i<iv.size();++i)
        //    std::cout << i << "\t" << iv[i] << "\t" << (iv[i] & (iv[i]+1)) << "\n";
        std::cout << "Rearranging columns\n";
        __gnu_parallel::sort(iv.begin(), iv.end(), std::bind(midless, std::placeholders::_1, std::placeholders::_2, maxsize));
        //for (size_t i=0;i<iv.size();++i)
        //    std::cout << i << "\t" << iv[i] << "\t" << (iv[i] & (iv[i]+1)) << "\n";
        std::cout << "Constructing rows\n";
        annot_t cannot;
        intvec_t iv_order;
        for (size_t i=0;i<maxsize;++i) {
            //TODO: replace this with copy_bits
            cannot=0;
            for (size_t j=0;j<iv.size();++j) {
                if (bit_test(iv[j], i)) {
                    bit_set(cannot, j);
                }
            }
            /*
            for (auto j = iv.begin();j!=iv.end();++j) {
                if (boost::multiprecision::bit_test(j->second, i)) {
                    boost::multiprecision::bit_set(cannot, j->first);
                }
            }
            */
            iv_order.push_back(cannot);
        }
        if (sort) {
            std::cout << "Sorting rows\n";
            __gnu_parallel::sort(iv_order.begin(), iv_order.end(), intlex);
        }
        //sanity check
        /*
        prefix_t prefix1 = longest_prefix(iv_order);
        prefix_t prefix2 = longest_prefix_transpose(iv, iv_order.size());
        if (std::get<0>(prefix1) != std::get<0>(prefix2) || std::get<1>(prefix1) != std::get<1>(prefix2)) {
            std::cout << iv_order.size() << "\t" << iv.size() << "\n";
            std::cerr << std::get<0>(prefix1) << "\t" << std::get<1>(prefix1) << "\n";
            std::cerr << std::get<0>(prefix2) << "\t" << std::get<1>(prefix2) << "\n";
            assert(false);
        }
        */
        return std::make_pair(iv_order, iv);
    }
    void Node::print(std::ostream &os, bool recursive) {
        os << this->alpha << ";" << this->beta << "\t";
        if (this->child[0] != NULL)
            os << "l";
        if (this->child[1] != NULL)
            os << "r";
        os << "\n";
        if (this->child[0] != NULL)
            this->child[0]->print(os);
        if (this->child[1] != NULL)
            this->child[1]->print(os);
    }

    void Node::print(bool recursive) {
        this->print(std::cout, recursive);
    }

    void WTR::print(std::ostream &os) {
        os << "Printing\n";
        this->root->print(os, true);
        os << "\n";
    }

    void WTR::print() {
        this->print(std::cout);
    }

    WTR::WTR() {
    }

}

#endif
