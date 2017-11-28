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
#define MAXNUM (std::numeric_limits<std::uint16_t>::max())
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
    typedef bool (*bittest_t)(const cpp_int&, uint32_t i);
    bittest_t bit_test = boost::multiprecision::bit_test;
    /*
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
    //typedef sdsl::sd_vector<> beta_t;
    typedef sdsl::rrr_vector<> beta_t;
    typedef beta_t::rank_1_type rank1_t;
    typedef beta_t::rank_0_type rank0_t;
    //typedef array_int::array_int alpha_t;
    typedef cpp_int alpha_t;
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

    prefix_t longest_prefix(intvec_t::iterator ivb, intvec_t::iterator ive, size_t begin=0, size_t end=MAXNUM, bool sorted=false, bool debug=false) {
        std::pair<alpha_t, size_t> prefix(1,0);
        bool allequal=true;
        intvec_t::iterator i = ivb;
        if ((ive - ivb) > 0) {
            //if (*ivb != 0 && boost::multiprecision::msb(*ivb) >= end) {
            //    return std::make_tuple(1, allequal, 0);
            //}
            //std::get<0>(prefix) = *ivb;
            std::get<0>(prefix) = array_int::toint(*ivb);
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
            size_t add(intvec_t *iv, size_t begin, size_t end, prefix_t prefix, bool sorted=false, bool debug=false);
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
            WTR(intvec_t *iv, bool transpose, size_t osize, size_t batch, bool sorted=false, bool debug=false);
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

    typedef std::tuple<Node*, intvec_t*, prefix_t, size_t> nodestate_t;

    bool process_node(nodestate_t &curnode) {
        while (std::get<1>(curnode)) {
            if (std::get<1>(curnode)->size()) {
                std::get<0>(curnode)->alpha = std::get<0>(std::get<2>(curnode));
                size_t length = msb(std::get<0>(curnode)->alpha);
                if (!std::get<1>(std::get<2>(curnode))) {
                    intvec_t* children[2] = {new intvec_t(), new intvec_t()};
                    std::get<0>(curnode)->child[0] = new Node();
                    std::get<0>(curnode)->child[1] = new Node();
                    bv_t bv(std::get<1>(curnode)->size());
                    annot_t prefices[2];
                    bit_set(prefices[0], MAXNUM);
                    bit_set(prefices[1], MAXNUM);
                    prefix_t cprefix[2] = {std::make_pair(0, true), std::make_pair(0, true)};
                    bool set[2] = {false, false};
                    size_t lengths[2] = {MAXNUM, MAXNUM};
                    for (size_t i=0;i<bv.size();++i) {
                        bv[i] = array_int::bit_test(std::get<1>(curnode)->operator[](i), length);
                        children[bv[i]]->push_back(std::get<1>(curnode)->operator[](i) >> (length+1));
                        if (children[bv[i]]->back() != prefices[bv[i]]) {
                            if (!set[bv[i]]) {
                                prefices[bv[i]] = array_int::toint(children[bv[i]]->back());
                                set[bv[i]]=true;
                            } else {
                                lengths[bv[i]] = std::min(lengths[bv[i]], (size_t)lsb((children[bv[i]]->back()) ^ prefices[bv[i]]));
                                clear_after(prefices[bv[i]], lengths[bv[i]]);
                            }
                        }
                    }
                    for (size_t i=0;i<2;++i) {
                        std::get<0>(cprefix[i]) = array_int::toint(prefices[i]);
                        if (lengths[i] == MAXNUM) {
                            lengths[i] = (std::get<0>(cprefix[i]) != 0) ? msb((alpha_t)std::get<0>(cprefix[i]))+1 : 0;
                        } else {
                            std::get<1>(cprefix[i]) = false;
                        }
                        //set the last bit to indicate the length of the prefix
                        bit_set(std::get<0>(cprefix[i]), lengths[i]);
                    }
                    beta_t beta(bv);
                    std::get<0>(curnode)->beta = beta;
                    sdsl::util::init_support(std::get<0>(curnode)->rank1, &(std::get<0>(curnode)->beta));
                    sdsl::util::init_support(std::get<0>(curnode)->rank0, &(std::get<0>(curnode)->beta));
                    assert(children[0]->size());
                    assert(children[1]->size());
                    delete std::get<1>(curnode);
                    nodestate_t leftnode = std::make_tuple(std::get<0>(curnode)->child[0], children[0], cprefix[0], std::get<3>(curnode));
                    nodestate_t rightnode = std::make_tuple(std::get<0>(curnode)->child[1], children[1], cprefix[1], std::get<3>(curnode));
                    if (lengths[0] > lengths[1]) {
                        #pragma omp task if (children[0]->size() > TASKMIN)
                        process_node(leftnode);
                        curnode = rightnode;
                    } else {
                        #pragma omp task if (children[1]->size() > TASKMIN)
                        process_node(rightnode);
                        curnode = leftnode;
                    }
                } else {
                    std::get<0>(curnode)->breakleaf(std::get<1>(curnode)->size());
                    delete std::get<1>(curnode);
                    std::get<1>(curnode) = NULL;
                    return true;
                }
            } else {
                std::get<0>(curnode)->breakleaf(std::get<1>(curnode)->size());
                delete std::get<1>(curnode);
                std::get<1>(curnode) = NULL;
                return true;
            }
            //delete std::get<1>(curnode);
            //std::get<1>(curnode) = NULL;
            //if (std::get<0>(return_value.first)) {
                /*
                #pragma omp critical
                {
                    nodes->push_back(return_value.first);
                    //nodes->push_back(return_value.second);
                }
                */
                //#pragma omp task
                //process_node(return_value.first);
            //}
            //curnode = return_value.second;
        }
        return true;
    }

    size_t Node::add(intvec_t *iv, size_t begin, size_t end, prefix_t prefix, bool sorted, bool debug) {
        //std::vector<nodestate_t> nodes;
        //nodes.push_back(std::make_tuple(this, iv, prefix, 0));
        //size_t i=0;
        nodestate_t start = std::make_tuple(this, iv, prefix, 0);
        //#pragma omp parallel
        //#pragma omp single nowait
        #pragma omp task
        process_node(start);
        /*
        {
            while (true) {
                nodestate_t curnode = nodes[i];
                size_t dependi = std::get<3>(curnode);
                //#pragma omp task depend(in:done[dependi]) depend(out:done[i])
                #pragma omp task
                {
                    process_node(&nodes, curnode);
                    //done[i] = true;
                }
                ++i;
                #pragma omp taskwait
                break;
                if (i == nodes.size()) {
                    //If the end has been reached, wait to make sure all pushes have been done
                    #pragma omp taskwait
                    {
                        //std::cerr << "Waiting\n";
                        if (i == nodes.size()) {
                            break;
                        }
                    }
                }
            }
        }
        */
        return 0;
    }

    size_t Node::add(intvec_t::iterator ivb, intvec_t::iterator ive, size_t begin, size_t end, prefix_t prefix, bool sorted, bool debug) {
        if (ive > ivb && end > begin) {
            if (debug)
                std::cout << ive - ivb << "\n";
            if (debug)
                std::cout << ive-ivb << "\t" << std::get<0>(prefix) << "\t" << std::get<1>(prefix) << "\n";
            this->alpha = std::get<0>(prefix);
            size_t length = msb(std::get<0>(prefix));
            if (!std::get<1>(prefix)) {
                intvec_t* children[2] = {new intvec_t(), new intvec_t()};
                //std::cout << ive-ivb << "\t" << std::get<1>(prefix) << "\n";
                if (debug) {
                    assert(this->child[0] == NULL);
                    assert(this->child[1] == NULL);
                }
                this->child[0] = new Node();
                this->child[1] = new Node();
                bv_t bv(ive-ivb);
                //alpha_t temp = 0;
                //bit_set(temp, end);
                //--temp;
                annot_t prefices[2];
                bit_set(prefices[0], end+1);
                bit_set(prefices[1], end+1);
                //temp.set_up_to(end);
                prefix_t cprefix[2] = {std::make_pair(0, true), std::make_pair(0, true)};
                bool set[2] = {false, false};
                size_t lengths[2] = {end+1, end+1};
                for (size_t i=0;i<bv.size();++i) {
                    bv[i] = array_int::bit_test(*(ivb+i), length);
                    children[bv[i]]->push_back((*(ivb+i)) >> (length+1));
                    if (children[bv[i]]->back() != prefices[bv[i]]) {
                    //if (children[bv[i]]->back() != std::get<0>(cprefix[bv[i]])) {
                        if (!set[bv[i]]) {
                            //std::get<0>(cprefix[bv[i]]) = array_int::toint(children[bv[i]]->back());
                            prefices[bv[i]] = array_int::toint(children[bv[i]]->back());
                            set[bv[i]]=true;
                        } else {
                            lengths[bv[i]] = std::min(lengths[bv[i]], (size_t)lsb((children[bv[i]]->back()) ^ prefices[bv[i]]));
                            //lengths[bv[i]] = std::min(lengths[bv[i]], (size_t)lsb((children[bv[i]]->back()) ^ std::get<0>(cprefix[bv[i]])));
                            //std::get<0>(cprefix[bv[i]]) &= (alpha_t(1) << lengths[bv[i]])-1;
                            //clear_after(std::get<0>(cprefix[bv[i]]), lengths[bv[i]]);
                            clear_after(prefices[bv[i]], lengths[bv[i]]);
                        }
                    }
                }
                for (size_t i=0;i<2;++i) {
                    std::get<0>(cprefix[i]) = array_int::toint(prefices[i]);
                    //if (lengths[i] == MAXNUM) {
                    if (lengths[i] == end + 1) {
                        lengths[i] = (std::get<0>(cprefix[i]) != 0) ? msb((alpha_t)std::get<0>(cprefix[i]))+1 : 0;
                    } else {
                        std::get<1>(cprefix[i]) = false;
                    }
                    //set the last bit to indicate the length of the prefix
                    bit_set(std::get<0>(cprefix[i]), lengths[i]);
                }
                beta_t beta(bv);
                this->beta = beta;
                sdsl::util::init_support(this->rank1, &(this->beta));
                sdsl::util::init_support(this->rank0, &(this->beta));
                assert(children[0]->size());
                assert(children[1]->size());
                #pragma omp task if (children[0]->size() > TASKMIN)
                this->child[0]->add(children[0]->begin(), children[0]->end(), begin+length+1, end, cprefix[0], sorted, debug);
                #pragma omp taskwait
                delete children[0];
                #pragma omp task if (children[1]->size() > TASKMIN)
                this->child[1]->add(children[1]->begin(), children[1]->end(), begin+length+1, end, cprefix[1], sorted, debug);
                #pragma omp taskwait
                delete children[1];
                #pragma omp taskwait
                {
                    //delete children[0];
                    //delete children[1];
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
        } else {
            this->breakleaf(ive-ivb);
        }
        return 0lu;
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

    //TODO: rewrite without recursion
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
    WTR::WTR(intvec_t *iv, bool transpose, size_t osize, size_t batch, bool sorted, bool debug) {
        //this->sorted = sorted;
        if (std::getenv("DEBUG") != NULL)
            debug=true;
        this->root = new Node();
        if (batch != 0 && iv->size() != 0) {
            size_t tsize=0;
            prefix_t prefix;
            //std::cout << "Computing alpha\n";
            prefix = longest_prefix(iv->begin(), std::min(iv->begin()+batch, iv->end()), 0, osize, sorted, debug); 
            //std::cout << "Computing wavelet trie\n";
            tsize = iv->size();
            #pragma omp parallel
            #pragma omp single nowait
            this->root->add(iv, 0, osize, prefix, sorted, debug); 
            //this->root->add(iv.begin(), std::min(iv.begin()+batch, iv.end()), 0, osize, prefix, sorted, debug); 
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
                //#pragma omp parallel
                //#pragma omp single nowait
                //other->add(iv.begin()+cursize, std::min(iv.begin()+cursize+batch, iv.end()), 0, osize, sorted, debug); 
                //std::cout << "Computing wavelet trie\n";
                prefix = longest_prefix(iv->begin()+cursize, std::min(iv->begin()+cursize+batch, iv->end()), 0, osize, sorted, debug); 
                #pragma omp parallel
                #pragma omp single nowait
                others[cursize/batch-1]->add(iv, 0, osize, prefix, sorted, debug); 
                //others[cursize/batch-1]->add(iv.begin()+cursize, std::min(iv.begin()+cursize+batch, iv.end()), 0, osize, prefix, sorted, debug); 
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

    std::pair<intvec_t*, intvec_t*> construct(std::ifstream &pfile, size_t begin=0, size_t end=MAXNUM, bool debug=false) {
        intvec_t* iv = new intvec_t();
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
            iv->push_back(cannot);
            if (iv->size() == begin)
                break;
            if (debug)
                std::cout << "\t" << iv->back() << "\n";
        }
        if (debug)
            std::cout << "\n";
        return std::make_pair(iv, new intvec_t(maxedi,0));
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
