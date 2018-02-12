#include "wavelet_trie.hpp"

namespace annotate {

    bool bit_test(const cpp_int &a, const size_t &col) {
        assert(col < -1llu);
        return mpz_tstbit(a.backend().data(), col);
    }

    void bit_set(mpz_t &a_d, const size_t &col) {
        assert(col < -1llu);
        mpz_setbit(a_d, col);
    }

    void bit_set(cpp_int &a, const size_t &col) {
        assert(col < -1llu);
        mpz_setbit(a.backend().data(), col);
    }

    void bit_unset(cpp_int &a, const size_t &col) {
        assert(col < -1llu);
        mpz_clrbit(a.backend().data(), col);
    }

    size_t next_bit(const cpp_int &a, const size_t &col) {
        assert(col < -1llu);
        return mpz_scan1(a.backend().data(), col);
    }

    void clear_after(mpz_t &a_d, const size_t &col) {
        assert(col < -1llu);
        if (!col) {
            mpz_clear(a_d);
            mpz_init(a_d);
        } else {
            mpz_tdiv_r_2exp(a_d, a_d, col);
        }
    }

    void clear_after(cpp_int &a, const size_t &col) {
        assert(col < -1llu);
        mpz_t& a_d = a.backend().data();
        clear_after(a_d, col);
    }

    size_t msb(const cpp_int &a) {
        assert(a != 0);
        const mpz_t& a_mpz = a.backend().data();
        size_t i = mpz_scan1(a_mpz, 0);
        do {
            size_t j = mpz_scan1(a_mpz, i + 1);
            if (j == -1llu)
                break;
            i = j;
        } while (true);
        return i;
    }

    size_t lsb(const cpp_int &a) {
        assert(a != 0);
        return next_bit(a, 0);
    }

    //TODO: align get_int to blocks in source
    template <typename Vector>
    bv_t insert_zeros(const Vector &target, const size_t count, const size_t i) {
        if (!count)
            return target;
        if (!target.size()) {
            return bv_t(count);
        }
        bv_t merged;
        size_t j;
        //assume limb size of 64
        if (i) {
            merged = target;
            merged.resize(target.size() + count);
            /*
            j = 0;
            for (; j + 64 <= count; j += 64) {
                merged.set_int(i + j, 0);
            }
            if (count - j)
                merged.set_int(i + j, 0, count - j);
            */
            j = ((i + 63) & -64llu);
            //WARNING: the region from merged.size() to j is not initialized
            merged.set_int(i, 0, std::min(j, merged.size()) - i);
            std::fill(merged.data() + (j >> 6), merged.data() + (merged.capacity() >> 6), 0);
        } else {
            merged = bv_t(count);
            merged.resize(target.size() + count);
        }
        j = i;
        for (; j + 64 <= target.size(); j += 64) {
            merged.set_int(count + j, target.get_int(j));
        }
        if (target.size() - j)
            merged.set_int(count + j, target.get_int(j, target.size() - j), target.size() - j);
        return merged;
    }

    template <typename Vector>
    bv_t insert_range(const Vector &target, const Vector &source, const size_t i) {
        if (!target.size()) {
            assert(i == 0);
            return source;
        }
        if (!source.size())
            return target;

        bv_t merged;
        size_t j;
        if (i) {
            merged = target;
            merged.resize(target.size() + source.size());
            j = 0;
            for (; j + 64 <= source.size(); j += 64) {
                merged.set_int(i + j, source.get_int(j));
            }
            if (source.size() - j)
                merged.set_int(i + j, source.get_int(j, source.size() - j), source.size() - j);
        } else {
            merged = source;
            merged.resize(target.size() + source.size());
        }
        j = i;
        for (; j + 64 <= target.size(); j += 64) {
            merged.set_int(source.size() + j, target.get_int(j));
        }
        if (target.size() - j)
            merged.set_int(source.size() + j, target.get_int(j, target.size() - j), target.size() - j);
        /*
        //TODO: move this to a unit test
        //super slow reference implementation
        std::string target_s = sdsl::util::to_string(target);
        std::string source_s = sdsl::util::to_string(source);
        std::string merged_s = target_s.substr(0, i) + source_s + target_s.substr(i);
        bv_t merged_test(target.size() + source.size());
        for (size_t j = 0; j < merged_s.length(); ++j) {
            if (merged_s[j] == '1')
                merged_test[j] = 1;
        }
        assert(merged == merged_test);
        */
        return merged;
    }

    WaveletTrie::WaveletTrie() : root(NULL) { }

    void WaveletTrie::serialize(std::ostream &out) const {
        root->serialize(out);
    }

    void WaveletTrie::Node::serialize(std::ostream &out) const {
        //serialize alpha
        size_t a;
        void *alpha_raw = mpz_export(NULL, &a, 1, 1, 0, 0, alpha_.backend().data());
        out.write((char*)&a, sizeof(a));
        out.write((char*)alpha_raw, a);
        free(alpha_raw);

        //beta
        rrr_t(beta_).serialize(out);

        std::string inds("01");
        if (child_[0]) {
            out.write(inds.data() + 1, 1);
            child_[0]->serialize(out);
        } else {
            out.write(inds.data(), 1);
        }
        if (child_[1]) {
            out.write(inds.data() + 1, 1);
            child_[1]->serialize(out);
        } else {
            out.write(inds.data(), 1);
        }
    }

    void WaveletTrie::print() {
        std::stack<std::pair<Node*,std::string>> nodestack;
        nodestack.emplace(root,"");
        while (nodestack.size()) {
            Node *curnode = nodestack.top().first;
            std::string parent = nodestack.top().second;
            nodestack.pop();
            std::cout << parent << ":\t" << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->all_zero << "\n";
            if (curnode->child_[0])
                nodestack.emplace(curnode->child_[0], parent + std::string("L"));
            if (curnode->child_[1])
                nodestack.emplace(curnode->child_[1], parent + std::string("R"));
        }
    }

    WaveletTrie::Node::Node(const WaveletTrie::Node &that)
        : alpha_(that.alpha_), beta_(that.beta_),
          rank1_(that.rank1_), rank0_(that.rank0_),
          all_zero(that.all_zero),
          popcount(that.popcount),
          support(that.support) {
        if (that.child_[0]) {
            child_[0] = new Node(*that.child_[0]);
        }
        if (that.child_[1]) {
            child_[1] = new Node(*that.child_[1]);
        }
    }

    template <class Iterator>
    WaveletTrie::Node::Node(const Iterator &row_begin, const Iterator &row_end,
            const size_t &col, Prefix prefix) {
        if (row_end > row_begin) {
            assert(prefix.col != -1llu);
            assert(!prefix.allequal);
            size_t col_end = prefix.col;

            //set alpha
            //alpha_ = (*row_begin >> col) % mask; //reference
            mpz_t& alph = alpha_.backend().data();
            mpz_tdiv_q_2exp(alph, row_begin->backend().data(), col);
            clear_after(alph, col_end - col);
            bit_set(alph, col_end - col);

            //set beta and compute common prefices
            bv_t beta;
            beta.resize(row_end - row_begin);
            Prefix prefices[2];
            Iterator split = row_begin;
            std::vector<cpp_int> right_children;
            sdsl::util::set_to_value(beta, 0);
            for (auto it = row_begin; it != row_end; ++it) {
                if (bit_test(*it, col_end)) {
                    beta[it - row_begin] = 1;
                    popcount++;
                    right_children.emplace_back(*it);
                    prefices[1].col = next_different_bit_(
                            right_children.begin(), right_children.end() - 1,
                            col_end + 1, prefices[1].col);
                } else {
                    *split = *it;
                    prefices[0].col = next_different_bit_(
                            row_begin, split,
                            col_end + 1, prefices[0].col);
                    split++;
                }
            }
            set_beta_(beta);
            assert(popcount == rank1(size()));
            //distribute to left and right children
            assert(split != row_begin && split != row_end);
            std::swap_ranges(right_children.begin(), right_children.end(), split);

            //TODO: copied code here
            if (prefices[0].col != -1llu) {
                prefices[0].allequal = false;
                child_[0] = new Node(row_begin, split, col_end + 1, prefices[0]);
                assert(child_[0]->size() == rank0(beta.size()));
            } else {
                cpp_int new_alpha = *row_begin;
                mpz_t& alph = new_alpha.backend().data();
                mpz_tdiv_q_2exp(alph, alph, col_end + 1);
                if (new_alpha) {
                    bit_set(new_alpha, msb(new_alpha) + 1);
                } else {
                    bit_set(new_alpha, 0);
                }
#pragma omp task
                child_[0] = new Node(new_alpha, split - row_begin);
                //assert(child_[0]->size() == rank0(beta.size()));
            }

            if (prefices[1].col != -1llu) {
                prefices[1].allequal = false;
                child_[1] = new Node(split, row_end, col_end + 1, prefices[1]);
                assert(child_[1]->size() == rank1(beta.size()));
            } else {
                cpp_int new_alpha = *split;
                mpz_t& alph = new_alpha.backend().data();
                mpz_tdiv_q_2exp(alph, alph, col_end + 1);
                if (new_alpha) {
                    bit_set(new_alpha, msb(new_alpha) + 1);
                } else {
                    bit_set(new_alpha, 0);
                }
#pragma omp task
                child_[1] = new Node(new_alpha, row_end - split);
                //assert(child_[1]->size() == rank1(beta.size()));
            }
        }
    }

    void WaveletTrie::Node::swap(WaveletTrie::Node&& that) {
        this->alpha_ = that.alpha_;
        this->beta_ = that.beta_;
        this->all_zero = that.all_zero;
        this->popcount = that.popcount;
        this->rank1_ = that.rank1_;
        this->rank0_ = that.rank0_;
        this->support = that.support;
        this->child_[0] = that.child_[0];
        this->child_[1] = that.child_[1];
        that.child_[0] = NULL;
        that.child_[1] = NULL;
    }

    WaveletTrie::Node::Node(const size_t count)
      : all_zero(true) {
        set_beta_(beta_t(count));
    }

    WaveletTrie::Node::Node(const cpp_int &alpha, const size_t count)
      : Node(count) {
        alpha_ = alpha;
    }

    template <class Iterator>
    WaveletTrie::WaveletTrie(Iterator row_begin, Iterator row_end) {
        if (row_end > row_begin) {
            Prefix prefix = WaveletTrie::Node::longest_common_prefix(row_begin, row_end, 0);
            if (prefix.allequal) {
                cpp_int alpha = *row_begin;
                if (alpha) {
                    bit_set(alpha, msb(alpha) + 1);
                } else {
                    bit_set(alpha, 0);
                }
                root = new Node(alpha, row_end - row_begin);
            } else {
#pragma omp parallel
#pragma omp single nowait
                root = new Node(row_begin, row_end, 0, prefix);
            }
        } else {
            root = NULL;
        }
#ifndef NPRINT
        print();
        std::cout << "\n";
#endif
    }

    WaveletTrie::~WaveletTrie() {
        delete root;
    }

    cpp_int WaveletTrie::at(size_t i) {
        assert(i < size());
        Node *node = root;
        size_t length = 0;
        cpp_int annot;
        while (!node->is_leaf()) {
            annot |= node->alpha_ << length;
            length += msb(node->alpha_) + 1;
            if (node->beta_[i]) {
                assert(node->child_[1]);
                i = node->rank1(i);
                node = node->child_[1];
            } else {
                bit_unset(annot, length - 1);
                assert(node->child_[0]);
                i = node->rank0(i);
                node = node->child_[0];
            }
        }
        annot |= node->alpha_ << length;
        bit_unset(annot, msb(annot));
        return annot;
    }

    size_t WaveletTrie::size() {
        return root ? root->size() : 0;
    }

    bool WaveletTrie::Node::is_leaf() {
        return all_zero;
    }

    bool WaveletTrie::Node::check(bool ind) {
        size_t rank = ind ? popcount : size() - popcount;
        assert(rank1(size()) == popcount);
        if (child_[ind]) {
            if (all_zero)
                return false;
            if (child_[ind]->size() != rank)
                return false;
            return child_[ind]->check(0) && child_[ind]->check(1);
        } else {
            if (!all_zero && (!child_[!ind] || rank > 0))
                return false;
        }
        return true;
    }

    void WaveletTrie::Node::fill_left(bool rightside) {
        size_t lrank = size() - popcount;
        assert(lrank == rank0(size()));
        Node *jnode = NULL;
        if (lrank) {
            jnode = child_[0];
            while (jnode && lrank > jnode->size()) {
                jnode->move_label_down_(lsb(jnode->alpha_));
                jnode->set_beta_(insert_zeros(jnode->beta_, lrank - jnode->size(), rightside ? jnode->size() : 0));
                assert(jnode->popcount == jnode->rank1(jnode->size()));
                lrank -= jnode->popcount;
                if (!jnode->child_[0] && lrank) {
                    jnode->all_zero = false;
                    jnode->child_[0] = new Node(lrank);
                    break;
                }
                jnode = jnode->child_[0];
            }
            assert(!child_[0] || child_[0]->check(0));
        }
    }

    void WaveletTrie::Node::fill_ancestors(Node *othnode, bool ind, const size_t i) {
        if (child_[ind]) {
            if (!othnode->child_[ind] && !ind) {
                //TODO: fix position when i != size()
                fill_left(true);
            }
        } else {
            assert(child_[!ind] || all_zero);
            if (othnode->child_[ind]) {
                std::swap(child_[ind], othnode->child_[ind]);
                all_zero = othnode->all_zero;
                if (!ind) {
                    //TODO: correct position when i != size() ?
                    fill_left(false);
                }
                assert(!all_zero);
            } else if (!ind) {
                size_t lrank = size() - popcount;
                assert(lrank == rank0(size()));
                if (lrank) {
                    all_zero = false;
                    child_[ind] = new Node(lrank);
                }
            }
        }
    }

    void WaveletTrie::insert(WaveletTrie&& wtr, size_t i) {
        if (!wtr.root) {
            return;
        }
        if (!root) {
            std::swap(root, wtr.root);
            return;
        }
        assert(root && wtr.root);
#pragma omp parallel
#pragma omp single nowait
        Node::merge(root, wtr.root, i);
    }

    void WaveletTrie::Node::merge(Node *curnode, Node *othnode, size_t i) {
        if (i == -1llu) {
            i = curnode->size();
        }


        assert(curnode->size());
        assert(othnode->size());
        assert(i <= curnode->size());

        while (curnode && othnode) {

            //assert(curnode && othnode);
            assert(curnode->size());
            assert(i <= curnode->size());
            assert(curnode->check(0));
            assert(curnode->check(1));
#ifndef NPRINT
            std::cout << path << "\t" << i << "\t"
                      << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->all_zero << "\t"
                      << othnode->alpha_ << ":" << othnode->beta_ << ";" << othnode->all_zero << "\t->\t";
#endif

            Node::overlap_prefix_(curnode, othnode);

            //update insertion point and merge betas
            size_t il = i == curnode->size()
                ? curnode->size() - curnode->popcount
                : curnode->rank0(i);
            size_t ir = i == curnode->size()
                ? curnode->popcount
                : curnode->rank1(i);
            Node::merge_beta_(curnode, othnode, i);

#ifndef NPRINT
            std::cout << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->all_zero << "\t"
                      << othnode->alpha_ << ":" << othnode->beta_ << ";" << othnode->all_zero << "\n";
#endif
            bool left  = curnode->child_[0] && othnode->child_[0];
            bool right = curnode->child_[1] && othnode->child_[1];
            curnode->fill_ancestors(othnode, 0, il);
            curnode->fill_ancestors(othnode, 1, ir);
            if (left) {
                assert(!curnode->all_zero);
                assert(il <= curnode->child_[0]->size());
                assert(curnode->size() - curnode->popcount
                        == curnode->child_[0]->size() + othnode->child_[0]->size()
                );
#pragma omp task if (curnode->child_[0]->size() > 32)
                Node::merge(curnode->child_[0], othnode->child_[0], il);
            }
            if (right) {
                assert(!curnode->all_zero);
                assert(ir <= curnode->child_[0]->size());
                assert(curnode->popcount == curnode->child_[1]->size() + othnode->child_[1]->size());
                curnode = curnode->child_[1];
                othnode = othnode->child_[1];
                i = ir;
            } else {
                curnode = NULL;
                othnode = NULL;
            }
        }
    }

    template <typename T>
    void WaveletTrie::insert(const T &a, size_t i) {
        Node *next = new Node(&a, &a + 1, 0);
        if (!i) {
            root = next;
        } else {
            WaveletTrie wtr;
            wtr.root = next;
            insert(wtr, i);
        }
    }

    WaveletTrie::Node::~Node() {
        if (child_[0]) {
            delete child_[0];
        }
        if (child_[1]) {
            delete child_[1];
        }
    }

    //return true if a change happened
    bool WaveletTrie::Node::overlap_prefix_(Node *curnode, Node *othnode) {
        assert(curnode && othnode);
        assert(curnode->alpha_ && othnode->alpha_);

        if (curnode->alpha_ == othnode->alpha_) {
            return false;
        }
        size_t common_pref = next_different_bit_alpha(curnode, othnode);
#ifndef NPRINT
        std::cout << common_pref << "\t";
#endif
        int cur = curnode->move_label_down_(common_pref);
#ifndef NPRINT
        std::cout << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->all_zero;
        if (cur > -1) {
            assert(!curnode->all_zero);
            std::cout << "," << curnode->child_[cur]->alpha_ << ":" << curnode->child_[cur]->beta_ << ";" << curnode->child_[cur]->all_zero;
        }
        std::cout << "\t";
#endif
        int oth = othnode->move_label_down_(common_pref);
#ifndef NPRINT
        std::cout << othnode->alpha_ << ":" << othnode->beta_ << ";" << othnode->all_zero;
        if (oth > -1) {
            assert(!othnode->all_zero);
            std::cout << "," << othnode->child_[oth]->alpha_ << ":" << othnode->child_[oth]->beta_ << ";" << othnode->child_[oth]->all_zero;
        }
        std::cout << "\t->\t";
#endif
        cur++;
        oth++;
        assert(curnode->alpha_ == othnode->alpha_);
        return true;
    }

    //find first different bit between two cpp_ints
    template <class Iterator>
    size_t WaveletTrie::Node::next_different_bit_(const Iterator &a, const Iterator &b,
            const size_t col, size_t next_col) {
        if (col == next_col)
            return next_col;
        size_t ranges[3] = {next_bit(*a, col), next_bit(*b, col), next_col};
        while (ranges[0] == ranges[1] && ranges[0] < ranges[2] && ranges[1] < ranges[2]) {
            ranges[0] = next_bit(*a, ranges[0] + 1);
            ranges[1] = next_bit(*b, ranges[1] + 1);
        }
        next_col = std::min(ranges[0], ranges[1]);
        if (ranges[0] != ranges[1] && next_col < ranges[2]) {
            return next_col;
        }
        return ranges[2];
    }

    size_t WaveletTrie::Node::next_different_bit_alpha(Node *curnode, Node *othnode) {
        assert(curnode->alpha_ != 0);
        assert(othnode->alpha_ != 0);
        mpz_t& cur = curnode->alpha_.backend().data();
        mpz_t& oth = othnode->alpha_.backend().data();
        size_t ranges[4] = {next_bit(cur, 0), next_bit(oth, 0), 0, 0};
        while (ranges[2] == ranges[3] && ranges[2] < -1llu && ranges[3] < -1llu) {
            ranges[2] = next_bit(ranges[2], ranges[0] + 1);
            ranges[3] = next_bit(ranges[3], ranges[1] + 1);
            if (ranges[2] < -1llu)
                ranges[0] = ranges[2];
            if (ranges[3] < -1llu)
                ranges[1] = ranges[3];
        }
        ranges[1] = std::min(ranges[0], ranges[1]);
        if (!curnode->all_zero || !othnode->all_zero) {
            return *(std::min_element(ranges + 1, ranges + 4));
        }
        return ranges[1];
    }

    template <class Iterator>
    Prefix WaveletTrie::Node::longest_common_prefix(const Iterator &row_begin, const Iterator &row_end, const size_t &col) {
        Prefix prefix;
        //empty prefix
        if (row_begin >= row_end) {
            prefix.col = col;
            prefix.allequal = true;
            return prefix;
        }
        //prefix.col = -1llu;
        for (auto it = row_begin + 1; it != row_end; ++it) {
            prefix.col = next_different_bit_(row_begin, it, col, prefix.col);
            if (prefix.col == col)
                break;
        }
        if (prefix.col == -1llu) {
            //all zeros or all equal
            if (*row_begin == 0) {
                prefix.col = col;
            } else {
                prefix.col = std::max(msb(*row_begin), col);
            }
            return prefix;
        }
        prefix.allequal = false;
        return prefix;
    }

    int WaveletTrie::Node::move_label_down_(size_t length) {
        size_t len = msb(alpha_);
        if (length > len) {
            bit_unset(alpha_, len);
            bit_set(alpha_, length);
        } else if (length < len) {
            //split alpha and compute new beta
            //TODO: clean up
            Node *child = new Node();
            mpz_t& child_alpha = child->alpha_.backend().data();
            mpz_tdiv_q_2exp(child_alpha, alpha_.backend().data(), length + 1);
            //TODO: replace with operator=
            child->set_beta_(beta_);
            child->child_[0] = child_[0];
            child->child_[1] = child_[1];
            child->all_zero = all_zero;
            child->popcount = popcount;
            all_zero = false;
            bool beta_bit = bit_test(alpha_, length);
            set_beta_(bv_t(size(), beta_bit));
            //only want length bits left
            clear_after(alpha_, length);
            bit_set(alpha_, length);
            if (beta_bit) {
                child_[0] = NULL;
                child_[1] = child;
                popcount = size();
                assert(beta_.size() == child_[1]->size());
            } else {
                child_[0] = child;
                child_[1] = NULL;
                popcount = 0;
                assert(beta_.size() == child_[0]->size());
            }
            assert(popcount == rank1(size()));
            assert(msb(alpha_) == length);
            assert((child_[0] == NULL) ^ (child_[1] == NULL));
        }
        assert(check(0));
        assert(check(1));
        assert((child_[0] && !all_zero)
            || (child_[1] && !all_zero)
            || (!child_[0] && !child_[1] && all_zero));
        if (length >= len)
            return -1;
        if (child_[0])
            return 0;
        return 1;
    }

    void WaveletTrie::Node::merge_beta_(Node *curnode, Node *othnode, size_t i) {
        assert(othnode);
        assert(curnode->alpha_ == othnode->alpha_);
        if (!othnode->size()) {
            return;
        }
        if (i == -1llu) {
            i = curnode->beta_.size();
        }
        assert(i <= curnode->beta_.size());
#ifndef NPRINT
        std::cout << i << "\t";
#endif
        auto beta_new = insert_range(curnode->beta_, othnode->beta_, i);
        assert(beta_new.size() == curnode->size() + othnode->size());
#ifndef NDEBUG
        //SANITY CHECK
        //TODO: move to unit test
        size_t j = 0;
        for (; j < i; ++j) {
            assert(beta_new[j] == curnode->beta_[j]);
        }
        for (; j - i < othnode->size(); ++j) {
            assert(beta_new[j] == othnode->beta_[j - i]);
        }
        for (; j < beta_new.size(); ++j) {
            assert(beta_new[j] == curnode->beta_[j - othnode->size()]);
        }
#endif
        curnode->set_beta_(beta_new);
        curnode->popcount += othnode->popcount;
        assert(curnode->popcount == curnode->rank1(curnode->size()));
    }

    template <class Vector>
    void WaveletTrie::Node::set_beta_(const Vector &bv) {
        assert(bv.size());
        beta_ = beta_t(bv);
    }

    size_t WaveletTrie::Node::rank0(const size_t i) {
        if (!support) {
            sdsl::util::init_support(rank1_, &beta_);
            sdsl::util::init_support(rank0_, &beta_);
            support = true;
        }
        return rank0_(i);
    }

    size_t WaveletTrie::Node::rank1(const size_t i) {
        if (!support) {
            sdsl::util::init_support(rank1_, &beta_);
            sdsl::util::init_support(rank0_, &beta_);
            support = true;
        }
        return rank1_(i);
    }

    template WaveletTrie::WaveletTrie(std::vector<cpp_int>::iterator, std::vector<cpp_int>::iterator);


};
