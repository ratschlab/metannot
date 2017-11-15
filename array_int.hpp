#ifndef __ARRAYINT__
#define __ARRAYINT__

#include <iostream>
#include <set>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

namespace array_int {

    typedef std::set<size_t> indset_type;
    class array_int {
        public:
            indset_type indices;
        public:
            array_int();
            array_int(size_t i);
            array_int(boost::multiprecision::cpp_int i);
            array_int(array_int &a);
            array_int(const array_int &a);
            template <typename T>
            array_int operator<<(T i) const;
            template <typename T>
            array_int operator>>(T i) const;
            template <typename T>
            array_int operator&(T i) const;
            template <typename T>
            array_int operator^(T i) const;
            template <typename T>
            bool operator>(T i) const;
            template <typename T>
            array_int& operator|=(T i);
            size_t popcount();
            bool operator==(const size_t& b) const;
            bool operator==(const array_int& b) const;
            template <typename T>
            bool operator!=(const T& b) const;
            array_int& operator&=(const array_int &that);
            bool operator<(const array_int &that) const;
            //template <typename T>
            //array_int operator+(T b) const;
            friend std::ostream& operator<<(std::ostream& os, const array_int& a);
            //array_int& operator--();
            indset_type& backend();
            void set_up_to(size_t i);
    };

    array_int::array_int() {
    }

    array_int::array_int(boost::multiprecision::cpp_int i) {
        size_t lsbit;
        if (i != 0) {
            while (i) {
                lsbit = boost::multiprecision::lsb(i);
                this->indices.insert(lsbit);
                boost::multiprecision::bit_unset(i, lsbit);
            }
        }
    }


    array_int::array_int(array_int &a) {
        this->indices.insert(a.indices.begin(), a.indices.end());
    }
    
    array_int::array_int(const array_int &a) {
        this->indices.insert(a.indices.begin(), a.indices.end());
    }

    void array_int::set_up_to(size_t i) {
        for (size_t j=0;j<i;++j)
            this->indices.insert(this->indices.end(), j);
    }
    
    boost::multiprecision::cpp_int& toint(boost::multiprecision::cpp_int &a) {
        return a;
    }
    boost::multiprecision::cpp_int toint(const array_int &a) {
        boost::multiprecision::cpp_int temp=0;
        for (auto it = a.indices.begin(); it != a.indices.end(); ++it) {
            boost::multiprecision::bit_set(temp, *it);
        }
        return temp;
    }


    indset_type& array_int::backend() {
        return this->indices;
    }
    
    std::ostream& operator<<(std::ostream& os, const array_int& a) {
        //TODO
        os << toint(a);
        return os;
    }

    array_int& array_int::operator&=(const array_int &that) {
        indset_type newinds;
        std::set_intersection(this->indices.begin(), this->indices.end(), that.indices.begin(), that.indices.end(), std::inserter(newinds, newinds.begin()));
        this->indices = newinds;
        return *this;
    }

    /*
    template <typename T>
    array_int array_int::operator+(T b) const {
        //TODO
        return array_int(toint(*this) + b);
    }
    */

    template <typename T>
    array_int array_int::operator&(T b) const {
        //TODO
        array_int a = *this;
        a &= array_int(b);
        return a;
    }

    template <typename T>
    array_int array_int::operator^(T b) const {
        //TODO
        return array_int(toint(*this) ^ toint(b));
    }

    template <typename T>
    array_int& array_int::operator|=(T b) {
        //TODO
        auto a = array_int(b);
        this->indices.insert(a.indices.begin(), a.indices.end());
        return *this;
    }

    template <typename T>
    bool array_int::operator>(T b) const {
        //TODO
        return toint(*this) > toint(b);
    }

    array_int::array_int(size_t i) {
        size_t ind = 0;
        while (i > 0lu) {
            if (i % 2lu) {
                this->indices.insert(ind);
            }
            ind++;
            i >>= 1lu;
        }
    }
    template <typename T>
    array_int array_int::operator<<(T i) const {
        if (i == 0)
            return *this;
        assert(i > 0);
        array_int newinds;
        for (auto it = this->indices.begin(); it != this->indices.end(); ++it)
            newinds.indices.insert(*it + i);
        return newinds;
    }
    template <typename T>
    array_int array_int::operator>>(T i) const {
        if (i == 0)
            return *this;
        assert(i > 0);
        array_int newinds;
        for (auto it = this->indices.begin(); it != this->indices.end(); ++it) {
            if (*it >= i)
                newinds.indices.insert(*it - i);
        }
        return newinds;
    }

    bool array_int::operator==(const size_t& b) const {
        if (b == 0) {
            return !(this->indices.size());
        }
        boost::multiprecision::cpp_int bigb(b);
        size_t lsb;
        auto it = this->indices.begin();
        while (bigb != 0 && it != this->indices.end()) {
            lsb = boost::multiprecision::lsb(bigb);
            if (lsb != *it) {
                return false;
            } else {
                boost::multiprecision::bit_unset(bigb, lsb);
                ++it;
            }
        }
        return bigb == 0 && it == this->indices.end();
    }

    bool array_int::operator==(const array_int &b) const {
        if (this->indices.size() != b.indices.size())
            return false;
        auto it = this->indices.begin();
        auto jt = b.indices.begin();
        while (it != this->indices.end() && jt != b.indices.end()) {
            if (*it != *jt)
                return false;
            ++it;
            ++jt;
        }
        return true;
    }

    template <typename T>
    bool array_int::operator!=(const T& b) const {
        return !(this->operator==(b));
    }

    bool array_int::operator<(const array_int &that) const {
        auto it = this->indices.rbegin();
        auto jt = that.indices.rbegin();
        while (it != this->indices.rend() && jt != that.indices.rend()) {
            if (*it != *jt) {
                return *it < *jt;
            }
            ++it;
            ++jt;
        }
        if (it == this->indices.rend()) {
            return jt != that.indices.rend();
        } else {
            return false;
        }
    }
   
    /*
    array_int& array_int::operator--() {
        //TODO
        auto a = --toint(*this);
        auto b = array_int(a);
        this->indices = b.indices;
        return *this;
    }
    */



bool bit_test(const array_int& i, size_t j) {
    return i.indices.find(j) != i.indices.end();
}

}

array_int::array_int bit_set(array_int::array_int& i, size_t j) {
    i.indices.insert(j);
    return i;
}

array_int::array_int bit_unset(array_int::array_int& i, size_t j) {
    i.indices.erase(j);
    return i;
}

size_t lsb(const array_int::array_int& a) {
    assert(a.indices.size());
    return *a.indices.begin();
}

size_t msb(const array_int::array_int& a) {
    assert(a.indices.size());
    return *a.indices.rbegin();
}

#endif
