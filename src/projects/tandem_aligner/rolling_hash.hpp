#pragma once
//
// Created by anton on 7/20/20.
//

#include <sequences/sequence.hpp>
#include <deque>

template<typename T, typename U>
T pow(T base, U p) {
    if (p == 0)
        return 1;
    T tmp = pow(base, p / 2);
    if (p %2 == 1)
        return base * tmp * tmp;
    else
        return tmp * tmp;
}

template<typename htype>
class RollingHash {
public:
    const size_t k;
    const htype hbase;
    const htype kpow;
    const htype inv;

    RollingHash(size_t _k, htype _hbase) : k(_k), hbase(_hbase),
            kpow(pow(hbase, k - 1)), inv(pow(hbase, (htype(1u) << (sizeof(htype) * 8u - 1u)) - 1u)){
        VERIFY(inv * hbase == htype(1));
    }

    RollingHash extensionHash() const {
        return RollingHash(k + 1, hbase);
    }

    htype hash(const Sequence &seq, size_t pos) const {
        htype hash = 0;
        for(size_t i = pos; i < pos + k; i++) {
            hash = hash * hbase + seq[i];
        }
        return hash;
    }

    htype extendRight(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
        return hash * hbase + c;
    }

    htype extendLeft(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
        return hash + c * kpow * hbase;
    }

    htype shiftRight(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
        return (hash - kpow * seq[pos]) * hbase + c;
    }

    htype shiftLeft(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
        return (hash - seq[pos + k - 1]) * inv + c * kpow;
    }

    htype next(const Sequence &seq, size_t pos, htype hash) const {
        return shiftRight(seq, pos, hash, seq[pos + k]);
    }

    htype prev(const Sequence &seq, size_t pos, htype hash) const {
        return shiftLeft(seq, pos, hash, seq[pos - 1]);
    }

    bool hasNext(const Sequence &seq, size_t pos) const {
        return pos + k < seq.size();
    }

    bool hasPrev(const Sequence &seq, size_t pos) const {
        return pos > 0;
    }
};

template<typename htype>
class KWH {
private:
    KWH(const RollingHash<htype> & _hasher, const Sequence &_seq, size_t _pos, htype _fhash, htype _rhash):
            hasher(_hasher), seq(_seq), pos(_pos), fhash(_fhash), rhash(_rhash) {
    }

    htype fhash;
    htype rhash;
    Sequence seq;
public:
    const RollingHash<htype> & hasher;
    size_t pos;

    KWH(const RollingHash<htype> & _hasher, const Sequence &_seq, size_t _pos):
            hasher(_hasher), seq(_seq), pos(_pos), fhash(_hasher.hash(_seq, _pos)),
            rhash(_hasher.hash(!_seq, _seq.size() - _pos - _hasher.k)) {
    }

    KWH(const KWH &other) = default;

    Sequence getSeq() const {
        return seq.Subseq(pos, pos + hasher.k);
    }

    KWH<htype> operator!() const {
        return KWH<htype>(hasher, !seq, seq.size() - pos - hasher.k, rhash, fhash);
    }

    htype hash() const {
        return std::min(fhash, rhash);
    }

    htype extendRight(unsigned char c) const {
        return std::min(hasher.extendRight(seq, pos, fhash, c), hasher.extendLeft(!seq, seq.size() - pos - hasher.k, rhash, c ^ 3u));
    }

    htype extendLeft(unsigned char c) const {
        return std::min(hasher.extendLeft(seq, pos, fhash, c), hasher.extendRight(!seq, seq.size() - pos - hasher.k, rhash, c ^ 3u));
    }

    KWH next() const {
        return {hasher, seq, pos + 1, hasher.next(seq, pos, fhash), hasher.prev(!seq, seq.size() - pos - hasher.k, rhash)};
    }

    KWH prev() const {
        return {hasher, seq, pos - 1, hasher.prev(seq, pos, fhash), hasher.next(!seq, seq.size() - pos - hasher.k, rhash)};
    }

    bool hasNext() const {
        return hasher.hasNext(seq, pos);
    }

    bool hasPrev() const {
        return hasher.hasPrev(seq, pos);
    }

    KWH &operator=(const KWH &other) {
        if(this == &other)
            return *this;
        seq = other.seq;
        pos = other.pos;
        fhash = other.fhash;
        rhash = other.rhash;
        return *this;
    }

    bool isCanonical() const {
        return fhash < rhash;
    }
};


template<typename htype>
class MinQueue {
    std::deque<KWH<htype>> q;
public:
    MinQueue() = default;

    void push(const KWH<htype> &kwh) {
        while(!q.empty() && q.back().hash() > kwh.hash()) {
            q.pop_back();
        }
        q.push_back(kwh);
    }

    void pop(size_t pos) {
        if(!q.empty() && q.front().pos < pos) {
            q.pop_front();
        }
    }

    bool empty() const {
        return q.empty();
    }

    KWH<htype> get() const {
        return q.front();
    }

    size_t size() const {
        return q.size();
    }
};

template<typename htype>
class MinimizerCalculator {
private:
    const Sequence seq;
    const size_t w;
    KWH<htype> kwh;
    size_t pos;
    MinQueue<htype> queue;
public:
    MinimizerCalculator(const Sequence& _seq, const RollingHash<htype> &_hasher, size_t _w) :
            seq(_seq), w(_w), kwh(_hasher, seq, 0), pos(-1) {
        VERIFY(w >= 2); //This code does not work for w = 1
        VERIFY(seq.size() >= _hasher.k + w - 1)
        queue.push(kwh);
        for(size_t i = 1; i < w; i++) {
            kwh = kwh.next();
            queue.push(kwh);
        }
    }

    KWH<htype> next() {
        pos += 1;
        queue.pop(pos);
        kwh = kwh.next();
        queue.push(kwh);
        return queue.get();
    }

    bool hasNext() const {
        return kwh.hasNext();
    }

    std::vector<htype> minimizerHashs() {
        std::vector<htype> res;
        res.push_back(queue.get().hash());
        while(hasNext()) {
            htype val = next().hash();
            if(val != res.back()) {
                res.push_back(val);
            }
        }
        return std::move(res);
    }

    std::vector<KWH<htype>> minimizers() {
        std::vector<KWH<htype>> res;
        res.push_back(next());
        while(hasNext()) {
            KWH<htype> val = next();
            if(val.pos != res.back().pos) {
                res.push_back(val);
            }
        }
        return std::move(res);
    }
};

