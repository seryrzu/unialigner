#pragma once

#include "sequence.hpp"
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include "verify.hpp"
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include <unordered_map>

using std::string;
using std::string;
namespace basic {
    inline string Reverse(const string &s){
        if (s[0] == '-')
            return s.substr(1);
        else
            return "-" + s;
    }
}

template<class T>
class Segment{
    T *contig_ptr;
public:
    size_t left;
    size_t right;
    Segment(T &contig_, size_t left_, size_t right_) : left(left_), right(right_), contig_ptr(&contig_){
        VERIFY(0 <= left and left <= right and right <= contig_ptr->size())
    }

    T &contig() const {
        return *contig_ptr;
    }

    size_t size() const {
        return right - left;
    }

    Sequence seq() const {
        return contig_ptr->seq.Subseq(left, right);
    }

    size_t dist(const Segment<T> &other) const {
        VERIFY(contig_ptr == other.contig_ptr);
        if (right <= other.left)
            return other.left - right;
        else if (other.right <= left)
            return left - other.right;
        else
            return 0;
    }

    Segment<T> RC() const {
        return Segment(contig_ptr->rc(), contig_ptr->size() - right, contig_ptr->size() - left);
    }

    bool inter(const Segment &other) const {
        return contig_ptr == other.contig_ptr and not (right <= other.left or left >= other.right);
    }

    int interSize(const Segment &other) const {
        if (not inter(other))
            return -1;
        else
            return int(std::min(right, other.right) - std::max(left, other.left));
    }

    bool operator<(const Segment<T> &other) const {
        return contig().id < other.contig().id ||
                    (contig().id == other.contig().id &&
                            (left < other.left || left == other.left && right < other.right));
    }
};

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Segment<T>& seg)
{
    os << seg.contig().id << "[" << seg.left << ":";
    if (seg.right > seg.contig().size() * 3 / 4)
        os << seg.contig().size() << "-" << (seg.contig().size() - seg.right);
    else
        os << seg.right;
    os << "]";
    return os;
}

template<class T>
class NamedSequence {
public:
    string id;
    Sequence seq;
//protected:
//    T * _rc;
public:
//    NamedSequence(const Sequence &_seq, string _id, T *_rc) : seq(_seq), id(std::move(_id)), _rc(_rc){
//    }

    NamedSequence(const Sequence &_seq, string _id) : seq(_seq), id(std::move(_id)){
//        _rc = new T(!seq, basic::Reverse(id), static_cast<T*>(this));
    }

    Segment<T> asSegment() const {
        return Segment<T>(*this, 0u, size());
    }

    Segment<T> segment(size_t left, size_t right) const {
        return Segment<T>(*(static_cast<const T*>(this)), left, right);
    }

    Segment<T> suffix(size_t pos) const {
        if (pos < 0)
            pos = size() + pos;
        if (pos < 0)
            pos = 0;
        if (pos > size())
            pos = size();
        return Segment<T>(*this, pos, size());
    }

    Segment<T> prefix(size_t len) const {
        len = min(len, size());
        return Segment<T>(*this, 0, len);
    }

    size_t size() const {
        return seq.size();
    }

//    NamedSequence &rc() const {
//        return *_rc;
//    }

    bool operator==(const NamedSequence &other) const {
        return id == other.id;
    }

    char operator[](size_t ind) const {
        return nucl(seq[ind]);
    }

    string str() const {
        return seq.str();
    }

    bool isNull() const {
        return seq.empty();
    }
};

class Contig: public NamedSequence<Contig> {
public:
    Contig(): NamedSequence(Sequence(), ""){
    }

    Contig(const Sequence &_seq, const string &_id): NamedSequence(_seq, _id) {
    }

//    Contig(const Sequence &_seq, const string &_id, Contig *_rc): NamedSequence(_seq, _id, _rc) {
//    }

    Contig(const string &_seq, const string &_id): NamedSequence(Sequence(_seq), _id) {
    }

    Contig RC() {
        if(id[0] == '-')
            return Contig(!seq, id.substr(1, id.size() - 1));
        else
            return Contig(!seq, "-" + id);
    }

//    Contig(const string &_seq, const string &_id, Contig *_rc): NamedSequence(Sequence(_seq), _id, _rc) {
//    }
};

class StringContig {
public:
    std::string id;
    std::string seq;
    static bool needs_compressing;

    StringContig() : id(""), seq("") {
    }

    StringContig(std::string && _seq, std::string &&_id) : id(_id), seq(_seq) {
    }

    StringContig(StringContig && other) = default;

    StringContig& operator=(StringContig && other) = default;

    void compress() {
        seq.erase(std::unique(seq.begin(), seq.end()), seq.end());
    }

    Contig makeContig() {
        if(needs_compressing)
            compress();
        return Contig(Sequence(seq), id);
    }

//    Contig makeCompressedContig() {
//        compress();
//        return makeContig();
//    }

    Sequence makeSequence() {
        if(needs_compressing)
            compress();
        return Sequence(seq);
    }

//    Sequence makeCompressedSequence() {
//        compress();
//        return makeSequence();
//    }

    bool isNull() const {
        return id.empty() && seq.empty();
    }

    size_t size() const {
        return seq.size();
    }
};



template <class T>
class SequenceCollection {
private:
    std::unordered_map<std::string, T *> items;
public:
    SequenceCollection() {
    }

    SequenceCollection(const std::vector<T*> & sequences) {
        for(T *item: sequences){
            items[item->id] = item;
        }
    }
};