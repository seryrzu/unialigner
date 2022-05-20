//
// Created by Andrey Bzikadze on 05/19/22.
//

#pragma once

#include <vector>
#include <fstream>

namespace minseq {

enum class CigarMode { M, I, D };

inline char cigar_mode2str(const CigarMode &fragment);

struct CigarFragment {
    size_t length{0};
    CigarMode mode{};
};

class Cigar {
    std::vector<CigarFragment> cigar_vec;

 public:
    Cigar() = default;

    [[nodiscard]] bool empty() const { return cigar_vec.empty(); }
    void extend(size_t length, CigarMode mode);

    [[nodiscard]] size_t QueryLength() const;
    [[nodiscard]] size_t TargetLength() const;

    [[nodiscard]] double GetMutability() const;

    void AssertValidity(const std::string &target,
                        const std::string &query) const;

    friend std::ostream &operator<<(std::ostream &os, const Cigar &cigar);

    decltype(cigar_vec)::iterator begin() { return cigar_vec.begin(); }
    decltype(cigar_vec)::iterator end() { return cigar_vec.end(); }
    [[nodiscard]] decltype(cigar_vec)::const_iterator begin() const {
        return cigar_vec.begin();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator end() const {
        return cigar_vec.end();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator cbegin() const {
        return cigar_vec.cbegin();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator cend() const {
        return cigar_vec.cend();
    }
};

std::ostream &operator<<(std::ostream &os, const Cigar &cigar);

} // namespace minseq