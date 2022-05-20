//
// Created by Andrey Bzikadze on 05/19/22.
//

#include "cigar.hpp"
#include "sequences/verify.hpp"

using namespace minseq;

inline char minseq::cigar_mode2str(const CigarMode &fragment) {
    if (fragment==CigarMode::M) {
        return 'M';
    } else if (fragment==CigarMode::I) {
        return 'I';
    } else {
        VERIFY(fragment==CigarMode::D);
        return 'D';
    }
}

void Cigar::extend(const size_t length, const CigarMode mode) {
    if (length == 0) {
        return;
    }
    if (empty() or cigar_vec.back().mode!=mode) {
        cigar_vec.push_back({length, mode});
    } else {
        cigar_vec.back().length += length;
    }
}


[[nodiscard]] size_t Cigar::QueryLength() const {
    size_t length{0};
    for (const CigarFragment& fragment : cigar_vec) {
        if (fragment.mode != CigarMode::D) {
            length += fragment.length;
        }
    }
    return length;
}

[[nodiscard]] size_t Cigar::TargetLength() const {
    size_t length{0};
    for (const CigarFragment& fragment : cigar_vec) {
        if (fragment.mode != CigarMode::I) {
            length += fragment.length;
        }
    }
    return length;
}

[[nodiscard]] double Cigar::GetMutability() const {
    double nmatches{0};
    for (const CigarFragment& fragment : cigar_vec) {
        if (fragment.mode == CigarMode::M) {
            nmatches += fragment.length;
        }
    }
    std::cout << 1 - nmatches / TargetLength() << " " << 1 - nmatches / QueryLength() << "\n";
    return 1 - (nmatches / TargetLength() + nmatches / QueryLength()) / 2;
}

void Cigar::AssertValidity(const std::string &target,
                           const std::string &query) const {
    VERIFY(target.size() == TargetLength());
    VERIFY(query.size() == QueryLength());

    int i{0}, j{0};
    for (const CigarFragment& fragment : cigar_vec) {
        if (fragment.mode == CigarMode::M) {
            VERIFY(target.substr(i, fragment.length) ==
                   query.substr(j, fragment.length));
            i += fragment.length;
            j += fragment.length;
        } else if (fragment.mode == CigarMode::I) {
            j += fragment.length;
        } else if (fragment.mode == CigarMode::D) {
            i += fragment.length;
        }
    }
}

std::ostream &minseq::operator<<(std::ostream &os, const Cigar &cigar) {
    for (const CigarFragment &fragment : cigar.cigar_vec) {
        os << fragment.length << cigar_mode2str(fragment.mode);
    }
    return os;
}