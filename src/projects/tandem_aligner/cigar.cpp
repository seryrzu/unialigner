//
// Created by Andrey Bzikadze on 05/19/22.
//

#include "cigar.hpp"
#include "sequences/verify.hpp"

using namespace tandem_aligner;

inline char tandem_aligner::cigar_mode2str(const CigarMode &fragment) {
    if (fragment==CigarMode::M) {
        return 'M';
    } else if (fragment==CigarMode::I) {
        return 'I';
    } else if (fragment==CigarMode::D) {
        return 'D';
    } else {
        VERIFY(fragment==CigarMode::X);
        return 'X';
    }
}

void Cigar::extend(const int64_t length, const CigarMode mode) {
    if (length==0) {
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
    for (const CigarFragment &fragment : cigar_vec) {
        if (fragment.mode!=CigarMode::D) {
            length += fragment.length;
        }
    }
    return length;
}

[[nodiscard]] size_t Cigar::TargetLength() const {
    size_t length{0};
    for (const CigarFragment &fragment : cigar_vec) {
        if (fragment.mode!=CigarMode::I) {
            length += fragment.length;
        }
    }
    return length;
}

[[nodiscard]] double Cigar::GetMutability() const {
    double nmatches{0};
    for (const CigarFragment &fragment : cigar_vec) {
        if (fragment.mode==CigarMode::M) {
            nmatches += fragment.length;
        }
    }
    std::cout << 1 - nmatches/TargetLength() << " "
              << 1 - nmatches/QueryLength() << "\n";
    return 1 - (nmatches/TargetLength() + nmatches/QueryLength())/2;
}

void Cigar::AssertValidity(const std::string &target,
                           const std::string &query) const {
    VERIFY(target.size()==TargetLength());
    VERIFY(query.size()==QueryLength());

    int i{0}, j{0};
    for (const CigarFragment &fragment : cigar_vec) {
        if (fragment.mode==CigarMode::M or fragment.mode==CigarMode::X) {
            if (fragment.mode==CigarMode::M) {
                VERIFY(target.substr(i, fragment.length)==
                    query.substr(j, fragment.length));
            } else {
                for (int64_t k = 0; k < fragment.length; ++k) {
                    VERIFY(target[i + k]!=query[j + k]);
                }
            }
            i += fragment.length;
            j += fragment.length;
        } else if (fragment.mode==CigarMode::I) {
            j += fragment.length;
        } else if (fragment.mode==CigarMode::D) {
            i += fragment.length;
        }
    }
    VERIFY(i==target.size());
    VERIFY(j==query.size());
}

void Cigar::Summary(logging::Logger &logger) const {
    std::vector<int64_t> conf_cnt(4), base_cnt(4);
    for (const CigarFragment &fragment : cigar_vec) {
        conf_cnt[(int) fragment.mode] += fragment.length;
        base_cnt[(int) fragment.mode]++;
    }
    logger.info() << "Counts of configurations in the optimal alignment\n";
    for (int i = 0; i < 4; ++i) {
        logger.info() << "\t" << cigar_mode2str((CigarMode) i) << " "
                      << conf_cnt[i] << "\n";
    }
    logger.info() << "Counts of bases in the optimal alignment\n";
    for (int i = 0; i < 4; ++i) {
        logger.info() << "\t" << cigar_mode2str((CigarMode) i) << " "
                      << base_cnt[i] << "\n";
    }
}

std::ostream &tandem_aligner::operator<<(std::ostream &os, const Cigar &cigar) {
    for (const CigarFragment &fragment : cigar.cigar_vec) {
        os << fragment.length << cigar_mode2str(fragment.mode);
    }
    return os;
}