//
// Created by Andrey Bzikadze on 03/09/22.
//

#pragma once

#include <cmath>
#include <optional>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "rolling_hash/rolling_hash.hpp"

namespace cubic_naive {

using RollingHash = RollingHash<unsigned __int128>;
using KWH = KWH<unsigned __int128>;
using SeqIndex = std::unordered_map<unsigned __int128, uint64_t>;

std::vector<RollingHash> GetRH(uint64_t len, int base = 239) {
    std::vector<RollingHash> rhs;
    for (uint64_t k = 1; k <= len; ++k) {
        rhs.emplace_back(k, base);
    }
    return rhs;
}

std::vector<SeqIndex> IndexSubSeqs(const Contig &contig,
                                   const std::vector<RollingHash> &hashers) {
    std::vector<SeqIndex> indexes;
    for (uint64_t k = 1; k <= hashers.size(); ++k) {
        SeqIndex &index = indexes.emplace_back();
        std::cout << k << "\n";
        const RollingHash &hasher = hashers.at(k - 1);
        KWH kwh(hasher, contig.seq, 0);
        while (true) {
            index[kwh.Fhash()] += 1;
            if (!kwh.hasNext()) {
                break;
            }
            kwh = kwh.next();
        }
    }
    return indexes;
}

void FilterIndexes(std::vector<SeqIndex> &indexes_first,
                   std::vector<SeqIndex> &indexes_second) {
    VERIFY(indexes_first.size()==indexes_second.size());
    for (uint64_t k = 1; k <= indexes_first.size(); ++k) {
        SeqIndex &index_first = indexes_first[k - 1];
        SeqIndex &index_second = indexes_second[k - 1];

        for (auto it = index_first.begin(); it!=index_first.end();) {
            if (not index_second.contains(it->first))
                it = index_first.erase(it);
            else
                ++it;
        }

        for (auto it = index_second.begin(); it!=index_second.end();) {
            if (not index_first.contains(it->first))
                it = index_second.erase(it);
            else
                ++it;
        }
    }
}

std::vector<std::vector<std::pair<uint64_t, uint64_t>>>
    Align(const Contig &first, const Contig &second,
          const std::vector<SeqIndex> &indexes_first,
          const std::vector<SeqIndex> &indexes_second,
          const double pindel = 0.3,
          const uint64_t base = 239) {
    std::vector<std::vector<double>> scores;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> backtrack;
    for (uint64_t i = 0; i <= first.size(); ++i) {
        scores.emplace_back(second.size() + 1);
        backtrack.emplace_back(second.size() + 1);
    }

    const double lpindel = std::log(pindel);
    const double lpm = std::log(1 - 2 * pindel);
    for (uint64_t i = 1; i <= first.size(); ++i) {
        scores[i][0] = scores[i - 1][0] + lpindel;
        backtrack[i][0] = std::make_pair(i - 1, 0);
    }
    for (uint64_t j = 1; j <= second.size(); ++j) {
        scores[0][j] = scores[0][j - 1] + lpindel;
        backtrack[0][j] = std::make_pair(0, j - 1);
    }

    const std::vector<RollingHash> hashers = GetRH(first.size(), base);
    for (uint64_t i = 1; i <= first.size(); ++i) {
        std::cout << i << "\n";
        for (uint64_t j = 1; j <= second.size(); ++j) {
            scores[i][j] = scores[i - 1][j] + lpindel;
            backtrack[i][j] = std::make_pair(i - 1, j);
            const double hscore = scores[i][j - 1] + lpindel;
            if (hscore > scores[i][j]) {
                scores[i][j] = hscore;
                backtrack[i][j] = std::make_pair(i, j - 1);
            }

            const double mism_score = scores[i - 1][j - 1] + lpm;
            if (first[i - 1]!=second[j - 1]) {
                if (mism_score > scores[i][j]) {
                    scores[i][j] = mism_score;
                    backtrack[i][j] = std::make_pair(i - 1, j - 1);
                }
                continue;
            }
            uint64_t k = 1;

            KWH kwh(hashers[0], first.seq, i - k);
            while (i >= k and j >= k and
                    first[i - k]==second[j - k] and
                    k <= indexes_first.size()) {
                // std::cout << i << " " << j << " " << k << " " << kwh.getSeq().str() << std::endl;

                const double diag_score =
                    scores[i - k][j - k]
                    - std::log(indexes_first[k - 1].at(kwh.Fhash()))
                    - std::log(indexes_second[k - 1].at(kwh.Fhash()))
                    + lpm;
                if (diag_score > scores[i][j]) {
                    scores[i][j] = diag_score;
                    backtrack[i][j] = std::make_pair(i - k, j - k);
                }
                ++k;
                if (i >= k) {
                    kwh = kwh.extendLeft(hashers[k-1]);
                }
            }
        }
    }
    return backtrack;
}

enum class CigarMode { M, I, D, S };

inline char cigar_mode2str(const CigarMode& fragment) {
    if (fragment == CigarMode::M) {
        return 'M';
    } else if (fragment == CigarMode::I) {
        return 'I';
    } else if (fragment == CigarMode::D) {
        return 'D';
    } else {
        VERIFY(fragment == CigarMode::S);
        return 'S';
    }
}

struct CigarFragment {
    size_t length{0};
    CigarMode mode{};
};

class Cigar {
    std::vector<CigarFragment> cigar_vec;

 public:
    Cigar() = default;
    Cigar(Cigar &) = default;
    Cigar(Cigar &&) = default;
    Cigar &operator=(const Cigar &) = default;
    Cigar &operator=(Cigar &&) = default;
    ~Cigar() = default;

    [[nodiscard]] bool empty() const { return cigar_vec.empty(); }

    void extend(const size_t length, const CigarMode mode) {
        if (empty() or cigar_vec.back().mode != mode) {
            cigar_vec.push_back({length, mode});
        } else {
            cigar_vec.back().length += length;
        }
    }

    void reverse() {
        std::reverse(cigar_vec.begin(), cigar_vec.end());
    }

    friend std::ostream& operator<<(std::ostream &os, const Cigar &cigar);
};

std::ostream& operator<<(std::ostream &os, const Cigar &cigar);

Cigar Backtrack(std::vector<std::vector<std::pair<uint64_t, uint64_t>>> &backtrack,
                const Contig &first, const Contig &second) {
    uint64_t i = first.size() - 1;
    uint64_t j = second.size() - 1;
    Cigar cigar;
    while (i != 0 or j != 0) {
        auto [new_i, new_j] = backtrack[i][j];
        if (new_i == i) {
            cigar.extend(1, CigarMode::I);
        } else if (new_j == j) {
            cigar.extend(1, CigarMode::D);
        } else if (i - new_i > 1) {
            cigar.extend(i - new_i, CigarMode::M);
        } else if (first[new_i] == second[new_j]) {
            cigar.extend(1, CigarMode::M);
        } else {
            cigar.extend(1, CigarMode::S);
        }
        i = new_i;
        j = new_j;
    }
    cigar.reverse();
    return cigar;
}

void CubicNaive(const std::experimental::filesystem::path &first_path,
                const std::experimental::filesystem::path &second_path,
                const std::experimental::filesystem::path &outdir,
                logging::Logger &logger) {
    io::SeqReader first_reader(first_path);
    std::vector<Contig> first_vec{first_reader.readAllContigs()};
    VERIFY(first_vec.size()==1);
    Contig first{std::move(first_vec[0])};
    logger.info() << "First length " << first.seq.size() << ", name "
                  << first.id
                  << " from " << first_path << std::endl;

    io::SeqReader second_reader(second_path);
    std::vector<Contig> second_vec{second_reader.readAllContigs()};
    VERIFY(second_vec.size()==1);
    Contig second{std::move(second_vec[0])};
    logger.info() << "Second length " << second.seq.size() << ", name "
                  << second.id
                  << " from " << second_path << std::endl;

    std::vector<RollingHash> hashers =
        GetRH(std::min(first.size(), second.size()));
    logger.info() << "Finish forming hashers\n";

    logger.info() << "Forming index for first sequence\n";
    std::vector<SeqIndex> indexes_first = IndexSubSeqs(first, hashers);
    logger.info() << "Finished forming index for first sequence\n";

    logger.info() << "Forming index for second sequence\n";
    std::vector<SeqIndex> indexes_second = IndexSubSeqs(second, hashers);
    logger.info() << "Finished forming index for second sequence\n";

    FilterIndexes(indexes_first, indexes_second);

    logger.info() << "Running alignment\n";
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> backtrack =
        Align(first, second, indexes_first, indexes_second);
    logger.info() << "Finished running alignment\n";

    logger.info() << "Running backtracking\n";
    Cigar cigar = Backtrack(backtrack, first, second);
    logger.info() << "Finished running backtracking\n";

    std::ofstream cigar_outfn(outdir / "cigar.tsv");
    cigar_outfn << cigar;
};

} // namespace cubic_naive
