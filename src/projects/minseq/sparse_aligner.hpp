//
// Created by Andrey Bzikadze on 05/18/22.
//

#pragma once

#include "cigar.hpp"

namespace minseq {

class MinFreqInterval {
    int len{0};
    int fst_freq{0}, snd_freq{0};
    int fst_coord{0}, snd_coord{0};

 public:
    MinFreqInterval(int len, int fst_freq, int snd_freq,
                    int fst_coord, int snd_coord) :
        len{len}, fst_freq{fst_freq}, snd_freq{snd_freq},
        fst_coord{fst_coord}, snd_coord{snd_coord} {}

    friend bool operator<(const MinFreqInterval &lhs,
                          const MinFreqInterval &rhs);
    friend class SparseAligner;

    [[nodiscard]] double Score() const {
        return double(len)/(fst_freq*snd_freq);
    }
};

bool operator<(const MinFreqInterval &lhs, const MinFreqInterval &rhs) {
    return lhs.fst_coord < rhs.fst_coord or
        lhs.fst_coord==rhs.fst_coord and lhs.snd_coord==rhs.snd_coord;
}

class SparseAligner {
    logging::Logger &logger;
    static std::vector<MinFreqInterval> Cols2Vec(const MinIntervalCollections &cols) {
        std::vector<MinFreqInterval> vec;
        for (const MinIntervalCollection &col : cols)
            for (const auto &[_, interval] : col)
                for (const int fst_coord : interval.GetFstCoords())
                    for (const int snd_coord : interval.GetSndCoords())
                        vec.emplace_back(interval.GetLen(),
                                         col.GetFstFreq(), col.GetSndFreq(),
                                         fst_coord, snd_coord);
        return vec;
    }

    std::vector<MinFreqInterval>
    GetAlignmentVec(const std::vector<MinFreqInterval> &vec) {
        std::vector<double> scores{vec.front().Score()};
        std::vector<int> backtracks{-1};

        double max_score{0};
        int argmax_score{0};
        for (int i = 1; i < vec.size(); ++i) {
            double &score = scores.emplace_back(vec[i].Score());
            int &backtrack = backtracks.emplace_back(-1);
            for (int j = i - 1; j >= 0; --j) {
                if (vec[i].fst_coord >= vec[j].fst_coord + vec[j].len and
                    vec[i].snd_coord >= vec[j].snd_coord + vec[j].len) {
                    double new_score = scores[j] + vec[i].Score();
                    if (score < new_score) {
                        score = new_score;
                        backtrack = j;
                    }
                }
            }
            if (score > max_score) {
                max_score = score;
                argmax_score = i;
            }
        }
        logger.info() << max_score << " " << argmax_score << "\n";

        std::vector<MinFreqInterval> alignment_vec;
        while (argmax_score!=-1) {
            alignment_vec.emplace_back(vec[argmax_score]);
            argmax_score = backtracks[argmax_score];
        }
        std::reverse(alignment_vec.begin(), alignment_vec.end());
        logger.info() << alignment_vec.size() << "\n";

        return alignment_vec;
    }

    Cigar AlignmentVec2Cigar(std::vector<MinFreqInterval> vec,
                             const std::string &fst, const std::string &snd) {
        vec.emplace_back(0, 0, 0, fst.size(), snd.size());
        Cigar cigar;
        int left_border_fst{0}, left_border_snd{0},
            right_border_fst{0}, right_border_snd{0};
        for (int i = 1; i < vec.size() - 1; ++i) {
            const MinFreqInterval &cur = vec[i], &next = vec[i + 1];
            right_border_fst = next.fst_coord;
            right_border_snd = next.snd_coord;
            int left_fst = cur.fst_coord, left_snd = cur.snd_coord;
            while (left_fst >= left_border_fst
                and left_snd >= left_border_snd
                and fst[left_fst]==snd[left_snd]) {
                left_fst--, left_snd--;
            }
            left_fst++, left_snd++;
            int right_fst = cur.fst_coord + cur.len,
                right_snd = cur.snd_coord + cur.len;
            while (right_fst < right_border_fst
                and right_snd < right_border_snd
                and fst[right_fst]==snd[right_snd]) {
                right_fst++, right_snd++;
            }
            cigar.extend(left_fst - left_border_fst, CigarMode::D);
            cigar.extend(left_snd - left_border_snd, CigarMode::I);
            cigar.extend(right_fst - left_fst, CigarMode::M);

            left_border_fst = right_fst;
            left_border_snd = right_snd;
        }
        cigar.extend(fst.size() - left_border_fst, CigarMode::D);
        cigar.extend(snd.size() - left_border_snd, CigarMode::I);

        cigar.AssertValidity(fst, snd);
        return cigar;
    }

 public:
    explicit SparseAligner(logging::Logger &logger) : logger{logger} {}

    Cigar Align(const MinIntervalCollections &cols,
                const std::string &fst, const std::string &snd) {
        std::vector<MinFreqInterval> vec = Cols2Vec(cols);
        logger.info() << vec.size() << "\n";
        std::sort(vec.begin(), vec.end());
        logger.info() << vec.size() << "\n";

        const std::vector<MinFreqInterval> alignment_vec = GetAlignmentVec(vec);

        return AlignmentVec2Cigar(std::move(alignment_vec), fst, snd);
    }
};

} // namespace minseq