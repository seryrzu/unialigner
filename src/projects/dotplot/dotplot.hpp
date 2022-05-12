#include <utility>

//
// Created by Andrey Bzikadze on 05/05/22.
//

#pragma once

#include <common/logging.hpp>
#include <experimental/iterator>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <suffix_array/suffix_array.hpp>
#include <min_queue/min_queue.hpp>

namespace dotplot {

using LCP = suffix_array::LCP<std::string>;
using SuffixArray = suffix_array::SuffixArray<std::string>;

class LCPInterval {
    const LCP &lcp;
    const SuffixArray &suf_arr;
    const int64_t fst_len{0};
    int64_t left{0}, right{0}; // [): coordinates in LCP array
    int64_t fst_freq{0}, snd_freq{0};
    min_queue::MinQueue<LCP> min_queue;

 public:
    LCPInterval(const LCP &lcp,
                const SuffixArray &suf_arr,
                const int64_t fst_len) :
        lcp{lcp}, suf_arr{suf_arr},
        fst_len{fst_len}, min_queue{lcp} {
        if (suf_arr.size()) {
            if (suf_arr[0] <= fst_len) {
                fst_freq++;
            } else {
                snd_freq++;
            }
        }
    }

    void ExtendRight() {
        VERIFY(right < lcp.size());
        if (suf_arr[right + 1] <= fst_len) {
            fst_freq++;
        } else {
            snd_freq++;
        }
        min_queue.PushBack();
        right++;
    }

    void TrimLeft() {
        VERIFY(left!=right);
        if (suf_arr[left] <= fst_len) {
            fst_freq--;
        } else {
            snd_freq--;
        }
        min_queue.PopFront();
        left++;
    }

    [[nodiscard]] int64_t GetMin() const { return min_queue.GetMin(); }
    [[nodiscard]] int64_t GetFstFreq() const { return fst_freq; }
    [[nodiscard]] int64_t GetSndFreq() const { return snd_freq; }
};

class Interval {
    int64_t start{0};
    int64_t end{0};

 public:
    Interval(const int64_t start, const int64_t end) : start{start}, end{end} {}
};

class ShortInterval {
    int64_t len{0};
    std::vector<int64_t> fst_coords;
    std::vector<int64_t> snd_coords;

 public:
    ShortInterval() = default;
    explicit ShortInterval(const int64_t len) : len{len} {}

    void PushBackFst(const int64_t start) { fst_coords.push_back(start); }
    void PushBackSnd(const int64_t start) { snd_coords.push_back(start); }

    [[nodiscard]] int64_t GetLen() const { return len; }
    [[nodiscard]] const std::vector<int64_t> &GetFstCoords() const { return fst_coords; }
    [[nodiscard]] const std::vector<int64_t> &GetSndCoords() const { return snd_coords; }
};

std::ostream &operator<<(std::ostream &os, const ShortInterval &interval) {
    os << interval.GetLen() << "\t";
    std::copy(interval.GetFstCoords().begin(),
              interval.GetFstCoords().end(),
              std::experimental::make_ostream_joiner(os, ","));
    os << "\t";
    std::copy(interval.GetSndCoords().begin(),
              interval.GetSndCoords().end(),
              std::experimental::make_ostream_joiner(os, ","));
    os << "\n";
    return os;
}

class ShortIntervalCollection {
    int fst_freq{1}, snd_freq{1};
    std::unordered_map<int64_t, ShortInterval> intervals;

 public:
    ShortIntervalCollection(const int fst_freq, const int snd_freq) :
        fst_freq{fst_freq}, snd_freq{snd_freq} {}

    ShortInterval &operator[](const int64_t &clas) { return intervals[clas]; }

    decltype(intervals)::iterator begin() { return intervals.begin(); }
    decltype(intervals)::iterator end() { return intervals.end(); }
    [[nodiscard]] decltype(intervals)::const_iterator begin() const {
        return intervals.begin();
    }
    [[nodiscard]] decltype(intervals)::const_iterator end() const {
        return intervals.end();
    }
    [[nodiscard]] decltype(intervals)::const_iterator cbegin() const {
        return intervals.cbegin();
    }
    [[nodiscard]] decltype(intervals)::const_iterator cend() const {
        return intervals.cend();
    }

    template<class... Args>
    std::pair<decltype(intervals)::iterator, bool> Emplace(Args &&... args) {
        return intervals.emplace(args...);
    }

    int GetFstFreq() const { return fst_freq; }
    int GetSndFreq() const { return snd_freq; }
};

std::ostream &operator<<(std::ostream &os, const ShortIntervalCollection &col) {
    for (const auto &[clas, interval] : col) {
        os << col.GetFstFreq() << "\t" << col.GetSndFreq() << "\t" << interval;
    }
    return os;
}

using ShortIntervalCollections = std::vector<ShortIntervalCollection>;

std::ostream &operator<<(std::ostream &os,
                         const ShortIntervalCollections &cols) {
    os << "FstFreq\tSndFreq\tLength\tFstStarts\tSndStarts\n";
    for (const ShortIntervalCollection &col : cols) {
        os << col;
    }
    return os;
}

class RareSharedSegmentFinder {
    const int max_freq{1};
    const int min_freq{1};

    struct ClassCoord {
        int64_t clas{-1};
        int64_t coord{-1};
    };

    struct ShortestRarePrefix {
        const int64_t fst_freq{0}, snd_freq{0};
        std::vector<ClassCoord> srp; // shortest rare prefix

        ShortestRarePrefix(const int64_t fst_freq, const int64_t snd_freq,
                           const int64_t length) :
            fst_freq{fst_freq}, snd_freq{snd_freq}, srp(length) {}
    };

    [[nodiscard]] std::vector<ShortestRarePrefix> GetSRP(const LCP &lcp,
                                                         const int64_t fst_len) const {
        VERIFY(lcp.PSufArr()!=nullptr);
        const SuffixArray suf_arr = *lcp.PSufArr();
        std::vector<int64_t> inv = suf_arr.Inverse();
        std::vector<ShortestRarePrefix> srp_vec;
        for (int w = 2; w <= 2*max_freq; w++) {
            for (int i = 1; i <= w - 1; ++i) {
                srp_vec.emplace_back(i, w - i, suf_arr.size());
            }
            LCPInterval lcp_interval(lcp, suf_arr, fst_len);
            int64_t prev{0}, cur{0}, next{0};
            int64_t fst_freq{0}, snd_freq{0};
            int64_t i{0}, j{0};
            for (; j < w - 1; ++j) {
                lcp_interval.ExtendRight();
            }
            while (j < suf_arr.size()) {
                cur = lcp_interval.GetMin();
                fst_freq = lcp_interval.GetFstFreq();
                snd_freq = lcp_interval.GetSndFreq();
                if (j < suf_arr.size() - 1) {
                    lcp_interval.ExtendRight();
                    next = lcp_interval.GetMin();
                } else {
                    next = 0;
                }
                if (cur > std::max(prev, next)
                    and std::min(fst_freq, snd_freq) >= min_freq and
                    std::max(fst_freq, snd_freq) <= max_freq) {
                    //std::cout << i << " " << suf_arr[i] << " " << j << " "
                    //          << suf_arr[j] << " "
                    //          << fst_freq << " " << snd_freq << " "
                    //          << prev << " " << cur << " " << next << "\n";
                    // std::cout << suf_arr[i] << ", " << suf_arr[j] << ", " << cur
                    //           << " " << fst_freq << " " << snd_freq << "\n";
                    int64_t srp_index = (w - 2)*(w - 1)/2 + fst_freq - 1;
                    VERIFY(srp_vec[srp_index].fst_freq==fst_freq);
                    VERIFY(srp_vec[srp_index].snd_freq==snd_freq);
                    for (int64_t k = i; k <= j; ++k) {
                        // VERIFY(srp_vec[srp_index].srp[suf_arr[k]]==0);
                        srp_vec[srp_index].srp[suf_arr[k]] =
                            {i, std::max(prev, next) + 1};
                    }
                }
                prev = next;
                lcp_interval.TrimLeft();
                i++, j++;
            }
        }
        return srp_vec;
    }

    [[nodiscard]] ShortIntervalCollections
    ShortenPrefixes(const std::vector<ShortestRarePrefix> &srp_vec,
                    const int64_t fst_len) const {
        ShortIntervalCollections cols;
        for (const ShortestRarePrefix &srp : srp_vec) {
            auto &col = cols.emplace_back(srp.fst_freq, srp.snd_freq);
            std::vector<ClassCoord> en2st(srp.srp.size());
            for (int64_t st = 0; st < srp.srp.size(); ++st) {
                if (srp.srp[st].coord==-1)
                    continue;
                int64_t en = st + srp.srp[st].coord;
                if (st > en2st[en].coord) {
                    en2st[en] = {srp.srp[st].clas, st};
                }
            }
            for (int64_t en = 0; en < srp.srp.size(); ++en) {
                int64_t clas{en2st[en].clas}, st{en2st[en].coord};
                if (clas!=-1) {
                    col.Emplace(clas, en - st);
                    if (st < fst_len) {
                        col[clas].PushBackFst(st);
                    } else {
                        col[clas].PushBackSnd(st - fst_len - 1);
                    }
                }
            }
        }
        return cols;
    }

 public:
    explicit RareSharedSegmentFinder(int max_freq) : max_freq{max_freq} {}

    ShortIntervalCollections Find(const suffix_array::LCP<std::string> &lcp,
                                  const int64_t fst_len) const {
        const std::vector<ShortestRarePrefix> srp_vec = GetSRP(lcp, fst_len);
        return ShortenPrefixes(srp_vec, fst_len);
    }

};

class DotPlotConstructor {
    logging::Logger &logger;
    const std::experimental::filesystem::path output_dir;

    [[nodiscard]] std::string ReadContig(const std::experimental::filesystem::path &path) const {
        io::SeqReader reader(path);
        std::vector<Contig> vec{reader.readAllContigs()};
        VERIFY(vec.size()==1);
        Contig contig{std::move(vec[0])};
        logger.info() << "Length " << contig.seq.size() << ", name "
                      << contig.id
                      << " from " << path << std::endl;
        return contig.str();
    }

    [[nodiscard]] std::string ConcatContigs(const std::string &first,
                                            const std::string &second) const {
        std::stringstream concat_stream;
        concat_stream << first << '$' << second << '#';
        return concat_stream.str();
    }

 public:
    DotPlotConstructor(logging::Logger &logger,
                       std::experimental::filesystem::path output_dir) :
        logger{logger}, output_dir{std::move(output_dir)} {}

    void Construct(const std::experimental::filesystem::path &first_path,
                   const std::experimental::filesystem::path &second_path) const {
        io::SeqReader first_reader(first_path);
        std::vector<Contig> first_vec{first_reader.readAllContigs()};
        VERIFY(first_vec.size()==1);
        const std::string first = ReadContig(first_path);
        const std::string second = ReadContig(second_path);
        const std::string concat = ConcatContigs(first, second);

        logger.info() << "Building suffix array...\n";
        const suffix_array::SuffixArray<std::string> suf_arr(concat);
        logger.info() << "Building LCP array...\n";
        const suffix_array::LCP<std::string> lcp(suf_arr);

        RareSharedSegmentFinder segment_finder(10);
        logger.info() << "Computing rare segments...\n";
        const ShortIntervalCollections
            int_col = segment_finder.Find(lcp, first.size());
        std::ofstream os(output_dir / "shortest_matches.tsv");
        os << int_col;
    }
};

} // namespace dotplot

