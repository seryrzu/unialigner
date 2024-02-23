//
// Created by Andrey Bzikadze on 05/18/22.
//

#pragma once

#include <utility>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <experimental/iterator>
#include <experimental/filesystem>
#include "lcp_interval.hpp"
#include "common/dir_utils.hpp"

namespace tandem_aligner {

class MinInterval {
    int len{0};
    std::vector<int> fst_coords;
    std::vector<int> snd_coords;

 public:
    MinInterval() = default;
    explicit MinInterval(const int len) : len{len} {}

    void PushBackFst(const int start) { fst_coords.push_back(start); }
    void PushBackSnd(const int start) { snd_coords.push_back(start); }

    [[nodiscard]] int GetLen() const { return len; }
    [[nodiscard]] const std::vector<int> &GetFstCoords() const { return fst_coords; }
    [[nodiscard]] const std::vector<int> &GetSndCoords() const { return snd_coords; }
};

std::ostream &operator<<(std::ostream &os, const MinInterval &interval);

class MaxDisjointIntervalCollection {
    int fst_freq{1}, snd_freq{1};
    std::unordered_map<int, MinInterval> intervals;

 public:
    MaxDisjointIntervalCollection(const int fst_freq, const int snd_freq) :
        fst_freq{fst_freq}, snd_freq{snd_freq} {}

    MinInterval &operator[](const int &clas) { return intervals[clas]; }

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

std::ostream &operator<<(std::ostream &os, const MaxDisjointIntervalCollection &col);

using MaxDisjointIntervalCollections = std::vector<MaxDisjointIntervalCollection>;

std::ostream &operator<<(std::ostream &os, const MaxDisjointIntervalCollections &cols);

class MaxDisjointIntervalFinder {
    const int max_freq{1};
    const int min_freq{1};
    bool force_highfreq_search{false};
    bool exprt{true};
    std::experimental::filesystem::path outdir;

    struct MinMaxRarePrefixArray {
        // represents shortest and longest prefixes present fst_freq (snd_freq) times in the first (second) seq
        struct MinMaxRarePrefix {
            int clas{-1};
            int min_len{std::numeric_limits<int>::max()}; // length of shortest (fst_freq,snd_freq)-prefix
            int max_len{std::numeric_limits<int>::max()}; // length of longest (fst_freq,snd_freq)-prefix

            [[nodiscard]] bool IsInit() const {
                return clas!=-1 and min_len!=std::numeric_limits<int>::max();
            }
        };

        const int fst_freq{0}, snd_freq{0};
        std::vector<MinMaxRarePrefix> mrp; // shortest rare prefix

        MinMaxRarePrefixArray(const int fst_freq, const int snd_freq,
                      const int length) :
            fst_freq{fst_freq}, snd_freq{snd_freq}, mrp(length) {}

        MinMaxRarePrefix &operator[](int i) { return mrp[i]; }
        const MinMaxRarePrefix &operator[](int i) const { return mrp[i]; }

        [[nodiscard]] int Size() const { return mrp.size(); }
    };

    [[nodiscard]] std::vector<MinMaxRarePrefixArray> GetMRP(const LCP &lcp,
                                                    const int fst_len) const;

    struct ClassCoord {
        int lcp_class{-1};
        int coord{-1};
    };

    [[nodiscard]] MaxDisjointIntervalCollections
    PrefixesToIntervals(const std::vector<MinMaxRarePrefixArray> &mrp_vec,
                        const int fst_len) const;

    void OutputMRP(const std::vector<MinMaxRarePrefixArray> &srp_vec,
                   std::experimental::filesystem::path outpath) const;

 public:
    MaxDisjointIntervalFinder(int max_freq,
                      bool force_highfreq_search,
                      bool exprt,
                      std::experimental::filesystem::path outdir) :
                      max_freq{max_freq}, exprt{exprt},
                      force_highfreq_search{force_highfreq_search},
                      outdir{std::move(outdir)} {
        ensure_dir_existance(this->outdir);
    }

    [[nodiscard]] MaxDisjointIntervalCollections Find(const suffix_array::LCP<std::string> &lcp,
                                              const int fst_len) const;
};

} // namespace MinInterval
