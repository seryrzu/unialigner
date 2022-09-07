//
// Created by Andrey Bzikadze on 05/17/22.
//

#pragma once

#include <min_queue/min_queue.hpp>
#include <suffix_array/suffix_array.hpp>

namespace tandem_aligner {

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
                const int64_t fst_len);

    void ExtendRight();
    void TrimLeft();

    [[nodiscard]] int64_t GetMin() const { return min_queue.GetMin(); }
    [[nodiscard]] int64_t GetFstFreq() const { return fst_freq; }
    [[nodiscard]] int64_t GetSndFreq() const { return snd_freq; }
};

} // End namespace tandem_aligner