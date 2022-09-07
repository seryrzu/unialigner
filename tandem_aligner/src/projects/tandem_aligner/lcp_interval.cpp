//
// Created by Andrey Bzikadze on 05/17/22.
//

#include "lcp_interval.hpp"

using namespace tandem_aligner;

LCPInterval::LCPInterval(const LCP &lcp,
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

void LCPInterval::ExtendRight() {
    VERIFY(right < lcp.size());
    if (suf_arr[right + 1] <= fst_len) {
        fst_freq++;
    } else {
        snd_freq++;
    }
    min_queue.PushBack();
    right++;
}

void LCPInterval::TrimLeft() {
    VERIFY(left!=right);
    if (suf_arr[left] <= fst_len) {
        fst_freq--;
    } else {
        snd_freq--;
    }
    min_queue.PopFront();
    left++;
}
