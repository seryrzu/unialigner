//
// Created by Andrey Bzikadze on 05/18/22.
//

#include "min_interval.hpp"

using namespace minseq;

std::ostream &minseq::operator<<(std::ostream &os,
                                 const MinInterval &interval) {
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

std::ostream &minseq::operator<<(std::ostream &os,
                                 const MinIntervalCollection &col) {
    for (const auto &[clas, interval] : col) {
        os << col.GetFstFreq() << "\t" << col.GetSndFreq() << "\t" << interval;
    }
    return os;
}

std::ostream &minseq::operator<<(std::ostream &os,
                                 const MinIntervalCollections &cols) {
    os << "FstFreq\tSndFreq\tLength\tFstStarts\t"
          "SndStarts\n";
    for (const MinIntervalCollection &col : cols) {
        os << col;
    }
    return os;
}

[[nodiscard]] std::vector<MinIntervalFinder::MinRarePrefix>
MinIntervalFinder::GetMRP(const LCP &lcp, const int fst_len) const {
    VERIFY(lcp.PSufArr()!=nullptr);
    const SuffixArray suf_arr = *lcp.PSufArr();
    std::vector<MinRarePrefix> mrp_vec;
    for (int w = 2; w <= 2*max_freq; w++) {
        for (int i = 1; i <= w - 1; ++i) {
            mrp_vec.emplace_back(i, w - i, suf_arr.size());
        }
        LCPInterval lcp_interval(lcp, suf_arr, fst_len);
        int prev{0}, cur{0}, next{0};
        int fst_freq{0}, snd_freq{0};
        int i{0}, j{0};
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
                int srp_index = (w - 2)*(w - 1)/2 + fst_freq - 1;
                VERIFY(mrp_vec[srp_index].fst_freq==fst_freq);
                VERIFY(mrp_vec[srp_index].snd_freq==snd_freq);
                for (int k = i; k <= j; ++k) {
                    // VERIFY(srp_vec[srp_index].srp[suf_arr[k]]==0);
                    // srp_vec[srp_index].srp[suf_arr[k]] =
                    //     {i, std::max(prev, next) + 1};
                    mrp_vec[srp_index][suf_arr[k]] =
                        {i, std::max(prev, next) + 1};
                }
            }
            prev = next;
            lcp_interval.TrimLeft();
            i++, j++;
        }
    }
    return mrp_vec;
}

[[nodiscard]] MinIntervalCollections
MinIntervalFinder::PrefixesToIntervals(const std::vector<MinRarePrefix> &mrp_vec,
                                       const int fst_len) const {
    MinIntervalCollections cols;
    for (const MinRarePrefix &mrp : mrp_vec) {
        std::vector<ClassCoord> en2st(mrp.Size() + 1);
        int cur_en{-1};
        for (int st = 0; st < mrp.Size(); ++st) {
            if (not mrp[st].IsInit()) {
                continue;
            }

            int en = st + mrp[st].len;
            int cur_st = en2st[en].coord;
            if (st > cur_st or cur_st==-1)
                en2st[en] = {mrp[st].clas, st};
        }

        auto &col = cols.emplace_back(mrp.fst_freq, mrp.snd_freq);
        for (int en = 0; en < mrp.Size(); ++en) {
            int clas{en2st[en].lcp_class}, st{en2st[en].coord};
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

void MinIntervalFinder::OutputMRP(const std::vector<MinRarePrefix> &srp_vec,
                                  std::experimental::filesystem::path outpath) const {
    ensure_dir_existance(outpath);
    for (const MinRarePrefix &srp : srp_vec) {
        std::stringstream fn_ss;
        fn_ss << "FstFreq" << srp.fst_freq << "_SndFreq" << srp.snd_freq
              << ".tsv";
        std::ofstream os(outpath/fn_ss.str());
        os << "Pos\tLength\n";
        for (auto it = srp.mrp.begin(); it!=srp.mrp.end(); ++it) {
            os << it - srp.mrp.begin() << "\t" << it->len << "\n";
        }
    }
}

[[nodiscard]] MinIntervalCollections
MinIntervalFinder::Find(const suffix_array::LCP<std::string> &lcp,
                        const int fst_len) const {
    const std::vector<MinRarePrefix> srp_vec = GetMRP(lcp, fst_len);
    if (exprt)
        OutputMRP(srp_vec, outdir/"min_rare_prefixes");
    return PrefixesToIntervals(srp_vec, fst_len);
}
