//
// Created by Andrey Bzikadze on 05/06/22.
//


#pragma once

#include <vector>

namespace suffix_array {

template<class T, class ValueType = int64_t>
class SuffixArray {
 public:
    using value_type = ValueType;
 private:
    const T *p_str{nullptr};
    std::vector<value_type> arr;

 public:

    SuffixArray(const T &str, std::vector<value_type> arr) : p_str{&str},
                                                             arr{std::move(arr)} {}

    explicit SuffixArray(const T &str) : p_str{&str}, arr(str.size()) {
        int64_t n = str.size();
        value_type alphabet = 1 + *std::max_element(str.begin(), str.end());
        std::vector<int64_t> cnt(std::max(n, alphabet)), classes(n);
        for (int64_t i = 0; i < n; ++i) {
            cnt[str[i]]++;
        }
        for (value_type i = 1; i < alphabet; i++) {
            cnt[i] += cnt[i - 1];
        }
        for (int64_t i = 0; i < n; ++i) {
            arr[--cnt[str[i]]] = i;
        }

        int64_t cl{0};
        for (int64_t i = 1; i < n; ++i) {
            if (str[arr[i]]!=str[arr[i - 1]]) {
                cl++;
            }
            classes[arr[i]] = cl;
        }

        std::vector<value_type> arr2(n), classes2(n);
        for (int64_t k = 1; cl + 1!=n; k *= 2) {
            for (int64_t i = 0; i < n; ++i) {
                arr2[i] = (arr[i] + n - k)%n;
            }

            std::fill(cnt.begin(), cnt.begin() + cl + 1, 0);
            for (int64_t i = 0; i < n; ++i) {
                cnt[classes[i]]++;
            }
            for (int64_t i = 1; i <= cl; i++) {
                cnt[i] += cnt[i - 1];
            }

            for (int64_t i = n - 1; i >= 0; --i) {
                arr[--cnt[classes[arr2[i]]]] = arr2[i];
            }

            cl = 0;
            classes2[arr[0]] = 0;
            for (int64_t i = 1; i < n; ++i) {
                if (classes[arr[i]]!=classes[arr[i - 1]] or
                    classes[arr[i] + k%n]!=classes[arr[i - 1] + k%n]) {
                    cl++;
                }
                classes2[arr[i]] = cl;
            }
            classes.swap(classes2);
        }
    }

    typename std::vector<value_type>::iterator begin() { return arr.begin(); }
    typename std::vector<value_type>::iterator end() { return arr.end(); }
    [[nodiscard]] typename std::vector<value_type>::const_iterator begin() const {
        return arr.begin();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator end() const {
        return arr.end();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator cbegin() const {
        return arr.cbegin();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator cend() const {
        return arr.cend();
    }

    value_type operator[](const int64_t i) { return arr[i]; }
    value_type operator[](const int64_t i) const { return arr[i]; }

    const T *PStr() const { return p_str; }
    [[nodiscard]] size_t size() const { return arr.size(); }

    std::vector<value_type> Inverse() const {
        std::vector<value_type> inverse(size());
        for (int64_t i = 0; i < size(); ++i) {
            inverse[arr[i]] = i;
        }
        return inverse;
    }
};

template<class T>
class LCP {
 public:
    using value_type = typename SuffixArray<T>::value_type;

 private:
    const SuffixArray<T> *p_suf_arr{nullptr};
    std::vector<value_type> lcp;

 public:
    explicit LCP(const SuffixArray<T> &suf_arr) : p_suf_arr{&suf_arr} {
        if (p_suf_arr==nullptr or p_suf_arr->PStr()==nullptr) {
            return;
        }
        const T str = *(suf_arr.PStr());
        const int64_t n = suf_arr.size();
        std::vector<value_type> inv(n);
        for (int i = 0; i < n; ++i) {
            inv[suf_arr[i]] = i;
        }
        lcp.resize(n - 1);
        value_type cur_lcp{0};
        for (int64_t i = 0; i < n; ++i) {
            if (inv[i]==n - 1) {
                cur_lcp = 0;
                continue;
            }
            value_type j = suf_arr[inv[i] + 1];

            while (std::max(i + cur_lcp, j + cur_lcp) < n
                and str[i + cur_lcp]==str[j + cur_lcp]) {
                cur_lcp++;
            }
            lcp[inv[i]] = cur_lcp;
            if (cur_lcp) {
                cur_lcp--;
            }
        }
    }

    typename std::vector<value_type>::iterator begin() { return lcp.begin(); }
    typename std::vector<value_type>::iterator end() { return lcp.end(); }
    [[nodiscard]] typename std::vector<value_type>::const_iterator begin() const {
        return lcp.begin();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator end() const {
        return lcp.end();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator cbegin() const {
        return lcp.cbegin();
    }
    [[nodiscard]] typename std::vector<value_type>::const_iterator cend() const {
        return lcp.cend();
    }

    value_type operator[](const int64_t i) { return lcp[i]; }
    value_type operator[] (const int64_t i) const { return lcp[i]; }

    const T *PStr() const {
        return p_suf_arr==nullptr ? nullptr : p_suf_arr->p_str;
    }

    const SuffixArray<T> *PSufArr() const {
        return p_suf_arr;
    }

    [[nodiscard]] size_t size() const { return lcp.size(); }

};

} // namespace suffix_array