//
// Created by anton on 02.04.2020.
//

#pragma once
#include <string>
#include <algorithm>

static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

static bool startsWith(const std::string& str, const std::string& prefix)
{
    return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}

static inline void ltrim_inplace(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
}

static inline void rtrim_inplace(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

static inline std::string trim(std::string s) {
    ltrim_inplace(s);
    rtrim_inplace(s);
    return s;
}

static inline std::string & compress_inplace(std::string &s) {
    s.erase(std::unique(s.begin(), s.end()), s.end());
    return s;
}

inline std::vector<std::string> split(const std::string &s, const std::string &delimiter) {
    std::vector<std::string> res;
    size_t cur = 0;
    while(cur < s.size()) {
        size_t next = s.find(delimiter, cur);
        if (next == size_t(-1)) {
            next = s.size();
        }
        if (next > cur) {
            res.push_back(s.substr(cur, next - cur));
        }
        cur = next + delimiter.size();
    }
    return res;
}

inline std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> res;
    size_t cur = 0;
    std::string bad = " \n\t";
    while(cur < s.size()) {
        size_t next = cur;
        while(next < s.size() && bad.find(s[next]) == size_t(-1)) {
//            std::cout << s[cur] << " " << size_t(s[next]) << " " << size_t('\t') << std::endl;
            next += 1;
        }
        if (next > cur) {
            res.push_back(s.substr(cur, next - cur));
        }
        cur = next + 1;
    }
    return res;
}