//
// Created by Andrey Bzikadze on 05/07/22.
//


#include <iostream>
#include <string>
#include <vector>
#include "suffix_array/suffix_array.hpp"


int main() {
    int n = 0;
    std::string str;
    std::cin >> n;
    std::cin >> str;
    std::vector<int64_t> arr(n);
    for (int i = 0; i < n; ++i)
        std::cin >> arr[i], arr[i]--;

    suffix_array::SuffixArray<std::string> sufarr(str, std::move(arr));
    suffix_array::LCP<std::string> lcp(sufarr);
    for (const auto x : lcp) {
        std::cout << x << " ";
    }
}