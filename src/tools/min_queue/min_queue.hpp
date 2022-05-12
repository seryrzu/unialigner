//
// Created by Andrey Bzikadze on 05/10/22.
//

#pragma once

#include <deque>
#include "sequences/verify.hpp"

namespace min_queue {

// MinQueue is a monotonic queue that saves indexes in the corresponding random access container
template<class T, class Compare = std::less<typename T::value_type>>
class MinQueue {
    const T &vec;
    // elements in vector with deque indexes are non-increasing
    std::deque<int64_t> deque;
    int next{0}, first{0};
    Compare value_compare;

 public:
    explicit MinQueue(const T &vec) : vec{vec} {}
    void PushBack() {
        VERIFY(next < vec.size());
        while (not deque.empty() and value_compare(vec[next], vec[deque.back()])) {
            deque.pop_back();
        }
        deque.push_back(next);
        next++;
    }

    [[nodiscard]] int64_t GetMinIndex() const {
        VERIFY(not deque.empty());
        return deque.front();
    }

    typename T::value_type GetMin() const {
        return vec[GetMinIndex()];
    }

    void PopFront() {
        if (not deque.empty() and first==deque.front()) {
            deque.pop_front();
        }
        first++;
    }

};
} // namespace min_queue