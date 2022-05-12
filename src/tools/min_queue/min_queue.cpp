//
// Created by Andrey Bzikadze on 05/10/22.
//

#include <iostream>
#include <vector>
#include "min_queue.hpp"

using namespace min_queue;
using std::cout, std::cin, std::vector;

int main() {
    int N{0}, k{0};
    cin >> N >> k;
    vector<int> vec;
    while (N--) {
        int x;
        cin >> x;
        vec.push_back(x);
    }

    MinQueue<decltype(vec)> queue(vec);
    for (int i = 0; i < vec.size(); ++i) {
        queue.PushBack();
        if (i >= k - 1) {
            cout << queue.GetMin() << " ";
            queue.PopFront();
        }
    }
}
