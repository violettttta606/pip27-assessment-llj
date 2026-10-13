#include "bipartite_matching.hpp"

#include <iostream>
#include <vector>

int main() {
    std::vector<Edge> edges{
        {0, 0, 8},
        {0, 1, 5},
        {1, 0, 5}
    };

    BipartiteMatcher matcher(2, 2, edges);
    const auto answer = matcher.maximumWeightMatching();

    std::cout << "Matching:";
    for (const auto right : answer) {
        std::cout << ' ' << right;
    }
    std::cout << '\n';

    if (answer == std::vector<std::int32_t>{1, 0}) {
        std::cout << "PASS: maximum total weight is 10\n";
        return 0;
    }

    std::cout << "FAIL: expected matching 1 0\n";
    return 1;
}
