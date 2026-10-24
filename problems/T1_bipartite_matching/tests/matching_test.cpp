#include "bipartite_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Score = std::int64_t;
constexpr Score unreachable = std::numeric_limits<Score>::min() / 4;
std::uint64_t case_count = 0;

// Independent subset DP: only the previous row's states can add an edge.
Score optimum(std::uint32_t left, std::uint32_t right,
              const std::vector<Edge>& edges) {
    if (edges.empty()) return 0;
    std::vector<Score> dp(std::size_t{1} << right, unreachable);
    dp[0] = 0;
    for (std::uint32_t row = 0; row < left; ++row) {
        auto next = dp;
        for (const auto& edge : edges) {
            if (edge.left != row) continue;
            const std::size_t bit = std::size_t{1} << edge.right;
            for (std::size_t mask = 0; mask < dp.size(); ++mask) {
                if ((mask & bit) == 0 && dp[mask] != unreachable) {
                    next[mask | bit] = std::max(next[mask | bit],
                                               dp[mask] + edge.weight);
                }
            }
        }
        dp.swap(next);
    }
    // A nonzero subset is required even when the best score is negative.
    return *std::max_element(dp.begin() + 1, dp.end());
}

void validate(const std::string& name, std::uint32_t left,
              std::uint32_t right, const std::vector<Edge>& edges,
              const std::vector<std::int32_t>& answer, Score expected) {
    const auto fail = [&](const std::string& reason) {
        std::cerr << "FAIL " << name << " L=" << left << " R=" << right
                  << " expected=" << expected << ": " << reason << '\n';
        for (const auto& edge : edges) {
            std::cerr << edge.left << ' ' << edge.right << ' '
                      << edge.weight << '\n';
        }
        throw std::runtime_error(reason);
    };
    if (answer.size() != left) fail("wrong answer length");
    std::vector<bool> used(right, false);
    Score score = 0;
    std::size_t matched = 0;
    for (std::uint32_t row = 0; row < left; ++row) {
        const auto column = answer[row];
        if (column == -1) continue;
        if (column < 0 || static_cast<std::uint32_t>(column) >= right)
            fail("right vertex out of range");
        if (used[column]) fail("duplicate right vertex");
        used[column] = true;
        const auto edge = std::find_if(edges.begin(), edges.end(),
            [&](const Edge& e) {
                return e.left == row && e.right == static_cast<std::uint32_t>(column);
            });
        if (edge == edges.end()) fail("returned a missing edge");
        score += edge->weight;
        ++matched;
    }
    if (!edges.empty() && matched == 0) fail("empty matching with input edges");
    if (score != expected) fail("suboptimal score " + std::to_string(score));
}

void check(const std::string& name, std::uint32_t left,
           std::uint32_t right, const std::vector<Edge>& edges,
           Score expected) {
    BipartiteMatcher matcher(left, right, edges);
    validate(name, left, right, edges, matcher.maximumWeightMatching(), expected);
    ++case_count;
}

void directed() {
    check("empty 0x0", 0, 0, {}, 0);
    check("empty left", 0, 5, {}, 0);
    check("empty right", 5, 0, {}, 0);
    check("empty edges", 4, 7, {}, 0);
    check("more left", 5, 2, {{0, 0, 8}, {3, 1, 9}, {4, 0, 7}}, 17);
    check("more right", 2, 5, {{0, 0, 8}, {1, 4, 9}, {1, 0, 7}}, 17);
    check("all negative", 5, 2, {{0, 0, -5}, {4, 1, -3}}, -3);
    check("zero best", 2, 5, {{0, 4, -5}, {1, 3, 0}}, 0);
    check("negative edge remains unused", 2, 3, {{0, 0, 8}, {1, 1, -7}}, 8);
    check("weight before cardinality", 2, 2, {{0, 0, 100}, {0, 1, 1}, {1, 0, 1}}, 100);
    check("augmenting reassignment", 2, 2, {{0, 0, 8}, {0, 1, 5}, {1, 0, 5}}, 10);
    check("64 bit total", 3, 4,
          {{0, 0, 1000000000}, {1, 1, 1000000000}, {2, 3, 1000000000}},
          Score{3000000000});
    check("20x1", 20, 1, {{0, 0, 1}, {19, 0, 1000000000}}, 1000000000);
    check("2000x1 negative", 2000, 1, {{0, 0, -1000000000}, {1999, 0, -1}}, -1);
    check("1x2000", 1, 2000, {{0, 0, -1000000000}, {0, 1999, 1000000000}}, 1000000000);
    const std::vector<Edge> original{{0, 0, 8}, {0, 1, 5}, {1, 0, 5}};
    auto input = original;
    BipartiteMatcher matcher(2, 2, input);
    input.clear();
    for (int call = 0; call < 3; ++call) {
        validate("copied input and repeated calls", 2, 2, original,
                 matcher.maximumWeightMatching(), 10);
        ++case_count;
    }
    std::cout << "PASS directed: " << case_count << " cases\n";
}

void exhaustive() {
    const auto start = case_count;
    for (std::uint32_t left = 0; left <= 3; ++left) {
        for (std::uint32_t right = 0; right <= 3; ++right) {
            const std::uint32_t possibilities = std::uint32_t{1} << (2 * left * right);
            for (std::uint32_t code = 0; code < possibilities; ++code) {
                auto digits = code;
                std::vector<Edge> edges;
                for (std::uint32_t i = 0; i < left; ++i) {
                    for (std::uint32_t j = 0; j < right; ++j) {
                        const auto state = digits & 3;
                        digits >>= 2;
                        if (state != 0) edges.push_back({i, j, Score{state} - 2});
                    }
                }
                check("exhaustive " + std::to_string(code), left, right, edges,
                      optimum(left, right, edges));
            }
        }
    }
    std::cout << "PASS exhaustive: " << case_count - start << " cases\n";
}

void random_cases(std::uint32_t count) {
    const auto start = case_count;
    std::mt19937 random(20261024);
    for (std::uint32_t trial = 0; trial < count; ++trial) {
        const auto left = static_cast<std::uint32_t>(random() % 9);
        const auto right = static_cast<std::uint32_t>(random() % 9);
        const auto density = random() % 101;
        std::vector<Edge> edges;
        for (std::uint32_t i = 0; i < left; ++i) {
            for (std::uint32_t j = 0; j < right; ++j) {
                if (random() % 100 >= density) continue;
                Score weight = static_cast<Score>(random() % 21) - 10;
                if (trial % 4 == 0) weight *= 100000000;
                if (trial % 4 == 1) weight = -static_cast<Score>(random() % 1000000001);
                if (trial % 4 == 2) weight = random() % 2 == 0 ? 0 : 1000000000;
                edges.push_back({i, j, weight});
            }
        }
        std::shuffle(edges.begin(), edges.end(), random);
        check("random " + std::to_string(trial), left, right, edges,
              optimum(left, right, edges));
    }
    std::cout << "PASS random: " << case_count - start << " cases\n";
}
}  // namespace

int main(int argc, char**) {
    try {
        directed();
        if (argc == 1) exhaustive();
        random_cases(argc == 1 ? 100000 : 10000);
        std::cout << "ALL PASSED: " << case_count << " cases\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
