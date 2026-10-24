#include "bipartite_matching.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct TestGraph {
    std::string name;
    std::uint32_t left;
    std::uint32_t right;
    std::vector<Edge> edges;
    std::int64_t expected;
};

// Diagonal edges certify a full matching. Row-constant weights give an
// independently known upper bound attained by that matching.
TestGraph planted(const std::string& name, std::uint32_t left,
                  std::uint32_t right, std::uint32_t edge_count,
                  bool row_weights) {
    TestGraph graph{name, left, right, {}, 0};
    std::vector<bool> present(static_cast<std::size_t>(left) * right, false);
    const auto weight = [row_weights](std::uint32_t row) -> std::int64_t {
        return row_weights ? 1000000000 - row : 1;
    };
    const auto add = [&](std::uint32_t row, std::uint32_t column) {
        const auto key = static_cast<std::size_t>(row) * right + column;
        if (!present[key]) {
            present[key] = true;
            graph.edges.push_back({row, column, weight(row)});
        }
    };
    for (std::uint32_t i = 0; i < std::min(left, right); ++i) {
        add(i, i);
        graph.expected += weight(i);
    }
    std::mt19937 random(20261024);
    while (graph.edges.size() < edge_count) add(random() % left, random() % right);
    std::shuffle(graph.edges.begin(), graph.edges.end(), random);
    return graph;
}

TestGraph congested() {
    TestGraph graph{"sparse-congested", 2000, 2000, {}, 25};
    // Only 25 right vertices can be used. Dense virtual-column assignment
    // must nevertheless process all 2000 left vertices.
    for (std::uint32_t i = 0; i < graph.left; ++i) {
        for (std::uint32_t j = 0; j < 25; ++j) graph.edges.push_back({i, j, 1});
    }
    return graph;
}

TestGraph varying_weights(bool mixed) {
    TestGraph graph{mixed ? "sparse-mixed" : "sparse-varied", 2000, 2000, {}, 0};
    std::vector<bool> present(4000000, false);
    // The sum of row maxima certifies the planted diagonal's optimality.
    for (std::uint32_t i = 0; i < 2000; ++i) {
        const std::int64_t weight = 500000000 + i * 200000;
        graph.edges.push_back({i, i, weight});
        present[i * 2000 + i] = true;
        graph.expected += weight;
    }
    std::mt19937 random(1729);
    while (graph.edges.size() < 50000) {
        const auto row = random() % 2000;
        const auto column = random() % 2000;
        const auto key = row * 2000 + column;
        if (present[key]) continue;
        present[key] = true;
        std::int64_t weight = random() % (500000000 + row * 200000);
        if (mixed && random() % 2) weight = -weight;
        graph.edges.push_back({static_cast<std::uint32_t>(row),
                               static_cast<std::uint32_t>(column), weight});
    }
    std::shuffle(graph.edges.begin(), graph.edges.end(), random);
    return graph;
}

void benchmark(const TestGraph& graph) {
    const auto start = std::chrono::steady_clock::now();
    BipartiteMatcher matcher(graph.left, graph.right, graph.edges);
    const auto answer = matcher.maximumWeightMatching();
    const auto elapsed = std::chrono::steady_clock::now() - start;
    std::vector<std::int64_t> weight(
        static_cast<std::size_t>(graph.left) * graph.right, -1);
    for (const auto& edge : graph.edges) weight[edge.left * graph.right + edge.right] = edge.weight;
    if (answer.size() != graph.left) throw std::runtime_error("wrong length");
    std::vector<bool> used(graph.right, false);
    std::int64_t score = 0;
    for (std::uint32_t row = 0; row < graph.left; ++row) {
        const auto column = answer[row];
        if (column == -1) continue;
        if (column < 0 || static_cast<std::uint32_t>(column) >= graph.right || used[column])
            throw std::runtime_error("invalid matching");
        used[column] = true;
        const auto value = weight[row * graph.right + column];
        if (value <= 0) throw std::runtime_error("nonpositive or missing edge selected");
        score += value;
    }
    if (score != graph.expected) throw std::runtime_error("wrong optimum for " + graph.name);
    std::cout << graph.name << " L=" << graph.left << " R=" << graph.right
              << " E=" << graph.edges.size() << " PASS " << std::fixed
              << std::setprecision(3)
              << std::chrono::duration<double, std::milli>(elapsed).count()
              << " ms\n" << std::flush;
}
}  // namespace

int main() {
    try {
        benchmark(planted("sparse-equal", 2000, 2000, 50000, false));
        benchmark(planted("sparse-row-weights", 2000, 2000, 50000, true));
        benchmark(congested());
        benchmark(varying_weights(false));
        benchmark(varying_weights(true));
        benchmark(planted("sparse-rectangular", 2000, 1700, 50000, false));
        benchmark(planted("dense-equal", 500, 500, 250000, false));
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
