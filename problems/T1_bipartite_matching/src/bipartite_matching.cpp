#include "bipartite_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

struct BipartiteMatcher::Impl {
    Impl(std::uint32_t left_count,
         std::uint32_t right_count,
         const std::vector<Edge>& edges)
        : left_count(left_count),
          right_count(right_count),
          edges(edges) {}

    std::uint32_t left_count;
    std::uint32_t right_count;
    std::vector<Edge> edges;

    std::vector<std::int32_t> solve() const {
        std::vector<std::int32_t> answer(left_count, -1);
        // 无边时没有可行的非空匹配。
        if (edges.empty()) {
            return answer;
        }
        // 记录最大单边，用于处理全部非正权的情况。
        const Edge* best_edge = &edges.front();
        for (const Edge& edge : edges) {
            if (edge.weight > best_edge->weight) {
                best_edge = &edge;
            }
        }
        // 非空约束下，全部非正权时选择最大单边。
        if (best_edge->weight <= 0) {
            answer[best_edge->left] =
                static_cast<std::int32_t>(best_edge->right);
            return answer;
        }
        // 虚拟列表示左侧点不匹配。
        const std::uint32_t row_count = left_count;
        const std::uint32_t real_column_count = right_count;
        const std::uint32_t column_count = right_count + left_count;
        // 缺失边使用大费用，真实边费用取权重的相反数。
        const std::int64_t FORBIDDEN = std::int64_t{1} << 50;
        const std::int64_t INF =
            std::numeric_limits<std::int64_t>::max() / 4;
        std::vector<std::vector<std::int64_t>> cost(
            row_count,
            std::vector<std::int64_t>(real_column_count, FORBIDDEN));
        for (const Edge& edge : edges) {
            cost[edge.left][edge.right] = -edge.weight;
        }
        // 行列势、列匹配以及增广搜索的辅助记录。
        std::vector<std::int64_t> row_potential(
            row_count + 1, 0);

        std::vector<std::int64_t> column_potential(
            column_count + 1, 0);
        // 值为 0 表示该列尚未匹配。
        std::vector<std::uint32_t> column_match(
            column_count + 1, 0);
        // 增广路径中的前一列。
        std::vector<std::uint32_t> previous_column(
            column_count + 1, 0);
        // 搜索树到各列的最小调整费用。
        std::vector<std::int64_t> slack(
            column_count + 1, INF);
        // 标记已经进入搜索树的列。
        std::vector<bool> used(
            column_count + 1, false);
        // 逐行增广，每轮为一个左侧点寻找匹配。
        for (std::uint32_t row = 1; row <= row_count; ++row) {
            column_match[0] = row;
            std::fill(slack.begin(), slack.end(), INF);
            std::fill(used.begin(), used.end(), false);
            std::fill(previous_column.begin(), previous_column.end(), 0);
            std::uint32_t current_column = 0;
            // 寻找通往空闲列的增广路径。
            do {
                used[current_column] = true;

                const std::uint32_t current_row =
                    column_match[current_column];

                std::int64_t delta = INF;
                std::uint32_t next_column = 0;

                for (std::uint32_t column = 1;
                     column <= column_count;
                     ++column) {

                    if (used[column]) {
                        continue;
                    }
                    std::int64_t edge_cost = 0;

                    if (column <= real_column_count) {
                        edge_cost =
                            cost[current_row - 1][column - 1];
                    }

                    const std::int64_t reduced_cost =
                        edge_cost
                        - row_potential[current_row]
                        - column_potential[column];
                    // 记录到达该列的更优路径。
                    if (reduced_cost < slack[column]) {
                        slack[column] = reduced_cost;
                        previous_column[column] =
                            current_column;
                    }
                    // 选择调整费用最小的未访问列。
                    if (slack[column] < delta) {
                        delta = slack[column];
                        next_column = column;
                    }
                }
                // 更新势与未访问列的调整费用。
                for (std::uint32_t column = 0;
                     column <= column_count;
                     ++column) {
                    if (used[column]) {
                        row_potential[column_match[column]] += delta;

                        column_potential[column] -= delta;
                    } else if (column != 0) {
                        slack[column] -= delta;
                    }
                }
                current_column = next_column;
            } while (column_match[current_column] != 0);
            //更新匹配
            do {
                const std::uint32_t preceding_column =
                    previous_column[current_column];
                column_match[current_column] =
                    column_match[preceding_column];
                current_column = preceding_column;
            } while (current_column != 0);
        }
        //提取真实列
        for (std::uint32_t column = 1;
             column <= real_column_count;
             ++column) {
            const std::uint32_t row =
                column_match[column];

            if (row != 0) {
                answer[row - 1] =
                    static_cast<std::int32_t>(column - 1);
            }
        }
        return answer;
    }
};

BipartiteMatcher::BipartiteMatcher(
    std::uint32_t left_count,
    std::uint32_t right_count,
    const std::vector<Edge>& edges)
    : impl_(std::make_unique<Impl>(
          left_count, right_count, edges)) {}

BipartiteMatcher::~BipartiteMatcher() = default;

std::vector<std::int32_t> BipartiteMatcher::maximumWeightMatching() const {
    return impl_->solve();
}