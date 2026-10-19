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
        //1
        if (edges.empty()) {
            return answer;
        }
        //2
        const Edge* best_edge = &edges.front();
        for (const Edge& edge : edges) {
            if (edge.weight > best_edge->weight) {
                best_edge = &edge;
            }
        }
        //3
        if (best_edge->weight <= 0) {
            answer[best_edge->left] =
                static_cast<std::int32_t>(best_edge->right);
            return answer;
        }
        //4
        const std::uint32_t row_count = left_count;
        const std::uint32_t real_column_count = right_count;
        const std::uint32_t column_count =
            right_count + left_count;
        //5
        const std::int64_t FORBIDDEN =
            std::int64_t{1} << 50;
        const std::int64_t INF =
            std::numeric_limits<std::int64_t>::max() / 4;
        std::vector<std::vector<std::int64_t>> cost(
            row_count,
            std::vector<std::int64_t>(
                real_column_count, FORBIDDEN));
        for (const Edge& edge : edges) {
            cost[edge.left][edge.right] = -edge.weight;
        }
        //6辅助数组
        std::vector<std::int64_t> row_potential(
            row_count + 1, 0);

        std::vector<std::int64_t> column_potential(
            column_count + 1, 0);
        // 值为 0 表示该列尚未匹配。
        std::vector<std::uint32_t> column_match(
            column_count + 1, 0);
        //previous前一列
        std::vector<std::uint32_t> previous_column(
            column_count + 1, 0);
        // slack当前搜索的最小调整费用
        std::vector<std::int64_t> slack(
            column_count + 1, INF);
        // used是否已经加入搜索
        std::vector<bool> used(
            column_count + 1, false);
        //7
        for (std::uint32_t row = 1;
             row <= row_count;
             ++row) {
            column_match[0] = row;
            std::fill(slack.begin(), slack.end(), INF);
            std::fill(used.begin(), used.end(), false);
            std::fill(
                previous_column.begin(),
                previous_column.end(),
                0);
            std::uint32_t current_column = 0;
            // 寻找
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
                    // 最佳
                    if (reduced_cost < slack[column]) {
                        slack[column] = reduced_cost;
                        previous_column[column] =
                            current_column;
                    }
                    // 未加入
                    if (slack[column] < delta) {
                        delta = slack[column];
                        next_column = column;
                    }
                }
                // 更新
                for (std::uint32_t column = 0;
                     column <= column_count;
                     ++column) {
                    if (used[column]) {
                        row_potential[
                            column_match[column]] += delta;

                        column_potential[column] -= delta;
                    } else if (column != 0) {
                        slack[column] -= delta;
                    }
                }
                current_column = next_column;
            } while (column_match[current_column] != 0);
            //空闲列,倒着更新匹配
            do {
                const std::uint32_t preceding_column =
                    previous_column[current_column];
                column_match[current_column] =
                    column_match[preceding_column];
                current_column = preceding_column;
            } while (current_column != 0);
        }
        //8忽略虚拟列。
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
std::vector<std::int32_t>
BipartiteMatcher::maximumWeightMatching() const {
    return impl_->solve();
}