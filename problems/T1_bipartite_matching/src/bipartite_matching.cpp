#include "bipartite_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

using namespace std;

struct BipartiteMatcher::Impl {
    Impl(uint32_t left_count,
         uint32_t right_count,
         const vector<Edge>& edges)
        : left_count(left_count),
          right_count(right_count),
          edges(edges) {}

    uint32_t left_count;
    uint32_t right_count;
    vector<Edge> edges;

    vector<int32_t> solve() const {
        vector<int32_t> answer(left_count, -1);
        // 无边
        if (edges.empty()) {
            return answer;
        }

        const Edge* best_edge = &edges.front();
        for (const Edge& edge : edges) {
            if (edge.weight > best_edge->weight) {
                best_edge = &edge;
            }
        }
        // 非正权时max
        if (best_edge->weight <= 0) {
            answer[best_edge->left] =
                static_cast<int32_t>(best_edge->right);
            return answer;
        }
        // 左点不匹配
        const uint32_t row_count = left_count;
        const uint32_t real_column_count = right_count;
        const uint32_t column_count = right_count + left_count;
        // 缺失边使用大费用，真实边费用取权重的相反数
        const int64_t FORBIDDEN = int64_t{1} << 50;
        const int64_t INF =
            numeric_limits<int64_t>::max() / 4;
        vector<vector<int64_t>> cost(
            row_count,
            vector<int64_t>(real_column_count, FORBIDDEN));
        for (const Edge& edge : edges) {
            cost[edge.left][edge.right] = -edge.weight;
        }

        vector<int64_t> row_potential(
            row_count + 1, 0);

        vector<int64_t> column_potential(
            column_count + 1, 0);
        // 值为 0 表示该列尚未匹配
        vector<uint32_t> column_match(
            column_count + 1, 0);
        // 增广路径中的前一列
        vector<uint32_t> previous_column(
            column_count + 1, 0);
        // 搜索树到各列的最小调整费用
        vector<int64_t> slack(
            column_count + 1, INF);
        // 标记已经进入搜索树的列
        vector<bool> used(
            column_count + 1, false);
        // 逐行增广
        for (uint32_t row = 1; row <= row_count; ++row) {
            column_match[0] = row;
            fill(slack.begin(), slack.end(), INF);
            fill(used.begin(), used.end(), false);
            fill(previous_column.begin(), previous_column.end(), 0);
            uint32_t current_column = 0;
            // 找通往空闲列的增广路径
            do {
                used[current_column] = true;

                const uint32_t current_row =
                    column_match[current_column];

                int64_t delta = INF;
                uint32_t next_column = 0;

                for (uint32_t column = 1;
                     column <= column_count;
                     ++column) {

                    if (used[column]) {
                        continue;
                    }
                    int64_t edge_cost = 0;

                    if (column <= real_column_count) {
                        edge_cost =
                            cost[current_row - 1][column - 1];
                    }

                    const int64_t reduced_cost =
                        edge_cost
                        - row_potential[current_row]
                        - column_potential[column];
                    // 记录到达该列的更优路径
                    if (reduced_cost < slack[column]) {
                        slack[column] = reduced_cost;
                        previous_column[column] =
                            current_column;
                    }
                    // 选择调整费用最小的未访问列
                    if (slack[column] < delta) {
                        delta = slack[column];
                        next_column = column;
                    }
                }
                // 更新势与未访问列的调整费用
                for (uint32_t column = 0;
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
                const uint32_t preceding_column =
                    previous_column[current_column];
                column_match[current_column] =
                    column_match[preceding_column];
                current_column = preceding_column;
            } while (current_column != 0);
        }
        //提取真实列
        for (uint32_t column = 1;
             column <= real_column_count;
             ++column) {
            const uint32_t row =
                column_match[column];

            if (row != 0) {
                answer[row - 1] =
                    static_cast<int32_t>(column - 1);
            }
        }
        return answer;
    }
};

BipartiteMatcher::BipartiteMatcher(
    uint32_t left_count,
    uint32_t right_count,
    const vector<Edge>& edges)
    : impl_(make_unique<Impl>(
          left_count, right_count, edges)) {}

BipartiteMatcher::~BipartiteMatcher() = default;

vector<int32_t> BipartiteMatcher::maximumWeightMatching() const {
    return impl_->solve();
}
