#include "bipartite_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <utility>
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

    struct FlowEdge {
        uint32_t to;
        uint32_t reverse;
        int64_t capacity;
        int64_t cost;
    };

    using FlowGraph = vector<vector<FlowEdge>>;

    static void addFlowEdge(FlowGraph& graph, uint32_t from, uint32_t to,
                            int64_t capacity, int64_t cost) {
        const uint32_t forward_index =
            static_cast<uint32_t>(graph[to].size());
        const uint32_t reverse_index =
            static_cast<uint32_t>(graph[from].size());
        graph[from].push_back(
            {to, forward_index, capacity, cost});
        graph[to].push_back(
            {from, reverse_index, 0, -cost});
    }

    vector<int32_t> solveSparse(
        const vector<Edge>& positive_edges) const {
        vector<int32_t> answer(left_count, -1);

        // Put the smaller partition on the source side. This bounds the
        // number of augmentations by min(left_count, right_count).
        const bool transpose = left_count > right_count;
        const uint32_t source_count = transpose ? right_count : left_count;
        const uint32_t target_count = transpose ? left_count : right_count;
        const uint32_t source = 0;
        const uint32_t source_base = 1;
        const uint32_t target_base = source_base + source_count;
        const uint32_t sink = target_base + target_count;
        const uint32_t node_count = sink + 1;
        FlowGraph graph(node_count);

        for (uint32_t source_vertex = 0;
             source_vertex < source_count; ++source_vertex) {
            addFlowEdge(graph, source, source_base + source_vertex,
                        1, 0);
        }
        for (uint32_t target_vertex = 0;
             target_vertex < target_count; ++target_vertex) {
            addFlowEdge(graph, target_base + target_vertex, sink, 1, 0);
        }

        struct MatchingEdge {
            uint32_t left;
            uint32_t right;
            uint32_t from;
            uint32_t index;
        };
        vector<MatchingEdge> matching_edges;
        matching_edges.reserve(positive_edges.size());
        vector<int64_t> potential(node_count, 0);

        // Initial shortest-path potentials make every residual reduced cost
        // nonnegative, despite negative costs on positive-weight edges.
        for (const Edge& edge : positive_edges) {
            const uint32_t source_vertex =
                transpose ? edge.right : edge.left;
            const uint32_t target_vertex =
                transpose ? edge.left : edge.right;
            potential[target_base + target_vertex] = std::min(
                potential[target_base + target_vertex], -edge.weight);
            const uint32_t from = source_base + source_vertex;
            const uint32_t index = static_cast<uint32_t>(graph[from].size());
            addFlowEdge(graph, from, target_base + target_vertex,
                        1, -edge.weight);
            matching_edges.push_back(
                {edge.left, edge.right, from, index});
        }
        for (uint32_t target_vertex = 0;
             target_vertex < target_count; ++target_vertex) {
            potential[sink] = std::min(
                potential[sink], potential[target_base + target_vertex]);
        }

        const int64_t infinity = numeric_limits<int64_t>::max() / 4;
        using QueueItem = pair<int64_t, uint32_t>;
        struct QueueCompare {
            bool operator()(const QueueItem& first,
                            const QueueItem& second) const {
                if (first.first != second.first) {
                    return first.first > second.first;
                }
                // The sink has the largest index: settle it first on ties.
                return first.second < second.second;
            }
        };
        priority_queue<QueueItem, vector<QueueItem>, QueueCompare> queue;
        vector<int64_t> distance(node_count);
        vector<uint32_t> previous_node(node_count);
        vector<uint32_t> previous_edge(node_count);

        while (true) {
            fill(distance.begin(), distance.end(), infinity);
            distance[source] = 0;
            queue = {};
            queue.emplace(0, source);
            while (!queue.empty()) {
                const auto [current_distance, node] = queue.top();
                queue.pop();
                if (current_distance != distance[node]) {
                    continue;
                }
                if (node == sink) {
                    break;
                }
                for (uint32_t index = 0;
                     index < graph[node].size(); ++index) {
                    const FlowEdge& edge = graph[node][index];
                    if (edge.capacity == 0) {
                        continue;
                    }
                    const int64_t reduced_cost =
                        edge.cost + potential[node] - potential[edge.to];
                    const int64_t next_distance =
                        current_distance + reduced_cost;
                    if (next_distance < distance[edge.to]) {
                        distance[edge.to] = next_distance;
                        previous_node[edge.to] = node;
                        previous_edge[edge.to] = index;
                        queue.emplace(next_distance, edge.to);
                    }
                }
            }
            if (distance[sink] == infinity) {
                break;
            }
            const int64_t marginal_cost =
                distance[sink] + potential[sink] - potential[source];
            if (marginal_cost >= 0) {
                break;
            }
            // Dijkstra stops at the sink. Capping all other distances at
            // that settled distance keeps every residual reduced cost
            // nonnegative, including edges from unvisited vertices.
            for (uint32_t node = 0; node < node_count; ++node) {
                potential[node] += std::min(distance[node], distance[sink]);
            }
            uint32_t node = sink;
            while (node != source) {
                FlowEdge& edge =
                    graph[previous_node[node]][previous_edge[node]];
                --edge.capacity;
                ++graph[node][edge.reverse].capacity;
                node = previous_node[node];
            }
        }

        for (const MatchingEdge& matching_edge : matching_edges) {
            if (graph[matching_edge.from][matching_edge.index].capacity == 0) {
                answer[matching_edge.left] =
                    static_cast<int32_t>(matching_edge.right);
            }
        }
        return answer;
    }

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

        vector<Edge> positive_edges;
        positive_edges.reserve(edges.size());
        for (const Edge& edge : edges) {
            if (edge.weight > 0) {
                positive_edges.push_back(edge);
            }
        }
        const uint64_t matrix_size =
            static_cast<uint64_t>(left_count) * right_count;
        // Sparse min-cost flow avoids scanning absent edges. Keep the dense
        // Hungarian path for the smaller dense test points.
        if (positive_edges.size() <= 100000 &&
            matrix_size > static_cast<uint64_t>(positive_edges.size()) * 4) {
            return solveSparse(positive_edges);
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
