// PC-GMSTP GRASP with adaptive memory: independent reproduction for author review.
// Bruno Jose da Silva Barros (UFF, 2025), chapter 4. Not the author's original code.
// Compile: c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic grasp.cpp -o grasp
// The algorithm, GTP reader, CLI and small independent checks are all in this file.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace timing {
class CpuLimit : public std::runtime_error {
public:
    CpuLimit() : std::runtime_error("Process CPU limit reached") {}
};
double processCpuSeconds() {
    const auto value = std::clock();
    if (value == static_cast<std::clock_t>(-1)) throw std::runtime_error("Cannot read process CPU clock");
    return static_cast<double>(value) / CLOCKS_PER_SEC;
}
void checkDeadline(double deadline) {
    if (deadline > 0 && processCpuSeconds() >= deadline) throw CpuLimit();
}
}

// 1. GTP input and objective: rounded Euclidean MST minus selected vertex prizes.


namespace pcgmstp
{
struct Instance
{
    std::vector<std::pair<std::int64_t, std::int64_t>> coordinates;
    std::vector<std::vector<int>> groups;
    std::vector<int> groupOf;
    std::vector<std::int64_t> prizes;
    static Instance read(std::istream &input, double deadline = 0);
    std::int64_t distance(int i, int j) const;
    void validateRepresentatives(const std::vector<int> &representatives) const;
};
} // namespace pcgmstp


namespace pcgmstp
{
namespace
{
void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}
}

Instance Instance::read(std::istream &input, double deadline)
{
    std::string raw((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    while (!raw.empty() && (std::isspace(static_cast<unsigned char>(raw.back())) || raw.back() == '\x1a'))
        raw.pop_back();
    std::istringstream stream(raw);
    auto integer = [&]() {
        timing::checkDeadline(deadline);
        std::string token;
        require(static_cast<bool>(stream >> token), "Truncated GTP instance");
        std::size_t consumed = 0;
        const long double value = std::stold(token, &consumed);
        require(consumed == token.size() && std::isfinite(value) && std::trunc(value) == value &&
                    std::abs(value) <= 1000000000.0L, "Invalid or unsupported GTP integer");
        return static_cast<std::int64_t>(value);
    };
    const auto n = integer();
    require(n >= 1 && n <= 10000, "Unsupported vertex count");
    Instance result;
    for (int i = 0; i < n; ++i)
    {
        const auto x = integer();
        const auto y = integer();
        result.coordinates.emplace_back(x, y);
    }
    const auto k = integer();
    require(k >= 1 && k <= n, "Invalid group count");
    result.groupOf.assign(n, -1);
    for (int g = 0; g < k; ++g)
    {
        const auto size = integer();
        require(size >= 1 && size <= n, "Invalid group size");
        std::vector<int> group;
        for (int t = 0; t < size; ++t)
        {
            const auto vertex = integer() - 1;
            require(vertex >= 0 && vertex < n && result.groupOf[vertex] == -1,
                    "GTP groups are not a disjoint partition");
            result.groupOf[vertex] = g;
            group.push_back(static_cast<int>(vertex));
        }
        result.groups.push_back(group);
    }
    require(std::find(result.groupOf.begin(), result.groupOf.end(), -1) == result.groupOf.end(),
            "GTP partition omits a vertex");
    for (int g = 0; g < k; ++g)
    {
        const auto center = integer() - 1;
        require(center >= 0 && center < n && result.groupOf[center] == g,
                "GTP center is outside its group");
    }
    require(integer() == -999, "Missing GTP first sentinel");
    for (int i = 0; i < n; ++i)
    {
        const auto prize = integer();
        require(prize >= 0, "Negative GTP prize");
        result.prizes.push_back(prize);
    }
    require(integer() == -999, "Missing GTP final sentinel");
    std::string extra;
    require(!(stream >> extra), "Unexpected GTP trailing data");
    return result;
}

std::int64_t Instance::distance(int i, int j) const
{
    const auto dx = coordinates.at(i).first - coordinates.at(j).first;
    const auto dy = coordinates.at(i).second - coordinates.at(j).second;
    const auto squared = dx * dx + dy * dy; // coordinates bounded so this sum fits int64.
    auto root = static_cast<std::int64_t>(std::sqrt(static_cast<long double>(squared)));
    while ((root + 1) * (root + 1) <= squared) ++root;
    while (root * root > squared) --root;
    return root + static_cast<std::int64_t>(squared - root * root > root);
}

void Instance::validateRepresentatives(const std::vector<int> &representatives) const
{
    require(representatives.size() == groups.size(), "Wrong representative cardinality");
    for (std::size_t g = 0; g < groups.size(); ++g)
        require(representatives[g] >= 0 && representatives[g] < static_cast<int>(groupOf.size()) &&
                    groupOf[representatives[g]] == static_cast<int>(g), "Invalid representative");
}

} // namespace pcgmstp

// 2. Theorem 11: delete edges sequentially, retaining a third-group witness.
// Descending cost, then ascending endpoints. Missing witness edges invalidate a group.
namespace pcgmstp {
class EdgeReduction {
    std::size_t vertices;
    std::vector<unsigned char> kept;
public:
    std::size_t originalEdges = 0, removedEdges = 0;
    explicit EdgeReduction(const Instance &data, double deadline = 0)
        : vertices(data.coordinates.size()), kept(vertices * vertices, 0) {
        const auto inf = std::numeric_limits<std::int64_t>::max() / 4;
        const std::size_t k = data.groups.size();
        std::vector<std::int64_t> maximum(vertices * k, 0);
        struct Edge { std::int64_t cost; int u, v; };
        std::vector<Edge> edges;
        for (std::size_t u = 0; u < vertices; ++u) {
            timing::checkDeadline(deadline);
            maximum[u * k + data.groupOf[u]] = inf;
            for (std::size_t v = u + 1; v < vertices; ++v) {
                if (data.groupOf[u] == data.groupOf[v]) continue;
                const auto cost = data.distance(static_cast<int>(u), static_cast<int>(v));
                edges.push_back({cost, static_cast<int>(u), static_cast<int>(v)});
                kept[u * vertices + v] = kept[v * vertices + u] = 1;
                maximum[u * k + data.groupOf[v]] = std::max(maximum[u * k + data.groupOf[v]], cost);
                maximum[v * k + data.groupOf[u]] = std::max(maximum[v * k + data.groupOf[u]], cost);
            }
        }
        originalEdges = edges.size();
        std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b) {
            if (a.cost != b.cost) return a.cost > b.cost;
            if (a.u != b.u) return a.u < b.u;
            return a.v < b.v;
        });
        for (const auto &e : edges) {
            timing::checkDeadline(deadline);
            for (std::size_t g = 0; g < k; ++g) {
                if (maximum[e.u * k + g] > e.cost || maximum[e.v * k + g] > e.cost) continue;
                kept[e.u * vertices + e.v] = kept[e.v * vertices + e.u] = 0;
                ++removedEdges;
                maximum[e.u * k + data.groupOf[e.v]] = inf;
                maximum[e.v * k + data.groupOf[e.u]] = inf;
                break;
            }
        }
        timing::checkDeadline(deadline);
    }
    bool contains(int u, int v) const { return kept.at(u * vertices + v) != 0; }
};
}

// 3. GRASP-MA: roulette construction, representative-swap local search, memory.
// Section 4.5 / algorithms 7-8; interpretation choices are listed in README.md.


namespace pcgmstp::grasp_ma
{
using CpuLimit = timing::CpuLimit;

struct Config
{
    std::uint64_t seed = 1;
    std::size_t iterations = 0;
    double cpuLimit = 0;
    double cpuStartedAt = 0;
    bool partialTree = true;
};

struct Solution
{
    std::vector<int> representatives;
    // Endpoints are group indices, as required by Instance::treeCost.
    std::vector<std::pair<int, int>> tree;
    std::int64_t objective = 0;
};

struct Statistics
{
    std::uint64_t candidateEvaluations = 0;
    std::uint64_t localSearchPasses = 0;
    std::uint64_t acceptedMoves = 0;
    std::uint64_t partialForests = 0;
    std::uint64_t disconnectedPartialForests = 0;
};

struct Result
{
    bool complete = false;
    bool hasBest = false;
    std::string terminationReason;
    std::size_t requestedIterations = 0;
    std::size_t completedIterations = 0;
    std::size_t startedIterations = 0;
    std::size_t discardedPartialIterations = 0;
    Solution best;
    std::vector<std::uint64_t> memoryWeightsDoubled;
    Statistics statistics;
};

class Search
{
  public:
    using IterationCallback = std::function<void(std::size_t, const Solution &, double)>;
    Search(const Instance &instance, const EdgeReduction *reduction, Config config);
    Result run(const IterationCallback &onIteration = {});
    Solution evaluate(const std::vector<int> &representatives) const;

  private:
    struct Edge
    {
        std::int64_t cost;
        int a, b, u, v;
    };
    const Instance &data;
    const EdgeReduction *reduction;
    Config config;
    std::mt19937_64 generator;
    std::vector<std::int64_t> distances;
    std::vector<std::vector<int>> orderedGroups;
    mutable Statistics statistics;

    void checkCpu() const;
    std::uint64_t bounded(std::uint64_t bound);
    bool retained(int u, int v) const;
    Edge edge(int a, int b, const std::vector<int> &representatives) const;
    static bool edgeLess(const Edge &a, const Edge &b);
    std::vector<Edge> forest(const std::vector<int> &representatives, int omitted) const;
    Solution fromEdges(const std::vector<int> &representatives, std::vector<Edge> edges) const;
    Solution localSearch(Solution solution);
};
} // namespace pcgmstp::grasp_ma

namespace pcgmstp::grasp_ma
{
namespace
{
class DisjointSets
{
    std::vector<int> parent;

  public:
    explicit DisjointSets(std::size_t size) : parent(size)
    {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int v)
    {
        while (parent[v] != v)
        {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    }
    bool join(int a, int b)
    {
        a = find(a);
        b = find(b);
        if (a == b)
            return false;
        parent[b] = a;
        return true;
    }
};
} // namespace

Search::Search(const Instance &instance, const EdgeReduction *reduction, Config config)
    : data(instance), reduction(reduction), config(config), generator(config.seed),
      orderedGroups(instance.groups)
{
    if (!(config.cpuLimit > 0))
        throw std::runtime_error("A positive CPU limit is required");
    if (this->config.iterations == 0)
        this->config.iterations = std::max<std::size_t>(300, 3 * data.groups.size());
    for (auto &group : orderedGroups)
        std::sort(group.begin(), group.end());
    const std::size_t n = data.coordinates.size();
    if (n == 0 || data.groups.empty() || n > std::numeric_limits<std::size_t>::max() / n)
        throw std::runtime_error("Distance matrix dimensions overflow");
    distances.resize(n * n);
    for (std::size_t u = 0; u < n; ++u)
    {
        checkCpu();
        for (std::size_t v = u + 1; v < n; ++v)
            distances[u * n + v] = distances[v * n + u] = data.distance(u, v);
    }
}

void Search::checkCpu() const
{
    const double elapsed = timing::processCpuSeconds() - config.cpuStartedAt;
    if (elapsed >= config.cpuLimit)
        throw CpuLimit();
}

std::uint64_t Search::bounded(std::uint64_t bound)
{
    if (bound == 0)
        throw std::runtime_error("Empty random range");
    // Rejection sampling avoids modulo bias and implementation-dependent distributions.
    const std::uint64_t threshold = -bound % bound;
    for (;;)
    {
        const auto value = generator();
        if (value >= threshold)
            return value % bound;
    }
}

bool Search::retained(int u, int v) const
{
    return !reduction || reduction->contains(u, v);
}

Search::Edge Search::edge(int a, int b, const std::vector<int> &representatives) const
{
    const int u = std::min(representatives[a], representatives[b]);
    const int v = std::max(representatives[a], representatives[b]);
    return {distances[u * data.coordinates.size() + v], a, b, u, v};
}

bool Search::edgeLess(const Edge &a, const Edge &b)
{
    if (a.cost != b.cost)
        return a.cost < b.cost;
    if (a.u != b.u)
        return a.u < b.u;
    return a.v < b.v;
}

std::vector<Search::Edge> Search::forest(const std::vector<int> &representatives, int omitted) const
{
    checkCpu();
    std::vector<Edge> edges;
    const int k = representatives.size();
    for (int a = 0; a < k; ++a)
    {
        if (a == omitted)
            continue;
        checkCpu();
        for (int b = a + 1; b < k; ++b)
            if (b != omitted && retained(representatives[a], representatives[b]))
                edges.push_back(edge(a, b, representatives));
    }
    std::sort(edges.begin(), edges.end(), edgeLess);
    checkCpu();
    DisjointSets components(k);
    std::vector<Edge> result;
    for (const auto &candidate : edges)
        if (components.join(candidate.a, candidate.b))
            result.push_back(candidate);
    return result;
}

Solution Search::fromEdges(const std::vector<int> &representatives, std::vector<Edge> edges) const
{
    checkCpu();
    if (!std::is_sorted(edges.begin(), edges.end(), edgeLess))
        throw std::runtime_error("Kruskal input is not sorted");
    DisjointSets components(representatives.size());
    Solution result;
    result.representatives = representatives;
    for (const int vertex : representatives)
        result.objective -= data.prizes[vertex];
    for (const auto &candidate : edges)
        if (components.join(candidate.a, candidate.b))
        {
            result.objective += candidate.cost;
            result.tree.emplace_back(candidate.a, candidate.b);
        }
    if (result.tree.size() + 1 != representatives.size())
        throw std::runtime_error("Reduced representative graph is disconnected");
    checkCpu();
    return result;
}

Solution Search::evaluate(const std::vector<int> &representatives) const
{
    data.validateRepresentatives(representatives);
    return fromEdges(representatives, forest(representatives, -1));
}

Solution Search::localSearch(Solution solution)
{
    std::vector<int> order(data.groups.size());
    std::iota(order.begin(), order.end(), 0);
    bool improved;
    do
    {
        checkCpu();
        ++statistics.localSearchPasses;
        // Revisit every group until no strict improvement, following section 4.5's text.
        for (std::size_t i = order.size(); i > 1; --i)
            std::swap(order[i - 1], order[bounded(i)]);
        improved = false;
        for (const int group : order)
        {
            checkCpu();
            const auto partial =
                config.partialTree ? forest(solution.representatives, group) : std::vector<Edge>();
            if (config.partialTree)
            {
                ++statistics.partialForests;
                const std::size_t remaining = data.groups.size() - 1;
                if (remaining > 0 && partial.size() + 1 != remaining)
                    ++statistics.disconnectedPartialForests;
            }
            Solution best = solution;
            for (const int vertex : orderedGroups[group])
            {
                if (vertex == solution.representatives[group])
                    continue;
                checkCpu();
                ++statistics.candidateEvaluations;
                auto representatives = solution.representatives;
                representatives[group] = vertex;
                Solution candidate;
                if (config.partialTree)
                {
                    auto edges = partial;
                    for (std::size_t other = 0; other < representatives.size(); ++other)
                        if (static_cast<int>(other) != group &&
                            retained(vertex, representatives[other]))
                        {
                            const auto incident = edge(group, other, representatives);
                            edges.insert(
                                std::lower_bound(edges.begin(), edges.end(), incident, edgeLess),
                                incident);
                        }
                    // An MSF of the other representatives is sufficient even when removing
                    // this group disconnects their reduced graph. Keep every incident edge.
                    candidate = fromEdges(representatives, std::move(edges));
                }
                else
                    candidate = evaluate(representatives);
                if (candidate.objective < best.objective)
                    best = std::move(candidate);
            }
            if (best.objective < solution.objective)
            {
                solution = std::move(best);
                ++statistics.acceptedMoves;
                improved = true;
            }
        }
    } while (improved);
    return solution;
}

Result Search::run(const IterationCallback &onIteration)
{
    Result result;
    result.requestedIterations = config.iterations;
    // Integer weights are doubled: reset to 2, best-selected |Vc|, other-selected +2.
    result.memoryWeightsDoubled.assign(data.coordinates.size(), 2);
    try
    {
        for (std::size_t iteration = 1;
             iteration <= config.iterations; ++iteration)
        {
            checkCpu();
            result.startedIterations = iteration;
            std::vector<int> representatives;
            for (const auto &group : orderedGroups)
            {
                std::uint64_t total = 0;
                for (const int vertex : group)
                    total += result.memoryWeightsDoubled[vertex];
                auto draw = bounded(total);
                int chosen = -1;
                for (const int vertex : group)
                {
                    if (draw < result.memoryWeightsDoubled[vertex])
                    {
                        chosen = vertex;
                        break;
                    }
                    draw -= result.memoryWeightsDoubled[vertex];
                }
                if (chosen < 0)
                    throw std::runtime_error("Roulette did not select a vertex");
                representatives.push_back(chosen);
            }
            auto solution = localSearch(evaluate(representatives));
            checkCpu();
            const bool newBest = !result.hasBest || solution.objective < result.best.objective;
            if (newBest)
            {
                result.best = solution;
                result.hasBest = true;
                std::fill(result.memoryWeightsDoubled.begin(), result.memoryWeightsDoubled.end(),
                          2);
                for (std::size_t g = 0; g < orderedGroups.size(); ++g)
                    result.memoryWeightsDoubled[solution.representatives[g]] =
                        orderedGroups[g].size();
            }
            else
                for (const int vertex : solution.representatives)
                    result.memoryWeightsDoubled[vertex] += 2;
            result.completedIterations = iteration;
            if (onIteration)
                onIteration(iteration, solution, timing::processCpuSeconds() - config.cpuStartedAt);
            // Completion refers to complete requested iterations. Trace/output work
            // after the last one cannot turn them into an incomplete search.
            if (iteration < config.iterations) checkCpu();
        }
        result.complete = true;
        result.terminationReason = "iteration_limit";
    }
    catch (const CpuLimit &)
    {
        result.terminationReason = "cpu_limit";
    }
    result.discardedPartialIterations = result.startedIterations - result.completedIterations;
    result.statistics = statistics;
    return result;
}
} // namespace pcgmstp::grasp_ma

// 4. CLI. Vertex IDs in input and output start at 1; internal IDs start at 0.
namespace app {
using pcgmstp::Instance;
using pcgmstp::EdgeReduction;
using pcgmstp::grasp_ma::Config;
using pcgmstp::grasp_ma::Result;
using pcgmstp::grasp_ma::Search;
using pcgmstp::grasp_ma::Solution;

template <class T> void array(std::ostream &out, const std::vector<T> &values, int offset = 0) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i] + offset;
    out << ']';
}
void solutionJson(std::ostream &out, const Solution &s) {
    out << "{\"objective\":" << s.objective << ",\"representatives\":";
    array(out, s.representatives, 1);
    out << ",\"edges\":[";
    for (std::size_t i = 0; i < s.tree.size(); ++i) {
        const auto [a, b] = s.tree[i];
        out << (i ? "," : "") << '[' << s.representatives[a] + 1 << ',' << s.representatives[b] + 1 << ']';
    }
    out << "]}";
}
constexpr const char *example =
    "4\n0 0\n1 1\n4 0\n3 4\n2\n2 1 2\n2 3 4\n1 3\n-999\n9 1 10 2\n-999\n";

void help() {
    std::cout <<
        "Independent PC-GMSTP GRASP-MA (Barros, UFF 2025, chapter 4)\n"
        "Usage: grasp INSTANCE.GTP [options]\n"
        "  --seed N          Unsigned 64-bit seed (default: 1)\n"
        "  --iterations N    Complete GRASP iterations (default: max(300,3*groups))\n"
        "  --cpu-limit S     Cooperative process-CPU safety limit, seconds (default: 600)\n"
        "  --no-reduction    Keep all inter-group edges instead of theorem 11\n"
        "  --full-mst        Recompute each neighbor MST instead of partial-forest acceleration\n"
        "  --trace           Emit a JSON line for every completed local optimum\n"
        "  --example         Print a small synthetic GTP instance, then exit\n"
        "  --self-test       Run small independent mathematical and input checks\n"
        "  --help            Show this help\n"
        "Output: JSON/JSONL on stdout; errors on stderr. Exit 0 = all iterations complete,\n"
        "2 = CPU-stopped search (best completed incumbent, if any), 1 = invalid input/error.\n";
}
std::uint64_t unsignedInteger(const std::string &s) {
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Expected a nonnegative decimal integer");
    std::size_t used = 0;
    const auto n = std::stoull(s, &used);
    if (used != s.size()) throw std::runtime_error("Invalid integer");
    return n;
}
double positiveReal(const std::string &s) {
    std::size_t used = 0;
    const auto n = std::stod(s, &used);
    if (used != s.size() || !std::isfinite(n) || n <= 0)
        throw std::runtime_error("CPU limit must be finite and positive");
    return n;
}

// 5. Small independent checks. The tree oracle enumerates Prüfer sequences:
// no production Kruskal, distance rounding, MST routine or adaptive weights
// are used to compute the expected objective.
namespace checks {
void expect(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error("Self-test: " + why);
}
std::int64_t oracleDistance(const Instance &d, int a, int b) {
    const long double x = static_cast<long double>(d.coordinates[a].first) - d.coordinates[b].first;
    const long double y = static_cast<long double>(d.coordinates[a].second) - d.coordinates[b].second;
    return static_cast<std::int64_t>(std::floor(std::hypot(x, y) + .5L));
}
std::int64_t oracleMst(const Instance &d, const std::vector<int> &r) {
    const std::size_t k = r.size();
    std::int64_t prize = 0;
    for (int v : r) prize += d.prizes[v];
    if (k == 1) return -prize;
    std::size_t sequences = 1;
    for (std::size_t i = 0; i < k - 2; ++i) sequences *= k;
    auto best = std::numeric_limits<std::int64_t>::max();
    for (std::size_t code = 0; code < sequences; ++code) {
        auto digits = code;
        std::vector<int> sequence(k - 2), degree(k, 1);
        for (int &v : sequence) { v = static_cast<int>(digits % k); digits /= k; ++degree[v]; }
        std::int64_t cost = 0;
        for (int v : sequence) {
            const auto leaf = std::find(degree.begin(), degree.end(), 1) - degree.begin();
            cost += oracleDistance(d, r[leaf], r[v]);
            --degree[leaf]; --degree[v];
        }
        const auto a = std::find(degree.begin(), degree.end(), 1) - degree.begin();
        const auto b = std::find(degree.begin() + a + 1, degree.end(), 1) - degree.begin();
        cost += oracleDistance(d, r[a], r[b]);
        best = std::min(best, cost);
    }
    return best - prize;
}
Instance fixture(std::vector<std::pair<std::int64_t, std::int64_t>> xy,
                 std::vector<std::vector<int>> groups, std::vector<std::int64_t> prizes) {
    Instance d;
    d.coordinates = std::move(xy); d.groups = std::move(groups); d.prizes = std::move(prizes);
    d.groupOf.assign(d.coordinates.size(), -1);
    for (std::size_t g = 0; g < d.groups.size(); ++g)
        for (int v : d.groups[g]) d.groupOf[v] = static_cast<int>(g);
    return d;
}
void verifySolution(const Instance &d, const Solution &s, const EdgeReduction *reduction) {
    expect(s.representatives.size() == d.groups.size(), "representative cardinality");
    for (std::size_t g = 0; g < d.groups.size(); ++g)
        expect(std::find(d.groups[g].begin(), d.groups[g].end(), s.representatives[g]) != d.groups[g].end(),
               "representative outside group");
    expect(s.tree.size() + 1 == d.groups.size(), "tree cardinality");
    std::vector<std::vector<int>> adjacency(d.groups.size());
    std::int64_t cost = 0;
    for (const auto &[a, b] : s.tree) {
        expect(a >= 0 && b >= 0 && static_cast<std::size_t>(a) < d.groups.size() &&
               static_cast<std::size_t>(b) < d.groups.size() && a != b, "tree endpoint");
        adjacency[a].push_back(b); adjacency[b].push_back(a);
        expect(!reduction || reduction->contains(s.representatives[a], s.representatives[b]), "removed edge");
        cost += oracleDistance(d, s.representatives[a], s.representatives[b]);
    }
    std::vector<bool> visited(d.groups.size(), false);
    std::vector<int> queue{0}; visited[0] = true;
    for (std::size_t i = 0; i < queue.size(); ++i)
        for (int v : adjacency[queue[i]]) if (!visited[v]) { visited[v] = true; queue.push_back(v); }
    expect(queue.size() == d.groups.size(), "tree connectivity");
    for (int v : s.representatives) cost -= d.prizes[v];
    expect(cost == s.objective && cost == oracleMst(d, s.representatives), "independent tree objective");
}
Config config(std::uint64_t seed, bool partial) {
    Config c; c.seed = seed; c.iterations = 12; c.partialTree = partial;
    c.cpuStartedAt = timing::processCpuSeconds(); c.cpuLimit = 60;
    return c;
}
void run() {
    std::istringstream input(example);
    const auto parsed = Instance::read(input);
    expect(parsed.groups == std::vector<std::vector<int>>{{0, 1}, {2, 3}} && parsed.prizes[2] == 10,
           "GTP IDs and prizes");
    expect(oracleMst(parsed, {0, 2}) == -15, "known negative objective");
    const std::vector<std::pair<std::string, std::string>> invalid{
        {"0", "Unsupported vertex count"},
        {"2 0 0 1 1 1 0 1 -999 0 0 -999", "Invalid group size"},
        {"2 0 0 1 1 2 1 1 1 1 1 1 -999 0 0 -999", "disjoint partition"},
        {"2 0 0 1 1 1 1 1 1 -999 0 0 -999", "omits a vertex"},
        {"2 0 0 1 1 1 2 1 3 1 -999 0 0 -999", "disjoint partition"},
        {"1 nan 0 1 1 1 1 -999 0 -999", "Invalid or unsupported GTP integer"},
        {"1 0 0 1 1 1 1 -999 -1 -999", "Negative GTP prize"},
        {"1 0 0 1 1 1 1 -999 0 -999 extra", "trailing data"},
        {"1 0 0 1 1 1 1 -998 0 -999", "first sentinel"},
        {"1 0 0 1 1 1 1 -999 0 -998", "final sentinel"},
        {"2 0 0 1 1 2 1 1 1 2 2 1 -999 0 0 -999", "center is outside"}
    };
    for (const auto &[raw, reason] : invalid) {
        bool rejected = false;
        try { std::istringstream in(raw); (void)Instance::read(in); }
        catch (const std::exception &e) { rejected = std::string(e.what()).find(reason) != std::string::npos; }
        expect(rejected, "malformed GTP was not rejected for: " + reason);
    }
    const auto distances = fixture({{0, 0}, {3, 4}, {1, 1}, {2, 2}, {-3, -4},
                                    {1000000000, 1000000000}, {-1000000000, -1000000000}},
                                   {{0, 1, 2, 3, 4, 5, 6}}, {0, 0, 0, 0, 0, 0, 0});
    expect(distances.distance(0, 1) == 5 && distances.distance(0, 2) == 1 &&
           distances.distance(0, 3) == 3 && distances.distance(0, 4) == 5 &&
           distances.distance(5, 6) == 2828427125LL, "integer EUC_2D rounding or range");
    const std::vector<Instance> fixtures{
        fixture({{0, 0}, {3, 4}}, {{0, 1}}, {3, 8}),
        parsed,
        fixture({{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
                {{0, 1}, {2, 3}, {4, 5}}, {0, 0, 0, 0, 0, 0}),
        fixture({{0, 0}, {10, 0}, {10, 1}, {20, 0}}, {{0}, {1, 2}, {3}}, {0, 0, 8, 0}),
        fixture({{-8, 2}, {-3, 3}, {2, -2}, {5, 6}, {10, -1}, {3, 0}, {7, 5}, {9, 9}},
                {{0}, {1, 2, 3}, {4, 5}, {6, 7}}, {1, 2, 8, 0, 4, 11, 1, 15}),
        fixture({{-9, 0}, {10, 0}, {-7, 5}, {3, 1}, {-1, -8}, {7, 5}, {0, 0}, {12, 4}},
                {{0, 1}, {2, 3}, {4, 5}, {6, 7}}, {35, 0, 4, 25, 31, 0, 2, 28}),
        fixture({{0, 0}, {1, 2}, {6, 0}, {-3, 3}, {8, 2}, {2, 5}, {4, -2}, {9, -3}, {0, 7}},
                {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}}, {0, 3, 2, 4, 1, 6, 2, 7, 0})
    };
    std::size_t selections = 0, iterations = 0;
    bool disconnected = false, increment = false;
    for (const auto &d : fixtures) {
        const EdgeReduction reduced(d);
        std::vector<std::vector<int>> selectionsHere;
        std::vector<int> r(d.groups.size());
        std::function<void(std::size_t)> enumerate = [&](std::size_t g) {
            if (g == d.groups.size()) { selectionsHere.push_back(r); return; }
            for (int v : d.groups[g]) { r[g] = v; enumerate(g + 1); }
        };
        enumerate(0);
        for (const EdgeReduction *graph : {static_cast<const EdgeReduction *>(nullptr), &reduced}) {
            Search evaluator(d, graph, config(1, true));
            for (const auto &selection : selectionsHere) {
                verifySolution(d, evaluator.evaluate(selection), graph); ++selections;
            }
            for (const auto seed : {1ULL, 19ULL}) {
                std::vector<std::pair<std::vector<int>, std::int64_t>> acceleratedTrace;
                for (const bool partial : {true, false}) {
                    std::vector<std::pair<std::vector<int>, std::int64_t>> trace;
                    std::vector<std::uint64_t> expectedWeights(d.coordinates.size(), 2);
                    auto best = std::numeric_limits<std::int64_t>::max();
                    Search search(d, graph, config(seed, partial));
                    const auto result = search.run([&](std::size_t n, const Solution &s, double) {
                        expect(n == trace.size() + 1, "completed iteration numbering");
                        verifySolution(d, s, graph);
                        for (std::size_t g = 0; g < d.groups.size(); ++g)
                            for (int v : d.groups[g]) {
                                auto neighbor = s.representatives; neighbor[g] = v;
                                expect(oracleMst(d, neighbor) >= s.objective, "unfinished one-group local search");
                            }
                        if (s.objective < best) {
                            best = s.objective;
                            std::fill(expectedWeights.begin(), expectedWeights.end(), 2);
                            for (std::size_t g = 0; g < d.groups.size(); ++g)
                                expectedWeights[s.representatives[g]] = d.groups[g].size();
                        } else {
                            for (int v : s.representatives) expectedWeights[v] += 2;
                            increment = true;
                        }
                        trace.emplace_back(s.representatives, s.objective); ++iterations;
                    });
                    expect(result.complete && result.completedIterations == 12 &&
                           result.memoryWeightsDoubled == expectedWeights, "iteration count or adaptive memory");
                    const auto firstBest = std::min_element(trace.begin(), trace.end(),
                        [](const auto &a, const auto &b) { return a.second < b.second; });
                    expect(result.best.representatives == firstBest->first && result.best.objective == firstBest->second,
                           "best must retain first strictly best completed candidate");
                    if (partial) acceleratedTrace = trace;
                    else expect(trace == acceleratedTrace, "partial/full MST changed fixed-seed trajectory");
                    disconnected = disconnected || result.statistics.disconnectedPartialForests > 0;
                }
            }
        }
    }
    expect(disconnected && increment, "disconnected partial forest and memory reinforcement not exercised");
    std::cout << "{\"self_test\":\"passed\",\"fixtures\":" << fixtures.size()
              << ",\"fixed_selections\":" << selections << ",\"completed_local_optima\":" << iterations
              << ",\"malformed_inputs_rejected\":" << invalid.size()
              << ",\"oracle\":\"exhaustive_Pruefer_trees\"}\n";
}
} // namespace checks

int run(int argc, char **argv) {
    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--help")) { help(); return 0; }
    if (argc == 2 && std::string(argv[1]) == "--example") { std::cout << example; return 0; }
    if (argc == 2 && std::string(argv[1]) == "--self-test") { checks::run(); return 0; }
    Config c; c.cpuStartedAt = timing::processCpuSeconds(); c.cpuLimit = 600;
    bool reduce = true, trace = false;
    std::string path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            return argv[i];
        };
        if (arg == "--seed") c.seed = unsignedInteger(value());
        else if (arg == "--iterations") {
            const auto n = unsignedInteger(value());
            if (!n || n > 1000000000ULL) throw std::runtime_error("Iterations must be in [1,1000000000]");
            c.iterations = static_cast<std::size_t>(n);
        } else if (arg == "--cpu-limit") c.cpuLimit = positiveReal(value());
        else if (arg == "--no-reduction") reduce = false;
        else if (arg == "--full-mst") c.partialTree = false;
        else if (arg == "--trace") trace = true;
        else if (arg.empty() || arg[0] == '-') throw std::runtime_error("Unknown option: " + arg);
        else if (path.empty()) path = arg;
        else throw std::runtime_error("Expected exactly one GTP input");
    }
    if (path.empty()) throw std::runtime_error("Missing GTP input");
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read GTP instance: " + path);
    const double deadline = c.cpuStartedAt + c.cpuLimit;
    const auto data = Instance::read(input, deadline);
    // No reduction certificates, solver, miner, elite or frequency exports.
    // The reduction changes only retained edges; the complete input stays intact.
    std::optional<EdgeReduction> reduced;
    if (reduce) reduced.emplace(data, deadline);
    const auto *graph = reduced ? &*reduced : nullptr;
    Search search(data, graph, c);
    const auto result = search.run([&](std::size_t n, const Solution &s, double cpu) {
        if (!trace) return;
        std::cout << "{\"type\":\"iteration\",\"iteration\":" << n << ",\"cpu_seconds\":" << cpu << ",\"solution\":";
        solutionJson(std::cout, s); std::cout << "}\n";
    });
    std::cout << std::boolalpha << std::setprecision(12)
              << "{\"type\":\"result\",\"algorithm\":\"GRASP-MA\",\"seed\":" << c.seed
              << ",\"vertices\":" << data.coordinates.size() << ",\"groups\":" << data.groups.size()
              << ",\"requested_iterations\":" << result.requestedIterations
              << ",\"completed_iterations\":" << result.completedIterations
              << ",\"discarded_partial_iterations\":" << result.discardedPartialIterations
              << ",\"complete\":" << result.complete << ",\"termination\":\"" << result.terminationReason
              << "\",\"cpu_seconds\":" << timing::processCpuSeconds() - c.cpuStartedAt
              << ",\"reduction\":" << reduce << ",\"removed_edges\":" << (reduced ? reduced->removedEdges : 0)
              << ",\"partial_forest\":" << c.partialTree << ",\"best\":";
    if (result.hasBest) solutionJson(std::cout, result.best); else std::cout << "null";
    std::cout << ",\"memory_weights_doubled\":"; array(std::cout, result.memoryWeightsDoubled);
    std::cout << "}\n";
    return result.complete ? 0 : 2;
}
} // namespace app

int main(int argc, char **argv) {
    try { return app::run(argc, argv); }
    catch (const timing::CpuLimit &) {
        std::cout << "{\"type\":\"result\",\"complete\":false,\"termination\":\"cpu_limit_before_search\",\"best\":null}\n";
        return 2;
    }
    catch (const std::exception &e) { std::cerr << "error: " << e.what() << '\n'; return 1; }
}
