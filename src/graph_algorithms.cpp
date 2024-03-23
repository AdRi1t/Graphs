#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <vector>
#include <set>
#include <algorithm>
#include <climits>
#include <stack>
#include <memory>
#include "graphs.hpp"
#include "graph_algorithms.hpp"

/*
D: array of lists of vertices with the same degree
    for all vertices that are not in L.
K: the degeneracy
link: https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)
*/
long computeDegeneracy(const Graph& graph)
{
    long k = 0;
    long n = graph.getVerticesCount();
    std::list<long> L;
    std::vector<std::list<long>> D(graph.getMaxDegree() + 1);
    std::vector<long> tmp_degree(n);

    // Initialize D.
    for (auto vertex : graph.getVerticesMap()) {
        tmp_degree[vertex.second.getID()] = vertex.second.getDegree();
        D[vertex.second.getDegree()].emplace_back(vertex.second.getID());
    }
    for (long i = 0; i < n; i++) {
        long j = 0;
        // Find j such that D[j] is non-empty.
        for (auto l : D) {
            if (!l.empty()) {
                break;
            }
            j += 1;
        }
        k = std::max(k, j);
        // Select a vertex v from D[j], add v to the front of L, and remove it from D[j].
        long v = D[j].back();
        L.push_front(v);
        D[j].pop_back();

        // For each neighbor w of v that is not already in L, subtract one from dw and move w to the cell in D corresponding to the new value of dw.
        for (const long& neighbor_id : graph.getNeighborsList(v)) {
            bool found = false;
            for (const long& id : L) {
                if (id == neighbor_id) {
                    found = true;
                    break;
                }
            }
            if (found == false) {
                D[tmp_degree[neighbor_id]].remove(neighbor_id);
                tmp_degree[neighbor_id] -= 1;
                D[tmp_degree[neighbor_id]].push_back(neighbor_id);
            }
        }
    }
    L.clear();
    D.clear();
    tmp_degree.clear();
    return k;
}

// Bron-Kerbosch algorithm for finding all maximal cliques
void bronKerbosch(const Graph& graph, std::vector<std::set<long>>& cliques)
{
    std::set<long> R;
    std::set<long> P;
    std::set<long> X;
    for (auto& v : graph.getVerticesMap()) {
        P.insert(v.second.getID());
    }
    bronKerboschRoutine(R, P, X, graph, cliques);
}

void bronKerboschRoutine(std::set<long>& R, std::set<long>& P, std::set<long>& X,
                         const Graph& graph, std::vector<std::set<long>>& cliques)
{
    using namespace std;
    if (P.empty() && X.empty()) {
        // R is a maximal clique
        if (R.size() > 2) {
            cliques.push_back(R);
        }
        return;
    }
    // Choose a pivot vertex u from P U X
    long u = *P.begin();
    if (!X.empty()) {
        u = *X.begin();
    }
    set<long> neighbor_u;
    for (const long& neighbor : graph.getNeighborsList(u)) {
        neighbor_u.insert(neighbor);
    }
    // P \ N(u)
    set<long> P_without_neighbors_u;
    set_difference(P.begin(), P.end(), neighbor_u.begin(), neighbor_u.end(),
                   inserter(P_without_neighbors_u, P_without_neighbors_u.begin()));

    for (auto v = P_without_neighbors_u.begin(); v != P_without_neighbors_u.end();) {
        set<long> R1 = R; // R u {v}
        R1.insert(*v);

        set<long> P1; // P ∩ n(v)
        set<long> X1; // X ∩ n(v)
        for (const auto& neighbor : graph.getNeighborsList(*v)) {
            if (P.find(neighbor) != P.end()) {
                P1.insert(neighbor);
            }
            if (X.find(neighbor) != X.end()) {
                X1.insert(neighbor);
            }
        }
        bronKerboschRoutine(R1, P1, X1, graph, cliques);
        X.insert(*v);
        P.erase(*v);
        v = P_without_neighbors_u.erase(v);
    }
}

// Tomita algorithm for finding all maximal cliques
void tomita(const Graph& graph, std::vector<std::set<long>>& cliques)
{
    std::set<long> Q;
    std::set<long> SUBG;
    std::set<long> CAND;
    std::unordered_map<long, std::set<long>> adjacencyList;
    for (const auto& v : graph.getVerticesMap()) {
        const long id = v.second.getID();
        SUBG.insert(id);
        CAND.insert(id);
        std::vector<long> neighborsVector = graph.getNeighborsList(id);
        adjacencyList[id] = std::set<long>(neighborsVector.begin(), neighborsVector.end());
    }
    tomitaRoutine(Q, SUBG, CAND, adjacencyList, cliques);
}

void tomitaRoutine(std::set<long>& Q, std::set<long>& SUBG, std::set<long>& CAND,
                   const std::unordered_map<long, std::set<long>>& adjacencyList,
                   std::vector<std::set<long>>& cliques)
{
    using namespace std;
    if (SUBG.empty()) {
        if (Q.size() > 2) {
            cliques.push_back(Q);
        }
        return;
    }

    // Choose pivot u ∈ SUBG that maximizes |CAND ∩ v(u)|
    long u = *SUBG.begin();
    set<long> max_intersection;
    for (const long& v : SUBG) {
        set<long> neighbors = adjacencyList.at(v);
        set<long> intersection;
        for (const long& neighbor_id : neighbors) {
            if (CAND.find(neighbor_id) != CAND.end()) {
                intersection.insert(neighbor_id);
            }
        }
        if (intersection.size() > max_intersection.size()) {
            u = v;
            max_intersection = intersection;
        }
    }
    set<long> neighbors_u = adjacencyList.at(u);
    set<long> CAND_minus_neighbors_u;
    for (const long& cand_id : CAND) {
        if (neighbors_u.find(cand_id) == neighbors_u.end()) {
            CAND_minus_neighbors_u.insert(cand_id);
        }
    }

    while (!CAND_minus_neighbors_u.empty()) {
        // q ∈ CAND − v(u)
        long q = *CAND_minus_neighbors_u.begin();
        // Q + q
        Q.insert(q);
        // SUBG ∩ v(q)
        set<long> neighbors_q = adjacencyList.at(q);
        set<long> SUBG_intersection;
        set<long> CAND_intersection;
        for (const long& neighbor_id : neighbors_q) {
            if (SUBG.find(neighbor_id) != SUBG.end()) {
                SUBG_intersection.insert(neighbor_id);
            }
            if (CAND.find(neighbor_id) != CAND.end()) {
                CAND_intersection.insert(neighbor_id);
            }
        }

        // Tomita(Q, SUBG ∩ v(q), CAND ∩ v(q), E)
        tomitaRoutine(Q, SUBG_intersection, CAND_intersection, adjacencyList, cliques);
        // Q - q
        Q.erase(q);
        // CAND - q
        CAND.erase(q);
        CAND_minus_neighbors_u.erase(q);
    }
}

// Compute the Gi graph
Graph getGiGraph(const Graph& graph, const long vertex_id)
{
    Graph GiGraph;
    GiGraph.addVertex(vertex_id);
    std::set<long> neighbors;
    std::set<long> neighbors2;

    for (const long& neighbor_id : graph.getNeighborsList(vertex_id)) {
        neighbors.insert(neighbor_id);
        GiGraph.addVertex(neighbor_id);
        for (const long& neighbor_id2 : graph.getNeighborsList(neighbor_id)) {
            if (neighbor_id2 != vertex_id && neighbor_id != neighbor_id2 &&
                !graph.areNeighbors(neighbor_id2, vertex_id)) {
                neighbors2.insert(neighbor_id2);
                GiGraph.addVertex(neighbor_id2);
            }
        }
    }
    for (const long& neighbor_id : neighbors) {
        for (const long& neighbor_id2 : neighbors2) {
            if (neighbor_id != vertex_id && neighbor_id2 != vertex_id &&
                !graph.areNeighbors(neighbor_id, neighbor_id2)) {
                // xy ∈ Gi , if x ∈ Ni (vi ), y ∈ Ni2 (vi ) and xy !∈ E
                GiGraph.addEdge(neighbor_id, neighbor_id2);
            }
        }
    }
    for (const auto& vertex_id_1 : neighbors) {
        for (const long& vertex_id_2 : neighbors) {
            if (vertex_id_1 != vertex_id && vertex_id_2 != vertex_id &&
                graph.areNeighbors(vertex_id_1, vertex_id_2)) {
                // xy ∈ Gi , if x ∈ Ni (vi ), y ∈ Ni (vi ) and xy ∈ E
                GiGraph.addEdge(vertex_id_1, vertex_id_2);
            }
        }
    }
    for (const auto& vertex_id_1 : neighbors2) {
        for (const long& vertex_id_2 : neighbors2) {
            if (vertex_id_1 != vertex_id && vertex_id_2 != vertex_id &&
                graph.areNeighbors(vertex_id_1, vertex_id_2)) {
                // xy ∈ Gi , if x ∈ Ni2 (vi ), y ∈ Ni2 (vi ) and xy ∈ E
                GiGraph.addEdge(vertex_id_1, vertex_id_2);
            }
        }
    }
    neighbors.clear();
    neighbors2.clear();
    return GiGraph;
}

std::vector<Graph> compute_Gik(const Graph& graph, long vertex_id, std::vector<long> degeneracy_order)
{
    std::vector<Graph> Gik;
    std::set<long> neighbors;
    std::set<long> neighbors2;
    for (const long& neighbor_id : graph.getNeighborsList(vertex_id)) {
        neighbors.insert(neighbor_id);
        for (const long& neighbor_id2 : graph.getNeighborsList(neighbor_id)) {
            if (neighbor_id2 != vertex_id && neighbor_id != neighbor_id2 &&
                !graph.areNeighbors(neighbor_id2, vertex_id)) {
                neighbors2.insert(neighbor_id2);
            }
        }
    }

    for (const long& x : neighbors)
    {
        Graph G_ik = getGiGraph(graph, vertex_id);
        for (const long& n_2 : neighbors2)
        {
            bool found = false;
            for(const long& neighbor_x : graph.getNeighborsList(x))
            {
                if(neighbor_x == n_2)
                {
                    found = true;
                    break;
                }
            }
            if (found == false)
            {
                G_ik.removeVertex(n_2);
            }
        }
        Gik.push_back(G_ik);
    }
    neighbors.clear();
    neighbors2.clear();
    return Gik;
}

std::vector<Graph> computeBiclique(const Graph& graph)
{
    std::vector<Graph> bicliques;
    std::unordered_set<std::string> T;
    for (const auto& vertex : graph.getVerticesMap()) {
        Graph GiGraph = getGiGraph(graph, vertex.second.getID()).getComplementGraph();
        std::vector<std::set<long>> cliques;
        bronKerbosch(GiGraph, cliques);
        for (const auto& I : cliques) {
            std::pair<std::vector<long>, std::vector<std::vector<long>>> possible_biclique =
                graph.getSubAdjacencyList(I);
            if (possible_biclique.second.at(0).size() > 2) {
                std::string ordering = toString(possible_biclique.first);
                if (T.find(ordering) == T.end()) {
                    T.emplace(ordering);
                    bicliques.emplace_back(graph.getSubGraph(I));
                }
            }
        }
    }
    return bicliques;
}

std::vector<Graph> computeBiclique_V2(const Graph& graph)
{
    std::vector<Graph> bicliques;
    std::unordered_set<std::string> T;
    for (const auto& vertex : graph.getVerticesMap()) {
        std::vector<Graph> Gik_graph = compute_Gik(graph, vertex.second.getID(), computeDegeneracyOrdering(graph));
        for (auto& Gik : Gik_graph)
        {   
            const Graph graph_test = Gik.getComplementGraph();
            std::vector<std::set<long>> cliques;
            bronKerbosch(graph_test, cliques);
            for (const auto& I : cliques) {
                std::pair<std::vector<long>, std::vector<std::vector<long>>> possible_biclique =
                    graph.getSubAdjacencyList(I);
                if (possible_biclique.second.at(0).size() > 2) {
                    std::string ordering = toString(possible_biclique.first);
                    if (T.find(ordering) == T.end()) {
                        T.emplace(ordering);
                        bicliques.emplace_back(graph.getSubGraph(I));
                    }
                }
            }
        }
    }
    return bicliques;
}

Graph createSpanningTree(const Graph& graph, long start_id)
{
    Graph spanning_tree;
    std::unordered_map<long, bool> visited;
    std::deque<long> to_visit;
    for (const auto& v : graph.getVerticesMap()) {
        visited[v.second.getID()] = false;
    }
    to_visit.push_back(start_id);
    visited.at(start_id) = true;
    while (!to_visit.empty()) {
        long v = to_visit.front();
        to_visit.pop_front();
        spanning_tree.addVertex(v);
        for (const long& neighbor_id : graph.getNeighborsList(v)) {
            if (!visited[neighbor_id]) {
                visited[neighbor_id] = true;
                to_visit.push_back(neighbor_id);
                spanning_tree.addVertex(neighbor_id);
                spanning_tree.addEdge(v, neighbor_id);
            }
        }
    }
    return spanning_tree;
}

std::vector<long> computeDegeneracyOrdering(const Graph& graph)
{
    std::list<long> L = computeDegeneracyOrder(graph);
    std::vector<long> order(L.size());
    size_t i = 0;
    for (const auto& l : L) {
        order[l] = i;
        i += 1;
    }
    return order;
}

// Degeneracy ordering of a graph
std::list<long> computeDegeneracyOrder(const Graph& graph)
{
    long k = 0;
    long n = graph.getVerticesCount();
    std::list<long> L;
    std::vector<std::list<long>> D(graph.getMaxDegree() + 1);
    std::vector<long> tmp_degree(n);

    for (auto vertex : graph.getVerticesMap()) {
        tmp_degree[vertex.second.getID()] = vertex.second.getDegree();
        D[vertex.second.getDegree()].emplace_back(vertex.second.getID());
    }
    for (long i = 0; i < n; i++) {
        long j = 0;
        for (auto l : D) {
            if (!l.empty()) {
                break;
            }
            j += 1;
        }
        k = std::max(k, j);
        long v = D[j].back();
        L.push_front(v);
        D[j].pop_back();

        for (const long& neighbor_id : graph.getNeighborsList(v)) {
            bool found = false;
            for (const long& id : L) {
                if (id == neighbor_id) {
                    found = true;
                    break;
                }
            }
            if (found == false) {
                D[tmp_degree[neighbor_id]].remove(neighbor_id);
                tmp_degree[neighbor_id] -= 1;
                D[tmp_degree[neighbor_id]].push_back(neighbor_id);
            }
        }
    }
    D.clear();
    tmp_degree.clear();
    return L;
}

std::string toString(const Graph& graph)
{
    std::vector<long> vertices;
    vertices.reserve(graph.getVerticesCount());
    for (const auto& v : graph.getVerticesMap()) {
        vertices.push_back(v.second.getID());
    }

    std::sort(vertices.begin(), vertices.end());

    std::string str;
    for (const auto& v : vertices) {
        str += std::to_string(v) + "-";
    }
    str.pop_back();
    str.push_back('#');
    return str;
}

// Transform a vector of vertices into a string
std::string toString(std::vector<long> vertices)
{
    std::sort(vertices.begin(), vertices.end());
    std::string str;
    for (const auto& v : vertices) {
        str += std::to_string(v) + '-';
    }
    str.pop_back();
    str.push_back('#');
    str.shrink_to_fit();
    // std::cout << str << "\n";
    return str;
}

// BFS returns the distances from a vertex
std::vector<std::pair<long, int>> distanceFromVertex(std::vector<std::vector<long>> graph, long id)
{
    std::vector<std::pair<long, int>> distance;
    distance.reserve(graph.size());
    for (size_t i = 0; i < graph.size(); i++) {
        distance.emplace_back(i, INT_MAX);
    }
    distance.at(id).second = 0;
    std::deque<long> to_visit;
    to_visit.push_back(id);
    while (!to_visit.empty()) {
        long v = to_visit.front();
        to_visit.pop_front();
        for (size_t i = 0; i < graph.at(v).size(); i++) {
            if (distance.at(graph.at(v).at(i)).second == INT_MAX) {
                distance.at(graph.at(v).at(i)).second = distance.at(v).second + 1;
                to_visit.push_back(graph.at(v).at(i));
            }
        }
    }
    return distance;
}

// BFS returns the distances from a vertex
std::unordered_map<long, int> bfs(const Graph& graph, long source_id)
{
    std::deque<const Vertex*> to_visit;
    std::unordered_map<long, int> distance;
    if (graph.getVerticesMap().find(source_id) == graph.getVerticesMap().end()) {
        return distance;
    }
    distance.reserve(graph.getVerticesCount());
    for (auto& v : graph.getVerticesMap()) {
        distance.emplace(v.second.getID(), INT_MAX);
    }
    distance.at(source_id) = 0;
    to_visit.push_back(&(graph.getVerticesMap().at(source_id)));
    while (!to_visit.empty()) {
        const Vertex* v = to_visit.front();
        to_visit.pop_front();
        for (auto& neighbor_id : graph.getNeighborsList(v->getID())) {
            if (distance.at(neighbor_id) == INT_MAX) {
                distance.at(neighbor_id) = distance.at(v->getID()) + 1;
                to_visit.push_back(&(graph.getVerticesMap().at(neighbor_id)));
            }
        }
    }
    return distance;
}

int computeDiameter(const Graph& graph)
{
    int diameter = 0;
    for (const auto& v : graph.getVerticesMap()) {
        std::unordered_map<long, int> distances = bfs(graph, v.second.getID());
        for (const auto& d : distances) {
            if (d.second != INT_MAX && d.second > diameter) {
                diameter = d.second;
            }
        }
    }
    return diameter;
}

int computeDiameterBfs2(const Graph& graph, long& start_id, long& end_id)
{
    int diameter = 0;
    std::vector<std::vector<long>> adjacency_matrix = graph.getAdjacencyMatrix();
    start_id = -1;
    end_id = -1;
    for (const auto& v : graph.getVerticesMap()) {
        std::vector<std::pair<long, int>> distances =
            distanceFromVertex(adjacency_matrix, v.second.getID());
        for (const auto& d : distances) {
            if (d.second != INT_MAX && d.second > diameter) {
                diameter = d.second;
                start_id = v.second.getID();
                end_id = d.first;
            }
        }
    }
    adjacency_matrix.clear();
    return diameter;
}

void dfs(const Graph& graph, long source_id)
{
    std::vector<std::pair<long, bool>> visited;
    for (auto& v : graph.getVerticesMap()) {
        visited.emplace_back(v.second.getID(), false);
    }
    std::stack<long> to_visit;
    to_visit.push(source_id);
    while (!to_visit.empty()) {
        long v = to_visit.top();
        to_visit.pop();
        if (visited.at(v).second == false) {
            visited.at(v).second = true;
            for (auto& neighbor_id : graph.getNeighborsList(v)) {
                if (visited.at(neighbor_id).second == false) {
                    to_visit.push(neighbor_id);
                }
            }
        }
    }
    visited.clear();
}

int computeConnectedComponents(const Graph& graph)
{
    int numComponents = 0;
    std::vector<bool> visited(graph.getVerticesCount(), false);

    for (const auto& v : graph.getVerticesMap()) {
        if (!visited[v.second.getID()]) {
            DFS(graph, v.second.getID(), visited);
            numComponents += 1;
        }
    }
    return numComponents;
}

void DFS(const Graph& graph, long v, std::vector<bool>& visited)
{
    visited[v] = true;
    for (const long& neighbor_id : graph.getNeighborsList(v)) {
        if (!visited[neighbor_id]) {
            DFS(graph, neighbor_id, visited);
        }
    }
}

// R = A U B
Graph unionGraphs(Graph A, Graph B)
{
    Graph R;
    // Add Vertex
    for (auto& v : A.getVerticesMap()) {
        R.addVertex(v.second.getID());
    }
    for (auto& v : B.getVerticesMap()) {
        // Don't add twice
        if (A.getVerticesMap().find(v.second.getID()) == A.getVerticesMap().end()) {
            R.addVertex(v.second.getID());
        }
    }
    // Add edge
    for (auto& v : R.getVerticesMap()) {
        // If v is in A
        if (A.getVerticesMap().find(v.second.getID()) != A.getVerticesMap().end()) {
            // For all neighbors of v that are in A
            for (auto& neighbor_A : A.getNeighborsList(v.second.getID())) {
                R.addEdge(v.second.getID(), neighbor_A);
                // std::cout << neighbor_A->getID() << " - " << v.second.getID() << "\n";
            }
        }
        if (B.getVerticesMap().find(v.second.getID()) != B.getVerticesMap().end()) {
            for (auto& neighbor_B : B.getNeighborsList(v.second.getID())) {
                R.addEdge(v.second.getID(), neighbor_B);
                // std::cout << neighbor_B->getID() << " - " << v.second.getID() << "\n";
            }
        }
    }
    return R;
}

// R = A ∩ n(v)
Graph intersectNeighbors(Graph A, const std::vector<long> neighbors)
{
    Graph R;
    // Vertex
    for (const auto& v : A.getVerticesMap()) {
        bool found = false;
        for (const long& n : neighbors) {
            if (v.second.getID() == n) {
                found = true;
                break;
            }
        }
        if (found == true) {
            R.addVertex(v.second.getID());
        }
    }
    for (const auto& v : R.getVerticesMap()) {
        // For all neighbors of v that are in A
        for (const long& neighbor_A : A.getNeighborsList(v.second.getID())) {
            if (R.getVerticesMap().find(neighbor_A) != R.getVerticesMap().end()) {
                R.addEdge(neighbor_A, v.second.getID());
                // std::cout << neighbor_A->getID() << " - " << v.second.getID() << "\n";
            }
        }
    }
    return R;
}

// R = A \ N(v)
Graph excludeNeighbors(Graph A, const std::vector<long> neighbors)
{
    Graph R(A);
    for (const long& neighbors_id : neighbors) {
        R.removeVertex(neighbors_id);
    }
    return R;
}
