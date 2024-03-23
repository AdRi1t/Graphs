#include <vector>
#include <set>
#include <memory>
#include <unordered_map>
#include "graphs.hpp"

// Computes the degeneracy of a graph
long computeDegeneracy(const Graph& graph);

// Computes the distance from each vertex to a given vertex
std::unordered_map<long, int> bfs(const Graph& graph, long source_id);
std::vector<std::pair<long, int>> distanceFromVertex(std::vector<std::vector<long>> graph, long id);

// Bron-Kerbosch algorithm
void bronKerbosch(const Graph& graph, std::vector<std::set<long>>& cliques);
void bronKerboschRoutine(std::set<long>& R, std::set<long>& P, std::set<long>& X,
             const Graph& graph, std::vector<std::set<long>>& cliques);

// Tomita algorithm
void tomita(const Graph& graph, std::vector<std::set<long>>& cliques);
void tomitaRoutine(std::set<long>& Q, std::set<long>& SUBG, std::set<long>& CAND,
           const std::unordered_map<long, std::set<long>>& adjacencyList,
           std::vector<std::set<long>>& cliques);

// Computes the G_i graph
Graph getGiGraph(const Graph& graph, const long vertex_id);

// Computes the Gik graphs
std::vector<Graph> compute_Gik(const Graph& graph, long vertex_id, std::vector<long> degeneracy_order);

// Computes the bicliques
std::vector<Graph> computeBiclique(const Graph& graph);
// Computes the bicliques improved version
std::vector<Graph> computeBiclique_V2(const Graph& graph);

// Computes the spanning tree of a graph
Graph createSpanningTree(const Graph& graph, long start_id);

// Get the degenaracy ordering of a graph 
std::list<long> computeDegeneracyOrder(const Graph& graph);
// Get the order rank for each vertex
std::vector<long> computeDegeneracyOrdering(const Graph& graph);

// Converts vertices to string
std::string toString(const Graph& graph);
std::string toString(std::vector<long> vertices);

// Computes the diameter of a graph
int computeDiameter(const Graph& graph);
int computeDiameterBfs2(const Graph& graph, long& start_id, long& end_id);

// Depth-first search
void dfs(const Graph& graph, long source_id);
void DFS(const Graph& graph, long v, std::vector<bool>& visited);

// Computes the number of connected components
int computeConnectedComponents(const Graph& graph);

// Union of two graphs
Graph unionGraphs(Graph A, Graph B);

// Intersection of a graph and the neighbors of a vertex
Graph intersectNeighbors(Graph A, const std::vector<long> neighbors);

// Returns a graph without the neighbors of a vertex
Graph excludeNeighbors(Graph A, std::vector<long> neighbors);
