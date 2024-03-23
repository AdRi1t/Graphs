#pragma once
#include <list>
#include <unordered_map>
#include <vector>
#include <set>
#include <map>

class Vertex
{
  private:
    long vertex_id;
    long degree;

  public:
    Vertex() : Vertex(-1) {}
    Vertex(long id) : vertex_id(id), degree(0) {}
    Vertex(const Vertex& vertex)
    {
        this->vertex_id = vertex.getID();
        this->degree = vertex.getDegree();
    }
    long getID() const;
    long getDegree() const;
    void addDegree(int d);
};

// Simple & unrdirected graph
class Graph
{
  private:
    long graph_id;
    long num_vertices;
    long num_edges;
    std::unordered_map<long, Vertex> vertices;
    std::unordered_map<long, std::vector<long>> adjacency_list;

  public:
    Graph()
    {
        this->graph_id = 0;
        this->num_edges = 0;
        this->num_vertices = 0;
        vertices = std::unordered_map<long, Vertex>();
        adjacency_list = std::unordered_map<long, std::vector<long>>();
    }
    ~Graph()
    {
        this->graph_id = -1;
        this->num_edges = 0;
        this->num_vertices = 0;
        this->adjacency_list.clear();
        this->vertices.clear();
    }

    Graph(const Graph& graph) :
        graph_id(graph.graph_id), num_vertices(graph.num_vertices), num_edges(graph.num_edges),
        vertices(graph.vertices), adjacency_list(graph.adjacency_list)
    {
    }

    Graph(Graph&& graph) noexcept :
        graph_id(graph.graph_id), num_vertices(graph.num_vertices), num_edges(graph.num_edges),
        vertices(std::move(graph.vertices)), adjacency_list(std::move(graph.adjacency_list))
    {
    }

    Graph& operator=(const Graph& other)
    {
        if (this != &other) {
            graph_id = other.graph_id;
            num_vertices = other.num_vertices;
            num_edges = other.num_edges;
            vertices = other.vertices;
            adjacency_list = other.adjacency_list;
        }
        return *this;
    }

    template <typename Container> Graph getSubGraph(const Container& subset_vertices) const
    {
        Graph sub_graph;
        for (const auto& vertex_id : subset_vertices) {
            sub_graph.addVertex(vertex_id);
        }
        for (const auto& vertex_id : subset_vertices) {
            for (const auto& neighbor_id : this->getNeighborsList(vertex_id)) {
                if (std::find(subset_vertices.begin(), subset_vertices.end(), neighbor_id) !=
                    subset_vertices.end()) {
                    sub_graph.addEdge(vertex_id, neighbor_id);
                }
            }
        }
        return sub_graph;
    }

    void addVertex(const long id);
    void addEdge(const long id_A, const long id_B);
    void removeVertex(long id);
    void removeEdge(long id_A, long id_B);
    bool areNeighbors(long id_A, long id_B) const;
    bool hasVertex(long id) const;

    const std::unordered_map<long, Vertex>& getVerticesMap() const;
    const std::unordered_map<long, std::vector<long>>& getAdjacencyList() const;
    const std::vector<long>& getNeighborsList(const long id) const;
    std::vector<std::vector<long>> getAdjacencyMatrix() const;
    std::vector<Vertex> getVerticesList() const;
    std::pair<std::vector<long>, std::vector<std::vector<long>>> getSubAdjacencyList(
        const std::set<long>& vertices) const;
    Graph getComplementGraph();
    Vertex getMaxDegreeVertex();
    long getVerticesCount() const;
    long getEdgesCount() const;
    long getMaxDegree() const;

    void initRGG(const long nb_vertices, const double distance_trashold);
    void initBA(const int m_0, const long nb_vertices);
    void initWithProbability(const long nb_vertices, const double edge_probability);
    void initFixedEdges(const long nb_vertices, const long nb_edges);
    void initComplete(const long nb_vertices);
    void initFromMtxFile(const std::string file_path);

    void printVerticesList() const;
    void printEdgesList() const;
    void printMetaData() const;
};
