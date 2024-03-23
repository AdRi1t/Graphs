#include <vector>
#include <set>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "graphs.hpp"

long Vertex::getID() const { return this->vertex_id; }

long Vertex::getDegree() const { return this->degree; }

void Vertex::addDegree(int d) { this->degree += d; }

const std::unordered_map<long, Vertex>& Graph::getVerticesMap() const { return this->vertices; }

std::vector<Vertex> Graph::getVerticesList() const
{
    std::vector<Vertex> vertices_list;
    vertices_list.reserve(this->getVerticesCount());
    for (const auto& v : this->getVerticesMap()) {
        vertices_list.push_back(v.second);
    }
    return vertices_list;
}

const std::unordered_map<long, std::vector<long>>& Graph::getAdjacencyList() const
{
    return this->adjacency_list;
}

std::vector<std::vector<long>> Graph::getAdjacencyMatrix() const
{
    std::vector<std::vector<long>> adjacency_matrix;
    adjacency_matrix.reserve(this->getVerticesCount());
    adjacency_matrix.resize(this->getVerticesCount());

    for (const auto& v : this->getVerticesMap()) {
        std::vector<long> neighbors;
        neighbors.reserve(this->getNeighborsList(v.second.getID()).size());
        for (const long& neighbor_id : this->getNeighborsList(v.second.getID())) {
            neighbors.push_back(neighbor_id);
        }
        adjacency_matrix.at(v.second.getID()) = neighbors;
    }
    return adjacency_matrix;
}

const std::vector<long>& Graph::getNeighborsList(const long id) const
{
    return this->adjacency_list.at(id);
}

std::pair<std::vector<long>, std::vector<std::vector<long>>> Graph::getSubAdjacencyList(
    const std::set<long>& vertices_set) const
{
    std::pair<std::vector<long>, std::vector<std::vector<long>>> sub_adjacency_list;
    sub_adjacency_list.first.reserve(vertices_set.size());
    sub_adjacency_list.second.reserve(vertices_set.size());
    for (const auto& vertex_id : vertices_set) {
        sub_adjacency_list.first.push_back(vertex_id);
        sub_adjacency_list.second.push_back(std::vector<long>());
        for (const auto& neighbor_id : this->getNeighborsList(vertex_id)) {
            if (vertices_set.find(neighbor_id) != vertices_set.end()) {
                sub_adjacency_list.second.back().push_back(neighbor_id);
            }
        }
    }
    return sub_adjacency_list;
}

Graph Graph::getComplementGraph()
{
    Graph complement_graph;
    for (const auto& vertex : this->getVerticesMap()) {
        complement_graph.addVertex(vertex.second.getID());
    }
    for (const auto& vertex_1 : this->getVerticesMap()) {
        for (const auto& vertex_2 : this->getVerticesMap()) {
            if (vertex_1.second.getID() != vertex_2.second.getID() &&
                this->areNeighbors(vertex_1.second.getID(), vertex_2.second.getID()) == false) {
                complement_graph.addEdge(vertex_1.second.getID(), vertex_2.second.getID());
            }
        }
    }
    return complement_graph;
}

long Graph::getVerticesCount() const { return this->num_vertices; }

long Graph::getEdgesCount() const { return this->num_edges; }

long Graph::getMaxDegree() const
{
    long max_degree = 0;
    for (const auto& vertex : vertices) {
        max_degree = std::max(vertex.second.getDegree(), max_degree);
    }
    return max_degree;
}

Vertex Graph::getMaxDegreeVertex()
{
    long max_degree = 0;
    long max_id = 0;
    for (const auto& vertex : vertices) {
        if (vertex.second.getDegree() > max_degree) {
            max_degree = vertex.second.getDegree();
            max_id = vertex.second.getID();
        }
    }
    if (max_id == 0) {
        return this->vertices.begin()->second;
    } else {
        return this->vertices.at(max_id);
    }
}

void Graph::addEdge(const long id_A, const long id_B)
{
    if (this->vertices.find(id_A) == this->vertices.end() ||
        this->vertices.find(id_B) == this->vertices.end()) {
        std::cerr << "Adding an invalid edge : " + std::to_string(id_A) + " - " +
                         std::to_string(id_B) + "\n";
        return;
    }
    if (areNeighbors(id_A, id_B)) {
        //std::cerr << "Adding an already existing edge : " + std::to_string(Id_A) +" - "+ std::to_string(Id_B) + "\n";
        return;
    }
    this->adjacency_list.at(id_A).push_back(id_B);
    this->adjacency_list.at(id_B).push_back(id_A);
    this->vertices.at(id_A).addDegree(1);
    this->vertices.at(id_B).addDegree(1);
    this->num_edges += 1;
}

bool Graph::areNeighbors(long id_A, long id_B) const
{
    for (const long& neighbor_id : this->adjacency_list.at(id_A)) {
        if (neighbor_id == id_B) {
            return true;
        }
    }
    return false;
}

bool Graph::hasVertex(long id) const { return (this->vertices.find(id) != this->vertices.end()); }

void Graph::removeVertex(long id)
{
    auto vertex_it = this->vertices.find(id);
    if (vertex_it == this->vertices.end()) {
        // std::cerr << "Deleting an invalid vertex " + std::to_string(id) + "\n";
        return;
    }
    // Remove all edge where the vertex is involved
    std::vector<long> n_list = this->getNeighborsList(id);
    for (const auto n : n_list) {
        this->removeEdge(n, id);
    }
    const auto edges_it = this->adjacency_list.find(id);
    this->adjacency_list.erase(edges_it);
    this->vertices.erase(vertex_it);
    this->num_vertices -= 1;
    return;
}

void Graph::removeEdge(long id_A, long id_B)
{
    if (!areNeighbors(id_A, id_B)) {
        std::cerr << "Deleting an invalid edge : " + std::to_string(id_A) + " - " +
                         std::to_string(id_B) + "\n";
        return;
    }
    this->adjacency_list.at(id_A).erase(std::remove(this->adjacency_list.at(id_A).begin(),
                                                    this->adjacency_list.at(id_A).end(), id_B),
                                        this->adjacency_list.at(id_A).end());
    this->adjacency_list.at(id_B).erase(std::remove(this->adjacency_list.at(id_B).begin(),
                                                    this->adjacency_list.at(id_B).end(), id_A),
                                        this->adjacency_list.at(id_B).end());
    this->vertices.at(id_A).addDegree(-1);
    this->vertices.at(id_B).addDegree(-1);
    this->num_edges -= 1;
    return;
}

void Graph::addVertex(const long id)
{
    auto vertex_it = this->vertices.find(id);
    if (vertex_it != this->vertices.end()) {
        // std::cerr << "Adding an already existing vertex : " + std::to_string(id) + "\n";
        return;
    }
    this->vertices.emplace(id, id);
    this->adjacency_list.emplace(id, std::vector<long>());
    this->num_vertices += 1;
    return;
}

// Random Geometric Graph
void Graph::initRGG(long nb_vertices, double distance_max)
{
    auto distance = [](std::pair<float, float> a, std::pair<float, float> b) {
        return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
    };

    this->num_vertices = 0;
    this->num_edges = 0;
    this->vertices.clear();
    this->adjacency_list.clear();
    this->vertices.reserve(nb_vertices);
    this->adjacency_list.reserve(nb_vertices);

    std::vector<std::pair<float, float>> vertices_pos;
    this->vertices.clear();
    vertices_pos.reserve(nb_vertices);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (long i = 0; i < nb_vertices; i++) {
        vertices_pos.emplace_back(dis(gen), dis(gen));
        this->addVertex(i);
    }

    for (long i = 0; i < nb_vertices; i++) {
        for (long j = i; j < nb_vertices; j++) {
            if (i != j && distance(vertices_pos[i], vertices_pos[j]) < distance_max) {
                this->addEdge(i, j);
            }
        }
    }
}

// Barabasi-Albert model
void Graph::initBA(const int m, const long nb_vertices)
{
    long m_0 = m + 1;
    initComplete((long)m_0);
    std::random_device rd;
    std::mt19937 gen1(rd());
    std::mt19937 gen2(rd());
    std::uniform_real_distribution<double> randR(0, 1);

    double k_j = double(2 * m_0 * (m_0 - 1)) / 2.0;
    for (long i = m_0; i < nb_vertices; i++) {
        this->addVertex(i);
        std::uniform_int_distribution<long> randI(0, i - 1);

        while (this->vertices[i].getDegree() < m) {
            long candidate = randI(gen2);
            double k_i = (double)this->vertices[candidate].getDegree();
            double p_i = k_i / k_j;
            if (randR(gen1) <= p_i && candidate != i) {
                this->addEdge(i, candidate);
                k_j += 2;
            }
        }
    }
}

// Complete graph
void Graph::initComplete(const long nb_vertices)
{
    this->num_vertices = 0;
    this->num_edges = 0;
    this->vertices.clear();
    this->adjacency_list.clear();
    this->vertices.reserve(nb_vertices);
    this->adjacency_list.reserve(nb_vertices);
    for (long i = 0; i < nb_vertices; i++) {
        this->addVertex(i);
    }
    for (long i = 0; i < nb_vertices; i++) {
        for (long j = 0; j < nb_vertices; j++) {
            if (j != i) {
                this->addEdge(i, j);
            }
        }
    }
}

// Read graph from a .mtx file
void Graph::initFromMtxFile(const std::string file_path)
{
    std::ifstream file(file_path);
    std::string line;
    long vertex_count = 0;
    if (!file.is_open()) {
        std::cerr << "Error opening file " << file_path << "\n";
        return;
    }
    while (std::getline(file, line)) {
        if (line[0] != '%') {
            std::istringstream iss(line);
            if (!(iss >> vertex_count)) {
                std::cerr << "Failed to read graphs dimensions\n";
                return;
            }
            break;
        }
    }
    long id_A, id_B;
    while (file >> id_A >> id_B) {
        id_A -= 1;
        id_B -= 1;
        this->addVertex(id_A);
        this->addVertex(id_B);
        this->addEdge(id_A, id_B);
    }
    file.close();
}

// Init graph with a given number of vertices and a probability
void Graph::initWithProbability(const long nb_vertices, const double edge_probability)
{
    this->num_vertices = 0;
    this->num_edges = 0;
    this->vertices.clear();
    this->adjacency_list.clear();
    this->vertices.reserve(nb_vertices);
    this->adjacency_list.reserve(nb_vertices);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (long i = 0; i < nb_vertices; i++) {
        this->addVertex(i);
    }

    for (long i = 0; i < nb_vertices; i++) {
        for (long j = i; j < nb_vertices; j++) {
            if (j != i && dis(gen) < edge_probability) {
                this->addEdge(i, j);
            }
        }
    }
}

// Init graph with a given number of vertices and a fixed number of edges
void Graph::initFixedEdges(const long nb_vertices, const long nb_edges)
{
    assert(nb_edges < (nb_vertices * (nb_vertices - 1) / 2));
    this->num_vertices = 0;
    this->num_edges = 0;
    this->vertices.clear();
    this->adjacency_list.clear();
    this->vertices.reserve(nb_vertices);
    this->adjacency_list.reserve(nb_vertices);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<long> dis(0, nb_vertices - 1);

    for (long i = 0; i < nb_vertices; i++) {
        this->addVertex(i);
    }
    while (this->num_edges < nb_edges) {
        long i = dis(gen);
        long j = dis(gen);
        if (i != j && !this->areNeighbors(i, j)) {
            this->addEdge(i, j);
        }
    }
}

void Graph::printMetaData() const
{
    using namespace std;
    cout << "Graph ID           : " << this->graph_id << "\n";
    cout << "Number of vertices : " << this->num_vertices << "\n";
    cout << "Number of edges    : " << this->num_edges << "\n";
    cout << "Maximum degree     : " << this->getMaxDegree() << "\n";
    cout << "Density            : "
         << (double)(2 * this->num_edges) / (double)(this->num_vertices * (this->num_vertices - 1))
         << "\n";
    cout << "Average degree     : " << (double)(2 * this->num_edges) / (double)this->num_vertices << "\n\n";
}

void Graph::printVerticesList() const
{
    std::cout << std::setw(8) << "ID" << std::setw(8) << "Degree"
              << "\n";
    for (auto vertex : this->getVerticesMap()) {
        std::cout << std::setw(8) << vertex.second.getID() << std::setw(8)
                  << vertex.second.getDegree() << "\n";
    }
}

void Graph::printEdgesList() const
{
    std::cout << "Nb Edges : " << this->getEdgesCount() << "\n";
    for (auto vertex : this->getVerticesMap()) {
        std::cout << std::setw(8) << vertex.second.getID() << " --- ";
        for (auto neighbor_id : this->getNeighborsList(vertex.second.getID())) {
            std::cout << neighbor_id << ", ";
        }
        std::cout << "\n";
    }
}
