#include <graphviz/cgraph.h>
#include <graphviz/gvc.h>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <array>
#include <algorithm>
#include "graphs.hpp"
#include "graph_algorithms.hpp"
#include "graph_viz.hpp"

GraphViz::GraphViz(Graph* graph_ptr)
{
    this->p_graph = graph_ptr;
    this->render_dir = std::string();
    setGraph(graph_ptr);
}

GraphViz::GraphViz(void)
{
    this->p_graph = nullptr;
    this->render_dir = std::string();
}

GraphViz::GraphViz(Graph* graph_ptr, std::string directory_name)
{
    this->p_graph = graph_ptr;
    setRenderDirectory(directory_name);
    setGraph(graph_ptr);
}

// Set the graph to be visualized
void GraphViz::setGraph(Graph* graph_ptr)
{
    this->p_graph = graph_ptr;
    this->graph_vizz = agopen((char*)"graph", Agundirected, 0);
    agattr(this->graph_vizz, 0, (char*)"fontname", (const char*)"Helvetica,Arial,sans-serif");
    for (auto vertex : this->p_graph->getVerticesMap()) {
        Agnode_t* v =
            agnode(this->graph_vizz, (char*)std::to_string(vertex.second.getID()).c_str(), 1);
        agsafeset(v, (char*)"shape", (const char*)"circle", "");
        agsafeset(v, (char*)"style", (const char*)"bold", "");
        for (auto neighbor_id : this->p_graph->getNeighborsList(vertex.second.getID())) {
            Agnode_t* n = agnode(this->graph_vizz, (char*)std::to_string(neighbor_id).c_str(), 1);
            agsafeset(n, (char*)"shape", (const char*)"circle", "");
            Agedge_t* e = agedge(this->graph_vizz, v, n, "", 1);
        }
    }
}

// Color the vertices according to their distance from a given vertex
void GraphViz::colorDistanceVertex(long id)
{
    std::vector<std::pair<long, int>> distance_from_vertex;
    distance_from_vertex = distanceFromVertex(this->p_graph->getAdjacencyMatrix(), id);
    std::array<std::string, 7> gradient_color = {"#000000", "#00ffff", "#00ff40", "#ffff00",
                                                 "#e6ac00", "#ff4000", "#cc0000"};

    for (auto vertex_distance : distance_from_vertex) {
        Agnode_t* v =
            agnode(this->graph_vizz, (char*)std::to_string(vertex_distance.first).c_str(), 1);
        if (vertex_distance.second > 5) {
            agsafeset(v, (char*)"color", (const char*)"#4d004d", "");
        } else {
            agsafeset(v, (char*)"color",
                      (const char*)gradient_color[vertex_distance.second].c_str(), "");
        }
    }
    distance_from_vertex.clear();
}

// Color the vertices according to their clique membership
void GraphViz::colorCliques(int max_cliques)
{
    using namespace std;
    array<string, 5> colors = {"#FF0000", "#0000FF", "#008000", "#FFFF00", "#800080"};
    vector<pair<long, long>> edges;
    vector<set<long>> cliques;
    int color_index = 0;

    bronKerbosch(*(this->p_graph), cliques);
    sort(cliques.begin(), cliques.end(),
         [](const set<long>& a, const set<long>& b) { return a.size() > b.size(); });
    cliques.resize(max_cliques);
    for (auto clique : cliques) {
        edges.clear();
        for (auto vertex_id : clique) {
            for (auto neighbor_id : clique) {
                if (vertex_id != neighbor_id) {
                    edges.push_back(make_pair(vertex_id, neighbor_id));
                }
            }
        }
        for (const auto& edge : edges) {
            Agnode_t* v = agnode(this->graph_vizz, (char*)std::to_string(edge.first).c_str(), 1);
            Agnode_t* n = agnode(this->graph_vizz, (char*)std::to_string(edge.second).c_str(), 1);
            Agedge_t* e = agedge(this->graph_vizz, v, n, "", 1);
            agsafeset(e, (char*)"color", (const char*)colors[color_index % 5].c_str(), "");
            agsafeset(v, (char*)"color", (const char*)colors[color_index % 5].c_str(), "");
            agsafeset(n, (char*)"color", (const char*)colors[color_index % 5].c_str(), "");
        }
        color_index++;
    }
}

// Color the vertices according to their biclique membership
void GraphViz::colorBicliques(int max_bicliques)
{
    using namespace std;
    vector<Graph> bicliques;
    bicliques = computeBiclique(*(this->p_graph));
    sort(bicliques.begin(), bicliques.end(),
         [](const Graph& a, const Graph& b) { return a.getEdgesCount() > b.getEdgesCount(); });
    bicliques.resize(max_bicliques);
    for (const auto& biclique : bicliques) {
        const Vertex ref_vertex = biclique.getVerticesMap().begin()->second;
        for (auto vertex : biclique.getVerticesMap()) {
            Agnode_t* v =
                agnode(this->graph_vizz, (char*)std::to_string(vertex.second.getID()).c_str(), 1);
            if (vertex.second.getID() == ref_vertex.getID() ||
                !biclique.areNeighbors(ref_vertex.getID(), vertex.second.getID())) {
                agsafeset(v, (char*)"color", (const char*)"#FF0000", "");
            } else {
                agsafeset(v, (char*)"color", (const char*)"#00FF00", "");
            }
            for (auto neighbor_id : biclique.getNeighborsList(vertex.second.getID())) {
                Agnode_t* n =
                    agnode(this->graph_vizz, (char*)std::to_string(neighbor_id).c_str(), 1);
                Agedge_t* e = agedge(this->graph_vizz, v, n, "", 1);
                agsafeset(e, (char*)"color", (const char*)"#0000FF", "");
            }
        }
    }
    bicliques.clear();
}
void GraphViz::colorSpanningTree(long start_id)
{
    Graph spanning_tree = createSpanningTree(*(this->p_graph), start_id);
    for (auto vertex : spanning_tree.getVerticesMap()) {
        Agnode_t* v =
            agnode(this->graph_vizz, (char*)std::to_string(vertex.second.getID()).c_str(), 1);
        for (auto neighbor_id : spanning_tree.getNeighborsList(vertex.second.getID())) {
            Agnode_t* n = agnode(this->graph_vizz, (char*)std::to_string(neighbor_id).c_str(), 1);
            Agedge_t* e = agedge(this->graph_vizz, v, n, "", 1);
            agsafeset(e, (char*)"color", (const char*)"#FF0000", "");
        }
    }
}

// Close the graphviz context
void GraphViz::closeGraphVizz() { agclose(this->graph_vizz); }

// Read a graph from a file in the dot format
bool GraphViz::readDotFile(std::string file_path)
{
    this->p_graph = new Graph();
    GVC_t* gvc = gvContext();
    FILE* fp;
    fp = fopen(file_path.c_str(), "r");
    if (fp == NULL) {
        return false;
    }

    Agraph_t* g = agread(fp, NULL);
    gvLayout(gvc, g, (const char*)"dot");

    Agnode_t* v;
    Agedge_t* e;

    for (v = agfstnode(g); v; v = agnxtnode(g, v)) {
        this->p_graph->addVertex(AGSEQ(v));
    }
    for (v = agfstnode(g); v; v = agnxtnode(g, v)) {
        for (e = agfstout(g, v); e; e = agnxtout(g, e)) {
            this->p_graph->addEdge(AGSEQ(aghead(e)), AGSEQ(agtail(e)));
        }
    }
    gvFreeLayout(gvc, g);
    gvFreeContext(gvc);
    agclose(g);
    fclose(fp);
    setGraph(this->p_graph);
    return true;
}

GraphViz::~GraphViz() {}

// Render the graph to a PNG file
bool GraphViz::renderToPNG(std::string file_name)
{
    auto workingDir = std::filesystem::current_path();
    auto oldDir = workingDir;
    auto resultDir = workingDir.append(this->render_dir);
    chdir(resultDir.c_str());

    GVC_t* gvc = gvContext();
    std::string title;
    title.assign("Vertices: ");
    title.append(std::to_string(p_graph->getVerticesCount()));
    title.append("\tEdges: ");
    title.append(std::to_string(p_graph->getEdgesCount()));

    agattr(this->graph_vizz, 0, (char*)"label", title.c_str());
    agattr(this->graph_vizz, 0, (char*)"K", (const char*)"0.25");
    gvLayout(gvc, this->graph_vizz, (const char*)"fdp");
    gvRenderFilename(gvc, this->graph_vizz, (const char*)"png", (file_name).c_str());

    gvFreeLayout(gvc, this->graph_vizz);
    gvFreeContext(gvc);
    chdir(oldDir.c_str());

    return true;
}

// Render the graph to a DOT file
bool GraphViz::renderToDOT(std::string file_name)
{
    auto workingDir = std::filesystem::current_path();
    auto oldDir = workingDir;
    auto resultDir = workingDir.append(this->render_dir);
    chdir(resultDir.c_str());

    GVC_t* gvc = gvContext();
    gvLayout(gvc, this->graph_vizz, (const char*)"dot");
    gvRenderFilename(gvc, this->graph_vizz, (const char*)"dot", (file_name).c_str());

    gvFreeLayout(gvc, this->graph_vizz);
    gvFreeContext(gvc);
    chdir(oldDir.c_str());

    return true;
}

// Set the directory where the graph will be rendered
void GraphViz::setRenderDirectory(std::string directory_name)
{
    auto render_path = std::filesystem::current_path().append(directory_name);
    if (std::filesystem::exists(render_path) && std::filesystem::is_directory(render_path)) {
        this->render_dir = directory_name;
    } else {
        try {
            std::filesystem::create_directory(render_path);
            this->render_dir = directory_name;
        } catch (const std::filesystem::filesystem_error& ex) {
            std::cerr << "Error creating directory: " << ex.what() << std::endl;
        }
    }
}

// Print the degree of each vertex to a file
void GraphViz::printDegree(std::string file_name)
{
    auto workingDir = std::filesystem::current_path();
    auto oldDir = workingDir;
    auto resultDir = workingDir.append(this->render_dir);
    chdir(resultDir.c_str());
    std::ofstream file(file_name);
    for (const auto& v : this->p_graph->getVerticesMap()) {
        file << v.second.getID() << "\t" << v.second.getDegree() << std::endl;
    }
    file.close();
    chdir(oldDir.c_str());
    return;
}

// Print the edges of the graph to a file
void GraphViz::printEdges(std::string file_name)
{
    auto workingDir = std::filesystem::current_path();
    auto oldDir = workingDir;
    auto resultDir = workingDir.append(this->render_dir);
    chdir(resultDir.c_str());
    std::ofstream file(file_name);
    file << "Nb Edges : " << this->p_graph->getEdgesCount() << "\n";
    for (auto vertex : this->p_graph->getVerticesMap()) {
        file << std::setw(8) << vertex.second.getID() << " --- ";
        for (auto neighbor_id : this->p_graph->getNeighborsList(vertex.second.getID())) {
            file << neighbor_id << ", ";
        }
        file << "\n";
    }
    file.close();
    chdir(oldDir.c_str());
    return;
}
