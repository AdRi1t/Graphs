#pragma once
#include <graphviz/cgraph.h>

class GraphViz
{
  private:
    Graph* p_graph;
    Agraph_t* graph_vizz;
    std::string render_dir;

  public:
    GraphViz();
    GraphViz(Graph* graph_ptr);
    GraphViz(Graph* graph_ptr, std::string directory_name);
    void setGraph(Graph* graph_ptr);
    void colorDistanceVertex(long id);
    void colorCliques(int max_cliques=3);
    void colorBicliques(int max_bicliques=1);
    void colorSpanningTree(long start_id);
    void closeGraphVizz();
    bool readDotFile(std::string file_path);
    bool renderToPNG(std::string file_name);
    bool renderToDOT(std::string file_name);
    void setRenderDirectory(std::string directory_name);
    void printDegree(std::string file_name);
    void printEdges(std::string file_name);
    ~GraphViz();
};
