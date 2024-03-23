#include <iostream>
#include <vector>
#include "plf_nanotimer.h"
#include "configuration.hpp"
#include "graphs.hpp"
#include "graph_algorithms.hpp"
#include "graph_viz.hpp"

int main(int argc, char* argv[])
{
    Configuration& config = Configuration::getInstance();
    config.parseCommand(argc, argv);
    if (config.isConfigValid() == false) {
        std::cerr << "\nInvalid configuration\n";
        exit(0);
    }
    if (config.getVerboseLevel() > 0) {
        config.printConfig();
    }
    Graph graph;
    if (config.hasGraphFile()) {
        graph.initFromMtxFile(config.getGraphSourceFileName());
    } else {
        switch (config.getGGAlgorithm()) {
        case GGAlgorithm::Random:
            if (config.getEdgeProbability() != 0.0) {
                graph.initWithProbability(config.getNumVertices(), config.getEdgeProbability());
            } else if (config.getNumEdges() != 0) {
                graph.initFixedEdges(config.getNumVertices(), config.getNumEdges());
            }
            break;
        case GGAlgorithm::BA:
            graph.initBA(config.getMValue(), config.getNumVertices());
            break;
        case GGAlgorithm::RGG:
            graph.initRGG(config.getNumVertices(), config.getDistanceThreshold());
        default:
            break;
        }
    }
    
    GraphViz viz1(&graph, config.getGraphDrawingDir());
    GraphViz viz2(&graph, config.getGraphDrawingDir());
    GraphViz viz3(&graph, config.getGraphDrawingDir());
    GraphViz viz4(&graph, config.getGraphDrawingDir());
    viz1.colorDistanceVertex(0);
    viz2.colorCliques(1);
    viz3.colorBicliques(1);
    viz4.colorSpanningTree(0);
    viz1.renderToPNG("BFS.png");
    viz2.renderToPNG("Cliques.png");
    viz3.renderToPNG("Bicliques.png");
    viz4.renderToPNG("SpanningTree.png");
    viz1.closeGraphVizz();
    viz2.closeGraphVizz();
    viz3.closeGraphVizz();  
    viz4.closeGraphVizz();

    graph.printMetaData();
    long start_id = 0;
    long end_id = 0;
    plf::nanotimer timer;
    timer.start();

    double t0 = timer.get_elapsed_ms();
    long degeneracy = computeDegeneracy(graph);
    double t1 = timer.get_elapsed_ms();
    int diameter = computeDiameterBfs2(graph, start_id, end_id);
    double t2 = timer.get_elapsed_ms();
    std::vector<std::set<long>> cliques;
    bronKerbosch(graph, cliques);
    double t3 = timer.get_elapsed_ms();
    std::vector<Graph> bicliques = computeBiclique(graph);
    double t4 = timer.get_elapsed_ms();
    
    std::cout << "degeneracy : " << degeneracy << "\n";
    std::cout << "Diameter   : " << diameter << " " << start_id << "---" << end_id << "\n";
    std::cout << "BiCliques  : " << bicliques.size() << "\n";
    std::cout << "Cliques    : " << cliques.size() << "\n";
    std::cout << "Time to compute degeneracy : " << t1 - t0 << " ms\n";
    std::cout << "Time to compute diameter   : " << t2 - t1 << " ms\n";
    std::cout << "Time to compute cliques    : " << t3 - t2 << " ms\n";
    std::cout << "Time to compute bicliques  : " << t4 - t3 << " ms\n";

    return 0;
}

