#include <iostream>
#include "configuration.hpp"
#include "graphs.hpp"
#include "graph_algorithms.hpp"
#include "graph_viz.hpp"

int test_configurations(int argc, char* argv[]);
int test_voidGraph();
int test_addVertex();
int test_addEdge();
int test_removeVertex();
int test_writeDebugInfo();
int test_loadGraph();
int test_drawGraph();
int test_writeDot();
int test_unionGraph();
int test_bronKerbosch();
using namespace std;

int main(int argc, char* argv[])
{
    // Test de la classe Configuration
    test_configurations(argc, argv);

    // Crée un Graph vide
    test_voidGraph();

    // Ajout d'un sommet
    test_addVertex();

    // Ajout d'une arête
    test_addEdge();

    // Supression d'un sommet
    test_removeVertex();
    
    // Écris les informations du graph dans un fichier texte
    test_writeDebugInfo();

    // Charge un graph depuis un fichier
    test_loadGraph();

    // Création d'une image à partir d'un graph
    test_drawGraph();

    // Donne le fichier Dot associer au graph
    test_writeDot();

    // !
    test_bronKerbosch();

    /*
    // Test l'union de de graph
    test_unionGraph();

    // Test de l'intersection avec les voisins d'un sommet
    test_interNeighborGraph();

    // Test la suppression des voisins d'un sommet
    test_excludingNeighbors();

    // Test l'union d'un sommet d'un autre graph dans un graph
    test_unionVertex();
    */

    return 0;
}

int test_configurations(int argc, char* argv[])
{
    Configuration& config = Configuration::getInstance();
    config.parseCommand(argc, argv);
    config.printConfig();
    config.isConfigValid();
    return 0;
}

int test_voidGraph()
{
    Graph g;
    cout << "\nTest Graphs vide \n";
    g.printEdgesList();
    cout << "\n";
    return 0;
}

int test_addVertex()
{
    Graph g;
    cout << "\nTest ajout des sommets [0,1,42]\n";
    g.addVertex(0);
    g.addVertex(1);
    g.addVertex(42);
    // Test Ajout d'un sommet de même ID.
    g.addVertex(42);
    g.printVerticesList();
    cout << "\n";
    return 0;
}

int test_addEdge()
{
    Graph g;
    cout << "\nTest ajout des arêtes [0 - 1, 0 - 42]\n";
    g.addVertex(0);
    g.addVertex(1);
    g.addVertex(42);
    g.addEdge(0, 1);
    g.addEdge(0, 42);
    // Test Ajout de la même arête
    g.addEdge(0, 42);
    // Test Ajout d'une arête impossible
    g.addEdge(2, 3);
    g.printVerticesList();
    g.printEdgesList();
    cout << "\n";
    return 0;
}

int test_removeVertex()
{
    Graph g;
    cout << "\nTest suppression du sommet 1 dans [0-1, 0-2, 1-2]\n";
    g.addVertex(0);
    g.addVertex(1);
    g.addVertex(2);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.printVerticesList();
    g.printEdgesList();
    cout << "Suppression du sommet 1\n";
    g.removeVertex(1);
    g.printVerticesList();
    g.printEdgesList();
    // Test suppression d'un sommet inexistant
    g.removeVertex(42);
    g.printMetaData();
    cout << "\n";
    return 0;
}

int test_bronKerbosch()
{
    Graph A;
    std::vector<std::set<long>> cliques;
    A.initRGG(10, 0.4);
    GraphViz viz;
    viz.setGraph(&A);
    viz.setRenderDirectory("test");
    viz.renderToPNG("Bron.png");
    viz.closeGraphVizz();

    bronKerbosch(A, cliques);

    for (auto clique : cliques) {
        std::cout << "Clique : ";
        for (const long& vertex_id : clique) {
            std::cout << vertex_id << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Size : " << cliques.size() << "\n";
    return 0;
}

int test_drawGraph()
{
    Graph g;
    GraphViz viz;
    g.initWithProbability(10, 0.5);
    viz.setGraph(&g);
    viz.setRenderDirectory("test");
    viz.renderToPNG("test_PNG_render.png");
    viz.closeGraphVizz();
    return 0;
}

int test_writeDot()
{
    Graph g;
    GraphViz viz;
    g.initWithProbability(10, 0.5);
    viz.setGraph(&g);
    viz.setRenderDirectory("test");
    viz.renderToDOT("test_Dot.dot");
    viz.closeGraphVizz();
    return 0;
}

int test_writeDebugInfo()
{ 
    Graph g;
    GraphViz viz;
    Configuration& config = Configuration::getInstance();
    g.initWithProbability(10, 0.5);
    viz.setGraph(&g);
    viz.setRenderDirectory(config.getOutputDir());
    viz.printEdges("test_edges.txt");
    viz.printDegree("test_degree.txt");
    viz.closeGraphVizz();
    return 0; 
}

int test_loadGraph()
{   
    Graph g;
    Configuration& config = Configuration::getInstance();
    g.initFromMtxFile(config.getGraphSourceFileName());
    g.printMetaData();
    return 0;
}

int test_unionGraph()
{
    Graph A;
    Graph B;
    A.addVertex(0);
    A.addVertex(1);
    A.addVertex(2);

    A.addEdge(0, 2);
    A.addEdge(0, 1);

    B.addVertex(1);
    B.addVertex(4);
    B.addEdge(1, 4);

    // C = A u B
    Graph C = unionGraphs(A, B);
    cout << "\nTest sur l'union de deux graphs A U B = C\n";
    cout << "A\n";
    A.printEdgesList();
    cout << "B\n";
    B.printEdgesList();
    cout << "C\n";
    C.printEdgesList();
    cout << "\n";
    return 0;
}
