#include <fstream>
#include <string>
#pragma once

// Balise pour le fichier de configuration xml
#define TAG_BASE          "graph-config"
#define TAG_TITLE         "title"
#define TAG_OUT_DIR       "out-dir"
#define TAG_DRAW_DIR      "draw-dir"
#define TAG_GRAPH_FILE    "graph-file"
#define TAG_VERTICES      "vertices"
#define TAG_EDGES         "edges"
#define TAG_P_EDGE        "p-edge"
#define TAG_M             "m"
#define TAG_DISTANCE      "distance-trheshold"
#define TAG_GENERATOR     "graph-generator"
#define TAG_BENCHMARK     "benchmark"
#define TAG_MEASURE_ITER  "measure-iter"
#define TAG_VERBOSE_LEVEL "verbose-level"

// Algorithme de génération de graphe
enum class GGAlgorithm { Random = 0, RGG, BA, Complete, None };

// Classe de configuration de l'exécution
class Configuration
{
  private:
    long num_vertices;
    long num_edges;
    double edge_probability;
    double distance_threshold;
    int m_value;
    int verbose_level;
    int measure_iterations;
    bool enable_dump;
    bool enable_benchmark;
    bool use_graph_src_file;
    bool enable_graph_drawing;
    std::string graph_src_file_name;
    std::string graph_drawing_dir;
    std::string title;
    std::string output_dir;
    GGAlgorithm gga;

    void parseFile(std::string file_name);

  public:
    static Configuration& getInstance()
    {
        static Configuration instance;
        return instance;
    }

    Configuration(const Configuration&) = delete;
    Configuration& operator=(const Configuration&) = delete;

    bool isConfigValid();
    void printConfig();
    void parseCommand(int argc, char* argv[]);
    bool shouldEnableBenchmark() const;
    bool hasGraphFile() const;
    long getNumVertices() const;
    long getNumEdges() const;
    double getEdgeProbability() const;
    int getMValue() const;
    double getDistanceThreshold() const;
    int getMeasureIterations() const;
    int getVerboseLevel() const;
    GGAlgorithm getGGAlgorithm() const;
    std::string getOutputDir() const;
    std::string getGraphDrawingDir() const;
    std::string getGraphSourceFileName() const;

  private:
    Configuration();
    ~Configuration();
};

void printUsage();
GGAlgorithm getGGAFromString(std::string algorithm_name);

