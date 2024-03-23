#include <array>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include "configuration.hpp"
#include "../lib/tinyxml2/tinyxml2.h"

Configuration::Configuration()
{
    num_vertices = 0;
    num_edges = 0;
    edge_probability = 0.0;
    m_value = 0;
    distance_threshold = 0.0;
    enable_dump = true;
    enable_benchmark = false;
    verbose_level = 0;
    measure_iterations = 10;
    gga = GGAlgorithm::None;
    title = std::string();
    output_dir = std::string("out");
    use_graph_src_file = false;
    enable_graph_drawing = false;
    graph_src_file_name = std::string();
    graph_drawing_dir = std::string("draws");
}

Configuration::~Configuration() {}

// Parse the command line
void Configuration::parseCommand(int argc, char* argv[])
{
    extern char* optarg;
    int opt;
    if (argc == 1) {
        printUsage();
        exit(EXIT_SUCCESS);
    } else if (argc >= 2) {
        if (std::filesystem::is_regular_file(argv[1])) {
            try {
                this->parseFile(argv[1]);
            } catch (const std::exception& e) {
                std::cerr << e.what() << '\n';
            }
        }
    }
    while ((opt = getopt(argc, argv, "hdf:e:v:p:t:m:")) != -1) {
        switch (opt) {
        case 'h':
            printUsage();
            exit(EXIT_SUCCESS);
            exit(0);
            break;
        case 'v':
            this->num_vertices = atol(optarg);
            break;
        case 'e':
            this->num_edges = atol(optarg);
            break;
        case 'p':
            this->edge_probability = atof(optarg);
            gga = GGAlgorithm::Random;
            break;
        case 't':
            this->distance_threshold = atof(optarg);
            gga = GGAlgorithm::RGG;
            break;
        case 'm':
            this->m_value = atoi(optarg);
            gga = GGAlgorithm::BA;
            break;
        case 'd':
            this->enable_graph_drawing = true;
            break;
        case 'f':
            this->use_graph_src_file = true;
            this->graph_src_file_name.assign(optarg);
            break;
        }
    }
    return;
}

// Parse the XML configuration file
void Configuration::parseFile(std::string file_name)
{
    tinyxml2::XMLError xml_error;
    tinyxml2::XMLDocument xml_document;
    xml_error = xml_document.LoadFile(file_name.c_str());
    if (xml_error != tinyxml2::XML_SUCCESS) {
        std::cout << "No file to parse" << std::endl;
        return;
    }

    tinyxml2::XMLElement* BASE_element = xml_document.FirstChildElement(TAG_BASE);

    if (BASE_element == nullptr) {
        std::string error("Missing root tag in XML FILE : ");
        error.append(TAG_BASE);
        throw std::runtime_error(error);
    }

    tinyxml2::XMLElement* out_dir_element = BASE_element->FirstChildElement(TAG_OUT_DIR);
    tinyxml2::XMLElement* title_element = BASE_element->FirstChildElement(TAG_TITLE);
    tinyxml2::XMLElement* graph_file_element = BASE_element->FirstChildElement(TAG_GRAPH_FILE);
    tinyxml2::XMLElement* make_benchmark_element = BASE_element->FirstChildElement(TAG_BENCHMARK);
    tinyxml2::XMLElement* measure_iter_element = BASE_element->FirstChildElement(TAG_MEASURE_ITER);
    tinyxml2::XMLElement* verbose_element = BASE_element->FirstChildElement(TAG_VERBOSE_LEVEL);
    tinyxml2::XMLElement* drawDir_element = BASE_element->FirstChildElement(TAG_DRAW_DIR);
    tinyxml2::XMLElement* generation_element = BASE_element->FirstChildElement(TAG_GENERATOR);
    tinyxml2::XMLElement* edges_element = nullptr;
    tinyxml2::XMLElement* vertices_element = nullptr;
    tinyxml2::XMLElement* p_link_element = nullptr;
    tinyxml2::XMLElement* distance_element = nullptr;
    tinyxml2::XMLElement* m_element = nullptr;

    if (generation_element != nullptr) {
        const tinyxml2::XMLAttribute* AG_attribute = generation_element->FirstAttribute();
        if (AG_attribute != nullptr) {
            this->gga = getGGAFromString(AG_attribute->Value());
        }
        vertices_element = generation_element->FirstChildElement(TAG_VERTICES);
        edges_element = generation_element->FirstChildElement(TAG_EDGES);
        p_link_element = generation_element->FirstChildElement(TAG_P_EDGE);
        distance_element = generation_element->FirstChildElement(TAG_DISTANCE);
        m_element = generation_element->FirstChildElement(TAG_M);
    }
    if (vertices_element != nullptr) {
        vertices_element->QueryInt64Text(&this->num_vertices);
    }
    if (edges_element != nullptr) {
        edges_element->QueryInt64Text(&this->num_edges);
    }
    if (p_link_element != nullptr) {
        p_link_element->QueryDoubleText(&this->edge_probability);
    }
    if (distance_element != nullptr) {
        distance_element->QueryDoubleText(&this->distance_threshold);
    }
    if (m_element != nullptr) {
        m_element->QueryIntText(&this->m_value);
    }
    if (out_dir_element != nullptr) {
        output_dir = std::string(out_dir_element->GetText());
    }
    if (title_element != nullptr) {
        title = std::string(title_element->GetText());
    }
    if (drawDir_element != nullptr) {
        this->graph_drawing_dir = std::string(drawDir_element->GetText());
        this->enable_graph_drawing = true;
    }
    if (graph_file_element != nullptr) {
        try {
            graph_src_file_name = std::string(graph_file_element->GetText());
            use_graph_src_file = true;
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
            use_graph_src_file = false;
        }
    }
    if (make_benchmark_element != nullptr) {
        std::string bench = std::string(make_benchmark_element->GetText());
        if (bench.compare("yes") == 0 || bench.compare("1") == 0) {
            enable_benchmark = true;
        }
    }
    if (measure_iter_element != nullptr) {
        measure_iter_element->QueryIntText(&this->measure_iterations);
    }
    if (verbose_element != nullptr) {
        verbose_element->QueryIntText(&this->verbose_level);
    }
}

// Check if the configuration is correct
bool Configuration::isConfigValid()
{
    using namespace std::filesystem;
    if (num_vertices <= 0) {
        std::cerr << "\nBad number of vertices : " << num_vertices << "\n";
        return false;
    }
    if (edge_probability < 0.0 || edge_probability > 1.0) {
        std::cerr << "\nBad edge probability : " << edge_probability << "\n";
        return false;
    }
    if (edge_probability == 0.0 && num_edges == 0 && m_value == 0 && distance_threshold == 0.0) {
        std::cerr << "\nBad generation parameters\n";
        return false;
    }
    if (m_value < 0) {
        std::cerr << "\nBad m value : " << m_value << "\n";
        return false;
    }
    if (distance_threshold < 0.0) {
        std::cerr << "\nBad distance : " << distance_threshold << "\n";
        return false;
    }
    if (num_edges < 0 || num_edges > (num_vertices * (num_vertices - 1)) / 2) {
        std::cerr << "\nBad number of edges : " << num_edges << "\n";
        return false;
    }
    if (use_graph_src_file &&
        (!exists(graph_src_file_name) || !is_regular_file(graph_src_file_name))) {
        std::cerr << "\nCannot open : " << graph_src_file_name << "\n";
        return false;
    }
    if (verbose_level < 0 || verbose_level > 2) {
        std::cerr << "\nBad verbose level must be 0,1 or 2\n";
        return false;
    }
    auto workingDir = current_path();
    auto outDir = workingDir.append(output_dir);
    if (!exists(outDir) || !is_directory(outDir)) {
        if (!create_directory(outDir)) {
            std::cerr << "Failed to create " << outDir << std::endl;
        }
    }
    workingDir = current_path();
    auto drawDir = workingDir.append(this->graph_drawing_dir);
    if (!is_directory(drawDir)) {
        if (!create_directory(drawDir)) {
            std::cerr << "Failed to create " << drawDir << std::endl;
        }
    }
    return true;
}

// Print the configuration
void Configuration::printConfig()
{
    using namespace std;
    cout << "\n";
    cout.width(24);
    cout << left << "Title : " << right << title << endl;
    cout.width(24);
    cout << left << "Number of vertices : " << right << num_vertices << endl;
    cout.width(24);
    cout << left << "Number of edges : " << right << num_edges << endl;
    cout.width(24);
    cout << left << "Edge probability : " << right << edge_probability << endl;
    cout.width(24);
    cout << left << "m value : " << right << m_value << endl;
    cout.width(24);
    cout << left << "Distance threshold : " << right << distance_threshold << endl;
    cout.width(24);
    cout << left << "Out data directory : " << right << output_dir << endl;
    cout.width(24);
    cout << left << "Draw directory : " << right << graph_drawing_dir << endl;
    cout.width(24);
    cout << left << "Source file : " << right << graph_src_file_name << endl;
    cout.width(24);
    cout << left << "Verbose level : " << right << verbose_level << "\n" << endl;
}

void printUsage()
{
    std::cout << "\nUsage :\n"
              << "\tgraphs [file.xml]\n"
              << "\t-v [nb vertices]\n"
              << "\t-e [nb edges]\n"
              << "\t-p [probability]\n"
              << "\t-t [distance]\n"
              << "\t-m [m value]\n"
              << "\t-d\n"
              << "\t-f [file.mtx]\n";
}

GGAlgorithm getGGAFromString(std::string algo_name)
{
    GGAlgorithm gga = GGAlgorithm::None;
    std::array<std::string, 3> rand_name = {"Random", "random", "RANDOM"};
    std::array<std::string, 3> complete_name = {"Complete", "complete", "COMPLETE"};
    std::array<std::string, 2> BA_name = {"BA", "ba"};
    std::array<std::string, 2> RGG_name = {"RGG", "rgg"};
    for (const auto& name : rand_name) {
        if (algo_name.compare(name) == 0) {
            gga = GGAlgorithm::Random;
        }
    }
    for (const auto& name : complete_name) {
        if (algo_name.compare(name) == 0) {
            gga = GGAlgorithm::Complete;
        }
    }
    for (const auto& name : BA_name) {
        if (algo_name.compare(name) == 0) {
            gga = GGAlgorithm::BA;
        }
    }
    for (const auto& name : RGG_name) {
        if (algo_name.compare(name) == 0) {
            gga = GGAlgorithm::RGG;
        }
    }
    return gga;
}

bool Configuration::shouldEnableBenchmark() const { return enable_benchmark; }

bool Configuration::hasGraphFile() const { return use_graph_src_file; }

long Configuration::getNumVertices() const { return num_vertices; }

long Configuration::getNumEdges() const { return num_edges; }

double Configuration::getEdgeProbability() const { return edge_probability; }

int Configuration::getMValue() const { return this->m_value; }

double Configuration::getDistanceThreshold() const { return this->distance_threshold; }

GGAlgorithm Configuration::getGGAlgorithm() const { return this->gga; }

int Configuration::getMeasureIterations() const { return measure_iterations; }

int Configuration::getVerboseLevel() const { return verbose_level; }

std::string Configuration::getOutputDir() const { return output_dir; }

std::string Configuration::getGraphDrawingDir() const { return graph_drawing_dir; }

std::string Configuration::getGraphSourceFileName() const { return graph_src_file_name; }
