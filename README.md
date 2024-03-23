# Graphs programme

Program for testing various algorithms on graphs. The project is managed by CMake, and a script is available to run the build command. The directory also contains configuration files (in XML format) for execution.

## Install

1. Clone the source file.  
`$ git clone https://github.com/AdRi1t/graphs.git`
2. Run the setup script, this will compile the whole program  
`$ ./setup.sh`
3. Try to execute  
`$ cd build`  
`$./graph ../example_1.xml`

## Libraries

- [TinyXML2](https://github.com/leethomason/tinyxml2) : read and write in XML file
- [GraphViz](https://graphviz.org/docs/library/) : Transforms graphs into png images

## Example

This section contains images showing the results of the algorithms.

### Breadth First Search

<img src="https://github.com/AdRi1t/Graphs/tree/main/git/BFS6.png" alt="BFS" width="400"/>

### Cliques

<img src="https://github.com/AdRi1t/Graphs/tree/main/git/Cliques.png" alt="Cliques" width="300"/>

### Bicliques

<img src="https://github.com/AdRi1t/Graphs/tree/main/git/graphs_8.png" alt="Bicliques" width="300"/>


### Spanning tree

<img src="https://github.com/AdRi1t/Graphs/tree/main/git/SpanningTree2.png" alt="SpanningTree" width="200"/>
