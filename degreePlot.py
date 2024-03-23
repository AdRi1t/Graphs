from matplotlib import pyplot
import numpy
import os
import sys

def main():
  if(len(sys.argv) < 2):
    sys.stderr.write("Miss file to plot \n")
    return 0
  data_file = open(sys.argv[1],"rt")
  ID = []
  degree = []
  with data_file as file:
    for line in file:
      n1, n2 = map(int, line.split())
      ID.append(n1)
      degree.append(n2)  
  
  setDegree = list(range(min(degree), max(degree)+1,1))
  degreePOP = [0] * (max(degree) - min(degree) + 1)
  for id in ID:
    degreePOP[degree[id] - min(degree)] += 1

  fig1, ax1 = pyplot.subplots()
  ax1.bar(setDegree, degreePOP, label="Degree plot")
  ax1.set_axisbelow(True)
  ax1.grid(True)
  ax1.set_xlabel("Degree") 
  ax1.set_ylabel("Number of vertices")  
  ax1.set_title("Degree distribution")
  fig1.tight_layout()
  fig1.savefig("degreePlot.png", dpi=300)
  
  return 0


if __name__ == "__main__":
  main()
  
  
