#include "graph.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits.h>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
map<string, vector<pair<string, int>>> adjacencyList;
vector<string> tester;

// constructor, empty graph
// directionalEdges defaults to true
Graph::Graph(bool directionalEdges) {
  this->directionalEdges = directionalEdges;
}
//----------------------------------------------------------------------------------------
// destructor
Graph::~Graph() {
  for (auto &vertex : adjacencyList) {
    vertex.second.clear();
  }
  // clear the entire adjacency list
  adjacencyList.clear();
}
//----------------------------------------------------------------------------------------
// @return total number of vertices
int Graph::verticesSize() { return static_cast<int>(adjacencyList.size()); }
//----------------------------------------------------------------------------------------
// @return total number of edges
int Graph::edgesSize() {
  int numOfEdges = 0;
  for (const auto &vertex : adjacencyList) {
    numOfEdges += static_cast<int>(vertex.second.size());
  }
  if (!directionalEdges) {
    numOfEdges /= 2;
  }
  return numOfEdges;
}
//----------------------------------------------------------------------------------------
// @return number of edges from given vertex, -1 if vertex not found
int Graph::vertexDegree(const string &label) {
  // this shows that the vertex has not been found given the label
  if (adjacencyList.find(label) == adjacencyList.end()) {
    return -1;
  }
  return adjacencyList.at(label).size();
}
//----------------------------------------------------------------------------------------
// @return true if vertex added, false if it already is in the graph
bool Graph::add(const string &label) {
  if (adjacencyList.find(label) != adjacencyList.end()) {
    return false;
  }
  // adding the vertex, and there are no edges connected to it yet
  adjacencyList[label] = {};
  return true;
}
//----------------------------------------------------------------------------------------
/** return true if vertex already in graph */
bool Graph::contains(const string &label) {
  return adjacencyList.find(label) != adjacencyList.end();
}
//----------------------------------------------------------------------------------------
// @return string representing edges and weights, "" if vertex not found
// A-3->B, A-5->C should return B(3),C(5)
string Graph::getEdgesAsString(const string &label) {
  if (adjacencyList.find(label) != adjacencyList.end()) {
    // converting the second part to string here
    // this is also to sort the edges in lexographical order
    vector<pair<string, int>> edges = adjacencyList.at(label);
    sort(edges.begin(), edges.end());
    string result = "";
    for (const auto &edge : edges) {
      result += edge.first + "(" + to_string(edge.second) + ")" + ",";
    }

    // this is for the trailing comma that we have at the end
    if (!result.empty()) {
      result.pop_back();
    }
    return result;
  }
  return "";
}
//----------------------------------------------------------------------------------------
// @return true if successfully connected
bool Graph::connect(const string &from, const string &to, int weight) const {
  // check if its in the list
  // if they are the same
  if (from == to) {
    return false;
  }
  if (!this->contains(from)) {
    adjacencyList[from] = {};
  }
  if (!this->contains(to)) {
    adjacencyList[to] = {};
  }

  // this is to check if this edge already exists
  for (const auto &checker : adjacencyList[from]) {
    if (checker.first == to) {
      return false;
    }
  }
  // this is when we create the edge
  if (!directionalEdges) {
    adjacencyList[from].push_back(make_pair(to, weight));
    adjacencyList[to].push_back(make_pair(from, weight));
  } else {
    adjacencyList[from].push_back(make_pair(to, weight));
  }

  return true;
}
//---------------------------------------------------------------------------------------------------------------
bool Graph::disconnect(const string &from, const string &to) {
  // check if its in the list
  if (!this->contains(from) || !this->contains(to)) {
    return false;
  }
  // checking from the "from" vertex
  for (const auto &checker : adjacencyList[from]) {
    if (checker.first == to) {
      adjacencyList[from].pop_back();
      // edge has been disconnected
      return true;
    }
  }
  // edge has not been disconnected.
  return false;
}
//---------------------------------------------------------------------------------------------------------------
// depth-first traversal starting from given startLabel
void Graph::dfs(const string &startLabel, void visit(const string &label)) {
  if (adjacencyList.find(startLabel) == adjacencyList.end()) {
    // start vertex is not in the graph, do nothing
    return;
  }
  unordered_set<string> visitedVertex;
  stack<string> stack;
  visitedVertex.insert(startLabel);
  stack.push(startLabel);
  while (!stack.empty()) {
    string currentLabel = stack.top();
    stack.pop();
    visit(currentLabel);
    vector<string> sortedList;
    for (const auto &list : adjacencyList[currentLabel]) {
      sortedList.push_back(list.first);
    }
    sort(sortedList.begin(), sortedList.end());
    for (const auto &neighbor : adjacencyList[currentLabel]) {
      if (visitedVertex.count(neighbor.first) == 0) {
        visitedVertex.insert(neighbor.first);
        stack.push(neighbor.first);
      }
    }
  }
}
//---------------------------------------------------------------------------------------------------------------
// breadth-first traversal starting from startLabel
// the exact same as dfs but uses queue instead of stack.
void Graph::bfs(const string &startLabel, void visit(const string &label)) {
  if (adjacencyList.find(startLabel) == adjacencyList.end()) {
    // start vertex is not in the graph, do nothing
    return;
  }
  unordered_set<string> visitedVertex;
  queue<string> queue;
  visitedVertex.insert(startLabel);
  auto itr = visitedVertex.begin();
  queue.push(startLabel);
  while (!queue.empty()) {
    string currentLabel = queue.front();
    queue.pop();
    visit(currentLabel);
    // sort the adjacencyList
    vector<string> sortedList;
    for (const auto &list : adjacencyList[currentLabel]) {
      sortedList.push_back(list.first);
    }
    sort(sortedList.begin(), sortedList.end());

    for (const auto &neighbor : sortedList) {
      if (visitedVertex.count(neighbor) == 0) {
        visitedVertex.insert(neighbor);
        queue.push(neighbor);
      }
    }
  }
}
//---------------------------------------------------------------------------------------------------------------
// store the weights in a map
// store the previous label in a map
pair<map<string, int>, map<string, string>> Graph::dijkstra(const string &startLabel) {
    map<string, int> weights;
    map<string, string> previous;

    // Initialization
    for (const auto &vertex : adjacencyList) {
        weights[vertex.first] = INT_MAX;
        previous[vertex.first] = "";
    }
    weights[startLabel] = 0;

    set<string> unvisited;
    for (const auto &vertex : adjacencyList) {
        unvisited.insert(vertex.first);
    }

    while (!unvisited.empty()) {
        // Find the vertex with the minimum weight
        string min_vertex = *unvisited.begin();
        for (const string &vertex : unvisited) {
            if (weights[vertex] < weights[min_vertex]) {
                min_vertex = vertex;
            }
        }

        // Remove the selected vertex from unvisited
        unvisited.erase(min_vertex);

        // Update the weights and previous vertices
        for (const auto &neighbor : adjacencyList[min_vertex]) {
            int new_weight = weights[min_vertex] + neighbor.second;
            if (new_weight < weights[neighbor.first]) {
                weights[neighbor.first] = new_weight;
                previous[neighbor.first] = min_vertex;
            }
        }
    }

    weights.erase(startLabel);
    previous.erase(startLabel);

    return make_pair(weights, previous);
}

//---------------------------------------------------------------------------------------------------------------
// minimum spanning tree using Prim's algorithm
int Graph::mstPrim(const string &startLabel,
                   void visit(const string &from, const string &to,
                              int weight)) const {
  int weight = 0;
  unordered_set<string> visited;
  priority_queue<pair<int, string>, vector<pair<int, string>>,
                 greater<pair<int, string>>>
      pq;
  pq.push(make_pair(0, startLabel));

  while (!pq.empty()) {
    pair<int, string> current = pq.top();
    pq.pop();

    // we will first check if that node has been visited, if it has, then we
    // continue to the next
    if (visited.find(current.second) != visited.end()) {
      continue;
    }

    visited.insert(current.second);
    weight += current.first;

    // for loop for going ot the neigbors and adding to priority queue.
    for (const auto &neighbor : adjacencyList.at(current.second)) {
      if (visited.find(neighbor.first) == visited.end()) {
        pq.push(make_pair(neighbor.second, neighbor.first));
        visit(current.second, neighbor.first, neighbor.second);
      }
    }
  }
  return weight;
}
//---------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// read a text file and create the graph
bool Graph::readFile(const string &filename) const {
  ifstream myfile(filename);
  if (!myfile.is_open()) {
    cerr << "Failed to open " << filename << endl;
    return false;
  }
  int edges = 0;
  int weight = 0;
  string fromVertex;
  string toVertex;
  myfile >> edges;

  for (int i = 0; i < edges; ++i) {
    myfile >> fromVertex >> toVertex >> weight;
    connect(fromVertex, toVertex, weight);
  }

  cout << "Adjacency List:" << endl;
  for (const auto &vertex : adjacencyList) {
    cout << vertex.first << ": ";
    for (const auto &edge : vertex.second) {
      cout << "(" << edge.first << ", " << edge.second << ") ";
    }
    cout << endl;
  }
  myfile.close();
  return true;
}