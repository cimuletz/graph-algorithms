

#ifndef GRAPH_ALGORITHMS_GRAPH_H
#define GRAPH_ALGORITHMS_GRAPH_H


#include <vector>
#include <stack>

struct edgeCost{
    int x, c;
};

class Graph{
    std::vector<std::vector<int>> Ad; //adjacency list
    std::vector<std::vector<edgeCost>> AdCost; //adjacency with cost list, also used for flow
    int nodes; // no of nodes
    int edges; // no of edges
    bool directed;
    const int inf = 100000000;
public:
    Graph();
    Graph(const std::vector<std::pair<int,int>> &ad, int nodes, int edges,int directed = 0);
    Graph(const std::vector<std::vector<edgeCost>> &ad, int nodes, int edges,int directed = 0);
    Graph(const std::vector<std::pair<int,int>> &ad, const std::vector<int> &cost, int nodes, int edges,int directed = 0);
    std::vector<int> getComponents(); //get connected components
    std::vector<int> getBFSPathLength(int node); //path length from argument to all other nodes
    std::stack<int> getTopSort(); //topological sort
    std::vector<std::vector<int>> getFloydWarshall();
    std::vector<std::vector<int>> getSCC(); //returns strongly connected components
    std::vector<std::vector<int>> getBCC(); //returns biconex components
    std::vector<int> solveDisjoint(); //list with result of find operation between two nodes
    std::vector<std::tuple<int, int, int>> getMST(); //returns edges which make up the mst
    std::vector<int> getBellmanDistance(int node); //returns distances from node to all other nodes using bellman ford
    std::vector<int> getDijkstraDistance(int node); //returns distances from node to all other nodes using dijkstra
    int getTreeSize(); // returns longest path between two leafs in a tree
    int getFlow(int src, int dest); // returns max flow between src and dest

private:
    void DFS(int node, std::vector<int> &vis, int nrComp);
    void BFS(int node, std::vector<int> &vis);
    void topSortDFS(int node, std::vector<int> &vis, std::stack<int> &s); //dfs for topological sort
    int findMST(int x, std::vector<int> &disjSets); //find function in union find for mst
    void unionMST(int x, int y, std::vector<int> &disjSets, std::vector<int> &sizes); //union function in union find
    void dfsSCC(int node,int &idInc, std::vector<int> &id, std::vector<int> &low, std::vector<bool> &onStack, std::stack<int> &s, std::vector< std::vector<int>> &scc); //dfs for findind strongly connected components
    void dfsBCC(int node, int parent, int &time, std::vector<int> &dq, std::vector<int> &low, std::vector<int> &disc, std::vector<std::vector<int>> &bcc);
    void bellman(bool &cycle, std::vector<int> &dist, int node); //bellman ford algorithm
    void dijkstra(std::vector<int> &dist, int node); //dijkstra algorithm
    void dfsTS(std::vector<int> &vis, int node, int parent, int &last, int &level); //dfs for tree size
    int bfsFlow(int src, int dest, std::vector<std::vector<int>> &edgeFlow, std::vector<int> &vis, std::vector<int> &bottleneck); // bfs for max flow
};

#endif //GRAPH_ALGORITHMS_GRAPH_H
