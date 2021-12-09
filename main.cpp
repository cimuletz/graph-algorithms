#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <algorithm>
#include <deque>
#include <tuple>
#include <cstring>
#include "Graph.h"


std::ifstream fin("input.in");
std::ofstream fout("output.out");

std::vector<std::pair<int,int>> readAdList(int n, int m, int directed){
    std::vector<std::pair<int,int>> Ad;
    int x,y;
    for(int i = 0; i<m; ++i){
        fin >> x >> y;
        Ad.push_back({x,y});
    }
    return Ad;
}

std::vector<std::vector<edgeCost>> readAdCostListFromMatrix(int n){
    std::vector<std::vector<edgeCost>> ad(n+1);
    edgeCost aux;
    for(int i = 1; i<=n; ++i)
        for(int j = 1; j<=n; ++j){
            fin >> aux.c;
            aux.x = j;
            ad[i].push_back(aux);
        }
    return ad;
}

std::vector<std::vector<edgeCost>> readAdCostList(int n, int m){
    std::vector<std::vector<edgeCost>> Ad(n+1);
    int x,y,c;
    for(int i = 0; i<m; ++i){
        fin >> x >> y >> c;
        edgeCost aux;
        aux.x = y;
        aux.c = c;
        Ad[x].push_back(aux);
    }
    return Ad;
}


void dfsHelper(){
    int n, m, ctComp = 0;
    std::vector<int> components;
    fin >> n >> m;
    Graph G(readAdList(n, m, 0), n, m);
    components = G.getComponents();
    for(int i = 1; i<=n; ++i)
        if(components[i] > ctComp)
            ctComp = components[i];
    fout << ctComp;
}

void bfsHelper(){
    std::vector<int> pathLength;
    int n, m, source;
    fin >> n >> m >> source;
    Graph G(readAdList(n, m, 0), n, m, 1);
    pathLength = G.getBFSPathLength(source);
    for(int i = 1; i<=n; ++i)
        fout << pathLength[i] << " ";
}

void topsortHelper(){
    std::stack<int> topSort;
    int n, m;
    fin >> n >> m;
    Graph G(readAdList(n, m, 0), n, m, 0);
    topSort = G.getTopSort();
    while(!topSort.empty()){
        fout << topSort.top() << " ";
        topSort.pop();
    }
}

void bccHelper(){
    std::vector<std::vector<int>> bcc;
    int n,m;
    fin >> n >> m;
    Graph G(readAdList(n, m, 0), n, m);
    bcc = G.getBCC();
    fout << bcc.size() << "\n";
    int nr = bcc.size();
    for(int i = 0; i < nr; ++i){
        int mr = bcc[i].size();
        for(int j = 0; j < mr; ++j)
            fout << bcc[i][j] << " ";
        fout << "\n";
    }

}

void sccHelper(){
    std::vector<std::vector<int>> scc;
    int n,m;
    fin >> n >> m;
    Graph G(readAdList(n, m, 0), n, m, 1);
    scc = G.getSCC();
    fout << scc.size() << "\n";
    for(int i = 0; i<scc.size(); ++i){
        for(int j = 0; j<scc[i].size(); ++j)
            fout << scc[i][j] << " ";
        fout << "\n";
    }
}

void floydWarshallHelper(){
    std::vector<std::vector<int>> dist;
    int n, m;
    fin >> n;
    m = n*n;
    Graph G(readAdCostListFromMatrix(n), n, m, 1);
    dist = G.getFloydWarshall();
    for(int i = 1; i<=n; ++i){
        for(int j = 1; j<=n; ++j)
            fout << dist[i][j] << " ";
        fout << "\n";
    }
}

void treeSizeHelper(){
    int depth;
    int n;
    fin >> n;
    Graph G(readAdList(n, n-1, 0), n, n-1);
    depth = G.getTreeSize();
    fout << depth;
}

void maxFlowHelper(){
    int n,m;
    fin >> n >> m;
    Graph G(readAdCostList(n,m), n, m, 1);
    fout << G.getFlow(1, n);
}

void dijkstraHelper(){
    int n,m;
    std::vector<int> dist;
    fin >> n >> m;
    Graph G(readAdCostList(n,m), n, m, 1);
    dist = G.getDijkstraDistance(1);
    for(int i = 2; i<=n; ++i)
        if( dist[i] == 100000000 )
            fout << 0 << " ";
        else fout << dist[i] << " ";
}

void mstHelper(){
    int n, m;
    int total = 0;
    std::vector<std::tuple<int, int, int>> result;
    fin >> n >> m;
    Graph G(readAdCostList(n,m), n, m, 0);
    result = G.getMST();
    for(int i = 0; i<result.size(); ++i)
        total += std::get<0>(result[i]);
    fout << total << "\n" << result.size() << "\n";
    for(int i = 0; i<result.size(); ++i)
        fout << std::get<1>(result[i]) << " " << std::get<2>(result[i]) << "\n";
}

int main(){
    dijkstraHelper();
    return 0;
}
