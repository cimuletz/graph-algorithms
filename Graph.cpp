

#include <queue>
#include <tuple>
#include <algorithm>
#include "Graph.h"
Graph::Graph(const std::vector<std::pair<int,int>> &ad, int nodes, int edges, int directed) : nodes(nodes), edges(edges), directed(directed) {
    Ad.resize(nodes+1);
    for(int i = 0; i<edges; ++i){
        Ad[ad[i].first].push_back(ad[i].second);
        if(directed == 0)
            Ad[ad[i].second].push_back(ad[i].first);
    }
}

Graph::Graph(const std::vector<std::vector<edgeCost>> &ad, int nodes, int edges,int directed) : AdCost(ad), nodes(nodes), edges(edges), directed(directed){}

Graph::Graph(const std::vector<std::pair<int,int>> &ad, const std::vector<int> &cost, int nodes, int edges,int directed) : nodes(nodes), edges(edges), directed(directed) {
    edgeCost aux;
    AdCost.resize(nodes+1);
    for(int i = 0; i<edges; ++i){
        aux.c = cost[i];
        aux.x = ad[i].second;
        AdCost[ad[i].first].push_back(aux);
        if(directed == 0){
            aux.x = ad[i].first;
            AdCost[ad[i].second].push_back(aux);
        }
    }
}

void Graph::DFS(int node, std::vector<int> &vis, int nrComp){
    vis[node] = nrComp;
    int n = Ad[node].size();
    for(int i = 0; i<n; ++i)
        if(vis[Ad[node][i]] == 0){
            DFS(Ad[node][i], vis, nrComp);
        }
}

std::vector<int> Graph::getComponents() {
    std::vector<int> vis(nodes + 1, 0);
    int ct = 1;
    for(int i = 1;i <= nodes; ++i)
        if(vis[i] == 0){
            DFS(i, vis, ct);
            ct++;
        }
    return vis;
}

void Graph::BFS(int node, std::vector<int> &vis){
    int parent;
    std::queue<int> Q;
    Q.push(node);
    vis[node] = 0;
    while(!Q.empty()){
        parent = Q.front();
        Q.pop();
        int n = Ad[node].size();

        for(int i = 0; i<n; ++i)
            if(vis[Ad[parent][i]] == -1){
                Q.push(Ad[parent][i]);
                vis[Ad[parent][i]] = vis[parent] + 1;
            }
    }
}

std::vector<int> Graph::getBFSPathLength(int node){
    std::vector<int> vis(nodes+1, -1);
    vis[node] = 0;
    BFS(node, vis);
    return vis;
}

void Graph::topSortDFS(int node, std::vector<int> &vis, std::stack<int> &s){
    vis[node] = 1;
    int n = Ad[node].size();
    for(int i = 0; i<n; ++i)
        if(vis[Ad[node][i]] == 0)
            topSortDFS(Ad[node][i], vis, s);
    s.push(node);
}

std::stack<int> Graph::getTopSort(){
    std::vector<int> vis(nodes+1, 0);
    std::stack<int> s;
    for(int i = 1; i<=nodes; ++i)
        if(!vis[i])
            topSortDFS(i, vis, s);
    return s;
}

int Graph::findMST(int x, std::vector<int> &disjSets){
    int root = x, aux;
    while(disjSets[root] != root)
        root = disjSets[root];
    while(disjSets[x] != root){
        aux = disjSets[x];
        disjSets[x] = root;
        x = aux;
    }
    return root;
}

void Graph::unionMST(int x, int y, std::vector<int> &disjSets, std::vector<int> &sizes){
    int rootx = findMST(x, disjSets);
    int rooty = findMST(y, disjSets);
    if(sizes[rootx] >= sizes[rooty]){
        sizes[rootx] += sizes[rooty];
        disjSets[rooty] = rootx;
    }
    else{
        sizes[rooty] += sizes[rootx];
        disjSets[rootx] = rooty;
    }
}

std::vector<std::tuple<int, int, int>> Graph::getMST(){
    std::vector<std::tuple<int, int, int>> result;
    std::vector< std::tuple<int, int, int> > edges;
    std::vector<int> disjSets, sizes; //vectors for disjoint sets and their respective sizes
    std::vector<int> sol;
    int solSize = 0;
    int n, x, y;
    disjSets.push_back(0);
    sizes.push_back(0);
    for(int i = 1; i<=nodes; ++i){
        disjSets.push_back(i);
        sizes.push_back(1);
        int n = AdCost[i].size();
        for(int j = 0; j<n; ++j)
            edges.push_back(std::make_tuple(AdCost[i][j].c, i, AdCost[i][j].x));
    }
    std::sort(edges.begin(), edges.end());
    int nr = edges.size();
    for(int i = 0; i < nr && solSize < n-1; ++i){
        x = std::get<1>(edges[i]);
        y = std::get<2>(edges[i]);
        if( findMST(x,disjSets) != findMST(y, disjSets) ){
            unionMST(x,y, disjSets, sizes);
            result.push_back(edges[i]);
            ++solSize;
        }
    }
    return result;
}

void Graph::dfsSCC(int node,int &idInc, std::vector<int> &id, std::vector<int> &low, std::vector<bool> &onStack, std::stack<int> &s, std::vector< std::vector<int>> &scc){
    s.push(node);
    onStack[node] = true;
    id[node] = low[node] = idInc++;
    int n = Ad[node].size();
    for(int i = 0; i<n; ++i){
        int next = Ad[node][i];
        if(id[next] == 0)
            dfsSCC(next, idInc, id, low, onStack, s, scc);
        if(onStack[next] == true)
            low[node] = std::min(low[node], low[next]);
    }
    if(id[node] == low[node]){
        std::vector<int> comp; ///next component in scc vector
        while(!s.empty()){
            int nod = s.top();
            s.pop();
            onStack[nod] = false;
            comp.push_back(nod);
            low[nod] = id[node];
            if(nod == node)
                break;
        }
        scc.push_back(comp);
    }
}

std::vector< std::vector<int>> Graph::getSCC(){
    int idInc = 1; //increment node id
    std::vector<int> id(nodes+1);
    std::vector<int> low(nodes+1);
    std::vector<bool> onStack(nodes+1, 0);
    std::vector< std::vector<int>> scc; //vector of components
    std::stack<int> s;
    for(int i = 1; i<=nodes; ++i)
        id[i] = low[i] = 0;
    for(int i = 1; i<=nodes; ++i)
        if(id[i] == 0)
            dfsSCC(i,idInc, id, low, onStack, s, scc);
    return scc;
}

void Graph::dfsBCC(int node, int parent, int &time, std::vector<int> &dq, std::vector<int> &low, std::vector<int> &disc, std::vector<std::vector<int>> &bcc){
    int w;
    disc[node] = low[node] = ++time;
    int n = Ad[node].size();
    for(int i = 0; i<n; ++i){

        w = Ad[node][i];

        if(w == parent)
            continue;

        if(disc[w] != 0)
            low[node] = std::min(low[node], disc[w]);

        else{
            dq.push_back(w);
            dfsBCC(w, node, time, dq, low, disc, bcc);

            low[node] = std::min(low[w], low[node]);

            if(low[w] >= disc[node]){
                dq.push_back(node);
                std::vector<int> nextComp;
                while(!dq.empty()){
                    int next = dq.back();
                    dq.pop_back();
                    nextComp.push_back(next);
                    if(next == w)
                        break;
                }
                bcc.push_back(nextComp);
            }
        }
    }

}

std::vector< std::vector<int>> Graph::getBCC(){
    int time = 0;
    std::vector<std::vector<int>> bcc;
    std::vector<int> dq; ///deque for components in next bcc
    std::vector<int> low(nodes), disc(nodes);
    dfsBCC(1,0,time, dq, low, disc, bcc);
    return bcc;
}

void Graph::bellman(bool &cycle, std::vector<int> &dist, int node){
    edgeCost aux;
    std::vector<int> countQ(nodes+1, 0); //number of times a node has been added to queue
    std::vector<bool> inQ(nodes+1, 0); //node in queue

    std::deque<int> q;
    int currentNode;
    int w,c;
    q.push_back(node);
    while( !q.empty() && !cycle ){
        currentNode = q.front();
        q.pop_front();
        inQ[currentNode] = 0;
        int n = Ad[currentNode].size();
        for(int i = 0; i<n; ++i){
            aux = AdCost[currentNode][i];
            w = aux.x;
            c = aux.c;
            if( dist[w] > dist[currentNode] + c ){
                dist[w] = dist[currentNode] + c;
                //cout << currentNode << " " << w << " " << dist[w] << " " << c << "\n";
                if( !inQ[w] ){
                    if( countQ[w] > nodes ){ //if node has updated the distance more than n times, this means it's part of a negative cycle
                        cycle = 1;
                        break;
                    }
                    else{
                        inQ[w] = 1;
                        countQ[w]++;
                        q.push_back(w);
                    }
                }
            }
        }
    }
}

std::vector<int> Graph::getBellmanDistance(int node){
    bool cycle = false;
    std::vector<int> dist; //distance vector;
    dist.resize(nodes+1, inf);
    dist[node] = 0;
    bellman(cycle, dist, node);
    if(cycle)
        dist[node] = inf;
    return dist;

}

void Graph::dijkstra(std::vector<int> &dist, int node){
    edgeCost aux;
    std::priority_queue< std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>>> heap;
    std::vector<bool>inH(nodes+1, 0);
    int currentNode, w, c;
    heap.push({0,node});

    while(!heap.empty()){
        currentNode = heap.top().second;
        heap.pop();
        if(inH[currentNode])
            continue;
        else inH[currentNode] = true;
        int n = AdCost[currentNode].size();
        for(int i = 0; i < n; ++i){
            aux = AdCost[currentNode][i];
            w = aux.x;
            c = aux.c;

            if( dist[currentNode] + c < dist[w] ){
                dist[w] = dist[currentNode] + c;
                heap.push({ dist[w], w });
            }
        }
    }
}

std::vector<int> Graph::getDijkstraDistance(int node){
    std::vector<int> dist; //distance vector;
    dist.resize(nodes+1, inf);
    dist[1] = 0;
    dijkstra(dist, node);
    return dist;
}


std::vector<std::vector<int>> Graph::getFloydWarshall(){
    std::vector<std::vector<int>> dist(nodes+1);
    for(int i = 1; i<=nodes; ++i){
        dist[i].resize(nodes+1);
        int n = AdCost[i].size();
        for(int j = 0; j < n; ++j){
            edgeCost aux = AdCost[i][j];
            dist[i][aux.x] = aux.c;
            if(aux.c == 0 && i!=aux.x)
                dist[i][aux.x] = inf;
        }
    }
    for(int k = 1; k<=nodes; ++k)
        for(int i = 1; i<=nodes; ++i)
            for(int j = 1; j<=nodes; ++j)
                if( dist[i][j] > dist[i][k] + dist[k][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
    return dist;
}

int Graph::bfsFlow(int src, int dest, std::vector<std::vector<int>> &edgeFlow, std::vector<int> &vis, std::vector<int> &bottleneck){
    std::queue<int> Q;
    vis[src] = -1;
    std::fill(vis.begin(), vis.end(), 0);
    Q.push(src);
    bottleneck[src] = inf;

    while( !Q.empty() ){
        int node = Q.front();
        Q.pop();
        int n = AdCost[node].size();
        for(int i = 0; i < n; ++i){
            int w = AdCost[node][i].x;
            int c = AdCost[node][i].c;

            if( vis[w] == 0 ){
                if( c - edgeFlow[node][w] > 0 ){
                    vis[w] = node;
                    bottleneck[w] = std::min(bottleneck[node], c - edgeFlow[node][w] );
                    if( w == dest )
                        return bottleneck[w];
                    Q.push(w);
                }
            }
        }
    }
    return 0;
}

int Graph::getFlow(int src, int dest){
    int totalFlow = 0, flow = -1;
    //flow passing through one edge
    std::vector<std::vector<int>> edgeFlow(nodes+1);
    //visited vector with parent
    std::vector<int> vis(nodes+1, 0);
    //vector for bottleneck of each node in the path
    std::vector<int> bottleneck(nodes+1, 0);

    for(int i = 1; i<=nodes; ++i){
        edgeFlow[i].resize(nodes+1);
        std::fill(edgeFlow[i].begin(), edgeFlow[i].end(), 0);
    }
    while( (flow = bfsFlow( src, dest, edgeFlow, vis, bottleneck))){
        //flow will be 0 when no more paths exist
        if(flow == 0)
            break;
        totalFlow += flow;
        int w = dest;
        while( w != src ){
            int z = vis[w];
            edgeFlow[z][w] += flow;
            edgeFlow[w][z] -= flow;
            w = z;
        }
    }
    return totalFlow;
}

void Graph::dfsTS(std::vector<int> &vis, int node, int parent, int &last, int &level){
    vis[node] = vis[parent] + 1;
    if(vis[node] > level){
        last = node;
        level = vis[node];
    }
    int n = Ad[node].size();
    for(int i = 0; i < n; ++i)
    {
        int w = Ad[node][i];
        if(!vis[w])
            dfsTS(vis, w, node, last, level);
    }
}

int Graph::getTreeSize(){
    int level = 0;
    int last;
    std::vector<int> vis(nodes+1, 0);
    dfsTS(vis, 1, 0, last, level);
    std::fill(vis.begin(), vis.end(), 0);

    dfsTS(vis, last, 0, last, level);
    return vis[last];

}