// PP_beadando.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream> 
#include <list> 
#include <limits.h> 
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <omp.h>
//#include <bits/stdc++.h> 
using namespace std;


class Graph
{
protected:
	int V;//number of nodes
	list<int> *adj; // pointer to adjency list
	bool isCyclicUtil(int v, bool visited[], int parent);
	bool isCyclicUtil_omp(int v, bool visited[], int parent);
public:
	// Constructor 
	Graph(int V);   
	// to add an edge to graph 
	void addEdge(int v, int w); 
	//Method to execute Bfs serch in graph
	bool BFS(int src, int dest, int v, vector<int>& pred, vector<int>& dist);
	bool BFS_omp(int src, int dest, int v, vector<int>& pred, vector<int>& dist);
	// returns true if graph is tree 
	bool isTree();   
	bool isTree_omp();
	// Method to check if this graph is Eulerian or not 
	int isEulerian();
	int isEulerian_omp();
	// Method to check if all non-zero degree vertices are connected 
	bool isConnected();
	bool isConnected_omp();
	//method to determine shortest path in graph between nodes s and dest
	void printShortestDistance(int s, int dest, int v);
	void printShortestDistance_omp(int s, int dest, int v);
	// Function to do DFS starting from v. Used in isConnected(); 
	void DFSUtil(int v, vector<bool>& visited);
	void DFSUtil_omp(int v, vector<bool>& visited);
	//Method to determine number of edges for a  node and total sum of them
	void countEdges();
	void countEdges_omp();
	// utility function to display edges
	void printEdges(list<int> nodes, int sum);
};
Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w); // Add w to v’s list. 
	adj[w].push_back(v); // Add v to w’s list. 
}
bool Graph::BFS( int src, int dest, int v,vector<int>& pred, vector<int>& dist)
{
	// a queue to maintain queue of vertices whose 
	// adjacency list is to be scanned as per normal 
	// DFS algorithm 
	list<int> queue;

	// boolean vector visited which stores the 
	// information whether ith vertex is reached 
	// at least once in the Breadth first search 
	vector<bool> visited;

	// initially all vertices are unvisited 
	// so v[i] for all i is false 
	// and as no path is yet constructed 
	// dist[i] for all i set to infinity 
	
	for (int i = 0; i < V;i++) {
		visited.push_back(false);
		dist.push_back(INT_MAX);
		pred.push_back(-1);
	}

	// now source is first to be visited and 
	// distance from source to itself should be 0 

	visited.at(src) = true;
	dist.at(src) = 0;
	queue.push_back(src);

	// standard BFS algorithm 
	while (!queue.empty()) {
		int u = queue.front();
		queue.pop_front();
		list<int>::iterator i;
		
		for (i = adj[u].begin(); i != adj[u].end();i++) {
			if (visited.at(*i) == false) {
				visited[*i] = true;
				dist.at(*i) = (dist.at(u) + 1);
				pred.at(*i) = u;
				queue.push_back(*i);

				// We stop BFS when we find destination. 
				if (*i == dest)
					return true;
			}

		}
	}

	return false;
}

bool Graph::BFS_omp(int src, int dest, int v, vector<int>& pred, vector<int>& dist)
{
	// a queue to maintain queue of vertices whose 
	// adjacency list is to be scanned as per normal 
	// DFS algorithm 
	list<int> queue;

	// boolean vector visited which stores the 
	// information whether ith vertex is reached 
	// at least once in the Breadth first search 
	vector<bool> visited;

	// initially all vertices are unvisited 
	// so v[i] for all i is false 
	// and as no path is yet constructed 
	// dist[i] for all i set to infinity 
	
//#pragma omp parallel for
	for (int i = 0; i < V; i++) {
		visited.push_back(false);
		dist.push_back(INT_MAX);
		pred.push_back(-1);
	}

	// now source is first to be visited and 
	// distance from source to itself should be 0 

	visited.at(src) = true;
	dist.at(src) = 0;
	queue.push_back(src);

	// standard BFS algorithm 
	while (!queue.empty()) {
		int u = queue.front();
		queue.pop_front();
		list<int>::iterator i;
//#pragma omp parallel for
		for (i = adj[u].begin(); i != adj[u].end(); i++) {
			if (visited.at(*i) == false) {
				visited[*i] = true;
				dist.at(*i) = (dist.at(u) + 1);
				pred.at(*i) = u;
				queue.push_back(*i);

				// We stop BFS when we find  destination. 
				if (*i == dest)
					return true;
			}

		}
	}

	return false;
}

void Graph::printShortestDistance( int s,int dest, int v)
{
	// predecessor[i] array stores predecessor of 
	// i and distance array stores distance of i 
	// from s 
	vector<int> pred, dist;
	
	if (BFS( s, dest, v, pred, dist) == false) {
		cout << "Given source and destination"
			<< " are not connected";
		return;
	}
	// vector path stores the shortest path 
	vector<int> path;
	//using crawl as the index of the node
	int crawl = dest;
	path.push_back(crawl);
	while (pred.at(crawl) != -1) {
		path.push_back(pred.at(crawl));
		crawl = pred.at(crawl);
	}

	// distance from source is in distance array 
	cout << "Shortest path length is : "
		<< dist.at(dest);

	// printing path from source to destination 
	cout << "\nPath is:";
	for (int i = path.size() - 1; i >= 0; i--)
	{
		cout << path[i] << " ";
	}
	cout << endl;
}
void Graph::printShortestDistance_omp(int s, int dest, int v)
{
	// predecessor[i] array stores predecessor of 
	// i and distance array stores distance of i 
	// from s 
	vector<int> pred, dist;

	if (BFS_omp(s, dest, v, pred, dist) == false) {
		cout << "Given source and destination"
			<< " are not connected";
		return;
	}
	// vector path stores the shortest path 
	vector<int> path;
	int crawl = dest;
	path.push_back(crawl);
	while (pred.at(crawl) != -1) {
		path.push_back(pred.at(crawl));
		crawl = pred.at(crawl);
	}

	// distance from source is in distance array 
	cout << "Shortest path length is : "
		<< dist.at(dest);

	// printing path from source to destination 
	cout << "\nPath is:";
	for (int i = path.size() - 1; i >= 0; i--)
	{
		cout << path[i] << " ";
	}
	cout << endl;
}

void Graph::DFSUtil(int v, vector<bool>& visited)
{
	// Mark the current node as visited 
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex 
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i])
			DFSUtil(*i, visited);
}

void Graph::DFSUtil_omp(int v, vector<bool>& visited)
{
	// Mark the current node as visited
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex 
	list<int>::iterator i;
//#pragma omp parallel for lehetne omp 3.0 tol
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
	{
		if (!visited[*i])
		{
			DFSUtil_omp(*i, visited);
		}
	}
}
// Method to check if all non-zero degree vertices are connected. 
// It mainly does DFS traversal starting from 
bool Graph::isConnected()
{
	// Mark all the vertices as not visited 
	
	vector<bool> visited(V);
	int i;
	for (i = 0; i < V; i++)
		visited[i] = false;

	// Find a vertex with non-zero degree 
	for (i = 0; i < V; i++)
		if (adj[i].size() != 0)
			break;

	// If there are no edges in the graph, return true 
	if (i == V)
		return true;

	// Start DFS traversal from a vertex with non-zero degree 
	DFSUtil(i, visited);

	// Check if all non-zero degree vertices are visited 
	for (i = 0; i < V; i++)
		if (visited[i] == false && adj[i].size() > 0)
			return false;

	return true;
}

bool Graph::isConnected_omp()
{
	// Mark all the vertices as not visited 
	
	vector<bool> visited(V);
	int i;
	for (i = 0; i < V; i++)
		visited[i] = false;

	// Find a vertex with non-zero degree 
#pragma omp parallel for
	for (i = 0; i < V; i++)
	{
		if (adj[i].size() != 0)
		{
			break;
		}
	}
	// If there are no edges in the graph, return true 
	if (i == V)
		return true;

	// Start DFS traversal from a vertex with non-zero degree 
	DFSUtil(i, visited);

	// Check if all non-zero degree vertices are visited
//#pragma omp parallel for lehetne omp 3.0 tol
	for (i = 0; i < V; i++)
		if (visited[i] == false && adj[i].size() > 0)
			return false;

	return true;
}
/* The function returns one of the following values
   0 --> If grpah is not Eulerian
   1 --> If graph has an Euler path (Semi-Eulerian)
   2 --> If graph has an Euler Circuit (Eulerian)  */
int Graph::isEulerian()
{
	// Check if all non-zero degree vertices are connected 
	if (isConnected() == false)
		return 0;

	// Count vertices with odd degree 
	int odd = 0;
	for (int i = 0; i < V; i++)
		if (adj[i].size() & 1)
			odd++;

	// If count is more than 2, then graph is not Eulerian 
	if (odd > 2)
		return 0;

	// If odd count is 2, then semi-eulerian. 
	// If odd count is 0, then eulerian 
	// Note that odd count can never be 1 for undirected graph 
	return (odd) ? 1 : 2;
}

int Graph::isEulerian_omp()
{
	// Check if all non-zero degree vertices are connected 
	if (isConnected() == false)
		return 0;

	// Count vertices with odd degree 
	int odd = 0;
#pragma omp parallel for
	for (int i = 0; i < V; i++)
		if (adj[i].size() & 1)
			odd++;

	// If count is more than 2, then graph is not Eulerian 
	if (odd > 2)
		return 0;

	// If odd count is 2, then semi-eulerian. 
	// If odd count is 0, then eulerian 
	// Note that odd count can never be 1 for undirected graph 
	return (odd) ? 1 : 2;
}
// Function to run test cases 
void testeuler(Graph &g)
{
	int res = g.isEulerian();
	if (res == 0)
		cout << "graph is not Eulerian\n";
	else if (res == 1)
		cout << "graph has a Euler path\n";
	else
		cout << "graph has a Euler cycle\n";
}

void testeuler_omp(Graph &g)
{
	int res = g.isEulerian_omp();
	if (res == 0)
		cout << "graph is not Eulerian\n";
	else if (res == 1)
		cout << "graph has a Euler path\n";
	else
		cout << "graph has a Euler cycle\n";
}


void testTree(Graph &g)
{
	int res = g.isTree();
	if (res == 1)
	{
		cout << "Graph is Tree " << endl;
	}
	if (res == 0)
	{
		cout << "Graph is not Tree " << endl;
	}

}

void testTree_omp(Graph &g)
{
	int res = g.isTree_omp();
	if (res == 1)
	{
		cout << "Graph is Tree " << endl;
	}
	if (res == 0)
	{
		cout << "Graph is not Tree " << endl;
	}

}


// A recursive function that uses visited[] and parent to 
// detect cycle in subgraph reachable from vertex v. 
bool Graph::isCyclicUtil(int v, bool visited[], int parent)
{
	// Mark the current node as visited 
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex 
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
	{
		// If an adjacent is not visited, then recur for  
		// that adjacent 
		if (!visited[*i])
		{
			if (isCyclicUtil(*i, visited, v))
				return true;
		}

		// If an adjacent is visited and not parent of current 
		// vertex, then there is a cycle. 
		else if (*i != parent)
			return true;
	}
	return false;
}

bool Graph::isCyclicUtil_omp(int v, bool visited[], int parent)
{
	// Mark the current node as visited 
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex 
	list<int>::iterator i;
//#pragma omp parallel for lehetne omp 3.0 tol
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
	{
		// If an adjacent is not visited, then recur for  
		// that adjacent 
		if (!visited[*i])
		{
			if (isCyclicUtil_omp(*i, visited, v))
				return true;
		}

		// If an adjacent is visited and not parent of current 
		// vertex, then there is a cycle. 
		else if (*i != parent)
			return true;
	}
	return false;
}

// Returns true if the graph is a tree, else false. 
bool Graph::isTree()
{
	// Mark all the vertices as not visited and not part of  
	// recursion stack 
	bool *visited = new bool[V];
	for (int i = 0; i < V; i++)
		visited[i] = false;

	// The call to isCyclicUtil serves multiple purposes. 
	// It returns true if graph reachable from vertex 0  
	// is cyclcic. It also marks all vertices reachable  
	// from 0. 
	if (isCyclicUtil(0, visited, -1))
		return false;

	// If we find a vertex which is not reachable from 0  
	// (not marked by isCyclicUtil(), then we return false 
	for (int u = 0; u < V; u++)
		if (!visited[u])
			return false;

	return true;
}

bool Graph::isTree_omp()
{
	// Mark all the vertices as not visited and not part of  
	// recursion stack 
	bool *visited = new bool[V];
#pragma omp parallel for
	for (int i = 0; i < V; i++)
		visited[i] = false;

	// The call to isCyclicUtil serves multiple purposes. 
	// It returns true if graph reachable from vertex 0  
	// is cyclcic. It also marks all vertices reachable  
	// from 0. 
	if (isCyclicUtil(0, visited, -1))
		return false;

	// If we find a vertex which is not reachable from 0  
	// (not marked by isCyclicUtil(), then we return false 
//#pragma omp parallel for
	for (int u = 0; u < V; u++)
		if (!visited[u])
			return false;

	return true;
}
// Returns count of edge in undirected graph 
void Graph::countEdges()
{
	int sum = 0;
	list<int> nodes;
	//traverse all vertex 
	for (int i = 0; i < V; i++)
	{
		// add all edge that are linked to the 
		// current vertex 
		// add edges to to current node
		sum += adj[i].size();
		nodes.push_back(adj[i].size());
	}
	// The count of edge is always even because in 
	// undirected graph every edge is connected 
	// twice between two vertices 
	printEdges(nodes, sum/2);
}

void Graph::countEdges_omp()
{
	int sum = 0;
	list<int> nodes;
	//traverse all vertex 
#pragma omp parallel for
	for (int i = 0; i < V; i++)
	{
		// add all edge that are linked to the 
		// current vertex 
		// add edges to to current node
		sum += adj[i].size();
		nodes.push_back(adj[i].size());
	}
	// The count of edge is always even because in 
	// undirected graph every edge is connected 
	// twice between two vertices 
	printEdges(nodes, sum / 2);
}

void Graph::printEdges(list<int> nodes,int sum)
{
	list<int>::iterator i;
	int node_num = 0;
	for (i = nodes.begin(); i != nodes.end(); i++)
	{
		cout << "Node "<<node_num <<" has "<< *i<<" edges"<<endl;
		node_num ++;
	}
	cout << "total number of edges: " << sum << endl;
}

int main()
{

	Graph g6(9);
	g6.addEdge(0, 1);
	g6.addEdge(0, 7);
	g6.addEdge(1, 2);
	g6.addEdge(1, 7);
	g6.addEdge(2, 3);
	g6.addEdge(2, 8);
	g6.addEdge(2, 5);
	g6.addEdge(3, 4);
	g6.addEdge(3, 5);
	g6.addEdge(4, 5);
	g6.addEdge(5, 6);
	g6.addEdge(6, 7);
	g6.addEdge(6, 8);
	g6.addEdge(7, 8);
	double time1=0;
	for (int i = 0; i < 10; i++)
	{
		double t1 = omp_get_wtime();
		g6.countEdges();
		testeuler(g6);
		testTree(g6);
		g6.printShortestDistance( 0, 5, 9);
		t1 = omp_get_wtime() - t1;
		time1 += t1;
	}
	cout << "sequential took time :" << time1/10 <<" s"<< endl;
	double time2 = 0;
	for (int i = 0; i < 10; i++)
	{
		double t2 = omp_get_wtime();
		g6.countEdges_omp();
		testeuler_omp(g6);
		testTree_omp(g6);
		g6.printShortestDistance_omp(0, 5, 9);
		t2 = omp_get_wtime() - t2;
		time2 += t2;
	}
	cout << "sequential took time :" << time1 / 10 << " s" << endl;
	cout << "parallel took time :" << time2 / 10 << " s" << endl;



	return 0;
}


