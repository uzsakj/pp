// PP_beadando.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream> 
#include <list> 
#include <limits.h> 
#include <vector>
using namespace std;

class Graph
{
protected:
	int V;//nodok száma
	list<int> *adj; // pointer a szomszédsági lista tombre
	bool isCyclicUtil(int v, bool visited[], int parent);
public:
	Graph(int V);   // Constructor 
	void addEdge(int v, int w);   // to add an edge to graph 
	bool isTree();   // returns true if graph is tree 
		// Method to check if this graph is Eulerian or not 
	int isEulerian();

	// Method to check if all non-zero degree vertices are connected 
	bool isConnected();

	// Function to do DFS starting from v. Used in isConnected(); 
	void DFSUtil(int v, vector<bool>& visited);

	int countEdges();
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
void Graph::DFSUtil(int v, vector<bool>& visited)
{
	// Mark the current node as visited and print it 
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex 
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i])
			DFSUtil(*i, visited);
}
// Method to check if all non-zero degree vertices are connected. 
// It mainly does DFS traversal starting from 
bool Graph::isConnected()
{
	// Mark all the vertices as not visited 
	//bool visited[V];
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
// Returns count of edge in undirected graph 
int Graph::countEdges()
{
	int sum = 0;

	//traverse all vertex 
	for (int i = 0; i < V; i++)

		// add all edge that are linked to the 
		// current vertex 
		sum += adj[i].size();


	// The count of edge is always even because in 
	// undirected graph every edge is connected 
	// twice between two vertices 
	return sum / 2;
}

int main()
{
	Graph g1(5);
	g1.addEdge(1, 0);
	g1.addEdge(0, 2);
	g1.addEdge(0, 3);
	g1.addEdge(3, 4);
	g1.isTree() ? cout << "Graph is Tree\n" :
		cout << "Graph is not Tree\n";

	Graph g2(5);
	g2.addEdge(1, 0);
	g2.addEdge(0, 2);
	g2.addEdge(2, 1);
	g2.addEdge(0, 3);
	g2.addEdge(3, 4);
	g2.isTree() ? cout << "Graph is Tree\n" :
		cout << "Graph is not Tree\n";

	Graph g3(5);
	g3.addEdge(1, 0);
	g3.addEdge(0, 2);
	g3.addEdge(2, 1);
	g3.addEdge(0, 3);
	g3.addEdge(3, 4);
	testeuler(g3);

	Graph g4(5);
	g4.addEdge(1, 0);
	g4.addEdge(0, 2);
	g4.addEdge(2, 1);
	g4.addEdge(0, 3);
	g4.addEdge(3, 4);
	g4.addEdge(4, 0);
	testeuler(g4);

	Graph g5(5);
	g5.addEdge(1, 0);
	g5.addEdge(0, 2);
	g5.addEdge(2, 1);
	g5.addEdge(0, 3);
	g5.addEdge(3, 4);
	g5.addEdge(1, 3);
	testeuler(g5);

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
	cout << g6.countEdges() << endl;

	return 0;
}


