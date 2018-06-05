#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>

template <typename T>
class weighted_graph {

	private:

	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.

	std::map<T,std::map<T,int>*> vertices;
	// A map of vertexes and edges. The edges are pointers to another map of vertexes are weights.
		
	int _num_vertices;
	int _num_edges;
	int _total_weight;
	
	int remove_edge_core(const T&, const T&); //Disconnects one side of an edge. Must be called twice. Returns the weight of the edge disconnected for recording.
	int set_edge_weight_core(const T&, const T&, const int&); // Sets weight of one side of an edge. Must be called twice. Returns the original weight of the edge for recording.
	void depth_first_core(const T&, std::map<T,bool>*, std::vector<T>*); // Recursive core function for dft
	
	//The graph_iterator class provides an iterator
	//over the vertices of the graph.
	//This is one of the harder parts, so if you're
	//not too comfortable with C++ leave this for last.
	//If you are, there are many ways of doing this,
	//as long as it passes the tests, it's okay.
	class graph_iterator {

		private:

		//You may need data members here.
		std::vector<T> vertices;
		size_t _pos;

		public:
			graph_iterator(const weighted_graph &);
			graph_iterator(const weighted_graph &, size_t);
			~graph_iterator();
			graph_iterator operator=(const graph_iterator&);
			bool operator==(const graph_iterator&) const;
			bool operator!=(const graph_iterator&) const;
			graph_iterator operator++();
			graph_iterator operator++(int);
			const T operator*();
			const T* operator->();
	};

	//The neighbour_iterator class provides an iterator
	//over the neighbours of a given vertex. This is
	//probably harder (conceptually) than the graph_iterator.
	//Unless you know how iterators work.
	class neighbour_iterator {

		private:

		//You may need data members here.
		std::vector<std::pair<const T,int>*> neighbours;
		size_t _pos;

		public:
			neighbour_iterator(const weighted_graph &, const T&);
			neighbour_iterator(const weighted_graph &, const T&, size_t);
			~neighbour_iterator();
			neighbour_iterator operator=(const neighbour_iterator& it);
			bool operator==(const neighbour_iterator&) const;
			bool operator!=(const neighbour_iterator&) const;
			neighbour_iterator operator++();
			neighbour_iterator operator++(int);			
			const std::pair<T, int> operator*();
			const std::pair<const T, int>* operator->();
	};

	public:


	weighted_graph(); //A constructor for weighted_graph. It should start empty.
	~weighted_graph(); //A destructor. Depending on how you do things, this may
					   //not be necessary.

	bool are_adjacent(const T&, const T&) const; //Returns true if the two vertices are
														   //adjacent, false otherwise.
	bool has_vertex(const T&) const; //Returns true if the passed in vertex is 
										  //a vertex of the graph, false otherwise.

	void add_vertex(const T&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const T&, const T&, const int&); //Adds an edge between the two vertices
															 //with the given weight (as an int).

	void remove_vertex(const T&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const T&, const T&); //Removes the edge between the two vertices, if it exists.
	void set_edge_weight(const T&, const T&, const int&); //Changes the edge weight between the two
																	//vertices to the new weight (the int).

	int get_edge_weight(const T&, const T&) const; //Returns the weight on the edge between the two vertices.
	int degree(const T&) const; //Returns the degree of the vertex.
	int weighted_degree(const T&); //Returns the sum of the weights on all the edges incident to the vertex.
	int num_vertices() const; //Returns the total number of vertices in the graph.
	int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight).
	int total_weight(); //Returns the sum of all the edge weights in the graph.

	std::vector<T> get_vertices() const; //Returns a vector containing all the vertices.
	std::vector<T> get_neighbours(const T&) const; //Returns a vector containing the neighbours of the given vertex.

	graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
	graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

	neighbour_iterator neighbours_begin(const T&); //Returns a neighbour_iterator pointing to the start
														//of the neighbour set for the given vertex.
	neighbour_iterator neighbours_end(const T&); //Returns a neighbour_iterator pointing to one-past-the-end
													  //of the neighbour set for the given vertex.

	std::vector<T> depth_first(const T&); //Returns the vertices of the graph in the order they
													//are visited in by a depth-first traversal starting at
													//the given vertex.
	std::vector<T> breadth_first(const T&); //Returns the vertices of the graph in the order they
													  //are visisted in by a breadth-first traversal starting
													  //at the given vertex.

	weighted_graph<T> mst(); //Returns a minimum spanning tree of the graph.

};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename T> 
weighted_graph<T>::graph_iterator::graph_iterator(const weighted_graph & g)
{
	vertices = g.get_vertices();
	_pos = 0;
}

template <typename T> 
weighted_graph<T>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos)
{
	vertices = g.get_vertices();
	_pos = start_pos;
}

template <typename T> weighted_graph<T>::graph_iterator::~graph_iterator(){}

template <typename T> 
typename weighted_graph<T>::graph_iterator weighted_graph<T>::graph_iterator::operator=(const graph_iterator& it)
{
	if (*this == it) return * this; // nothing to do
	this->vertices = it.vertices;
	this->_pos = it._pos;
	return *this;
}

template <typename T> bool weighted_graph<T>::graph_iterator::operator==(const graph_iterator& it) const {return this->vertices == it.vertices && this->_pos == it._pos;}
template <typename T> bool weighted_graph<T>::graph_iterator::operator!=(const graph_iterator& it) const { return !(*this == it); }

template <typename T> 
typename weighted_graph<T>::graph_iterator weighted_graph<T>::graph_iterator::operator++()
{	
	_pos++; // increment before returning
	return * this;
}

template <typename T> 
typename weighted_graph<T>::graph_iterator weighted_graph<T>::graph_iterator::operator++(int)
{
	auto g = * this; // return before incrementing
	_pos++;
	return g; 
}

template <typename T> const T weighted_graph<T>::graph_iterator::operator*() {return vertices.at(_pos);}
template <typename T> const T* weighted_graph<T>::graph_iterator::operator->() { return &vertices.at(_pos); }

template <typename T> 
weighted_graph<T>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const T& u)
{
	for (T neighbour : g.get_neighbours(u))
	{
		neighbours.push_back(new std::pair<const T,int>(neighbour, g.get_edge_weight(u, neighbour)));
	}
	_pos = 0;	
}

template <typename T>
weighted_graph<T>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const T& u, size_t start_pos) 
{
	for (T neighbour : g.get_neighbours(u))
	{
		neighbours.push_back(new std::pair<const T,int>(neighbour, g.get_edge_weight(u, neighbour)));
	}
	_pos = start_pos;
}

template <typename T> weighted_graph<T>::neighbour_iterator::~neighbour_iterator()
{
	// Cannot destroy properly due to not being able to use use smart or unique pointers. Left in as a proof of intent
	// for (auto it = neighbours.begin(); it != neighbours.end(); ++it)
	// {
	// 	delete (*it);
	// }
	// neighbours.clear();
}

template <typename T> 
typename weighted_graph<T>::neighbour_iterator weighted_graph<T>::neighbour_iterator::operator=(const neighbour_iterator& it)
{
	if (*this == it) return * this; // nothing to do
	this->neighbours = it.neighbours;
	this->_pos = it._pos;
	return *this;
}

template <typename T> bool weighted_graph<T>::neighbour_iterator::operator==(const neighbour_iterator& it) const 
{
	if (this->_pos != it._pos || this->neighbours.size() != it.neighbours.size()) return false;
	for (int i = 0; i < this->neighbours.size(); i++) // comparing vectors of pointers directly doesn't work, so do it manually
	{
		if (*(this->neighbours.at(i)) != *(it.neighbours.at(i))) return false;
	}
	return true;
}
template <typename T> bool weighted_graph<T>::neighbour_iterator::operator!=(const neighbour_iterator& it) const { return !(*this == it); }

template <typename T> 
typename weighted_graph<T>::neighbour_iterator weighted_graph<T>::neighbour_iterator::operator++() 
{
	_pos++; // increment before returning
	return * this;
}

template <typename T> 
typename weighted_graph<T>::neighbour_iterator weighted_graph<T>::neighbour_iterator::operator++(int)
{
	auto n = * this; // return before incrementing
	_pos++;
	return n;
}

template <typename T> const std::pair<T, int> weighted_graph<T>::neighbour_iterator::operator*() {return *neighbours.at(_pos);}
template <typename T> const std::pair<const T, int>* weighted_graph<T>::neighbour_iterator::operator->() {return neighbours.at(_pos);}

template <typename T>	
typename weighted_graph<T>::graph_iterator weighted_graph<T>::begin()
{
	return graph_iterator(*this);
}

template <typename T>	
typename weighted_graph<T>::graph_iterator weighted_graph<T>::end() {
	return graph_iterator(*this, _num_vertices);
}

template <typename T>	
typename weighted_graph<T>::neighbour_iterator weighted_graph<T>::neighbours_begin(const T& u)
{
	return neighbour_iterator(*this, u);
}

template <typename T>	
typename weighted_graph<T>::neighbour_iterator weighted_graph<T>::neighbours_end(const T& u) 
{
	return neighbour_iterator(*this, u, degree(u));
}

template <typename T> 
weighted_graph<T>::weighted_graph()
{
	_num_vertices = 0;
	_num_edges = 0;
	_total_weight = 0;
}

template <typename T> 
weighted_graph<T>::~weighted_graph()
{
	for (T vertex : get_vertices())
	{
		remove_vertex(vertex);
	}
}

template <typename T> 
bool weighted_graph<T>::has_vertex(const T& u) const 
{
	for (auto vertex : vertices)
	{
		if (vertex.first == u) return true;
	}
	return false;
}

template <typename T> 
bool weighted_graph<T>::are_adjacent(const T& u, const T& v) const 
{
	// if v is a neighbour of u
	for (T vertex : get_neighbours(u))
	{
		if (vertex == v) return true;
	}
	return false;
}

template <typename T> 
void weighted_graph<T>::add_vertex(const T& v) 
{
	//add vertex pair with empty vector of edges
	if (!has_vertex(v))
	{
		vertices.insert(std::pair<T, std::map<T,int>*> (v, new std::map<T,int>));
		_num_vertices++;
	}
}

template <typename T> 
void weighted_graph<T>::add_edge(const T& u, const T& v, const int& weight)
{
	// ensure both u and v exist
	if (!has_vertex(u)) add_vertex(u);
	if (!has_vertex(v)) add_vertex(v);
	
	vertices[u]->insert(std::pair<T,int> (v, weight));
	
	vertices[v]->insert(std::pair<T,int> (u, weight));
	
	_num_edges++;
	_total_weight += weight;
}

template <typename T> 
void weighted_graph<T>::remove_vertex(const T& u) // If this is done correctly this can be used by the destructor
{
	if (!has_vertex(u)) return; // Nothing to do here
	
	auto edges = vertices[u];
	for (T edge : get_neighbours(u))
	{
		remove_edge(u, edge); // disconnect all edges
	}
	
	delete edges;
	vertices.erase(u);
	_num_vertices--;
}

template <typename T>
int weighted_graph<T>::remove_edge_core(const T& u, const T& v)
{
	if (!are_adjacent(u, v)) return 0; // Nothing was removed so the weight returned should be 0 (For recording)
	
	auto edges = vertices[u];
	
	int weight = (*edges)[v];
	edges->erase(v);
	return weight;
}

template <typename T> 
void weighted_graph<T>::remove_edge(const T& u, const T& v) // If this is done correctly this can be used by the desctructor
{
	// if either u or v doesn't exist, there is nothing to remove
	if (!has_vertex(u) || !has_vertex(v)) return;

	// disconnect v from u
	int weight = remove_edge_core(u, v);
	
	// disconnect u from v
	weight = remove_edge_core(v, u);
	
	_num_edges--;
	_total_weight -= weight;
}

template <typename T>
int weighted_graph<T>::set_edge_weight_core(const T& u, const T& v, const int& newWeight)
{
	if (!are_adjacent(u, v)) return 0; // Nothing was changed so the weight returned should be 0 (For recording)
	
	int oldWeight = vertices[u]->at(v);
	vertices[u]->at(v) = newWeight;
	return oldWeight;
}

template <typename T> 
void weighted_graph<T>::set_edge_weight(const T& u, const T& v, const int& newWeight) 
{
	// if either u or v doesn't exist, there is nothing to set
	if (!has_vertex(u) || !has_vertex(v)) return;
	
	// set weight of connection from u to v
	int oldWeight = set_edge_weight_core(u, v, newWeight);
	
	// set weight of connection from v to u
	oldWeight = set_edge_weight_core(v, u, newWeight);
	
	_total_weight -= oldWeight;
	_total_weight += newWeight;
}

template <typename T> 
int weighted_graph<T>::get_edge_weight(const T& u, const T& v) const 
{
	if (!are_adjacent(u, v)) return 0; // Nothing to see here so the weight returned should be 0
	
	auto edges = *(vertices.at(u)); // operator[] is not const for map. use at() instead
	return edges[v];
}

template <typename T> 
int weighted_graph<T>::degree(const T& u) const 
{
	auto edges = vertices.at(u); // operator[] is not const for map. use at() instead.
	return edges->size();
}

template <typename T> 
int weighted_graph<T>::weighted_degree(const T& u) 
{
	auto edges = vertices[u];
	
	int _weighted_degree = 0;
	for (auto edge : *edges)
	{
		_weighted_degree += edge.second; // tally up weight of all connecting edges
	}
	
	return _weighted_degree;
}

template <typename T> int weighted_graph<T>::num_vertices() const { return _num_vertices; }

template <typename T> int weighted_graph<T>::num_edges() const { return _num_edges; }

template <typename T> int weighted_graph<T>::total_weight() { return _total_weight; }

template <typename T>	
std::vector<T> weighted_graph<T>::get_vertices() const
{
	std::vector<T> vertices_vector;
	
	for (auto vertex : vertices)
	{
		vertices_vector.push_back(vertex.first);
	}
	
	return vertices_vector;
}

template <typename T>	
std::vector<T> weighted_graph<T>::get_neighbours(const T& u) const 
{
	std::vector<T> neighbours;
	
	// if u doesn't exist, return empty vector
	if (!has_vertex(u)) return neighbours;
	
	auto edges = *(vertices.at(u)); // operator[] is not const for map. use at() instead
	for (auto edge : edges)
	{
		neighbours.push_back(edge.first);
	}

	return neighbours;
}

template <typename T>
void weighted_graph<T>::depth_first_core(const T& current, std::map<T,bool> *visited, std::vector<T> *visit_order)
{
	visited->at(current) = true;
	visit_order->push_back(current);
	
	for (T neighbour : get_neighbours(current))
	{
		if (!visited->at(neighbour))
		{
			depth_first_core(neighbour, visited, visit_order);
		}
	}
}

template <typename T> 
std::vector<T> weighted_graph<T>::depth_first(const T& start_vertex)
{
	std::vector<T> *visit_order = new std::vector<T>();
	std::map<T,bool> *visited = new std::map<T,bool>;
	// prepare visited map with blanks
	for (T vertex : get_vertices())
	{
		visited->insert(std::pair<T,bool> (vertex, false));
	}
	
	// start recursion
	depth_first_core(start_vertex, visited, visit_order);
	
	return *visit_order;
}

template <typename T>
std::vector<T> weighted_graph<T>::breadth_first(const T& start_vertex)
{
	std::vector<T> visit_order;
	std::queue<T> unprocessed;
	unprocessed.push(start_vertex);
	
	std::map<T,bool> visited;
	// prepare visited map with blanks
	for (T vertex : get_vertices())
	{
		visited.insert(std::pair<T,bool> (vertex, false));
	}
	
	while (!unprocessed.empty())
	{
		T current = unprocessed.front();
		unprocessed.pop();
		
		if (!visited.at(current))
		{
			visit_order.push_back(current);
			visited.at(current) = true;
			
			for (T neighbour : get_neighbours(current))
			{
				unprocessed.push(neighbour);
			}
		}
	}
	
	return visit_order;
}

template <typename T>	
weighted_graph<T> weighted_graph<T>::mst() 
{
	weighted_graph<T> * tree = new weighted_graph<T>();
	std::vector<std::pair<T, std::pair<T,int> * > *> * seen_edges = new std::vector<std::pair<T, std::pair<T,int> * > *>(); 	
	// This is a vector of edges spotted along the way. 
	// The next smallest edges will be chosen from this.
	// T refers to the "to edge", and the pair is the source and weight of the edge. 
	// This is the opposite of how the vertices map works
	
	// Prim's algorithm
	auto current = (*(vertices.begin())).first; // get arbitrary first element
	tree->add_vertex(current);
	
	while (true)
	{
		// record new edges seen
		for (T neighbour : get_neighbours(current))
		{
			seen_edges->push_back(new std::pair<T,std::pair<T,int>*>(neighbour, new std::pair<T,int>(current, get_edge_weight(current, neighbour))));
		}

		// get next smallest edge
		std::pair<T,std::pair<T,int>*> * smallest_edge = nullptr;
		for (auto edge : *seen_edges)
		{
			if (!tree->has_vertex(edge->first) // don't add edge twice
				&& (smallest_edge == nullptr // if it's the first edge
			  		|| smallest_edge->second->second > edge->second->second))
			{
				smallest_edge = edge; // become the new smallest
			}
		}
		
		if (smallest_edge == nullptr) break; // no edges left to add, we're done.
		
		tree->add_edge(smallest_edge->first, smallest_edge->second->first, smallest_edge->second->second);
		current = smallest_edge->first;
	}
	
	// Ensure proper destruction of seen_edges
	if (!seen_edges->empty())
	{
		for (auto edge : *seen_edges)
		{
			delete edge->second;
			delete edge;
		}
	}
	delete seen_edges;
		
	return * tree;
}


#endif
