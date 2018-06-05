#ifndef GRAPH_ALGS
#define GRAPH_ALGS

#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <algorithm>
#include "weighted_graph.hpp"
#include "easy_weighted_graph_algorithms.cpp"

template <typename T>
std::vector<T> get_vertices(const weighted_graph<T>& g)
{
	std::vector<T> vertices;
	
	for (auto it = g.begin(); it != g.end(); it++)
	{
		vertices.push_back(*it);
	}
	return vertices;
}

template <typename T>
std::map<T,int> get_neighbours_weighted(const weighted_graph<T>& g, T v)
{
	std::map<T,int> neighbours;
	
	// if v doesn't exist, return an empty map
	if (!g.has_vertex(v)) return neighbours;
	
	for (T u : get_vertices(g))
	{
		if (g.are_adjacent(u, v))
		{
			neighbours.insert(std::pair<T,int>(u, g.get_edge_weight(u, v)));
		}
	}
	
	return neighbours;
}

template <typename T>
std::vector<T> get_neighbours(const weighted_graph<T>& g, T vertex)
{
	std::vector<T> neighbours;

	for (auto edge : get_neighbours_weighted(g, vertex))
	{
		neighbours.push_back(edge.first);
	}
	
	return neighbours;
}

//Returns true if the graph is connected, false otherwise.
template <typename T>
bool is_connected(const weighted_graph<T>& g)
{
	if (g.num_vertices() == 0)
	{
		return true; // save us from exceptions
	}
	
	auto vertices = get_vertices(g);
	// Pick an arbitrary start point
	T start = *(vertices.begin());
	auto connected_vertices = depth_first(g, start);
	
	for (T v : vertices)
	{
		bool found = false;
		for (T u : connected_vertices)
		{
			if (v == u)
			{
				found = true;
				break;
			}
		}
		
		if (!found)
		{
			return false;
		}
	}
	
	return true;
}

//Returns a vector of weighted graphs, where each weighted graph is a connected
//component of the input graph.
template <typename T>
std::vector<weighted_graph<T>> connected_components(const weighted_graph<T>& g)
{
	auto components = std::vector<weighted_graph<T>>();
	std::map<T,bool> visited;
	// prepare visited map with blanks
	for (T vertex : get_vertices(g))
	{
		visited.insert(std::pair<T,bool>(vertex, false));
	}
	
	// Pick an arbitrary start point
	auto start = visited.begin();
	while (start != visited.end())
	{
		weighted_graph<T> component = weighted_graph<T>();
		auto connected_vertices = depth_first(g, (*start).first);
		
		// mark all seen and build the current component
		for (T v : connected_vertices)
		{
			visited.at(v) = true;
			component.add_vertex(v);
		}
		
		for (T v : connected_vertices)
		{
			auto neighbours = get_neighbours_weighted(g, v);
			for (auto edge : neighbours)
			{
				if (!component.are_adjacent(v, edge.first))
				{
					component.add_edge(v, edge.first, edge.second);
				}
			}
		}
		
		components.push_back(component);
		
		// pick a new starting point
		bool starting_point_found = false;
		for (auto it = visited.begin(); it != visited.end(); it++)
		{
			if (!(*it).second)
			{
				start = it;
				starting_point_found = true;
			}
		}
		
		if (!starting_point_found) break;
	}
	
	return components;
}

//Returns a map of the vertices of the weighted graph g and their distances from
//the given starting vertex v.
template <typename T> 
std::map<T, int> dijkstras(const weighted_graph<T>& g, const T& start)
{
	std::map<T,bool> visited;
	// perpare visited map with blanks
	for (T vertex : get_vertices(g))
	{
		visited.insert(std::pair<T,bool>(vertex, false));
	}
	
	std::map<T,int> distances;
	// perpare distances map with blanks
	for (T vertex : get_vertices(g))
	{
		distances.insert(std::pair<T,int>(vertex, -1));
	}

	if (distances.size() == 0) return distances; // save us from exceptions
	distances.at(start) = 0;
	
	T current = start;
	
	while (true)
	{
		// For each unvisited neighbour
		for (auto neighbour : get_neighbours_weighted(g, current))
		{
			if (!visited.at(neighbour.first) // unvisisted
				&& (distances.at(neighbour.first) == -1 // doesn't have a tentative value
				|| distances.at(neighbour.first) > distances.at(current) + neighbour.second)) // or tentaive value > current path
			{
				distances.at(neighbour.first) = distances.at(current) + neighbour.second;
			}
		}
		
		// Mark current vertex as visited
		visited.at(current) = true;
		
		// While there are any unvisited vertices get the vertex with the smallest tentative distances
		std::pair<T,int> next_smallest_distance = std::pair<T,int>(T(),-1);
		for (auto vertex : visited)
		{
			if (!vertex.second // must be unvisited
				&& distances.at(vertex.first) != -1 // and have a tentative distance
				&& (next_smallest_distance.second == -1 // and either be the first time round
				|| next_smallest_distance.second > distances.at(vertex.first))) // or be smaller than the current smallest
			{
				next_smallest_distance = std::pair<T,int>(vertex.first, distances.at(vertex.first)); // become the new smallest
			}
		}

		if (next_smallest_distance.second == -1) break; // no unvisited vertices left
		current = next_smallest_distance.first;
	}
	
	return distances;
}

template <typename T>
void articulation_points_core(const weighted_graph<T>& g, T vertex, int depth, std::map<T,bool> * visited, std::map<T,int> * low, std::map<T,int> * depths, std::map<T,T> * parent, std::vector<T> * ap)
{
	visited->at(vertex) = true;
	(*depths)[vertex] = depth;
	(*low)[vertex] = depth;
	int child_count = 0;
	bool is_articulation = false;
	
	for (T neighbour : get_neighbours(g, vertex))
	{
		if (!visited->at(neighbour))
		{
			(*parent)[neighbour] = vertex;
			articulation_points_core(g, neighbour, depth+1, visited, low, depths, parent, ap);
			child_count++;
			if (low->at(neighbour) >= depths->at(vertex))
			{
				is_articulation = true;
			}
			(*low)[vertex] = std::min(low->at(vertex), low->at(neighbour));
		}
		else if (parent->find(vertex) != parent->end() && neighbour != parent->at(vertex))
		{
			(*low)[vertex] = std::min(low->at(vertex), depths->at(neighbour));
		}
	}
	
	if ((parent->find(vertex) != parent->end() && is_articulation)
	   || (parent->find(vertex) == parent->end() && child_count > 1))
	{
		ap->push_back(vertex);
	}
}

//Returns a vector containing all the articulation points of the
//input weighted graph g.
template <typename T>
std::vector<T> articulation_points(const weighted_graph<T>& g)
{
	auto vertices = get_vertices(g);
	std::map<T,bool> * visited = new std::map<T,bool>();
	std::map<T,int> * low = new std::map<T,int>();
	std::map<T,int> * depths = new std::map<T,int>();
	std::map<T,T> * parent = new std::map<T,T>();
	std::vector<T> * ap = new std::vector<T>();
	
	for (T vertex : vertices)
	{
		visited->insert(std::pair<T,bool>(vertex, false));
	}
	
	for (T vertex : vertices)
	{
		if (!visited->at(vertex))
		{
			articulation_points_core(g, vertex, 0, visited, low, depths, parent, ap);
		}
	}
	
	// Clean up after yourself
	delete visited;
	delete low;
	delete depths;
	delete parent;
	return * ap;
}

#endif
