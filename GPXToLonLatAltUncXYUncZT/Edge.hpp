#pragma once
#include "Vertex.hpp"
#include "Line.hpp"
#include <list>
struct Edge {
private:
	Vertex vertex1, vertex2;	//first, second endpoint of edge
	double curveStart, curveEnd; // contains start, end point of white interval on curve
	double edgeStart, edgeEnd; // contains corresponding points on that edge
	Line line; // line from vertex v1 to v2
	bool done; // done is marked as true when an edge is done being compared with all segments of a curve.
	int curveStartIndex, curveEndIndex; // contains indices of start and end index curves for white intervals corresponding to the edge 
	std::list<double> edgeSplitPositions;   //Contains a sorted list (in descending order) of positions indicating where to split this edge.
	std::list<int> edgeSplitVertices;	//For each split position the corresponding new vertex is saved in this list.
public:
};