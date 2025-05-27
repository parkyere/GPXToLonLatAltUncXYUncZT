#pragma once
#include "Vertex.hpp"
#include "Line.hpp"
#include <list>
#include <vector>
#include <limits>
#include <compare>
#include <memory>
struct Edge {
private:
	VertexPtr vertex1, vertex2;	//first, second endpoint of edge
	double curveStart, curveEnd; // contains start, end point of white interval on curve
	double edgeStart, edgeEnd; // contains corresponding points on that edge
	LinePtr line; // line from vertex v1 to v2
	bool done; // done is marked as true when an edge is done being compared with all segments of a curve.
	int curveStartIndex, curveEndIndex; // contains indices of start and end index curves for white intervals corresponding to the edge 
	std::shared_ptr<std::vector<double>> edgeSplitPositions;   //Contains a sorted list (in descending order) of positions indicating where to split this edge.
	std::shared_ptr<std::vector<int>> edgeSplitVertices;	//For each split position the corresponding new vertex is saved in this list.
public:
	void reset() {
		curveStart = std::numeric_limits<double>::max();
		curveEnd = -1.0;
		edgeStart = std::numeric_limits<double>::max();
		edgeEnd = -1.0;
		line = new Line(vertex1, vertex2);
		done = false;
	}
	Edge(VertexPtr v1, VertexPtr v2) : vertex1{ v1 }, vertex2{ v2 },
		edgeSplitPositions{ new std::vector<double>()},
		edgeSplitVertices{ new std::vector<int>() } {
		reset();
	}
	void set(Edge edge) {
		curveStart = edge.curveStart;
		curveEnd = edge.curveEnd;
		edgeStart = edge.edgeStart;
		edgeEnd = edge.edgeEnd;
		line = Line(vertex1, vertex2);
		done = edge.done;
	}
	VertexPtr getVertex1() { return vertex1; }
	VertexPtr getVertex2() { return vertex2; }
	LinePtr getLine() { return line; }
	bool getDone() { return done; }
	void setDone(bool done) { this->done = done; }
	double getCurveStart() { return curveStart; }
	void setCurveStart(double curveStart) { this->curveStart = curveStart; }
	double getCurveEnd() { return curveEnd; }
	void setCurveEnd(double curveEnd) { this->curveEnd = curveEnd; }
	double getEdgeStart() { return edgeStart; }
	void setEdgeStart(double edgeStart) { this->edgeStart = edgeStart; }
	double getEdgeEnd() { return edgeEnd; }
	void setEdgeEnd(double edgeEnd) { this->edgeEnd = edgeEnd; }
	int getCurveStartIndex() { return curveStartIndex; }
	void setCurveStartIndex(int startIndex) {
		if (startIndex >= 0) {
			this->curveStartIndex = startIndex;
		}
	}
	int getCurveEndIndex() { return curveEndIndex; }
	void setCurveEndIndex(int endIndex) {
		curveEndIndex = endIndex;
	}
	std::vector<int>& getEdgeSplitVertices() { return *edgeSplitVertices; }
	std::vector<double>& getEdgeSplitPositions() { return *edgeSplitPositions; }
	// Inserts a new split position if the list doesn't have it, otherwise return.
	// position: indicvate the new split position
	// vertex : indicates the vertex which should be inserted in this edge.
	void addSplit(double position, int vertex) {
		int i = 0;
		for (i = 0; i < edgeSplitPositions.size(); i++) {
			if (edgeSplitPositions[i] == position) {
				return;
			}
			else if (edgeSplitPositions[i] > position) {
				auto it1 = edgeSplitPositions.begin() + i;
				auto it2 = edgeSplitVertices.begin() + i;
				edgeSplitPositions.insert(it1, position);
				edgeSplitVertices.insert(it2, vertex);
				return;
			}
		}
		edgeSplitPositions.push_back(position);
		edgeSplitVertices.push_back(vertex);
	}
	std::string toString() const {
		return vertex1.toString() + vertex2.toString();
	}
	std::strong_ordering operator<=>(const Edge& other) const {
		if (curveStart < other.curveStart) return std::strong_ordering::less;
		if (curveStart > other.curveStart) return std::strong_ordering::greater;
		if (curveEnd < other.curveEnd) return std::strong_ordering::less;
		if (curveEnd > other.curveEnd) return std::strong_ordering::greater;
		return std::strong_ordering::equal;
	}
};
using EdgePtr = std::shared_ptr<Edge>;