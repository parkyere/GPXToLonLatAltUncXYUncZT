#include "Edge.hpp"
#include "Line.hpp"
#include "Vertex.hpp"

inline void Edge::reset() {
	curveStart = std::numeric_limits<double>::max();
	curveEnd = -1.0;
	edgeStart = std::numeric_limits<double>::max();
	edgeEnd = -1.0;
	line = LinePtr{ new Line(vertex1, vertex2) };
	done = false;
}

inline Edge::Edge(VertexPtr v1, VertexPtr v2) : vertex1{ v1 }, vertex2{ v2 },
edgeSplitPositions{ new std::vector<double>() },
edgeSplitVertices{ new std::vector<int>() } {
	reset();
}
inline void Edge::set(Edge& edge) {
	curveStart = edge.curveStart;
	curveEnd = edge.curveEnd;
	edgeStart = edge.edgeStart;
	edgeEnd = edge.edgeEnd;
	line = LinePtr{ new Line(vertex1, vertex2) };
	done = edge.done;
}

inline VertexPtr Edge::getVertex1() const { return vertex1; }

inline VertexPtr Edge::getVertex2() const { return vertex2; }

inline LinePtr Edge::getLine() { return line; }

inline bool Edge::getDone() { return done; }

inline void Edge::setDone(bool done) { this->done = done; }

inline double Edge::getCurveStart() { return curveStart; }

inline void Edge::setCurveStart(double curveStart) { this->curveStart = curveStart; }

inline double Edge::getCurveEnd() { return curveEnd; }

inline void Edge::setCurveEnd(double curveEnd) { this->curveEnd = curveEnd; }

inline double Edge::getEdgeStart() { return edgeStart; }

inline void Edge::setEdgeStart(double edgeStart) { this->edgeStart = edgeStart; }

inline double Edge::getEdgeEnd() { return edgeEnd; }

inline void Edge::setEdgeEnd(double edgeEnd) { this->edgeEnd = edgeEnd; }

inline int Edge::getCurveStartIndex() { return curveStartIndex; }

inline void Edge::setCurveStartIndex(int startIndex) {
	if (startIndex >= 0) {
		this->curveStartIndex = startIndex;
	}
	else {
		std::cerr << "Invalid assignment of Edge.startIndex\n";
	}
}

inline int Edge::getCurveEndIndex() { return curveEndIndex; }

inline void Edge::setCurveEndIndex(int endIndex) {
	curveEndIndex = endIndex;
}


// Inserts a new split position if the list doesn't have it, otherwise return.
// position: indicvate the new split position
// vertex : indicates the vertex which should be inserted in this edge.
inline void Edge::addSplit(double position, int vertex) {
	int i = 0;
	for (i = 0; i < edgeSplitPositions->size(); i++) {
		if (isEqual((*edgeSplitPositions)[i], position)) {
			return;
		}
		else if ((*edgeSplitPositions)[i] > position) {
			auto it1 = edgeSplitPositions->begin() + i;
			auto it2 = edgeSplitVertices->begin() + i;
			edgeSplitPositions->insert(it1, position);
			edgeSplitVertices->insert(it2, vertex);
			return;
		}
	}
	edgeSplitPositions->push_back(position);
	edgeSplitVertices->push_back(vertex);
}

inline std::string Edge::toString() const {
	return vertex1->toString() + vertex2->toString();
}
