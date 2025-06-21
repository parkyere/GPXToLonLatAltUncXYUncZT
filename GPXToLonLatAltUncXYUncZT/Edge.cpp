#include "Edge.hpp"
#include "Line.hpp"
//#include "Vertex.hpp"

 void Edge::reset() {
	curveStart = std::numeric_limits<double>::max();
	curveEnd = -1.0;
	edgeStart = std::numeric_limits<double>::max();
	edgeEnd = -1.0;
	line = LinePtr{ new Line(vertex1, vertex2) };
	done = false;
}

 Edge::Edge(VertexPtr v1, VertexPtr v2) : vertex1{ v1 }, vertex2{ v2 },
edgeSplitPositions{ new std::vector<double>() },
edgeSplitVertices{ new std::vector<int>() } {
	reset();
}
 void Edge::set(Edge& edge) {
	curveStart = edge.curveStart;
	curveEnd = edge.curveEnd;
	edgeStart = edge.edgeStart;
	edgeEnd = edge.edgeEnd;
	line = LinePtr{ new Line(vertex1, vertex2) };
	done = edge.done;
}

 VertexPtr Edge::getVertex1() const { return vertex1; }

 VertexPtr Edge::getVertex2() const { return vertex2; }

 LinePtr Edge::getLine() { return line; }

 bool Edge::getDone() { return done; }

 void Edge::setDone(bool done) { this->done = done; }

 double Edge::getCurveStart() { return curveStart; }

 void Edge::setCurveStart(double curveStart) { this->curveStart = curveStart; }

 double Edge::getCurveEnd() { return curveEnd; }

 void Edge::setCurveEnd(double curveEnd) { this->curveEnd = curveEnd; }

 double Edge::getEdgeStart() { return edgeStart; }

 void Edge::setEdgeStart(double edgeStart) { this->edgeStart = edgeStart; }

 double Edge::getEdgeEnd() { return edgeEnd; }

 void Edge::setEdgeEnd(double edgeEnd) { this->edgeEnd = edgeEnd; }

 int Edge::getCurveStartIndex() { return curveStartIndex; }

 void Edge::setCurveStartIndex(int startIndex) {
	if (startIndex >= 0) {
		this->curveStartIndex = startIndex;
	}
	else {
		std::cerr << "Invalid assignment of Edge.startIndex\n";
	}
}

 int Edge::getCurveEndIndex() { return curveEndIndex; }

 void Edge::setCurveEndIndex(int endIndex) {
	curveEndIndex = endIndex;
}


// Inserts a new split position if the list doesn't have it, otherwise return.
// position: indicvate the new split position
// vertex : indicates the vertex which should be inserted in this edge.
 void Edge::addSplit(double position, int vertex) {
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

 std::string Edge::toString() const {
	return vertex1->toString() + vertex2->toString();
}
