#include "Common.hpp"
//#include "Line.hpp"
#include "Vertex.hpp"
struct Line;
using LinePtr = std::shared_ptr<Line>;
struct Edge;
using EdgePtr = std::shared_ptr<Edge>;
struct Vertex;
using VertexPtr = std::shared_ptr<Vertex>;

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
	void reset();
	Edge(VertexPtr v1, VertexPtr v2);
	void set(Edge& edge);
	VertexPtr getVertex1() const;
	VertexPtr getVertex2() const;
	LinePtr getLine();
	bool getDone();
	void setDone(bool done);
	double getCurveStart();
	void setCurveStart(double curveStart);
	double getCurveEnd();
	void setCurveEnd(double curveEnd);
	double getEdgeStart();
	void setEdgeStart(double edgeStart);
	double getEdgeEnd();
	void setEdgeEnd(double edgeEnd);
	int getCurveStartIndex();
	void setCurveStartIndex(int startIndex);
	int getCurveEndIndex();
	void setCurveEndIndex(int endIndex);
	std::vector<int>& getEdgeSplitVertices() { return *edgeSplitVertices; }
	std::vector<double>& getEdgeSplitPositions() { return *edgeSplitPositions; }
	// Inserts a new split position if the list doesn't have it, otherwise return.
	// position: indicvate the new split position
	// vertex : indicates the vertex which should be inserted in this edge.
	void addSplit(double position, int vertex);
	std::string toString() const;
	std::strong_ordering operator<=>(const Edge& other) const {
		if (curveStart < other.curveStart) return std::strong_ordering::less;
		if (curveStart > other.curveStart) return std::strong_ordering::greater;
		if (curveEnd < other.curveEnd) return std::strong_ordering::less;
		if (curveEnd > other.curveEnd) return std::strong_ordering::greater;
		return std::strong_ordering::equal;
	}
};

using EdgePtr = std::shared_ptr<Edge>;