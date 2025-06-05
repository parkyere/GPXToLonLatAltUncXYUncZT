#include "Common.hpp"

struct Vertex {
private:
	double x, y, z; // ECEF coordinates
	double lng, lat, alt; // Geodetic coordinates
	bool done = false;
	std::shared_ptr<std::vector<int>> adjacencyList; // Contains the indices of adjacent vertices
	double timestamp = -1; //Timestamp in milliseconds
public:
	Vertex();
	Vertex(double x, double y, double z);
	Vertex(double x, double y, double z, double timestamp);
	Vertex(double lat, double lng, double alt, double x, double y, double z);
	Vertex(double lat, double lng, double alt, double x, double y, double z, double timestamp);
	double getX() const;
	double getY() const;
	double getZ() const;
	double getLat() const;
	double getLng() const;
	double getAlt() const;
	double norm()   const;
	static double dotProd(const Vertex& vector1, const Vertex& vector2);
	int getDegree() const;
	bool getDone() const;
	double getTimestamp() const;
	void setDone(bool done);
	std::vector<int>& getAdjacencyList();
	void addElementAdjList(int v);
	int getIndexAdjacent(int v);
	int getAdjacentElementAt(int index);
	void setAdjacentElementAt(int index, int value);
	double dist(Vertex& v2) const;
	void reset(); //resets a vertex's processing state
	std::string toString() const;
	Vertex deepCopy();
};
using VertexPtr = std::shared_ptr<Vertex>;