#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <format>
#include <memory>
#include "rts_smoother.h"

struct Vertex {
private:
	double x, y, z; // ECEF coordinates
	double lng, lat, alt; // Geodetic coordinates
	bool done = false;
	std::shared_ptr<std::vector<int>> adjacencyList; // Contains the indices of adjacent vertices
	double timestamp = -1; //Timestamp in milliseconds
public:
	Vertex() : adjacencyList{ new std::vector<int>() } {}
	Vertex(double x, double y, double z) :	Vertex() {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Vertex(double x, double y, double z, double timestamp) : Vertex(x, y, z) {
		this->timestamp = timestamp;
	}
	Vertex(double lat, double lng, double alt, double x, double y, double z) : Vertex(x, y, z) {
		this->lat = lat;
		this->lng = lng;
		this->alt - alt;
	}
	Vertex(double lat, double lng, double alt, double x, double y, double z, double timestamp) : Vertex(lat, lng, alt, x, y, z) {
		this->timestamp = timestamp;
	}
	double getX() const { return x; }
	double getY() const { return y; }
	double getZ() const { return z; }
	double getLng() const { return lng; }
	double getLat() const { return lat; }
	double getAlt() const { return alt; }
	double norm()   const { return std::sqrt( x*x + y*y + z*z ); }
	static double dotProd(const Vertex& vector1, const Vertex& vector2) {
		return vector1.getX() * vector2.getX()
			+ vector1.getY() * vector2.getY()
			+ vector1.getZ() * vector2.getZ();
	}
	int getDegree() const { return adjacencyList->size(); }
	bool getDone() const { return done; }
	double getTimestamp() const { return timestamp; }
	void setDone(bool done) { this->done = done; }
	std::vector<int>& getAdjacencyList() { return *adjacencyList; }
	void addElementAdjList(int v) {
		for (int i = 0; i < this->getDegree(); i++) {
			if ((*(this->adjacencyList))[i] == v) {
				return;
			}
		}
		adjacencyList->push_back(v);
	}
	int getIndexAdjacent(int v) {
		auto it = std::find(adjacencyList->begin(), adjacencyList->end(), v);
		if (it != adjacencyList->end())
			return std::distance(adjacencyList->begin(), it);
		else
			return -1; // Not found
	}
	int getAdjacentElementAt(int index) {
		if (index >= 0 && index < adjacencyList->size()) {
			return (*adjacencyList)[index];
		}
		return -1; // Invalid index
	}
	void setAdjacentElementAt(int index, int value) {
		if (index >= 0 && index < adjacencyList->size()) {
			(*adjacencyList)[index] = value;
		}
	}
	double dist(Vertex& v2) const {
		double dx = x - v2.getX();
		double dy = y - v2.getY();
		return std::sqrt(dx * dx + dy * dy);
	}
	void reset() { done = false; }
	std::string toString() const {
		return std::format("{0:.8f} {1:.8f} {2:.8f}",x,y,z);
	}
	Vertex deepCopy() {
		Vertex vertex(lat, lng, alt, x, y, z, timestamp);
		vertex.done = done;
		for (int i = 0; i < adjacencyList->size(); i++) {
			vertex.adjacencyList->push_back((*adjacencyList)[i]);
		}
		return vertex;
	}
};
using VertexPtr = std::shared_ptr<Vertex>;