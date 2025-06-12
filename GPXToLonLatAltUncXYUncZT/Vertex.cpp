#pragma once
#include "Vertex.hpp"
bool isEqual(double a, double b) { return std::abs(a - b) < EqualityCriterion; }
Vertex::Vertex() : adjacencyList{ new std::vector<int>() } {}
Vertex::Vertex(double x, double y, double z) : Vertex() {
	this->x = x;
	this->y = y;
	this->z = z;
}
Vertex::Vertex(double x, double y, double z, double timestamp) : Vertex(x, y, z) {
	this->timestamp = timestamp;
}
Vertex::Vertex(double lat, double lng, double alt, double x, double y, double z) : Vertex(x, y, z) {
	this->lat = lat;
	this->lng = lng;
	this->alt = alt;
}
Vertex::Vertex(double lat, double lng, double alt, double x, double y, double z, double timestamp) : Vertex(lat, lng, alt, x, y, z) {
	this->timestamp = timestamp;
}
double Vertex::getX() const { return x; }
double Vertex::getY() const { return y; }
double Vertex::getZ() const { return z; }

double Vertex::getLat() const { return lat; }

double Vertex::getLng() const { return lng; }

double Vertex::getAlt() const { return alt; }

double Vertex::norm() const { return std::sqrt(x * x + y * y + z * z); }

double Vertex::dotProd(const Vertex& vector1, const Vertex& vector2) {
	return vector1.getX() * vector2.getX()
		+ vector1.getY() * vector2.getY()
		+ vector1.getZ() * vector2.getZ();
}

int Vertex::getDegree() const { return adjacencyList->size(); }

bool Vertex::getDone() const { return done; }

double Vertex::getTimestamp() const { return timestamp; }

void Vertex::setDone(bool done) { this->done = done; }

std::vector<int>& Vertex::getAdjacencyList() { return *adjacencyList; }

void Vertex::addElementAdjList(int v) {
	for (int i = 0; i < this->getDegree(); i++) {
		if ((*(this->adjacencyList))[i] == v) {
			return;
		}
	}
	adjacencyList->push_back(v);
}

int Vertex::getIndexAdjacent(int v) {
	auto it = std::find(adjacencyList->begin(), adjacencyList->end(), v);
	if (it != adjacencyList->end())
		return std::distance(adjacencyList->begin(), it);
	else
		return -1; // Not found
}

int Vertex::getAdjacentElementAt(int index) {
	if (index >= 0 && index < adjacencyList->size()) {
		return (*adjacencyList)[index];
	}
	return -1; // Invalid index
}

void Vertex::setAdjacentElementAt(int index, int value) {
	if (index >= 0 && index < adjacencyList->size()) {
		(*adjacencyList)[index] = value;
	}
}

double Vertex::dist(Vertex& v2) const {
	double dx = x - v2.getX();
	double dy = y - v2.getY();
	return std::sqrt(dx * dx + dy * dy);
}

void Vertex::reset() { done = false; }

//resets a vertex's processing state
std::string Vertex::toString() const {
	return std::format("{0:.16f} {1:.16f} {2:.16f}", x, y, z);
}

Vertex Vertex::deepCopy() {
	Vertex vertex(lat, lng, alt, x, y, z, timestamp);
	vertex.done = done;
	for (int i = 0; i < adjacencyList->size(); i++) {
		vertex.adjacencyList->push_back((*adjacencyList)[i]);
	}
	return vertex;
}
