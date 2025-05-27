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
        struct Data {
                double x{0}, y{0}, z{0}; // ECEF coordinates
                double lng{0}, lat{0}, alt{0}; // Geodetic coordinates
                bool done{false};
                std::vector<int> adjacencyList; // Contains the indices of adjacent vertices
                double timestamp{-1}; //Timestamp in milliseconds
        };
        std::shared_ptr<Data> data;
public:
        Vertex() : data(std::make_shared<Data>()) {}
        Vertex(double x, double y, double z) : data(std::make_shared<Data>()) {
                data->x = x; data->y = y; data->z = z;
        }
        Vertex(double x, double y, double z, double timestamp) : Vertex(x, y, z) {
                data->timestamp = timestamp;
        }
        Vertex(double lat, double lng, double alt, double x, double y, double z) : Vertex(x, y, z) {
                data->lat = lat;
                data->lng = lng;
                data->alt = alt;
        }
        Vertex(double lat, double lng, double alt, double x, double y, double z, double timestamp) : Vertex(lat, lng, alt, x, y, z) {
                data->timestamp = timestamp;
        }
        double getX() const { return data->x; }
        double getY() const { return data->y; }
        double getZ() const { return data->z; }
        double getLng() const { return data->lng; }
        double getLat() const { return data->lat; }
        double getAlt() const { return data->alt; }
        double norm()   const { return std::sqrt( data->x*data->x + data->y*data->y + data->z*data->z ); }
	static double dotProd(const Vertex& vector1, const Vertex& vector2) {
		return vector1.getX() * vector2.getX()
			+ vector1.getY() * vector2.getY()
			+ vector1.getZ() * vector2.getZ();
	}
        int getDegree() { return data->adjacencyList.size(); }
        bool getDone() { return data->done; }
        double getTimestamp() { return data->timestamp; }
        void setDone(bool done) { data->done = done; }
        std::vector<int>& getAdjacencyList() { return data->adjacencyList; }
        void addElementAdjList(int v) {
                for (int i = 0; i < this->getDegree(); i++) {
                        if (data->adjacencyList[i] == v) {
                                return;
                        }
                }
                data->adjacencyList.push_back(v);
        }
        int getIndexAdjacent(int v) {
                auto it = std::find(data->adjacencyList.begin(), data->adjacencyList.end(), v);
                if (it != data->adjacencyList.end())
                        return std::distance(data->adjacencyList.begin(), it);
                else
                        return -1; // Not found
        }
        int getAdjacentElementAt(int index) {
                if (index >= 0 && index < data->adjacencyList.size()) {
                        return data->adjacencyList[index];
                }
                return -1; // Invalid index
        }
        void setAdjacentElementAt(int index, int value) {
                if (index >= 0 && index < data->adjacencyList.size()) {
                        data->adjacencyList[index] = value;
                }
        }
        double dist(Vertex& v2) {
                double dx = data->x - v2.getX();
                double dy = data->y - v2.getY();
                return std::sqrt(dx * dx + dy * dy);
        }
        void reset() { data->done = false; }
        std::string toString() const {
                return std::format("{0:.8f} {1:.8f} {2:.8f}",data->x,data->y,data->z);
        }
        Vertex deepCopy() {
                Vertex vertex(data->lat, data->lng, data->alt, data->x, data->y, data->z, data->timestamp);
                vertex.data->done = data->done;
                for (int i = 0; i < data->adjacencyList.size(); i++) {
                        vertex.data->adjacencyList.push_back(data->adjacencyList[i]);
                }
                return vertex;
        }
};
using VertexPtr = std::shared_ptr<Vertex>;
