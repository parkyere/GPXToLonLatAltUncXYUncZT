/**
 * Frechet-based map construction 2.0 Copyright 2013 Mahmuda Ahmed and Carola Wenk
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  ------------------------------------------------------------------------
 *
 *  This software is based on the following article. Please cite this article when using this code
 * as part of a research publication:
 *
 *  Mahmuda Ahmed and Carola Wenk, "Constructing Street Networks from GPS Trajectories", European
 * Symposium on Algorithms (ESA): 60-71, Ljubljana, Slovenia, 2012
 *
 *  ------------------------------------------------------------------------
 *
 * Author: Mahmuda Ahmed Filename: MapConstruction.java
 *
 */
#pragma once
#include <string>
#include <vector>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include "Vertex.hpp"
#include "Edge.hpp"
#include "rts_smoother.h"

struct PoseFile {
	std::string fileName;
	std::vector<Vertex> curve;
	PoseFile() {
		fileName = "";
	}
	PoseFile(std::string fileName, std::vector<Vertex> curve) 
		: fileName{ fileName }, curve{ std::move(curve) } {}
	std::string getFileName() const { return fileName; }
	std::vector<Vertex>& getPose() { return curve; }
	double getLength() {
		double length = 0.0;
		for (size_t i = 1; i < curve.size(); ++i) {
			length += curve[i - 1].dist(curve[i]);
		}
		return length;
	}
	static PoseFile readFile(const std::string& fileName, const std::vector<StepData>& SmoothedTrajectory, bool hasAltitude = true) {
		PoseFile poseFile;
		poseFile.fileName = fileName;
		for (auto& v : SmoothedTrajectory) {
			Vertex vertex(v.filtered_state.pos[0], v.filtered_state.pos[1],
				hasAltitude ? v.filtered_state.pos[2] : 0.0,
				v.timestamp);
			poseFile.curve.push_back(vertex);
		}
		return poseFile;
	}
};
/**
 * An object that takes a set of poses as input, construct graph and write two
 * files one for vertices and one for edges.
 */
struct MapConstruction {
	static int curveid;
	static std::string curveName;
/**
 * Writes the constructed map into files.
 */

	static void writeToFile(std::vector<Vertex> vList, std::string fileName) {
		try {
			int count = 0;
			std::ostringstream bwedges;
			std::ostringstream bvertex;
			if (vList.empty()) {
				std::cerr << "No vertices to write." << std::endl;
				return;
			}
			for (int i = 0; i < vList.size(); i++) {
				Vertex& v = vList[i];
				auto llh = ecefToGeodetic(v.getX(), v.getY(), v.getZ());
				bvertex << std::format("{0}, {1:.8f}, {2:.8f}, {3:.8f}\n", i, llh[0], llh[1], llh[2]);
				for (int j = 0; j < v.getDegree(); j++) {
					if (i != v.getAdjacentElementAt(j)) {
						bwedges << std::format("{0}, {1}, {2}\n", count, i, v.getAdjacentElementAt(j));
	
						count++;
					}
				}
			}
			std::ofstream ofsEdges(fileName + "edges.txt");
			if (ofsEdges.is_open()) {
				ofsEdges << bwedges.str();
				ofsEdges.close();
			}
			else {
				std::cerr << "Error opening edges file for writing: " << fileName + "edges.txt" << std::endl;
			}
			std::ofstream ofsVertices(fileName + "vertices.txt");
			if (ofsVertices.is_open()) {
				ofsVertices << bvertex.str();
				ofsVertices.close();
			}
			else {
				std::cerr << "Error opening vertices file for writing: " << fileName + "vertices.txt" << std::endl;
			}
		}
		catch (std::exception& ex) {
			std::cerr << ex.what();
		}
	}
/**
 * Computes interval on edge e for a line segment consists of
 * (currentIndex-1)-th and currentIndex-th vertices of pose and return true
 * if edge e has a part of white interval else false.
 */
	bool isWhiteInterval(Edge& edge, std::vector<Vertex>& pose, int currentIndex, double eps, double altEps) {
		Line line(pose[currentIndex - 1], pose[currentIndex]);
		if (std::abs(line.avgAltitude() - edge.getLine().avgAltitude()) <= altEps) {
			return line.pIntersection(edge, eps);
		}
		else {
			return false;
		}
	}
/**
 * Sets corresponding interval endpoints on Edge.
 */
	void setEndPointsOnEdge(Edge& edge, int startIndex, int endIndex,
		double cstart, double vstart) {
		edge.setCurveStartIndex(startIndex);
		edge.setCurveStart(startIndex + cstart);
		edge.setEdgeStart(vstart);

		edge.setCurveEnd(endIndex - 1 + edge.getCurveEnd());
		edge.setCurveEndIndex(endIndex);
	}
	/**
	 * Scans for next white interval on an Edge starting from index newstart of
	 * pose.
	 */
	void computeNextInterval(Edge& edge, std::vector<Vertex>& pose, int newstart,
		double eps, double altEps) {
		// Compute next white interval on edge.
		bool first = true;
		bool debug = false;

		int startIndex = 0;
		double cstart = 0, vstart = 0;
		if (newstart >= pose.size()) {
			edge.setCurveEndIndex(pose.size());
			edge.setDone(true);
			return;
		}
		for (int i = newstart; i < pose.size(); i++) {
			bool result = isWhiteInterval(edge, pose, i, eps, altEps);

			// first = true means we are still looking for our first interval
			// starting from newstart.
			// !result indicate Line(pose.get(i), pose.get(i+1)) doesn't contain
			// white interval.
			// we can just ignore if(first && !result).

			if (first && result) {
				// first segment on the white interval
				first = false;
				startIndex = i - 1;
				cstart = edge.getCurveStart();
				vstart = edge.getEdgeStart();

				// if the white interval ends within the same segment
				if (edge.getCurveEnd() < 1) {
					setEndPointsOnEdge(edge, startIndex, i, cstart, vstart);
					return;
				}
			}
			else if (!first && result) {
				// not the first segment on the white interval
				if (edge.getCurveEnd() < 1) {
					// if the white interval ends within that segment
					setEndPointsOnEdge(edge, startIndex, i, cstart, vstart);
					return;
				}
			}
			else if (!first && !result) {
				// the white interval ends at 1.0 of previous segment
				setEndPointsOnEdge(edge, startIndex, i, cstart, vstart);
				return;
			}
		}

		if (first) {
			// if the last segment on the curve is the first segment of that
			// interval
			edge.setCurveEndIndex(pose.size());
			edge.setDone(true);
		}
		else {
			edge.setCurveStartIndex(startIndex);
			edge.setCurveStart(startIndex + cstart);
			edge.setEdgeStart(vstart);
			edge.setCurveEnd(pose.size() - 2 + edge.getCurveEnd());
			edge.setCurveEndIndex(pose.size() - 2);
		}
		return;
	}
/**
 * Updates constructedMap by adding an Edge. Detail description of the
 * algorithm is in the publication.
 */
	void updateMap(std::vector<Vertex>& constructedMap,
		std::map<std::string, int>& map, Edge& edge) {

		// update the map by adding a new edge
		Vertex v;
		int parent = -1;
		int child = -1;

		std::string keyParent = edge.getVertex1().toString();
		std::string keyChild = edge.getVertex2().toString();
		// find the index of parent node
		if (map.contains(keyParent)) {
			parent = map[keyParent];
		}
		else {
			v = edge.getVertex1();
			constructedMap.push_back(v);
			parent = std::find(constructedMap.begin(),constructedMap.end(),v)-constructedMap.begin();
			map[keyParent] = parent;
		}
		// find the index of child node
		if (map.contains(keyChild)) {
			child = map[keyChild];
		}
		else {
			v = edge.getVertex2();
			constructedMap.push_back(v);
			child = std::find(constructedMap.begin(),constructedMap.end(),v)-constructedMap.begin();
			map[keyChild]= child;
		}
		// update the map
		if (parent == -1 || child == -1) {
			std::cerr << "inconsistent graph child, parent :" << child << ", " << parent << std::endl;
		}
		else if (parent != child) {
			constructedMap[parent].addElementAdjList(child);
			constructedMap[child].addElementAdjList(parent);
			//logger.log(Level.FINEST, "child, parent :" + child + ", " + parent);
			//logger.log(Level.FINEST, "child, parent :" + parent + ", " + child);
		}
	}
/**
 * Adds a split point on an Edge.
 *
 * @param newVertexPosition
 *            represents position of a new Vertex
 */
	void edgeSplit(std::vector<Vertex>& constructedMap,
		std::map<std::string, int>& map, Edge& edge, double newVertexPosition) {
		Vertex v1 = edge.getVertex1();
		Vertex v2 = edge.getVertex2();

		auto key1 = v1.toString();
		auto key2 = v2.toString();

		// call of this method always after updateMap which ensures
		// map.containsKey(key1) is
		// always true.
		int index1 = map[key1];
		int index2 = map[key2];
		Vertex v = edge.getLine().getVertex(newVertexPosition);

		// splitting an edge on split point vertex v
		auto key = v.toString();
		int index = map[key];

		if (index == index1 || index == index2) {
			return;
		}
		//logger.log(Level.FINER, "Index = " + index1 + " " + index2 + " "
		//	+ index);
		edge.addSplit(newVertexPosition, index);
	}
};