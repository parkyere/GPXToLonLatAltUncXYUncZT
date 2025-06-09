//#include "Vertex.hpp"
//#include "Line.hpp"
#pragma once
#ifndef COMMON_HPP
#define COMMON_HPP
#include <list>
#include <vector>
#include <limits>
#include <compare>
#include <memory>
#include <iostream>
#include <cmath>
#include <optional>
#include <utility>
#include <stdexcept>
#include <array>
#include <numbers>
#include <algorithm>
#include <string>
#include <format>
#include <filesystem>
#include <sstream>
#include <map>
#include <queue>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <cctype>
//#include "rts_smoother.h"
double deg2rad(double deg);
double rad2deg(double rad);
constexpr double EqualityCriterion = 1e-36; // Criterion for equality of two doubles
bool isEqual(double a, double b);
#endif
//struct Vertex;
//struct Edge;
//struct Line;
//using VertexPtr = std::shared_ptr<Vertex>;
//using EdgePtr = std::shared_ptr<Edge>;
//using LinePtr = std::shared_ptr<Line>;
