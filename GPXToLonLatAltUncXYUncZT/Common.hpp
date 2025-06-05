//#include "Vertex.hpp"
//#include "Line.hpp"
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
#include "rts_smoother.h"
constexpr double EqualityCriterion = 1e-36; // Criterion for equality of two doubles
inline bool isEqual(double a, double b) { return std::abs(a - b) < EqualityCriterion; }
#endif
//struct Vertex;
//struct Edge;
//struct Line;
//using VertexPtr = std::shared_ptr<Vertex>;
//using EdgePtr = std::shared_ptr<Edge>;
//using LinePtr = std::shared_ptr<Line>;
