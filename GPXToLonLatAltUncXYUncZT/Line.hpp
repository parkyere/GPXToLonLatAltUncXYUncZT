#include "Common.hpp"
#include "Vertex.hpp"

struct Line;
using LinePtr = std::shared_ptr<Line>;
struct Edge;
using EdgePtr = std::shared_ptr<Edge>;
struct Vertex;
using VertexPtr = std::shared_ptr<Vertex>;

struct Line {
private:
    VertexPtr p1, p2;
    double xdiff, ydiff, zdiff;
    double c{ 0 }; //y-intersect of the line equation
    double m{ 0 }; // slope of the line
    double theta; //smallest angle between x-axis and this Line object
public:
	Line() = default;
    Line(VertexPtr p1, VertexPtr p2);
    VertexPtr getP1();
    VertexPtr getP2();
    double getXdiff() const;
    double getYdiff() const;
    double getZdiff() const;
    double getM() const;
    double getC() const;
    double getTheta() const;
    void setM(double m);
    void setC(double c);
    void setTheta(double theta);
    std::string toString() const;
    double avgAltitude() const;
    /**
 * Compute intersection of eps-disc around p and line
 * @param p the center of the disc
 * @param eps the radius of the disc
 * @return a optional<array<double,2>> containing two intersection points or null when they don't intersect.
 */
    std::optional<std::array<double,2>> pIntersection(const Vertex& p, double eps, bool precisionCompromise) const;
    /**
 * Returns intersection points of this line and two lines in Line[] (neighborhood). For detail see
 * {@link package-info}
 *
 * @return an array of four double values, first two values represents interval start and end
 *         points on this line and next two corresponds to intersections with lines[0] and
 *         lines[1] respectively.
 */
    std::array<double, 4> getIntervals(std::array<LinePtr,2>& lines);
 /**
 * Get a Vertex on this line with parameter t.
 * @param t the parameter
 * @return a vertex on this line with parameter t.
 */
    VertexPtr getVertex(double t) const;

    // sets curveStart, curveEnd, edgeStart, edgeEnd on edge.
    void setIntervalOnEdge(Edge& e,
        double eps,
        double cIntervalStart,
        double cIntervalEnd,
        double vstart,
        double vend);
	static std::array<LinePtr, 2> getEpsilonNeighborhood(Line& vline, double eps);
    /**
   * Computes intersections between eps-region around this line and a line segment when two lines
   * are parallel.
   * @param line is a line parallel to this line
   * @return a optional<array<double,2>> containing two intersection points or empty optional when they don't intersect
   */
    std::optional<std::array<double,2>> getTParallel(Line& line, double eps);
 /**
 * Compute intersection of eps-region around edge e and this line.
 * @param e is the edge around which we would consider eps-region
 * @param eps the radius of the disc
 * @return true when they have intersection or false when they don't intersect.
 */
 // TODO(mahmuda): Add unit test for this method.
    bool pIntersection(Edge& e, double eps);
/**
 * Computes the distance between this line and a point.
 * @param p the vertex from which we will compute distance
 * @return a double value containing distance
 */
    double distance(const Vertex& p);
/**
 * Computes t value on this line for Vertex v, t = 0 at p1 and t = 1 at p2.
 * @return a double value.
 */
    double tValueOnLine(const Vertex& v) const;
/**
 * Check if the vertex, v lies on this line.
 * @return boolean true, if the vertex lies on this line or false otherwise.
 */
    bool onLine(const Vertex& v) const;
};
//using LinePtr = std::shared_ptr<Line>;