#pragma once
#include "Line.hpp"
#include "Edge.hpp"
//#include "Vertex.hpp"

Line::Line(VertexPtr p1, VertexPtr p2)
    : p1{ p1 }, p2{ p2 },
    xdiff{ p2->getX() - p1->getX() },
    ydiff{ p2->getY() - p1->getY() },
    zdiff{ p2->getZ() - p1->getZ() } {
    /*
     * y = mx + c this equation was used unless this line is parallel to y axis.
     *
     * y1 = mx1 + c y2 = mx2 + c
     *
     * c = ((y1+y2)-m(x1+x2))/2 m = (y2-y1)/(x2-x1)
     */

    if (!isEqual(xdiff, 0.0)) {
        m = ydiff / xdiff;
        //VertexPtr vector{ new Vertex(p2->getX() - p1->getX(), p2->getY() - p1->getY(), 0) };
        //double angle1 = Vertex::dotProd(*vector, Vertex(1.0, 0.0, 0.0))
        //    / vector->norm();
		//double angle1 = vector->getX() / vector->norm();
        //double cosAngle1 = std::acos(angle1);
        double cosAngle1 = std::acos(xdiff/(std::sqrt(xdiff*xdiff + ydiff*ydiff)));
        if (ydiff >= 0) {
            theta = rad2deg(cosAngle1);
        }
        else {
            theta = -1 * rad2deg(cosAngle1);
        }
        c = ((p1->getY() + p2->getY()) - m * (p1->getX() + p2->getX())) / 2.0;
    }
    else if (ydiff > 0.0) {
        theta = rad2deg(std::numbers::pi / 2.0);
    }
    else {
        theta = -rad2deg(std::numbers::pi / 2.0);
    }
}
 VertexPtr Line::getP1() { return p1; }
 VertexPtr Line::getP2() { return p2; }
 double Line::getXdiff() const { return xdiff; }
 double Line::getYdiff() const { return ydiff; }
 double Line::getZdiff() const { return zdiff; }
 double Line::getM() const { return m; }
 double Line::getC() const { return c; }
 double Line::getTheta() const { return theta; }
 void Line::setM(double m) { this->m = m; }
 void Line::setC(double c) { this->c = c; }
 void Line::setTheta(double theta) { this->theta = theta; }
 std::string Line::toString() const {
    return std::format("[{0:.16f} {1:.16f} {2:.16f}; {3:.16f} {4:.16f} {5:.16f}]",
        p1->getX(), p1->getY(), p1->getZ(),
        p2->getX(), p2->getY(), p2->getZ());
}
 double Line::avgAltitude() const {
    return (p1->getZ() + p2->getZ()) / 2.0;
}

/**
* Compute intersection of eps-disc around p and line
* @param p the center of the disc
* @param eps the radius of the disc
* @return a optional<array<double,2>> containing two intersection points or null when they don't intersect.
*/
 std::optional<std::array<double, 2>> Line::pIntersection(const Vertex& p, double eps, bool precisionCompromise) const {
    /*
    * Line 1: x1+(x2-x1)*t = x; x1+xdiff*t = x y1+(y2-y1)*t = y; y1+ydiff*t = y
    *
    * Equation of disc around vertex p: (x-a)^2+(y-b)^2=eps^2
    * (x1+xdiff*t-a)^2+(y1+ydiff*t-b)^2=eps^2 (xdiff^2+ydiff^2)t^2 + (2(x1-a)xdiff+2(y1-b)ydiff)t +
    * (x1-a)^2+(y1-b)^2-eps^2=0
    *
    * Quadratic solution for t, gives us intersection of disc points with line
    *
    * t = (-(2(x1-a)xdiff+2(y1-b)ydiff)) +- sqrt((2(x1-a)xdiff+2(y1-b)ydiff))^2 -
    * 4(xdiff^2+ydiff^2)(x1^2+y1^2-eps^2)))/(2*(xdiff^2+ydiff^2))
    *
    * a*t^2 + b*t + c = 0
    */
    double t[2];
    double b = 2 * ((p1->getX() - p.getX()) * xdiff + (p1->getY() - p.getY()) * ydiff);
    double c = (p1->getX() - p.getX()) * (p1->getX() - p.getX()) + (p1->getY() - p.getY())
        * (p1->getY() - p.getY()) - eps * eps;
    double a = xdiff * xdiff + ydiff * ydiff;
    if (isEqual(a, 0.0)) {
        return {};
    }
    double determinant = b * b - 4 * a * c;
    double dist = std::sqrt(a);
    if (determinant >= 0) {
        t[1] = (-b + std::sqrt(determinant)) / (2 * a);
        t[0] = (-b - std::sqrt(determinant)) / (2 * a);
    }
    else if (precisionCompromise) {
        double newEps = eps + 0.1;
        double newC = (p1->getX() - p.getX()) * (p1->getX() - p.getX()) + (p1->getY() - p.getY())
            * (p1->getY() - p.getY()) - newEps * newEps;
        double newDeterminant = b * b - 4 * a * newC;
        if (newDeterminant >= 0) {
            t[1] = (-b + std::sqrt(newDeterminant)) / (2 * a);
            t[0] = (-b - std::sqrt(newDeterminant)) / (2 * a);
        }
        else {
            return {};
        }
    }
    else {
        return {};
    }
    double min = std::min(t[0], t[1]);
    double max = std::max(t[0], t[1]);

    t[0] = min;
    t[1] = max;

    return { { t[0],t[1] } };
}

/**
* Returns intersection points of this line and two lines in Line[] (neighborhood). For detail see
* {@link package-info}
*
* @return an array of four double values, first two values represents interval start and end
*         points on this line and next two corresponds to intersections with lines[0] and
*         lines[1] respectively.
*/
 std::array<double, 4> Line::getIntervals(std::array<LinePtr, 2>& lines) {
    double cIntervalStart;
    double cIntervalEnd;
    double neighborhoodStart;
    double neighborhoodEnd;

    if (isEqual(lines[0]->getXdiff(), 0.0)) {
        // when the edge is a vertical line
        cIntervalStart = m * lines[0]->getP1()->getX() + c;
        cIntervalEnd = m * lines[1]->getP1()->getX() + c;

        neighborhoodStart = (cIntervalStart - lines[0]->getP1()->getY()) / lines[0]->getYdiff();
        neighborhoodEnd = (cIntervalEnd - lines[1]->getP1()->getY()) / lines[1]->getYdiff();
        if (!isEqual(ydiff, 0.0)) {
            cIntervalStart = (cIntervalStart - getP1()->getY()) / ydiff;
            cIntervalEnd = (cIntervalEnd - getP1()->getY()) / ydiff;
        }
        else {
            cIntervalStart = (lines[0]->getP1()->getX() - p1->getX()) / xdiff;
            cIntervalEnd = (lines[1]->getP1()->getX() - p1->getX()) / xdiff;
        }
    }
    else if (isEqual(xdiff, 0.0)) {
        // when the line is a vertical line
        cIntervalStart = lines[0]->getM() * getP1()->getX() + lines[0]->getC();
        cIntervalEnd = lines[1]->getM() * getP1()->getX() + lines[1]->getC();

        if (!isEqual(lines[0]->getYdiff(), 0.0)) {
            neighborhoodStart = (cIntervalStart - lines[0]->getP1()->getY()) / lines[0]->getYdiff();
            neighborhoodEnd = (cIntervalEnd - lines[1]->getP1()->getY()) / lines[1]->getYdiff();
        }
        else {
            neighborhoodStart = (p1->getX() - lines[0]->getP1()->getX()) / lines[0]->getXdiff();
            neighborhoodEnd = (p1->getX() - lines[1]->getP1()->getX()) / lines[1]->getXdiff();
        }
        cIntervalStart = (cIntervalStart - getP1()->getY()) / ydiff;
        cIntervalEnd = (cIntervalEnd - getP1()->getY()) / ydiff;
    }
    else {

        cIntervalStart = -(c - lines[0]->getC()) / (m - lines[0]->getM());
        cIntervalEnd = -(c - lines[1]->getC()) / (m - lines[1]->getM());

        neighborhoodStart = (cIntervalStart - lines[0]->getP1()->getX()) / lines[0]->getXdiff();
        neighborhoodEnd = (cIntervalEnd - lines[1]->getP1()->getX()) / lines[1]->getXdiff();

        cIntervalStart = (cIntervalStart - getP1()->getX()) / xdiff;
        cIntervalEnd = (cIntervalEnd - getP1()->getX()) / xdiff;
    }
    return { cIntervalStart , cIntervalEnd , neighborhoodStart, neighborhoodEnd };
}

/**
* Get a Vertex on this line with parameter t.
* @param t the parameter
* @return a vertex on this line with parameter t.
*/
 VertexPtr Line::getVertex(double t) const {
    return VertexPtr{ new Vertex(p1->getX() + xdiff * t, p1->getY() + ydiff * t,
        p1->getZ() + zdiff * t) };
}
// sets curveStart, curveEnd, edgeStart, edgeEnd on edge.
 void Line::setIntervalOnEdge(Edge& e, double eps, double cIntervalStart, double cIntervalEnd, double vstart, double vend) {
	double interval[2];
	interval[0] = std::max(0.0, cIntervalStart);
	if (isEqual(vstart, -1)) {
		auto in1 = e.getLine()->pIntersection(*(getVertex(interval[0])), eps, true);
		if (!(in1.has_value())) {
			//logger.log(Level.SEVERE, "Problem computing Line intersection: in1.");
			throw std::runtime_error("Problem computing Line intersection: in1.");
		}

		if ((*in1)[0] >= 0 && (*in1)[0] <= 1 && (*in1)[1] >= 0 && (*in1)[1] <= 1) {
			vstart = ((*in1)[0] + (*in1)[1]) / 2;
		}
		else if ((*in1)[0] >= 0 && (*in1)[0] <= 1) {
			vstart = (*in1)[0];
		}
		else if ((*in1)[1] >= 0 && (*in1)[1] <= 1) {
			vstart = (*in1)[1];
		}
		else {
			vstart = 0;
		}
	}
	interval[1] = std::min(1.0, cIntervalEnd);
	if (isEqual(vend, -1)) {
		auto in2 = e.getLine()->pIntersection(*(getVertex(interval[1])), eps, true);
		if (!(in2.has_value())) {
			//logger.log(Level.SEVERE, "Problem computing Line intersection: in2.");
			throw std::runtime_error("Problem computing Line intersection: in2.");
		}
		if ((*in2)[0] >= 0 && (*in2)[0] <= 1 && (*in2)[1] >= 0 && (*in2)[1] <= 1) {
			vend = ((*in2)[0] + (*in2)[1]) / 2;
		}
		else if ((*in2)[0] >= 0 && (*in2)[0] <= 1) {
			vend = (*in2)[0];
		}
		else if ((*in2)[1] >= 0 && (*in2)[1] <= 1) {
			vend = (*in2)[1];
		}
		else {
			vend = 1;
		}
	}
	e.setCurveStart(interval[0]);
	e.setCurveEnd(interval[1]);
	e.setEdgeStart(vstart);
	e.setEdgeEnd(vend);
}

 std::array<LinePtr, 2> Line::getEpsilonNeighborhood(Line& vline, double eps) {
    // compute the equations of boundaries of eps-region around the line
    std::array<LinePtr, 2> lines;

    double dTheta;
    if (!isEqual(vline.getXdiff(), 0.0)) {
        dTheta = std::atan(vline.getM()) + std::numbers::pi / 2.0;
    }
    else if (vline.ydiff > 0.0) {
        dTheta = std::numbers::pi / 2.0;
    }
    else {
        dTheta = -std::numbers::pi / 2.0;
    }
    double dx, dy;
    dx = eps * std::cos(dTheta);
    dy = eps * std::sin(dTheta);

    lines[0] = LinePtr{ new Line(
        VertexPtr{ new Vertex(vline.getP1()->getX() - dx, vline.getP1()->getY() - dy, vline.getP1()->getZ()) },
        VertexPtr{ new Vertex(vline.getP2()->getX() - dx, vline.getP2()->getY() - dy, vline.getP2()->getZ()) }) };
    lines[1] = LinePtr{ new Line(
        VertexPtr{ new Vertex(vline.getP1()->getX() + dx, vline.getP1()->getY() + dy, vline.getP1()->getZ()) },
        VertexPtr{ new Vertex(vline.getP2()->getX() + dx, vline.getP2()->getY() + dy, vline.getP2()->getZ()) }) };

    if (!isEqual(lines[0]->getM(), lines[1]->getM())) {
        lines[0]->setM(lines[1]->getM());
        lines[0]->setTheta(lines[1]->getTheta());
    }
    return lines;
}


/**
* Computes intersections between eps-region around this line and a line segment when two lines
* are parallel.
* @param line is a line parallel to this line
* @return a optional<array<double,2>> containing two intersection points or empty optional when they don't intersect
*/
 std::optional<std::array<double, 2>> Line::getTParallel(Line& line, double eps) {
    double t[2];
    double newm;
    double x1, y1, x2, y2;
    if (isEqual(std::abs(line.getTheta()), std::numbers::pi / 2)) {
        newm = 0;
        x1 = p1->getX();
        y1 = line.getP1()->getY();
        x2 = p2->getX();
        y2 = line.getP2()->getY();
    }
    else if (isEqual(line.getTheta(), 0.0)) {
        newm = 0;
        x1 = line.getP1()->getX();
        y1 = p1->getY();
        x2 = line.getP2()->getX();
        y2 = p2->getY();
    }
    else {
        newm = 1 / line.getM();
        double c1 = line.getP1()->getY() + newm * line.getP1()->getX();
        double c2 = line.getP2()->getY() + newm * line.getP2()->getX();

        x1 = (c1 - c) / (m + newm);
        y1 = m * x1 + c;

        x2 = (c2 - c) / (m + newm);
        y2 = m * x2 + c;
    }

    if (std::sqrt((line.getP1()->getX() - x1) * (line.getP1()->getX() - x1)
        + (line.getP1()->getY() - y1) * (line.getP1()->getY() - y1)) > eps) {
        return {};
    }
    double intersection1;
    double intersection2;

    if (!isEqual(xdiff, 0.0)) {
        intersection1 = (x1 - p1->getX()) / xdiff;
        intersection2 = (x2 - p1->getX()) / xdiff;
    }
    else {
        intersection1 = (y1 - p1->getY()) / ydiff;
        intersection2 = (y2 - p1->getY()) / ydiff;
    }
    t[0] = std::min(intersection1, intersection2);
    t[1] = std::max(intersection1, intersection2);

    if (t[1] < 0 || t[0] > 1) {
        return {};
    }

    t[0] = std::max(t[0], 0.0);
    t[1] = std::min(t[1], 1.0);
    return { { t[0],t[1] } };
}

/**
* Compute intersection of eps-region around edge e and this line.
* @param e is the edge around which we would consider eps-region
* @param eps the radius of the disc
* @return true when they have intersection or false when they don't intersect.
*/
// TODO(mahmuda): Add unit test for this method.
 bool Line::pIntersection(Edge& e, double eps) {
    /*
    * Line 1: x1+(x2-x1)*t = x y1+(y2-y1)*t = y
    * Line 2: y = mx + c
    * y1+(y2-y1)*t = (x1+(x2-x1)*t)*m + c (y2-y1)*t - (x2-x1)*t*m = x1*m + c- y1
    * t = (x1*m + c - y1)/((y2-y1)-(x2-x1)*m)
    */
    LinePtr vline = e.getLine();

    auto lines = Line::getEpsilonNeighborhood(*vline, eps);
    double cIntervalStart;
    double cIntervalEnd;
    double neighborhoodStart;
    double neighborhoodEnd;

    double vstart = -1;
    double vend = -1;

    if (isEqual(theta, vline->getTheta())
        || (isEqual(std::abs(theta - vline->getTheta()), 180.0))) {// For parallel lines

        auto t = getTParallel(*vline, eps);

        if (!(t.has_value())) {
            return false;
        }

        cIntervalStart = (*t)[0];
        cIntervalEnd = (*t)[1];

        t = vline->pIntersection(*(getVertex(cIntervalStart)), eps, true);
        if (!(t.has_value())) {
            return false;
        }
        neighborhoodStart = ((*t)[0] + (*t)[1]) / 2.0;
        t = vline->pIntersection(*(getVertex(cIntervalEnd)), eps, true);
        if (!(t.has_value())) {
            return false;
        }
        neighborhoodEnd = ((*t)[0] + (*t)[1]) / 2.0;

    }
    else {
        auto temp = getIntervals(lines);
        cIntervalStart = temp[0];
        cIntervalEnd = temp[1];
        neighborhoodStart = temp[2];
        neighborhoodEnd = temp[3];
    }


    if (cIntervalStart > cIntervalEnd) {
        double temp = cIntervalStart;
        cIntervalStart = cIntervalEnd;
        cIntervalEnd = temp;

        temp = neighborhoodStart;
        neighborhoodStart = neighborhoodEnd;
        neighborhoodEnd = temp;
    }

    // computing intersection with endpoint p1
    auto interval1 = pIntersection(*(vline->getP1()), eps, false);
    // computing intersection with endpoint p2
    auto interval2 = pIntersection(*(vline->getP2()), eps, false);

    double minInterval1 = 0;
    double minInterval2 = 0;
    double maxInterval1 = 1;
    double maxInterval2 = 1;
    // line doesn't intersect either of the eps-disc at end points
    if ((!(interval1.has_value())) && (!(interval2.has_value()))) {
        // intersection of line and eps-neighborhood of e is non-empty
        if (cIntervalStart > 1 || cIntervalEnd < 0) {
            return false;
        }
        // intersection of line and eps-neighborhood of e is empty
        if ((neighborhoodStart > 1 && neighborhoodEnd > 1)
            || (neighborhoodStart < 0 && neighborhoodEnd < 0)) {
            return false;
        }
        // intersection needs to be in the eps-region
        if (neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
            && neighborhoodEnd <= 1) {
            cIntervalStart = std::max(0.0, cIntervalStart);
            cIntervalEnd = std::min(1.0, cIntervalEnd);
        }
        else {
            return false;
        }
    }
    // line doesn't intersect with eps-disc of first endpoint

    if (interval1.has_value()) {
        minInterval1 = std::min((*interval1)[0], (*interval1)[1]);
        maxInterval1 = std::max((*interval1)[0], (*interval1)[1]);

        if (((neighborhoodStart > 1 && neighborhoodEnd > 1)
            || (neighborhoodStart < 0 && neighborhoodEnd < 0))
            && (minInterval1 > 1 || maxInterval1 < 0)) {
            return false;
        }
        if (neighborhoodStart < 0) {
            if (minInterval1 <= 1) {
                cIntervalStart = std::max(cIntervalStart, minInterval1);
            }
            if (isEqual(cIntervalStart, minInterval1)) {
                vstart = 0;
            }
        }
        if (neighborhoodEnd < 0) {
            if (maxInterval1 >= 0) {
                cIntervalEnd = std::min(cIntervalEnd, maxInterval1);
            }
            if (isEqual(cIntervalEnd, maxInterval1)) {
                vend = 0;
            }
        }
    }
    // line doesn't intersect with eps-disc of second end point
    if (interval2.has_value()) {
        minInterval2 = std::min((*interval2)[0], (*interval2)[1]);
        maxInterval2 = std::max((*interval2)[0], (*interval2)[1]);

        if (((neighborhoodStart > 1 && neighborhoodEnd > 1)
            || (neighborhoodStart < 0 && neighborhoodEnd < 0))
            && (minInterval2 > 1 || maxInterval2 < 0)) {
            return false;
        }

        if (neighborhoodStart > 1) {
            if (minInterval2 <= 1) {
                cIntervalStart = std::max(cIntervalStart, minInterval2);
            }
            if (isEqual(cIntervalStart, minInterval2)) {
                vstart = 1;
            }
        }

        if (neighborhoodEnd > 1) {
            if (maxInterval2 >= 0) {
                cIntervalEnd = std::min(cIntervalEnd, maxInterval2);
            }
            if (isEqual(cIntervalEnd, maxInterval2)) {
                vend = 1;
            }
        }
    }

    if ((!interval1.has_value()) && (interval2.has_value())) {
        if (!(neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
            && neighborhoodEnd <= 1)) {
            if (neighborhoodStart >= 0 && neighborhoodStart <= 1) {
                if (cIntervalStart > 1 && minInterval2 > 1) {
                    return false;
                }
                if (cIntervalStart < 0 && maxInterval2 < 0) {
                    return false;
                }
            }
            if (neighborhoodEnd >= 0 && neighborhoodEnd <= 1) {
                if (cIntervalEnd > 1 && minInterval2 > 1) {
                    return false;
                }
                if (cIntervalEnd < 0 && maxInterval2 < 0) {
                    return false;
                }
            }
        }
    }
    if ((interval1.has_value()) && (!interval2.has_value())) {
        if (!(neighborhoodStart >= 0 && neighborhoodStart <= 1 && neighborhoodEnd >= 0
            && neighborhoodEnd <= 1)) {

            if (neighborhoodStart >= 0 && neighborhoodStart <= 1) {
                if (cIntervalStart > 1 && minInterval1 > 1) {
                    return false;
                }
                if (cIntervalStart < 0 && maxInterval1 < 0) {
                    return false;
                }
            }
            if (neighborhoodEnd >= 0 && neighborhoodEnd <= 1) {
                if (cIntervalEnd > 1 && minInterval1 > 1) {
                    return false;
                }
                if (cIntervalEnd < 0 && maxInterval1 < 0) {
                    return false;
                }
            }
        }
    }

    if (cIntervalStart > cIntervalEnd) {
        double temp = cIntervalStart;
        cIntervalStart = cIntervalEnd;
        cIntervalEnd = temp;

        temp = vend;
        vend = vstart;
        vstart = temp;
    }

    if (cIntervalStart > 1 || cIntervalEnd < 0) {
        return false;
    }

    if (std::max(maxInterval1, maxInterval2) < 0 || std::min(minInterval1, minInterval2) > 1) {
        return false;
    }
    setIntervalOnEdge(e, eps, cIntervalStart, cIntervalEnd, vstart, vend);
    return true;
}



/**
* Computes the distance between this line and a point.
* @param p the vertex from which we will compute distance
* @return a double value containing distance
*/
 double Line::distance(const Vertex& p) {
    double distance = 0;
    if (!isEqual(xdiff, 0)) {
        distance =
            std::abs(-m * p.getX() + p.getY() + c) / std::sqrt(std::pow(m, 2) + 1);
    }
    else {
        distance = std::abs(p1->getX() - p.getX());
    }

    auto t = pIntersection(p, distance, false);

    if ((!t.has_value()) || (*t)[0] > 1 || (*t)[0] < 0) {
        return std::min(std::sqrt((p1->getX() - p.getX()) * (p1->getX() - p.getX())
            + (p1->getY() - p.getY()) * (p1->getY() - p.getY())), std::sqrt((p2->getX() - p.getX())
                * (p2->getX() - p.getX()) + (p2->getY() - p.getY()) * (p2->getY() - p.getY())));
    }
    else {
        return distance;
    }
}

/**
* Computes t value on this line for Vertex v, t = 0 at p1 and t = 1 at p2.
* @return a double value.
*/
 double Line::tValueOnLine(const Vertex& v) const {
    if (isEqual(xdiff, 0.0)) {
        return (v.getY() - p1->getY()) / ydiff;
    }
    else {
        return (v.getX() - p1->getX()) / xdiff;
    }
}

/**
* Check if the vertex, v lies on this line.
* @return boolean true, if the vertex lies on this line or false otherwise.
*/
 bool Line::onLine(const Vertex& v) const {

    if (isEqual(xdiff, 0.0) && isEqual(v.getX(), p1->getX())) {
        return true;
    }
    else {
        return isEqual((v.getY() - m * v.getX() - c), 0.0);
    }
}
