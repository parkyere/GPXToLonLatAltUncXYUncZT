#pragma once
#include "Vertex.hpp"
#include <cmath>
#include <optional>
#include <utility>
struct Line {
private:
	Vertex p1, p2;
	double xdiff, ydiff, zdiff;
	double c; //y-intersect of the line equation
	double m; // slope of the line
	double theta; //smallest angle between x-axis and this Line object
public:
	Line(Vertex p1, Vertex p2) 
		: p1{ p1 }, p2{ p2 },
		xdiff{ p2.getX() - p1.getX() },
		ydiff{ p2.getY() - p1.getY() },
		zdiff{ p2.getZ() - p1.getZ() } {
/*
 * y = mx + c this equation was used unless this line is parallel to y axis.
 *
 * y1 = mx1 + c y2 = mx2 + c
 *
 * c = ((y1+y2)-m(x1+x2))/2 m = (y2-y1)/(x2-x1)
 */

        if (xdiff != 0.0) {
            m = ydiff / xdiff;
            Vertex vector = Vertex(p2.getX() - p1.getX(), p2.getY() - p1.getY(), 0);
            double angle1 = Vertex::dotProd(vector, Vertex(1.0, 0.0, 0.0))
                / vector.norm();
            double cosAngle1 = std::acos(angle1);
            if (ydiff >= 0) {
                theta = rad2deg(cosAngle1);
            }
            else {
                theta = -1 * rad2deg(cosAngle1);
            }

            c = ((p1.getY() + p2.getY()) - m * (p1.getX() + p2.getX())) / 2.0;

        }
        else if (ydiff > 0.0) {
            theta = rad2deg(std::numbers::pi / 2.0);
        }
        else {
            theta = -rad2deg(std::numbers::pi / 2.0);
        }
	}
	Vertex& getP1() { return p1; }
	Vertex& getP2() { return p2; }
	double getXdiff() const { return xdiff; }
    double getYdiff() const { return ydiff; }
    double getZdiff() const { return zdiff; }
	double getM() const { return m; }
	double getC() const { return c; }
	double getTheta() const { return theta; }
	void setM(double m) { this->m = m; }
	void setC(double c) { this->c = c; }
	void setTheta(double theta) { this->theta = theta; }
	std::string toString() const {
		return std::format("[{0:.8f} {1:.8f} {2:.8f}; {3:.8f} {4:.8f} {5:.8f}]",
			p1.getX(), p1.getY(), p1.getZ(),
			p2.getX(), p2.getY(), p2.getZ());
	}
	double avgAltitude() const {
		return (p1.getZ() + p2.getZ()) / 2.0;
	}
    /**
 * Compute intersection of eps-disc around p and line
 * @param p the center of the disc
 * @param eps the radius of the disc
 * @return a optional<pair<double,double>> containing two intersection points or null when they don't intersect.
 */
	std::optional<std::pair<double, double>> pIntersection(Vertex& p, double eps, bool precisionCompromise) {
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
        double b = 2 * ((p1.getX() - p.getX()) * xdiff + (p1.getY() - p.getY()) * ydiff);
        double c = (p1.getX() - p.getX()) * (p1.getX() - p.getX()) + (p1.getY() - p.getY())
            * (p1.getY() - p.getY()) - eps * eps;
        double a = xdiff * xdiff + ydiff * ydiff;
        if (a < 1e-32) {
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
            double newC = (p1.getX() - p.getX()) * (p1.getX() - p.getX()) + (p1.getY() - p.getY())
                * (p1.getY() - p.getY()) - newEps * newEps;
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

        return { {t[0],t[1]} };
	}
    /**
 * Returns intersection points of this line and two lines in Line[] (neighborhood). For detail see
 * {@link package-info}
 *
 * @return an array of four double values, first two values represents interval start and end
 *         points on this line and next two corresponds to intersections with lines[0] and
 *         lines[1] respectively.
 */
     std::array<double,4> getIntervals(Line lines[2]) {
        double cIntervalStart;
        double cIntervalEnd;
        double neighborhoodStart;
        double neighborhoodEnd;

        if (lines[0].getXdiff() == 0.0) {
            // when the edge is a vertical line
            cIntervalStart = m * lines[0].getP1().getX() + c;
            cIntervalEnd = m * lines[1].getP1().getX() + c;

            neighborhoodStart = (cIntervalStart - lines[0].getP1().getY()) / lines[0].getYdiff();
            neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getY()) / lines[1].getYdiff();
            if (ydiff != 0.0) {
                cIntervalStart = (cIntervalStart - getP1().getY()) / ydiff;
                cIntervalEnd = (cIntervalEnd - getP1().getY()) / ydiff;
            }
            else {
                cIntervalStart = (lines[0].getP1().getX() - p1.getX()) / xdiff;
                cIntervalEnd = (lines[1].getP1().getX() - p1.getX()) / xdiff;
            }

        }
        else if (xdiff == 0.0) {
            // when the line is a vertical line
            cIntervalStart = lines[0].getM() * getP1().getX() + lines[0].getC();
            cIntervalEnd = lines[1].getM() * getP1().getX() + lines[1].getC();

            if (lines[0].getYdiff() != 0.0) {
                neighborhoodStart = (cIntervalStart - lines[0].getP1().getY()) / lines[0].getYdiff();
                neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getY()) / lines[1].getYdiff();
            }
            else {
                neighborhoodStart = (p1.getX() - lines[0].getP1().getX()) / lines[0].getXdiff();
                neighborhoodEnd = (p1.getX() - lines[1].getP1().getX()) / lines[1].getXdiff();
            }
            cIntervalStart = (cIntervalStart - getP1().getY()) / ydiff;
            cIntervalEnd = (cIntervalEnd - getP1().getY()) / ydiff;
        }
        else {

            cIntervalStart = -(c - lines[0].getC()) / (m - lines[0].getM());
            cIntervalEnd = -(c - lines[1].getC()) / (m - lines[1].getM());

            neighborhoodStart = (cIntervalStart - lines[0].getP1().getX()) / lines[0].getXdiff();
            neighborhoodEnd = (cIntervalEnd - lines[1].getP1().getX()) / lines[1].getXdiff();

            cIntervalStart = (cIntervalStart - getP1().getX()) / xdiff;
            cIntervalEnd = (cIntervalEnd - getP1().getX()) / xdiff;
        }
        auto temp = std::array<double,4>();
        temp[0] = cIntervalStart;
        temp[1] = cIntervalEnd;
        temp[2] = neighborhoodStart;
        temp[3] = neighborhoodEnd;
        return temp;
    }
};