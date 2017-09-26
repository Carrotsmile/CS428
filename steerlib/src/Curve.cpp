//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	// Robustness: make sure there is at least two control point: start and end points
	if (checkRobust()) {
		float f = 0.0f;

		Point * p0 = new Point();
		Point * p1 = new Point();

		calculatePoint(*p0, f);
		f += window;

		// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
		while (f < controlPoints[controlPoints.size() - 1].time) {
			calculatePoint(*p1, f);
			DrawLib::drawLine(*p0, *p1, curveColor, curveThickness);
			Point * temp = p0;
			p0 = p1;
			p1 = temp;
			f += window;
		}

		// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve
		calculatePoint(*p1, f);
		DrawLib::drawLine(*p0, *p1, curveColor, curveThickness);

		delete p0;
		delete p1;
	}
	else {
		return;
	}
	
	return;
#endif
}

//comparison function necessary to use standard library sort function
bool static cmpCurvePoints(CurvePoint &c1, CurvePoint &c2) {
	return c1.time < c2.time;
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::sort(controlPoints.begin(), controlPoints.end(), cmpCurvePoints);
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if ((type == hermiteCurve && controlPoints.size() >= 2) ||
		(type == catmullCurve && controlPoints.size() >= 3)) {
		return true;
	}
	else {
		return false;
	}
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	int index = 0;
	for (std::vector<CurvePoint>::iterator it = controlPoints.begin(); it != controlPoints.end(); it++) {
		if (it->time > time) {
			nextPoint = index;
			return true;
		}
		index++;
	}
	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	//this doesn't have to be like "Point newPosition = new Point(...)"?, won't the data get lost?
	Point newPosition;
	float normalTime, intervalTime;

	unsigned int startPoint = nextPoint - 1;

	//changed time so that we could interpolate on an arbitrary interval
	//src for reference:
	//https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_a_single_interval

	intervalTime = (controlPoints[nextPoint].time - controlPoints[startPoint].time);
	normalTime = (time - controlPoints[startPoint].time) / intervalTime;
	float timeCubed = std::pow(normalTime, 3);
	float timeSquared = std::pow(normalTime, 2);

	// Calculate position at t = time on Hermite curve
	newPosition = (2 * timeCubed - 3 * timeSquared + 1) * controlPoints[startPoint].position
		+ (timeCubed - 2 * timeSquared + time) * controlPoints[startPoint].tangent
		+ (0.0 - 2 * timeCubed + 3 * timeSquared) * controlPoints[nextPoint].position
		+ (timeCubed - timeSquared) * controlPoints[nextPoint].tangent;

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	Point p0 = controlPoints[nextPoint - 1].position;
	Point p1 = controlPoints[nextPoint].position;
	Point p2 = controlPoints[nextPoint + 1].position;
	Point p3 = controlPoints[nextPoint + 2].position;

	float time_Cubed = std::pow(time, 3);
	float time_Squared = std::pow(time, 2);

	// Calculate position at t = time on Catmull-Rom curve
	newPosition = 0.5 * ((2 * p1) +
		(-p0 + p2) * time +
		(2 * p0 - 5 * p1 + 4 * p2 - p3) * time_Squared +
		(-p0 + 3 * p1 - 3 * p2 + p3) * time_Cubed);
	
	// Return result
	return newPosition;
}