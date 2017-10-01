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
	
	if (checkRobust())
	{
		Point p0;
		Point p1;

		float delTime = 0;
		//draw curve between each pair of points
		for (int i = 0; i < controlPoints.size() - 1; i++)
		{
			p1 = controlPoints[i].position;
			//draw line segments using window as time stepper
			for (delTime = 0.0f; delTime < controlPoints[i + 1].time - controlPoints[i].time; delTime += window)
			{
				p0 = p1;
				calculatePoint(p1, controlPoints[i].time + delTime);
				DrawLib::drawLine(p0, p1, curveColor, curveThickness);
			}
			DrawLib::drawLine(p1, controlPoints[i+1].position, curveColor, curveThickness);
		}
	}
	else {
		return;
	}
	

#endif
}

bool static eqCurvePoints(const CurvePoint &c1, const CurvePoint &c2) {
	return c1.time == c2.time;
}
//comparison function necessary to use standard library sort function
bool static cmpCurvePoints(const CurvePoint &c1, const CurvePoint &c2) {
	return c1.time < c2.time;
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	//sorts control points according to time
	std::sort(controlPoints.begin(), controlPoints.end(), cmpCurvePoints);

	//removes duplicates (control points with equal time) from the list
	controlPoints.erase(unique(controlPoints.begin(), controlPoints.end(), eqCurvePoints), controlPoints.end());
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
	if (controlPoints.size() >= 2 && type == hermiteCurve) {
		return true;
	} 
	else if (controlPoints.size() >= 3 && type == catmullCurve) {
		return true;
	}
	else {
		return false;
	}
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	int index = 1;
	//next point is the first one with time more than the current time
	for (index = 1; index < controlPoints.size(); index++) {
		if (controlPoints[index].time > time) {
			nextPoint = index;
			return true;
		}
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
	float timeCubed = normalTime * normalTime * normalTime;
	float timeSquared = normalTime * normalTime;
	
	
	
	// Calculate position at t = time on Hermite curve
	newPosition = (2 * timeCubed - 3 * timeSquared + 1) * controlPoints[startPoint].position
		+ (timeCubed - 2 * timeSquared + normalTime) * (intervalTime * controlPoints[startPoint].tangent)
		+ (0.0 - 2 * timeCubed + 3 * timeSquared) * controlPoints[nextPoint].position
		+ (timeCubed - timeSquared) * (intervalTime * controlPoints[nextPoint].tangent);
	
	// Return result
	return newPosition;
	
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	//same source for reference as the hermite curve, only difference being how the tangents are calculated
	Point newPosition;
	CurvePoint p0 = controlPoints[nextPoint - 1];
	CurvePoint p1 = controlPoints[nextPoint];
	Vector v0, v1;

	// Calculate position at t = time on Catmull-Rom curve
	float normalTime, intervalTime;
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - p0.time) / intervalTime;
	float timeSquared = normalTime * normalTime;
	float timeCubed = normalTime * normalTime * normalTime;

	//select values for v0 and v1 based on if nextPoint is last, first, or neither
	//first point
	if (nextPoint == 1) {
		CurvePoint p2 = controlPoints[nextPoint + 1];
		v0 = (p1.position - p0.position);
		v1 = (p2.position - p0.position) / 2;
	} 
	//last point
	else if (nextPoint == controlPoints.size() - 1) {
		CurvePoint p2 = controlPoints[nextPoint - 2];
		v0 = (p1.position - p2.position)/2;
		v1 = (p1.position - p0.position);
	}
	//neithe
	else {
		CurvePoint p2 = controlPoints[nextPoint - 2];
		CurvePoint p3 = controlPoints[nextPoint + 1];
		v0 = (p1.position - p2.position) / 2;
		v1 = (p3.position - p0.position) / 2;
	}

	//wow
	newPosition = (2 * timeCubed - 3 * timeSquared + 1) * p0.position +
		(-2 * timeCubed + 3 * timeSquared) * p1.position +
		(timeCubed - 2 * timeSquared + normalTime) * v0 +
		(timeCubed - timeSquared) * v1;
	// Return result
	return newPosition;
}