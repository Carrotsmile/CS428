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

		if (!checkRobust())
		{
			return;
		}
		
		Point * p0 = new Point();
		Point * p1 = new Point(); 
		float f = 0.0f;
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

		/*
		if (!checkRobust())
		{
			return;
		}

		Point start;
		Point end;

		// draw curve from controlPoints[i] to controlPoints[i+1], use 'window' as number of segments.
		for (int i = 0; i < controlPoints.size() - 1; i++)
		{
			end = controlPoints[i].position;
			for (int j = 0; j < window - 1; j++)
			{
				start = end;
				calculatePoint(end, controlPoints[i].time + (float)(j + 1) / (float)(window)* (controlPoints[i + 1].time - controlPoints[i].time));
				DrawLib::drawLine(start, end, curveColor, curveThickness);
			}
			DrawLib::drawLine(end, controlPoints[i + 1].position, curveColor, curveThickness);
		}
		return;
		*/
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
		std::cout << time << std::endl;
		if (it->time <= time && index < controlPoints.size()) {
			index++;
		}
	}
	nextPoint = index;
	if (index == controlPoints.size()) {
		return false;
	} else {
		return true;
	}
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
		+ (timeCubed - 2 * timeSquared + normalTime) * controlPoints[startPoint].tangent
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

	// Calculate position at t = time on Catmull-Rom curve
	
	// Return result
	return newPosition;
}