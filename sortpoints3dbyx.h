#ifndef SortPoints3DByX_H
#define SortPoints3DByX_H

#include "Point3D.h"

//Sort points by X coordinate
struct SortPoints3DByX
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.x < p2.x;
	}
};

#endif