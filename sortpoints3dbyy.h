#ifndef SortPoints3DByY_H
#define SortPoints3DByY_H

#include "Point3D.h"

//Sort points by Y coordinate
struct SortPoints3DByY
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.y < p2.y;
	}
};

#endif
