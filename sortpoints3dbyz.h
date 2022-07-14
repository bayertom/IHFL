#ifndef SortPoints3DByZ_H
#define SortPoints3DByZ_H

#include "Point3D.h"

//Sort points by Z coordinate
struct SortPoints3DByZ
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.z < p2.z;
	}
};

#endif
