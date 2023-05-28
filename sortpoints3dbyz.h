// Sort 3D points according to Z coordinate

// Copyright (c) 2021 - 2023
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


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
