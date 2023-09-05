// Description: Spatial Indexing using 3D grid

// Copyright (c) 2010 - 2021
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

#include <iostream>
#include "gridindexing.h"
#include "sortpoints3dbyx.h"
#include "sortpoints3dbyy.h"
#include "sortpoints3dbyz.h"

void GridIndexing::initializeIndex(const TVector <Point3D>& points)
{
	//Create 3D index according to grid
	const Point3D p_xmax = *std::max_element(points.begin(), points.end(), SortPoints3DByX());
	const Point3D p_xmin = *std::min_element(points.begin(), points.end(), SortPoints3DByX());
	const Point3D p_ymax = *std::max_element(points.begin(), points.end(), SortPoints3DByY());
	const Point3D p_ymin = *std::min_element(points.begin(), points.end(), SortPoints3DByY());
	const Point3D p_zmax = *std::max_element(points.begin(), points.end(), SortPoints3DByZ());
	const Point3D p_zmin = *std::min_element(points.begin(), points.end(), SortPoints3DByZ());

	//Get extreme coordinates
	xmin = p_xmin.x;
	const double xmax = p_xmax.x;
	ymin = p_ymin.y;
	const double ymax = p_ymax.y;
	zmin = p_zmin.z;
	const double zmax = p_zmax.z;

	//Compute edge sizes
	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;

	//Compute amount of bins
	if (nx == -1)
	{
		nx = ceil(std::max(dx / bx, 1.0));
		ny = ceil(std::max(dy / by, 1.0));
		nz = ceil(std::max(dz / bz, 1.0));
	}

	//Compute bin size
	else if (bx == -1)
	{
		bx = dx / nx;
		by = dy / ny;
		bz = dz / nz;
	}
}


std::tuple<int, int, int> GridIndexing::get3DIndex(const Point3D& q) const
{
	//Get 3D index of a query point inside the grid
	const double k = 0.99;

	//Reduce coordinates
	const double xr = (dx == 0 ? 0 : (q.x - xmin) / dx);
	const double yr = (dy == 0 ? 0 : (q.y - ymin) / dy);
	const double zr = (dz == 0 ? 0 : (q.z - zmin) / dz);

	//Point outside the min-max box
	if ((xr < 0 || xr > 1) || (yr < 0 || yr > 1) || (zr < 0 || zr > 1))
		return { -1, -1, -1 };

	//Compute indices
	const int ix = int(k * nx * xr);
	const int iy = int(k * ny * yr);
	const int iz = int(k * nz * zr);

	return { ix, iy, iz };
}


int GridIndexing::get1DIndex(const Point3D& q) const 
{
	//Get 3D index of a query point inside the grid
	auto [ix, iy, iz] = get3DIndex(q);

	//Point outside the min-max box
	if ((ix < 0) || (iy < 0) || (iz < 0))
		return -1;

	//Convert 3D index to 1D index
	return index3DTo1D(ix, iy, iz);
}


TVector <int> GridIndexing::getAdjacentIndices(const Point3D& q) const
{
	//Get indices of the adjacent 27 bins to the query point
	TVector <int> idxs;

	//Get 3D index of a query point
	auto [ix, iy, iz] = get3DIndex(q);

	//Point outside the min-max box
	if ((ix < 0) || (iy < 0) || (iz < 0))
		return idxs;

	//Get indices of adjacent bins
	for (int i = std::max(ix - 1, 0); i <= std::min(ix + 1, nx - 1); i++)
	{
		for (int j = std::max(iy - 1, 0); j <= std::min(iy + 1, ny - 1); j++)
		{
			for (int k = std::max(iz - 1, 0); k <= std::min(iz + 1, nz - 1); k++)
			{
				//std::cout << i << " " << j << " " << k << '\n';
				idxs.push_back(index3DTo1D(i, j, k));
				//auto i1 = index3DTo1D(i, j, k);
				//std::cout << i1 << '\n';
			}
		}
	}

	return idxs;
}

