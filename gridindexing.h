// Description: Spatial Indexing using 3D grid  (formed by bins)

// Copyright (c) 2022 - 2023
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


#ifndef GridIndexing_H
#define GridIndexing_H

 
#include "point3d.h"
#include "tvector.h"

//Spatial indexing according to 3D grid (formed by bins)
class GridIndexing
{
	public:

		double bx;			//Bin size in X direction
		double by;			//Bin size in Y direction
		double bz;			//Bin size in Z direction
		int nx;				//Amount of bins in X direction
		int ny;				//Amount of bins in Y direction
		int nz;				//Amount of bins in Z direction
		double dx;			//Grid size in X direction
		double dy;			//Grid size in Y direction
		double dz;			//Grid size in Z direction
		double xmin;			//Minimum X coordinate of the grid point
		double ymin;			//Minimum Y coordinate of the grid point
		double zmin;			//Minimum Z coordinate of the grid point

	public:
		GridIndexing(const double bx_, const double by_, const double bz_) : bx(bx_), by(by_), bz(bz_), nx(-1), ny(-1), nz(-1),
			dx(0.0), dy(0.0), dz(0.0), xmin(0.0), ymin(0.0), zmin(0.0) {}
		
		GridIndexing(const int nx_, const int ny_, const int nz_) : bx(-1.0), by(-1.0), bz(-1.0), nx(nx_), ny(ny_), nz(nz_),
			dx(0.0), dy(0.0), dz(0.0), xmin(0.0), ymin(0.0), zmin(0.0) {}

		void initializeIndex(const TVector <Point3D>& points);

		std::tuple<int, int, int> get3DIndex(const Point3D &q) const;
		int get1DIndex(const Point3D& q) const;
		TVector <int> getAdjacentIndices (const Point3D& q) const;
		int index3DTo1D(const int ix, const int iy, const int iz) const { return ix + iy * nx + iz * nx * ny;}
};

#endif
