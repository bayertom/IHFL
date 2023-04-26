// Sort points according to 3D bins (Sloan)

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


#ifndef SortPointsTo3DBins_H
#define SortPointsTo3DBins_H


//Sort points to 3D bins according to their X,Y,Z
struct SortPointsTo3DBins
{
	unsigned int n;
	double xmin, ymin, zmin, xmax, ymax, zmax;

	SortPointsTo3DBins(const double xmin_, const double ymin_, const double zmin_, const double xmax_, const double ymax_, const double zmax_, const int n_) :
		xmin(xmin_), ymin(ymin_), zmin(zmin_), xmax(xmax_), ymax(ymax_), zmax(zmax_), n(n_) {}


	unsigned int getIndexLR1(unsigned int i, unsigned int j, unsigned int k, unsigned int n) const
	{
		//Return index for even i: from left to right
		//Measured from bottom-left element of the grid
		return i * n + j + 1 + k * n * n;
	}

	unsigned int getIndexRL1(unsigned int i, unsigned int j, unsigned int k, unsigned int n) const
	{
		//Return index for even i: from right to left
		//Measured from bottom-left element of the grid
		return (i + 1) * n - j + k * n * n;
	}

	unsigned int getIndexLR2(unsigned int i, unsigned int j, unsigned int k, unsigned int n) const
	{
		//Return index for even i: from left to right
		//Measured from top-right element of the grid
		return n * (n - i - 1) + 1 + j + k * n * n;
	}

	unsigned int getIndexRL2(unsigned int i, unsigned int j, unsigned int k, unsigned int n) const
	{
		//Return index for even i: from right to left
		//Measured from top-right element of the grid
		return n * (n - i) - j + k * n * n;
	}

	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		//Sort points to bins (Analogy to Sloan, 1993)
		const double dmax = std::max(std::max(xmax - xmin, ymax - ymin), zmax - zmin);
		const double x_maxr = (xmax - xmin) / dmax;
		const double y_maxr = (ymax - ymin) / dmax;
		const double z_maxr = (zmax - zmin) / dmax;

		//Reduce coordinates
		const double x1r = (p1.x - xmin) / dmax;
		const double y1r = (p1.y - ymin) / dmax;
		const double z1r = (p1.z - zmin) / dmax;
		const double x2r = (p2.x - xmin) / dmax;
		const double y2r = (p2.y - ymin) / dmax;
		const double z2r = (p2.z - zmin) / dmax;

		//Indices
		unsigned int i1 = 0, j1 = 0, k1 = 0, i2 = 0, j2 = 0, k2 = 0;

		if (x_maxr != 0 && y_maxr != 0 && z_maxr != 0)
		{
			//First point index
			i1 = (unsigned int)(0.99 * n * y1r / y_maxr);
			j1 = (unsigned int)(0.99 * n * x1r / x_maxr);
			k1 = (unsigned int)(0.99 * n * z1r / z_maxr);

			//Second point index
			i2 = (unsigned int)(0.99 * n * y2r / y_maxr);
			j2 = (unsigned int)(0.99 * n * x2r / x_maxr);
			k2 = (unsigned int)(0.99 * n * z2r / z_maxr);
		}

		//Is index "i1" odd or even?
		unsigned int b1, b2;

		//Even row, point p1
		if (k1 % 2 == 0)
		{
			//Left - right direction
			if (i1 % 2 == 0)
			{
				b1 = getIndexLR1(i1, j1, k1, n);
			}

			//Right - left direction
			else
			{
				b1 = getIndexRL1(i1, j1, k1, n);
			}
		}

		//Odd row, point p1
		else
		{
			//Right - left direction
			if (i1 % 2 == 0)
			{
				b1 = getIndexRL2(i1, j1, k1, n);
			}

			//Left - right direction
			else
			{
				b1 = getIndexLR2(i1, j1, k1, n);
			}
		}

		//Even row, point p2
		if (k2 % 2 == 0)
		{
			//Left - right direction
			if (i2 % 2 == 0)
			{
				b2 = getIndexLR1(i2, j2, k2, n);
			}

			//Right - left direction
			else
			{
				b2 = getIndexRL1(i2, j2, k2, n);
			}
		}

		//Odd row, point p2
		else
		{
			//Right - left direction
			if (i2 % 2 == 0)
			{
				b2 = getIndexRL2(i2, j2, k2, n);
			}

			//Left - right direction
			else
			{
				b2 = getIndexLR2(i2, j2, k2, n);
			}
		}

		return b1 < b2;
	}
};

#endif