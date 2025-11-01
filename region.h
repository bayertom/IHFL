// Description: Region definition, used for region growing

// Copyright (c) 2025 - 2025
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

#ifndef Region_H
#define Region_H

#include "tvector.h"

//Region definition
struct Region 
{
		double curv_aver;				//Average curvature of the region
		double tx_aver, ty_aver, tz_aver;		//Average tangent vector the region (direction)
		double nx_aver, ny_aver, nz_aver;		//Average normal vector the region
		TVector	<int> points_idx;			//Indices of points forming region

		Region() : curv_aver(0.0), tx_aver(0.0), ty_aver(0.0), tz_aver(0.0), nx_aver(0.0), ny_aver(0.0), nz_aver(0.0) {}	//Constructor

		void add(const int point_idx, const double curv, const double tx, const double ty, const double tz, const double nx, const double ny, const double nz)
		{		
			//Add point to the region and update region properties
			
			//Curvature
			curv_aver = (curv_aver * points_idx.size() + curv) / (points_idx.size() + 1);

			//Tangent vectors (directions)
			tx_aver = (tx_aver * points_idx.size() + tx) / (points_idx.size() + 1);
			ty_aver = (ty_aver * points_idx.size() + ty) / (points_idx.size() + 1);
			tz_aver = (tz_aver * points_idx.size() + tz) / (points_idx.size() + 1);

			//Normal vectors
			nx_aver = (nx_aver * points_idx.size() + nx) / (points_idx.size() + 1);
			ny_aver = (ny_aver * points_idx.size() + ny) / (points_idx.size() + 1);
			nz_aver = (nz_aver * points_idx.size() + nz) / (points_idx.size() + 1);

			//Add point to the region
			points_idx.push_back(point_idx);
		}

		void clear()
		{
			//Clear region and reset its properties
			curv_aver = 0.0;
			tx_aver = 0.0;
			ty_aver = 0.0;
			tz_aver = 0.0;
			nx_aver = 0.0;
			ny_aver = 0.0;
			nz_aver = 0.0;

			points_idx.clear();
		}
};

#endif
