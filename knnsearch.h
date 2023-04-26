// Description: KNN search using the nanoflann library

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


#ifndef KNNSearch_H
#define KNNSearch_H

#include "point3d.h"
#include "tvector.h"
#include "tvector2d.h"

//KNN search using Nano-Flann library
class KNNSearch
{
	private:
		TVector <Point3D> points;

	public:
		KNNSearch(const TVector <Point3D>& points_) : points(points_) {}
		void findAllKNN(const TVector <Point3D>& qpoints, const TVector<int>& K, TVector2D <size_t>& knn_id, TVector2D <float>& knn_dist);
};

#endif