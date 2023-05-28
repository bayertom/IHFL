// Description: Split dataset using k-D tree

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


#ifndef SplitDataset_H
#define SplitDataset_H

#include <string>
#include "point3d.h"
#include "tvector.h"

//Split 3D dataset using KD-Tree
class SplitDataset
{
	public:
		static void createKDPointTiles(const std::string& file_name, const int n_max, const double fc, TVector <std::string>& file_names);
		
		static void loadKDPointTileFileNames(const std::string& file_name, TVector<std::string>& file_name_point_tiles);
		static void saveKDPointTileFileNames(const std::string& file_name, const TVector<std::string>& file_name_point_tiles);

	private:
		static void splitPointCloudMedian(const TVector <Point3D>& U, const double median, const int dir, TVector <Point3D>& PL, TVector <Point3D>& PU);

};

#endif