// Description: Region growing algorithm

// Copyright (c) 2021 - 2025
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

#include <string>

#include "tvector.h"
#include "point3d.h"
#include "io.h"
#include "rg.h"


void mainrg()
{
	//Input data
	std::string file_name = "data\\Cone\\cone_10000.txt";

	//Load points
	TVector <Point3D> points;
	IO::loadPointCloud(file_name, points);

	//Perform region growing
	TVector2D <Point3D> regions;
	//RG::regionGrow(points, 0.5, 0.1, 30, regions);

}

