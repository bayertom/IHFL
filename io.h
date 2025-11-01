// Description: Read / write the input / output point cloud

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


#ifndef IO_H
#define IO_H

#include <string>
#include <map>
#include "point3d.h"
#include "tvector.h"
#include "tvector2D.h"
#include "facility.h"

//Input - output functions
class IO
{
	public:
		static void loadPointCloud(const std::string& file_name, TVector <Point3D>& U, const double fc = 1.0, const double multiplier = 1.0, const bool non_uniform_cl  = false, const bool point_id_enabled = false);
		static void savePointCloud(const std::string& file_name, const TVector <Point3D>& U);
		static void saveFacilitesIDXs(const std::string& file_name, const TVector2D <int>& idxs_all);
		static void saveClientsToFacilitesIDXs(const std::string& file_name, const std::map <int, int>& idxs_all);
		static void saveFacilitiesNoStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& ID);
		static void saveFacilitiesAndStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& ID, const TVector <int>& NC, const TVector <double>& RAD, const TVector <double>& ABN, const TVector <double>& DFP, const TVector <double>& ASP, const TVector <int>& DIM, const TVector <int>& OVER, const TVector <double>& SLO);

};

#endif
