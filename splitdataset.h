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