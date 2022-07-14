#include <iostream>
#include <fstream>
#include <stack>
#include <algorithm>

#include "splitdataset.h"

#include "sortpoints3dbyx.h"
#include "sortpoints3dbyy.h"
#include "sortpoints3dbyz.h"
#include "kdpointtile.h"
#include "io.h"


void SplitDataset::splitPointCloudMedian(const TVector <Point3D>& P, const double median, const int dir, TVector <Point3D>& PL, TVector <Point3D>& PU)
{
	//Split set of points in median according to the direction
	for (const auto &p : P)
	{
		//Split according to X
		if (dir == 0)
		{
			if (p.x <= median)
				PL.push_back(p);
			else
				PU.push_back(p);
		}

		//Split according to Y
		else if (dir == 1)
		{
			if (p.y <= median)
				PL.push_back(p);
			else
				PU.push_back(p);
		}

		//Split according to Z
		else
		{
			if (p.z <= median)
				PL.push_back(p);
			else
				PU.push_back(p);
		}
	}
}



void SplitDataset::createKDPointTiles(const std::string& file_name, const int n_max, const double fc, TVector <std::string>& file_names)
{
	//Create subsets from the input point cloud divided in median
	TVector <Point3D> U;
	
	//Load entire point cloud
	std::cout << ">> Loading point cloud: ";
	IO::loadPointCloud(file_name, U, fc);
	std::cout << "OK \n";

	std::cout << ">> Partition to subsets: ";
	KDPointTile kd_tile(0, 0, U);
	
	//Create stack
	std::stack <KDPointTile> S;
	S.push(kd_tile);
	
	//Divide in median until S is empty
	while (!S.empty())
	{
		//Get element on the top of the stack
		KDPointTile kd_tile_top = S.top();
		
		//Subset too large, split
		double median;
		if (kd_tile_top.points.size() > n_max)
		{
			//Find median in x
			if (kd_tile_top.dir == 0)
			{
				std::nth_element(kd_tile_top.points.begin(), kd_tile_top.points.begin() + kd_tile_top.points.size() / 2, kd_tile_top.points.end(), SortPoints3DByX());
				median = kd_tile_top.points[kd_tile_top.points.size() / 2].x;
			}
			
			//Find median in y
			else if (kd_tile_top.dir == 1)
			{
				std::nth_element(kd_tile_top.points.begin(), kd_tile_top.points.begin() + kd_tile_top.points.size() / 2, kd_tile_top.points.end(), SortPoints3DByY());
				median = kd_tile_top.points[kd_tile_top.points.size() / 2].y;
			}

			//Find median in z
			else
			{
				std::nth_element(kd_tile_top.points.begin(), kd_tile_top.points.begin() + kd_tile_top.points.size() / 2, kd_tile_top.points.end(), SortPoints3DByZ());
				median = kd_tile_top.points[kd_tile_top.points.size() / 2].z;
			}
			
			//Split in two data sets
			TVector <Point3D> PL, PU;
			splitPointCloudMedian(kd_tile_top.points, median, kd_tile_top.dir, PL, PU);

			//Change direction
			const int kd_dir = (kd_tile_top.dir + 1) % 3;

			//Increment depth
			const int kd_depth = kd_tile_top.depth + 1;

			//Create 2 KD subsets
			KDPointTile kd_tile_l(kd_depth, kd_dir, PL);
			KDPointTile kd_tile_u(kd_depth, kd_dir, PU);

			//Remove kd_ point tile
			S.pop();

			//Add new tiles to the stack
			S.push(kd_tile_l);
			S.push(kd_tile_u);
			
			std::cout << ".";
		}

		//Point cloud of the acceptable size: store to the file
		else
		{
			//Remove kd_tile
			S.pop();

			//Create file name
			std::string file_name_kdss = file_name + "_" + std::to_string(kd_tile_top.kd_tile_id) + "_" + std::to_string(kd_tile_top.depth) + "_" + std::to_string(kd_tile_top.dir);

			//Save point tile
			IO::savePointCloud(file_name_kdss, kd_tile_top.points);

			//Add file name to the list
			file_names.push_back(file_name_kdss);
		}
		
	}

	std::cout << "OK \n\n";
	
}


void SplitDataset::loadKDPointTileFileNames(const std::string& file_name, TVector<std::string>& file_name_point_tiles)
{
	//Load file names of point tiles from the list
	std::string line;
	std::ifstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Read line by line
		while (std::getline(file, line))
		{
			//Skip empty line
			if (line.size() > 0)
				file_name_point_tiles.push_back(line);
		}
	}

	//Throw exception
	catch (std::ifstream::failure e)
	{
		std::cerr << "Can not open file> :" << file_name << '\n';
	}
}




void SplitDataset::saveKDPointTileFileNames(const std::string& file_name, const TVector<std::string>& file_name_point_tiles)
{
	//Save names of point tiles to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Process all kd point tiles
		for (auto f : file_name_point_tiles)
			file << f << '\n';
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}