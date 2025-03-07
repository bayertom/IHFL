// Description: Incremental heuristic facility location clustering of the point cloud
// Main cpp file: read data, perform clusterization, write data

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


#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include <iostream>
#include <filesystem>
#include <format>

#include "pfnorm.h"
#include "ihfl.h"
#include "io.h"
#include "splitdataset.h"
#include "regressionplane.h"
#include "sortpoints3dbyx.h"
#include "sortpoints3dbyy.h"
#include "sortpoints3dbyz.h"
#include "sortpointsto3dbins.h"
#include "dxfexport.h"
//#include "rg.h"

int main(int argc, char* argv[])
{
	/*
	//Input data
	std::string file_namer = "data\\Cone\\cone_10000.txt";

	//Load points
	TVector <Point3D> points;
	IO::loadPointCloud(file_namer, points);

	//Perform region growing
	TVector2D <Point3D> regions;
	rg::regionGrow(points, 0.1, 0.1, 30, regions);

	//Edges
	for (int i = 0; i < regions.size(); i++)
	{
		std::string file_name_output = "data\\Cone\\edge" + std::to_string(i) + ".txt";
		IO::savePointCloud(file_name_output, regions[i]);
	}
	*/

	std::cout << "*** FACILITY LOCATION CLUSTERING WITH PSEUDO-METRICS *** \n";
	std::cout << "***   (ver. 2.6, 02/2025, bayertom@natur.cuni.cz)    *** \n";
	//Testing data
	//std::string file_name = "data\\test_pseudometrics.txt";
	//std::string file_name = "data\\ETH\\eth_mid.txt";
	//std::string file_name = "data\\Cone\\cone_10000.txt";
	//std::string file_name = "data\\Ferrata\\via_ferrata_xyz_rgb_small.txt"; 
	//std::string file_name = "data2\\test_pseudometrics.txt";
	//std::string file_name = "data\\Sima\\Canopy\\canopy_small.txt";
	//std::string file_name = "data\\MP\\boulder_small.txt";
	//std::string file_name = "data\\Plane\\points1000_1.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\Boulders\\data_boulder\\boulder4\\sub\\boulder4.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\IHFL\\data\\IK\\tree14.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\LBruha\\IHFL\\2000m\\test2.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\LBruha\\IHFLN\\10000m\\spopnz_non_unif.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\LS\\NU\\1000\\SPOPNZ_LK2212_A_C.txt";
	//std::string file_name = "e:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\LS\\NU\\250\\SPOPNZ_PocOB_AllAdr_2212_A_C_gt_0.txt";

	//Parameters of the clusterization algorithm
	bool non_uniform_cl = false, export_dxf = false, recompute_fac_cost = false, cluster_statistics = false;
	int knn = 50, ns = 100000, l = 1;
	double fc = 1.0, mul = 1.0, lambda = 1.0, bin = 1.0, mju = 0.95;
	double x_scale = 1.0, y_scale = 1.0, z_scale = 1.0;
	pfnorm fnorm = &IHFL::nL2;
	//pfnorm fnorm = &IHFL::nABN;
	//std::string /*/file_name,*/ fnorm_text = "abn";
	std::string file_name, fnorm_text = "l2";

	/* Testing data
	IHFL clust(non_uniform_cl, cluster_statistics, knn, lambda, mju, l, fnorm);
	TVector <Point3D> U;
	clust.generateCone(10, 5, 1000, U);
	IO::savePointCloud("E:\\Tomas\\CPP\\IHFL\\cone2.txt", U);
	*/

	//Not enough command line parameters
	//if (argc < 2)
	//	return 0;

	//Process command line parameters
	while (--argc > 0)
	{
		//Get - (A parameter follows)
		if (*argv[argc] == '-')
		{
			//Process parameter after -
			for (char* parameter = argv[argc]; ; )
			{
				switch (*(++parameter))
				{
					//Export to DXF
					case 'e':
					{
						export_dxf = true;
						break;
					}

					//Non-uniform clusterization, recompute facility costs according to the normals
					case 'n':
					{
						non_uniform_cl = true;
						break;
					}


					//Recompute cost of clusters according to normals
					case 'r':
					{
						recompute_fac_cost = true;
						break;
					}

					//Compute parameters of the clusters
					case 's':
					{
						cluster_statistics = true;
						break;
					}

					//Terminate character \0 of the argument
					case '\0':
						break;

					//Throw exception
					default:
						throw std::exception("Exception: Invalid parameter in command line!");
				}

				//Stop processing of the argument
				break;
			}
		}

		//Set values
		else if (*argv[argc] == '+')
		{
			//Get new command line parameter
			char* attribute = const_cast <char*> (argv[argc] + 1), * value = NULL;
			char* command = attribute;

			//Find splitter: =
			for (; (*command != '=') && (*command != '\0'); command++);

			//We found splitter, trim command and copy to value
			if ((*command == '=') && (*command != '\0'))
			{
				*command = '\0';
				value = command + 1;
			}

			//Throw exception
			if (attribute == NULL || value == NULL)
				throw std::exception("Exception: Invalid value in command line!");

			//Set clusterization norm
			if (!strcmp("norm", attribute))
			{
				//DIS (G1 pseudonorm)
				if (!strcmp("dis", value))
				{
					fnorm = &IHFL::nDIS;
					fnorm_text = "dis";

					recompute_fac_cost = true;
				}

				//ABN (G2 pseudonorm)
				else if (!strcmp("abn", value))
				{
					fnorm = &IHFL::nABN;
					fnorm_text = "abn";

					recompute_fac_cost = true;
				}

				//ABLP (G3 pseudonorm)
				else if (!strcmp("ablp", value))
				{
					fnorm = &IHFL::nABLP;
					fnorm_text = "ablp";

					recompute_fac_cost = true;
				}

				//DFP (G4 pseudonorm)
				else if (!strcmp("dfp", value))
				{
					fnorm = &IHFL::nDFP;
					fnorm_text = "dfp";

					recompute_fac_cost = true;
				}

				//L2 norm
				else if (!strcmp("l2", value))
				{
					fnorm = &IHFL::nL2;
					fnorm_text = "l2";
				}

				//L1 norm
				else if (!strcmp("l1", value))
				{
					fnorm = &IHFL::nL1;
					fnorm_text = "l1";
				}

				//L22 norm
				else if (!strcmp("l22", value))
				{
					fnorm = &IHFL::nL22;
					fnorm_text = "l22";
				}

				//Elliptic norm
				else if (!strcmp("ell", value))
				{
					fnorm = &IHFL::nEll;
					fnorm_text = "ell";
				}

				//Geographic norm
				else if (!strcmp("geo", value))
				{
					fnorm = &IHFL::nGeo;
					fnorm_text = "geo";
				}

				//Linearity pseudonorm
				else if (!strcmp("lin", value))
				{
					fnorm = &IHFL::nLIN;
					fnorm_text = "lin";
				}

				//Planarity pseudonorm
				else if (!strcmp("pla", value))
				{
					fnorm = &IHFL::nPLA;
					fnorm_text = "pla";
				}

				//Sphericity pseudonorm
				else if (!strcmp("sph", value))
				{
					fnorm = &IHFL::nSPH;
					fnorm_text = "sph";
				}

				//Omnivariance pseudonorm
				else if (!strcmp("omn", value))
				{
					fnorm = &IHFL::nOMN;
					fnorm_text = "omn";
				}

				//Anisotropy pseudonorm
				else if (!strcmp("ani", value))
				{
					fnorm = &IHFL::nANI;
					fnorm_text = "ani";
				}

				//Curvature change pseudonorm
				else if (!strcmp("cur", value))
				{
					fnorm = &IHFL::nCUR;
					fnorm_text = "cur";
				}

				//Eigen entropy pseudonorm
				else if (!strcmp("ent", value))
				{
					fnorm = &IHFL::nENT;
					fnorm_text = "ent";
				}

				//Verticality pseudonorm
				else if (!strcmp("ver", value))
				{
					fnorm = &IHFL::nVER;
					fnorm_text = "ver";
				}

				//Horizontality pseudonorm
				else if (!strcmp("hor", value))
				{
					fnorm = &IHFL::nHOR;
					fnorm_text = "hor";
				}

				//Exception
				else
					throw std::exception("Exception: Invalid clusterization norm type in command line!");
			}

			//Set the facility cost
			else if (!strcmp("fc", attribute))
			{
				fc = std::max(std::min(atof(value), 10000000.0), 0.0001);
			}

			//Set the facility cost multiplier
			else if (!strcmp("mul", attribute))
			{
				mul = std::max(std::min(atof(value), 1000000.0), 0.0001);
			}

			//Set maximum radius ball
			else if (!strcmp("lambda", attribute))
			{
				lambda = std::max(std::min(atof(value), 100000.0), 0.01);
				std::cout << "lam:" << lambda;
			}
			
			//Set size of the bin of the 3D grid (spatial indexing)
			else if (!strcmp("bin", attribute))
			{
				bin = std::max(std::min(atof(value), 100000.0), 0.01);
				//bin = std::max(std::min(atof(value), 10.0 * lambda), 0.25 * lambda);
			}

			//Split cloud to subsets (point tiles)
			else if (!strcmp("ns", attribute))
			{
				ns = std::max(std::min(atoi(value), 10000000), 100);
			}

			//Amount of knn
			else if (!strcmp("knn", attribute))
			{
				knn = std::max(std::min(atoi(value), 500), 5);
			}

			//Set isotropic ratio
			else if (!strcmp("mju", attribute))
			{
				mju = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set impact of penalization (power)
			else if (!strcmp("l", attribute))
			{
				l = std::max(std::min(atof(value), 3.0), 0.5);
			}

			//Set scale factor in the X direction
			else if (!strcmp("a", attribute))
			{
				x_scale = std::max(std::min(atof(value), 100.0), 0.01);
			}

			//Set scale factor in the Y direction
			else if (!strcmp("b", attribute))
			{
				y_scale = std::max(std::min(atof(value), 100.0), 0.01);
			}

			//Set scale factor in the Z direction
			else if (!strcmp("c", attribute))
			{
				z_scale = std::max(std::min(atof(value), 100.0), 0.01);
			}

			//Bad argument
			else
			{
				//std::cout << attribute << '\n';
				throw std::exception("Exception: Invalid attribute in command line!");
			}
		}

		//Process file
		else
		{
			file_name = argv[argc];
		}
	}
	
	//List of files
	std::string file_list = file_name + ".list";

	//Print parameters
	std::cout << "\nParameters: unif = " << !non_uniform_cl << ", rec_cost = " << recompute_fac_cost << ", norm = " << fnorm_text << ", f_cost = " << fc << ", mul = " << mul << ", lambda = " << lambda << ", bin = " << bin << ", l = " << l << ", mju = " << mju << ", knn = " << knn << ", n_split = " << ns << ".\n\n";

	//Load list of kd point tiles, if they exist
	TVector<std::string> file_name_point_tiles;
	SplitDataset::loadKDPointTileFileNames(file_list, file_name_point_tiles);

	//Otherwise, create list of point tiles
	if (file_name_point_tiles.size() == 0)
	{
		//Split point cloud to point tiles
		SplitDataset::createKDPointTiles(file_name, ns, fc, mul, non_uniform_cl, file_name_point_tiles);
		
		//Save file names of tiles
		SplitDataset::saveKDPointTileFileNames(file_list, file_name_point_tiles);
	}

	//Point tiles have been loaded
	else
	{
		std::cout << ">> KD subsets loaded (" << file_name_point_tiles.size() << ") ...";
	}

	//Create file name part unique for all files
	std::string file_name_part = "_" + fnorm_text + "_l" + std::to_string(l) + "_unif" + std::to_string(!non_uniform_cl) + "_fc" + std::format("{:.2f}", fc) + "_lam" + std::format("{:.2f}", lambda) + "_bin" + std::format("{:.2f}", bin) + "_mju" + std::format("{:.2f}", mju) + "_knn" + std::to_string(knn);

	//Process point tiles one by one
	unsigned int i = 0, n_subsets = file_name_point_tiles.size();
	TVector <int> NC_all, OVER_all, DIM_all;
	TVector <double> RAD_all, ABN_all, DFP_all, ASP_all, SLO_all;
	TVector <Point3D> output_facilities;
	for (const auto &f_name : file_name_point_tiles)
	{
		//Load KD-point tiles, do not change facility costs
		TVector <Point3D> kd_point_tile, output_facilities_tile;
		IO::loadPointCloud(f_name, kd_point_tile, fc, 1, true);

		//Empty file
		if (kd_point_tile.size() == 0)
			continue;

		std::cout << "\nFile: " << f_name << " (" << ++i << "/" << n_subsets << ", " << kd_point_tile.size() << " points)" << '\n';

		//Output files
            	std::string facil_file_subset = f_name + file_name_part + "_facil.txt";
		std::string facil_file_subset_stat = f_name + file_name_part + "_facil_stat.txt";
		std::string clients_to_facil_file_subset = f_name + file_name_part + "_facil2.txt";
		std::string dxf_file_subset = f_name + file_name_part + ".dxf";
		std::string result_file_subset = f_name + file_name_part + "_results.log";

		//Output file with facilities does not exist
		if (!std::filesystem::exists(facil_file_subset))
		{
			//Find extreme points
			const Point3D p_xmax = *std::max_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByX());
			const Point3D p_xmin = *std::min_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByX());
			const Point3D p_ymax = *std::max_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByY());
			const Point3D p_ymin = *std::min_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByY());
			const Point3D p_zmax = *std::max_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByZ());
			const Point3D p_zmin = *std::min_element(kd_point_tile.begin(), kd_point_tile.end(), SortPoints3DByZ());

			//Get first point ID
			unsigned int ids = kd_point_tile.front().id;

			//Sort points to 3D bins
			sort(kd_point_tile.begin(), kd_point_tile.end(), SortPointsTo3DBins(p_xmin.x, p_ymin.y, p_zmin.z, p_xmax.x, p_ymax.y, p_zmax.z, pow(kd_point_tile.size(), 1.0 / 6)));

			//Change point IDS
			for (int i = 0; i < kd_point_tile.size(); i++)
				kd_point_tile[i].id = ids + i;

			//Start measure time
			const clock_t begin_time = clock();

			//Initialize the grid index
			GridIndexing gi(bin, bin, bin);
			gi.initializeIndex(kd_point_tile);

			//Auxilliary structures
			TVector <RegressionPlane> RP;
			TVector2D <Facility> FG(gi.nx* gi.ny* gi.nz);

			//Incremental heuristic location (IHFL)
			IHFL clust(non_uniform_cl, recompute_fac_cost, knn, lambda, bin, mju, l, x_scale, y_scale, z_scale, fnorm);
			clust.clusterizeIHFL(kd_point_tile, fc, gi, FG, RP);

			//Convert indexed grid of facilities (2D vector) to 1D vector
			TVector <Facility> F;
			for (auto fg : FG)
				F.insert(F.end(), std::make_move_iterator(fg.begin()), std::make_move_iterator(fg.end()));

			//Statistics
			const double time = (clock() - begin_time) / (CLOCKS_PER_SEC);
			std::cout << "(time: " << time << "s, ";
			std::cout << F.size() << " facilities). ";
			std::cout << "OK \n";

			//Export facilities to the list
			TVector <int> clients_to_facilities_tile(kd_point_tile.size());

			for (int i = 0; i < F.size(); i++)
			{
				output_facilities_tile.push_back(kd_point_tile[abs(F[i].p_idx) - 1]);

				//Browse all connected clients
				for (int c_id : F[i].U_idxs)
				{
					clients_to_facilities_tile[abs(c_id) - 1] = abs(F[i].p_idx) - 1;
				}
			}

			//Save facilities to txt file
			if (cluster_statistics)
			{
				TVector <int> NC, OVER, DIM;
				TVector <double> RAD, ABN, DFP, ASP, SLO;
				
				std::cout << ">> Compute statistics: ";

				//Compute parameters of clusters
				clust.clusterStatistics(kd_point_tile, FG, gi, RP, NC, RAD, ABN, DFP, ASP, DIM, OVER, SLO);

				//Add results to the list
				NC_all.insert(NC_all.end(), NC.begin(), NC.end());
				OVER_all.insert(OVER_all.end(), OVER.begin(), OVER.end());
				DIM_all.insert(DIM_all.end(), DIM.begin(), DIM.end());
				RAD_all.insert(RAD_all.end(), RAD.begin(), RAD.end());
				ABN_all.insert(ABN_all.end(), ABN.begin(), ABN.end());
				DFP_all.insert(DFP_all.end(), DFP.begin(), DFP.end());
				ASP_all.insert(ASP_all.end(), ASP.begin(), ASP.end());
				SLO_all.insert(SLO_all.end(), SLO.begin(), SLO.end());

				//Store facilities and their statistics for the tile
				IO::savePointCloudAndStatistics(facil_file_subset_stat, output_facilities_tile, NC, RAD, ABN, DFP, ASP, DIM, OVER, SLO);
			
				std::cout << "OK \n";
			}
			else
				IO::savePointCloud(facil_file_subset, output_facilities_tile);


			//Export clusters to DXF
			if (export_dxf)
			{
				std::cout << ">> Exporting to DXF: ";

				DXFExport::exportClustersToDXF(dxf_file_subset, F, kd_point_tile, RP);

				std::cout << "OK \n";
			}


			//Save facilities with assigned point
			IO::saveClientsToFacilites(clients_to_facil_file_subset, clients_to_facilities_tile);

			//Compute total cost
			const double total_cost = clust.computeCost(F, kd_point_tile, RP);
			std::cout << "Total cost: " << total_cost << '\n';
		}

		//Output file with facilities exists, load point cloud
		else
		{
			IO::loadPointCloud(facil_file_subset, output_facilities_tile, fc);
		}

		//Add all output facilities to the list
		output_facilities.insert(output_facilities.end(), output_facilities_tile.begin(), output_facilities_tile.end());
	}


	//Save all facilities and their statistics into file
	std::string facil_file = file_name + file_name_part + "_facil_all.txt";
	
	if (cluster_statistics)
		IO::savePointCloudAndStatistics(facil_file, output_facilities, NC_all, RAD_all, ABN_all, DFP_all, ASP_all, DIM_all, OVER_all, SLO_all);

	//Save all facilities into file
	else
		IO::savePointCloud(facil_file, output_facilities);
	
	std::cout << "Finished... \n";

	return 0;
}
