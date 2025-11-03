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
#include <map>

#include "pfclustmetric.h"
#include "pfregionmetric.h"
#include "ihfl.h"
#include "io.h"
#include "splitdataset.h"
#include "regressionplane.h"
#include "sortpoints3dbyx.h"
#include "sortpoints3dbyy.h"
#include "sortpoints3dbyz.h"
#include "sortpointsto3dbins.h"
#include "dxfexport.h"
#include "region.h"

//#include "rg.h"



int main(int argc, char* argv[])
{
	/*
	//Input data
	//std::string file_namer = "data\\Cone\\cone_10000.txt";
	std::string file_namer = "data\\Cube\\cube_10000.txt";

	//Load points
	TVector <Point3D> points;
	IO::loadPointCloud(file_namer, points);

	//Define thresholds
	TVector <double> edge_thresholds({ 0.05, 0.05, 0.50, 0.50, 0.05 }); //ABN + ABT

	//Detect edge points
	auto [edge_points, edge_points_features] = RG::findEdgePointsCurv(points, 50, edge_thresholds[0]);

	//Perform region growing strategy
	pfRegionMetric pf_edge_metric = &RG::inSameRegion;
	TVector2D <Point3D> regions = RG::regionGrowDyn(edge_points, edge_points_features, pf_edge_metric, edge_thresholds, 50);
	
	//Edges
	for (int i = 0; i < regions.size(); i++)
	{
		std::string file_name_output = "data\\Cube\\edge" + std::to_string(i) + ".txt";
		IO::savePointCloud(file_name_output, regions[i]);
	}
	
	//Export to dxf
	DXFExport::exportPolylineToDXF("polyline_cube2_p1.dxf", regions[0], 1);

	return 0;
	*/
	std::cout << "*** FACILITY LOCATION CLUSTERING WITH PSEUDO-METRICS *** \n";
	std::cout << "***   (ver. 2.7, 11/2025, bayertom@natur.cuni.cz)    *** \n";
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
	//std::string file_name = "e:\\Tomas\\CPP\\IHFL\\data\\Cube\\cube_10000.txt";
	//std::string file_name = "j:\\Tomas\\CPP\\Clustering\Clustering\Tests\IHFLN\Prýmek\Tile\18_5_2.txt";
	//std::string file_name = "j:\\Tomas\\CPP\\Clustering\\Clustering\\Tests\\IHFLN\\Cone\\Test3\\V1\\cone_10000.txt";
	
	//Parameters of the clusterization algorithm
	bool non_uniform_cl = false, export_dxf = false, recompute_fac_cost = false, cluster_statistics = false;
	int knn = 30, ns = 100000, l = 1;
	double fc = 0.1, mul = 1.0, lambda = 1.0, bin = 1.0, mju = 0.95;
	double x_scale = 1.0, y_scale = 1.0, z_scale = 1.0;
	pfClustMetric pf_clust_metric = &IHFL::mL2;
	std::string file_name, pf_clust_metric_text = "l2";

	//pfClustMetric pf_clust_metric = &IHFL::pmABN;
	//std::string file_name, pf_clust_metric_text = "abn";

	//Point3D p(27, 28, 29, 30);
	//std::vector<Point3D> points;
	//points.push_back(p);

	/* Testing data
	IHFL clust(non_uniform_cl, cluster_statistics, knn, lambda, mju, l, pf_clust_metric);
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

			//Set clusterization (pseudo)metric
			if (!strcmp("met", attribute) || !strcmp("norm", attribute))
			{
				//DIS (G1 pseudometric)
				if (!strcmp("dis", value))
				{
					pf_clust_metric = &IHFL::pmDIS;
					pf_clust_metric_text = "dis";

					recompute_fac_cost = true;
				}

				//ABN (G2 pseudometric)
				else if (!strcmp("abn", value))
				{
					pf_clust_metric = &IHFL::pmABN;
					pf_clust_metric_text = "abn";

					recompute_fac_cost = true;
				}

				//ABLP (G3 pseudometric)
				else if (!strcmp("ablp", value))
				{
					pf_clust_metric = &IHFL::pmABLP;
					pf_clust_metric_text = "ablp";

					recompute_fac_cost = true;
				}

				//DFP (G4 pseudometric)
				else if (!strcmp("dfp", value))
				{
					pf_clust_metric = &IHFL::pmDFP;
					pf_clust_metric_text = "dfp";

					recompute_fac_cost = true;
				}

				//L2 metric
				else if (!strcmp("l2", value))
				{
					pf_clust_metric = &IHFL::mL2;
					pf_clust_metric_text = "l2";
				}

				//L1 metric
				else if (!strcmp("l1", value))
				{
					pf_clust_metric = &IHFL::mL1;
					pf_clust_metric_text = "l1";
				}

				//L22 metric
				else if (!strcmp("l22", value))
				{
					pf_clust_metric = &IHFL::mL22;
					pf_clust_metric_text = "l22";
				}

				//Elliptic metric
				else if (!strcmp("ell", value))
				{
					pf_clust_metric = &IHFL::mEll;
					pf_clust_metric_text = "ell";
				}

				//Geographic pseudometric
				else if (!strcmp("geo", value))
				{
					pf_clust_metric = &IHFL::pmGeo;
					pf_clust_metric_text = "geo";
				}

				//Linearity pseudometric
				else if (!strcmp("lin", value))
				{
					pf_clust_metric = &IHFL::pmLIN;
					pf_clust_metric_text = "lin";
				}

				//Planarity pseudometric
				else if (!strcmp("pla", value))
				{
					pf_clust_metric = &IHFL::pmPLA;
					pf_clust_metric_text = "pla";
				}

				//Sphericity pseudometric
				else if (!strcmp("sph", value))
				{
					pf_clust_metric = &IHFL::pmSPH;
					pf_clust_metric_text = "sph";
				}

				//Omnivariance pseudometric
				else if (!strcmp("omn", value))
				{
					pf_clust_metric = &IHFL::pmOMN;
					pf_clust_metric_text = "omn";
				}

				//Anisotropy pseudometric
				else if (!strcmp("ani", value))
				{
					pf_clust_metric = &IHFL::pmANI;
					pf_clust_metric_text = "ani";
				}

				//Curvature change pseudometric
				else if (!strcmp("cur", value))
				{
					pf_clust_metric = &IHFL::pmCUR;
					pf_clust_metric_text = "cur";
				}

				//Eigen entropy pseudometric
				else if (!strcmp("ent", value))
				{
					pf_clust_metric = &IHFL::pmENT;
					pf_clust_metric_text = "ent";
				}

				//Verticality pseudometric
				else if (!strcmp("ver", value))
				{
					pf_clust_metric = &IHFL::pmVER;
					pf_clust_metric_text = "ver";
				}

				//Horizontality pseudometric
				else if (!strcmp("hor", value))
				{
					pf_clust_metric = &IHFL::pmHOR;
					pf_clust_metric_text = "hor";
				}

				//Exception
				else
					throw std::exception("Exception: Invalid clusterization metric/norm type in command line!");
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
	std::cout << "\nParameters: unif = " << !non_uniform_cl << ", rec_cost = " << recompute_fac_cost << ", met = " << pf_clust_metric_text << ", f_cost = " << fc << ", mul = " << mul << ", lambda = " << lambda << ", bin = " << bin << ", l = " << l << ", mju = " << mju << ", knn = " << knn << ", n_split = " << ns << ".\n\n";

	//Load list of kd point tiles, if exists
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
	std::string file_name_part = "_" + pf_clust_metric_text + "_l" + std::to_string(l) + "_unif" + std::to_string(!non_uniform_cl) + "_fc" + std::format("{:.2f}", fc) + "_lam" +
		std::format("{:.2f}", lambda) + "_bin" + std::format("{:.2f}", bin) + "_mju" + std::format("{:.2f}", mju) + "_knn" + std::to_string(knn) + "_mult" + std::to_string(mul);

	//Process point tiles one by one
	unsigned int i = 0, n_subsets = file_name_point_tiles.size();
	TVector <int> ID_all, NC_all, OVER_all, DIM_all;
	std::map <int, int> clients_to_facilities_cloud_all;
	TVector <double> RAD_all, ABN_all, DFP_all, ASP_all, SLO_all;
	TVector <Point3D> output_facilities_all;
	TVector2D <int> cloud_idxs_all;

	for (const auto &f_name : file_name_point_tiles)
	{
		//Load KD-point tiles, do not change facility costs
		TVector <Point3D> kd_points_tile, output_facilities_tile;
		IO::loadPointCloud(f_name, kd_points_tile, fc, 1, true, true);

		//Extract file name from the path
		std::string base_f_name = f_name.substr(f_name.find_last_of("/\\") + 1);

		//Empty file
		if (kd_points_tile.size() == 0)
			continue;

		std::cout << "\nFile: " << f_name << " (" << ++i << "/" << n_subsets << ", " << kd_points_tile.size() << " points)" << '\n';

		//Output files
            	std::string facil_file_subset = f_name + file_name_part + "_facil.txt";
		std::string facil_file_subset_stat = f_name + file_name_part + "_facil_stat.txt";
		std::string clients_to_facil_file_subset = f_name + file_name_part + "_facil2.txt";
		std::string facil_and_clients_file_subset = f_name + file_name_part + "_facil3.txt";
		std::string dxf_file_subset = f_name + file_name_part + ".dxf";
		std::string result_file_subset = f_name + file_name_part + "_results.log";

		//Output file with facilities does not exist
		if (!std::filesystem::exists(facil_file_subset))
		{
			//Find extreme points
			const Point3D p_xmax = *std::max_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByX());
			const Point3D p_xmin = *std::min_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByX());
			const Point3D p_ymax = *std::max_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByY());
			const Point3D p_ymin = *std::min_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByY());
			const Point3D p_zmax = *std::max_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByZ());
			const Point3D p_zmin = *std::min_element(kd_points_tile.begin(), kd_points_tile.end(), SortPoints3DByZ());

			//Get first point ID
			const unsigned int id0 = kd_points_tile.front().id;

			//Sort points to 3D bins
			sort(kd_points_tile.begin(), kd_points_tile.end(), SortPointsTo3DBins(p_xmin.x, p_ymin.y, p_zmin.z, p_xmax.x, p_ymax.y, p_zmax.z, pow(kd_points_tile.size(), 1.0 / 6)));
			
			//Start measure time
			const clock_t begin_time = clock();

			//Initialize the grid index
			GridIndexing gi(bin, bin, bin);
			gi.initializeIndex(kd_points_tile);

			//Auxilliary structures
			TVector <RegressionPlane> RP;
			TVector2D <Facility> FG(gi.nx * gi.ny * gi.nz);

			//Incremental heuristic location (IHFL)
			IHFL clust(non_uniform_cl, recompute_fac_cost, knn, lambda, bin, mju, l, x_scale, y_scale, z_scale, pf_clust_metric);
			clust.clusterizeIHFL(kd_points_tile, gi, FG, RP);

			//Convert indexed grid of facilities (2D vector) to 1D vector
			TVector <Facility> F;
			for (auto fg : FG)
				F.insert(F.end(), std::make_move_iterator(fg.begin()), std::make_move_iterator(fg.end()));

			//Create output facilities
			for (int j = 0; j < F.size(); j++)
				output_facilities_tile.push_back(kd_points_tile[abs(F[j].p_idx) - 1]);

			//Statistics
			const double time = (clock() - begin_time) / (CLOCKS_PER_SEC);
			std::cout << "(time: " << time << "s, ";
			std::cout << F.size() << " facilities). ";
			std::cout << "OK \n";

			//Convert facilities to indices
			TVector <int> ID = clust.outputFaciliesToIDXs(kd_points_tile, FG);

			//Save facilities to txt file
			if (cluster_statistics)
			{
				TVector <int> NC, OVER, DIM;
				TVector <double> RAD, ABN, DFP, ASP, SLO;
				
				std::cout << ">> Compute statistics: ";

				//Compute parameters of clusters
				clust.clusterStatistics(kd_points_tile, FG, gi, RP, NC, RAD, ABN, DFP, ASP, DIM, OVER, SLO);

				//Add results to the lists
				NC_all.insert(NC_all.end(), NC.begin(), NC.end());
				OVER_all.insert(OVER_all.end(), OVER.begin(), OVER.end());
				DIM_all.insert(DIM_all.end(), DIM.begin(), DIM.end());
				RAD_all.insert(RAD_all.end(), RAD.begin(), RAD.end());
				ABN_all.insert(ABN_all.end(), ABN.begin(), ABN.end());
				DFP_all.insert(DFP_all.end(), DFP.begin(), DFP.end());
				ASP_all.insert(ASP_all.end(), ASP.begin(), ASP.end());
				SLO_all.insert(SLO_all.end(), SLO.begin(), SLO.end());

				//Store facilities and their statistics for the tile
				IO::saveFacilitiesAndStatistics(facil_file_subset_stat, output_facilities_tile, ID, NC, RAD, ABN, DFP, ASP, DIM, OVER, SLO);
			
				std::cout << "OK \n";
			}

			//Save facilities without statistics
			else
			{
				IO::saveFacilitiesNoStatistics(facil_file_subset, output_facilities_tile, ID);
			}
			
			//Export clusters to DXF
			if (export_dxf)
			{
				std::cout << ">> Exporting to DXF file: ";

				DXFExport::exportClustersToDXF(dxf_file_subset, F, kd_points_tile, RP);

				std::cout << "OK \n";
			}

			//Export facilities and clients
			TVector2D <int> cloud_idxs = clust.facilitiesToIDXs(kd_points_tile, F);

			//Save facilities and assigned points
			IO::saveFacilitesIDXs(facil_and_clients_file_subset, cloud_idxs);

			//Export facilities to the list
			//Format: row = index, number = centroid _id;
			std::map <int, int> clients_to_facilities_cloud = clust.clientsToFacilitiesIDXs(kd_points_tile, F);

			//Add results to the list
			ID_all.insert(ID_all.end(), ID.begin(), ID.end());
			cloud_idxs_all.insert(cloud_idxs_all.end(), cloud_idxs.begin(), cloud_idxs.end());
			clients_to_facilities_cloud_all.insert(clients_to_facilities_cloud.begin(), clients_to_facilities_cloud.end());
			output_facilities_all.insert(output_facilities_all.end(), output_facilities_tile.begin(), output_facilities_tile.end());

			//Compute total cost
			const double total_cost = clust.computeCost(F, kd_points_tile, RP);
			std::cout << "Total cost: " << total_cost << '\n';
		}

		//Output file with resulted facilities exists, load it as the point cloud
		else
		{
			IO::loadPointCloud(facil_file_subset, output_facilities_tile, fc);
		}
	}


	//Save all facilities and their statistics into file
	std::string facil_stat_file_all = file_name + file_name_part + "_facil_stat_all.txt";
	
	//Save facilities + statistics
	if (cluster_statistics)
		IO::saveFacilitiesAndStatistics(facil_stat_file_all, output_facilities_all, ID_all,  NC_all, RAD_all, ABN_all, DFP_all, ASP_all, DIM_all, OVER_all, SLO_all);

	//Save all facilities without statistics
	else
		IO::saveFacilitiesNoStatistics(facil_stat_file_all, output_facilities_all, ID_all);

	//Save points with assigned facilities
	std::string clients_to_facil_all_file = file_name + file_name_part + "_facil2_all.txt";
	IO::saveClientsToFacilitesIDXs(clients_to_facil_all_file, clients_to_facilities_cloud_all);

	//Save all facilities
	std::string facil_and_clients_all_file = file_name + file_name_part + "_facil3_all.txt";
	IO::saveFacilitesIDXs(facil_and_clients_all_file, cloud_idxs_all);
	
	std::cout << "Finished... \n";

	return 0;
}
