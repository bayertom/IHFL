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


int main(int argc, char* argv[])
{
	std::cout << "*** HYBRID CLUSTERIZATION OF THE POINT CLOUD *** \n";

	//Testing data
	//std::string file_name = "data\\eth_mid.txt";
	//std::string file_name = "data\\test_pseudometrics.txt";
	//std::string file_name = "data\\eth_mid.txt";
	//std::string file_name = "data\\Cone\\cone_10000.txt";
	//std::string file_name = "data\\Ferrata\\via_ferrata_xyz_rgb_small.txt"; 
	//std::string file_name = "data2\\test_pseudometrics.txt";

	//Parameters of the clusterization algorithm
	bool non_uniform_cl = false, export_dxf = false, cluster_statistics = false;
	int knn = 50, ns = 100000, l = 1;
	double fc = 0.01, lambda = 0.25, mju = 0.95;
	pfnorm fnorm = &IHFL::nDFP;
	std::string file_name, fnorm_text = "dfp", method_text = "ihfl";
	/*
	IHFL clust(non_uniform_cl, cluster_statistics, knn, lambda, mju, l, fnorm);
	TVector <Point3D> U;
	clust.generateCone(10, 5, 1000, U);
	IO::savePointCloud("E:\\Tomas\\CPP\\IHFL\\cone2.txt", U);
	*/
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
				}

				//ABN (G2 pseudonorm)
				else if (!strcmp("abn", value))
				{
					fnorm = &IHFL::nABN;
					fnorm_text = "abn";
				}

				//ABLP (G3 pseudonorm)
				else if (!strcmp("ablp", value))
				{
					fnorm = &IHFL::nABLP;
					fnorm_text = "ablp";
				}

				//DFP (G4 pseudonorm)
				else if (!strcmp("dfp", value))
				{
					fnorm = &IHFL::nDFP;
					fnorm_text = "dfp";
				}

				//L2
				else if (!strcmp("l2", value))
				{
					fnorm = &IHFL::nL2;
					fnorm_text = "l2";
				}

				else
					throw std::exception("Exception: Invalid clusterization norm type in command line!");
			}

			//Set the facility cost
			else if (!strcmp("fc", attribute))
			{
				fc = std::max(std::min(atof(value), 1000.0), 0.0001);
			}

			//Set maximum radius ball
			else if (!strcmp("lambda", attribute))
			{
				lambda = std::max(std::min(atof(value), 10.0), 0.01);
			}

			//Split cloud to subsets (point tiles)
			else if (!strcmp("ns", attribute))
			{
				ns = std::max(std::min(atoi(value), 500000), 100);
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
	std::cout << "\nParameters: norm = " << fnorm_text << ", method = " << method_text << ", f_cost = " << fc << ", lambda = " << lambda << ", l = " << l << ", mju = " << mju << ", max_clust = " << lambda << ", knn = " << knn << ", n_split = " << ns << ".\n\n";

	//Load list of kd point tiles, if they exist
	TVector<std::string> file_name_point_tiles;
	SplitDataset::loadKDPointTileFileNames(file_list, file_name_point_tiles);

	//Otherwise, create list of point tiles
	if (file_name_point_tiles.size() == 0)
	{
		//Split point cloud to point tiles
		SplitDataset::createKDPointTiles(file_name, ns, fc, file_name_point_tiles);
		
		//Save file names of tiles
		SplitDataset::saveKDPointTileFileNames(file_list, file_name_point_tiles);
	}

	//Point tiles have been loaded
	else
	{
		std::cout << ">> KD subsets loaded (" << file_name_point_tiles.size() << ") ...";
	}

	//Create file name part unique for all files
	std::string file_name_part = "_" + fnorm_text + "_l" + std::to_string(l) + "_unif" + std::to_string(!non_uniform_cl) + "_f" + std::format("{:.2f}", fc) + "_dc" + std::format("{:.2f}", lambda);

	//Process point tiles one by one
	unsigned int i = 0, n_subsets = file_name_point_tiles.size();
	TVector <Point3D> output_facilities;
	for (const auto &f_name : file_name_point_tiles)
	{
		//Load KD-point tiles
		TVector <Point3D> kd_point_tile, output_facilities_tile;
		IO::loadPointCloud(f_name, kd_point_tile, fc);

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

			//Apply facility location clusterization
			TVector <Facility> F;
			TVector <RegressionPlane> AN;

			//Auxilliary variables
			const clock_t begin_time = clock();

			//Incremental heuristic location (IFL)
			IHFL clust(non_uniform_cl, knn, lambda, mju, l, fnorm);
			clust.clusterizeIHFL(kd_point_tile, fc, F, AN);

			//Statistics
			const double time = (clock() - begin_time) / (CLOCKS_PER_SEC);
			std::cout << "(time: " << time << "s, ";
			std::cout << F.size() << " facilities). ";
			std::cout << "OK \n\n";

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

			//Export clusters to DXF
			if (export_dxf)
				DXFExport::exportClustersToDXF(dxf_file_subset, F, kd_point_tile, AN);

			//Save facilities to txt file
			if (cluster_statistics)
			{
				TVector <int> NC, OVER, DIM;
				TVector <double> RAD, ABN, DFP, ASP;

				//Compute parameters of clusters
				clust.clusterStatistics(kd_point_tile, F, AN, NC, RAD, ABN, DFP, ASP, DIM, OVER);

				IO::savePointCloudAndStatistics(facil_file_subset_stat, output_facilities_tile, NC, RAD, ABN, DFP, ASP, DIM, OVER);
			}

			else
				IO::savePointCloud(facil_file_subset, output_facilities_tile);

			//Save facilities with assigned point
			IO::saveClientsToFacilites(clients_to_facil_file_subset, clients_to_facilities_tile);


			//Compute total cost
			const double total_cost = clust.computeCost(F, kd_point_tile, AN);
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


	//Save output points
	std::string facil_file = file_name + file_name_part + "_facil_all.txt";
	IO::savePointCloud(facil_file, output_facilities);
	
	std::cout << "Finished... \n";

	return 0;
}
