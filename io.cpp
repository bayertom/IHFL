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


#include "io.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "tvector.h"

void IO::loadPointCloud(const std::string& file_name, TVector <Point3D>& U, const double fc, const double multiplier, const bool non_uniform_cl, const bool point_id_enabled)
{
	//Load pointcloud from the file
	int index = 0;
	std::string line;
	std::ifstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Read line by line
		while (std::getline(file, line))
		{
			char* cline = &line[0u];
			const char* item = strtok(cline, " \t");

			//Delimit the row
			TVector <std::string> row;
			while (item != NULL)
			{
				row.push_back(item);
				item = strtok(NULL, " \t");
			}

			//Add point P = [x, y, z] to the list
			if (row.size() == 3)
			{
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), multiplier * fc));
			}

			//Add point P = [x, y, z, fc] to the list
			else if (row.size() == 4)
			{
				//Add point P = [id, x, y, z] to the list
				if (point_id_enabled)
				{ 
					U.push_back(Point3D(std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), multiplier * fc));
				}

				//Add point P = [x, y, z, fc] to the list
				else
				{
					const double fac_cost = multiplier * (non_uniform_cl ? std::stod(row[3]) : fc);
					U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fac_cost));
				}
			}

			//Add point P = [x, y, z, r, g, b] to the list
			else if (row.size() == 6)
			{
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), multiplier * fc, std::stoi(row[3]), std::stoi(row[4]), std::stoi(row[5])));
			}

			//Add point P = [x, y, z, fc, r, g, b] to the list
			else if (row.size() == 7)
			{
				//Add point P = [id, x, y, z, r, g, b] to the list
				if (point_id_enabled)
				{
					U.push_back(Point3D(std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), multiplier * fc, std::stoi(row[4]), std::stoi(row[5]), std::stoi(row[6])));

				}

				//Add point P = [x, y, z, fc, r, g, b] to the list
				else
				{
					const double fac_cost = multiplier * (non_uniform_cl ? std::stod(row[3]) : fc);
					U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fac_cost, std::stoi(row[4]), std::stoi(row[5]), std::stoi(row[6])));
				}
			}

			//Add point P = [id, x, y, z, fc, r, g, b] to the list
			else if (row.size() == 8)
			{
				const double fac_cost = multiplier * (non_uniform_cl ? std::stod(row[4]) : fc);
				U.push_back(Point3D(std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), fac_cost, std::stoi(row[5]), std::stoi(row[6]), std::stoi(row[7])));

			}
			//Increment index
			index++;
		}
	}

	//Throw exception
	catch (std::ifstream::failure e)
	{
		std::cerr << "Can not open file> :" << file_name << '\n';
	}
}


void IO::savePointCloud(const std::string& file_name, const TVector <Point3D>& U)
{
	//Save point cloud to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all points
		for (const auto &point : U)
		{
			//Write coordinates
			file << std::setprecision(8);
			file << point.id << "  " << point.x << "  " << point.y << "  " << point.z << "  " << point.fc << "  ";

			//Write r, g, b components
			file << std::setprecision(0);
			file << point.r << "  " << point.g << "  " << point.b << '\n';
		}

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void IO::saveFacilitiesNoStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& ID)
{
	//Save pfacilities to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all points
		for (int i = 0; i < U.size(); i++)
		{
			//Write point ID
			file << std::setprecision(0);
			file << ID[i] << "  ";

			//Write coordinates
			file << std::setprecision(8);
			file << U[i].x << "  " << U[i].y << "  " << U[i].z << "  " << U[i].fc << "  ";

			//Write r, g, b components
			file << std::setprecision(0);
			file << U[i].r << "  " << U[i].g << "  " << U[i].b << " ";
		}

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void IO::saveFacilitiesAndStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& ID, const TVector <int>& NC, const TVector <double>& RAD, const TVector <double>& ABN, const TVector <double>& DFP, const TVector <double>& ASP, const TVector <int>& DIM, const TVector <int>& OVER, const TVector <double>& SLO)
{
	//Save facilities and statistics to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all points
		for (int i = 0; i < U.size(); i++)
		{
			//Write point ID
			file << std::setprecision(0);
			file << ID[i] << "  ";

			//Write coordinates
			file << std::setprecision(8);
			file << U[i].x << "  " << U[i].y << "  " << U[i].z << "  " << U[i].fc << "  ";

			//Write r, g, b components
			file << std::setprecision(0);
			file << U[i].r << "  " << U[i].g << "  " << U[i].b << " ";

			//Write statistics
			file << std::setprecision(0);
			file << NC[i] << "  ";

			file << std::setprecision(8);
			file << RAD[i] << "  " << ABN[i] << "  " << DFP[i] << "  " << ASP[i] << "  ";

			file << std::setprecision(0);
			file << DIM[i] << "  " << OVER[i] << "  ";

			//Write slope
			file << std::setprecision(6);
			file << SLO[i] << '\n';
		}

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}



void IO::saveFacilitesIDXs(const std::string& file_name, const TVector2D <int> &idxs_all)
{
	//Save facilities and assigned clients
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all facilities
		for (auto idxs_clust : idxs_all)
		{
			//Write facility
			file << idxs_clust[0] << ":\n";

			//Are clients available?
			if (idxs_clust.size() > 1)
			{
				//Write its clients
				for (int i = 1; i < idxs_clust.size(); i++)
					file << idxs_clust[i] << " ";

				file << '\n';
			}

			file << '\n';
		}

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void IO::saveClientsToFacilitesIDXs(const std::string& file_name, const std::map <int, int>& idxs_all)
{
	//Assign client to a given facility
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all items
		for (const auto idx : idxs_all)
			file << idx.second << "\n";

		file.close();	
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}