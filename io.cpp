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

void IO::loadPointCloud(const std::string& file_name, TVector <Point3D>& U, const double fc, const double mul)
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
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fc * mul));

			//Add point P = [x, y, z, fc] to the list
			else if (row.size() == 4)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]) * mul));

			//Add point P = [x, y, z, r, g, b] to the list
			else if (row.size() == 6)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fc * mul, std::stod(row[3]) * fc, std::stoi(row[4]), std::stoi(row[5])));

			//Add point P = [x, y, z, fc, r, g, b] to the list
			else if (row.size() == 7)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), std::stoi(row[4]), std::stoi(row[5]), std::stoi(row[6])));

			//Unknown structure
			//Add point P = [x, y, z, r, g, b] to the list
			else if (row.size() > 7)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fc, std::stoi(row[3]), std::stoi(row[4]), std::stoi(row[5])));

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
			file << point.x << "  " << point.y << "  " << point.z << "  " << point.fc << "  ";

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


void IO::savePointCloudAndStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& NC, const TVector <double>& RAD, const TVector <double>& ABN, const TVector <double>& DFP, const TVector <double>& ASP, const TVector <int>& DIM, const TVector <int>& OVER, const TVector <double>& SLO)
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
		for (int i = 0; i < U.size(); i++)
		{
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


void IO::saveClientsToFacilites(const std::string& file_name, const TVector <int>& CL)
{
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all items
		for (auto id : CL)
			file << id << "\n";

		file.close();
			
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}