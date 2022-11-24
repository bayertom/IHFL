#include "io.h"

#include <iostream>
#include <fstream>
#include <iomanip>

//#include "kdpointtile.h"
#include "tvector.h"

void IO::loadPointCloud(const std::string& file_name, TVector <Point3D>& U, const double fc)
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
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fc));

			//Add point P = [x, y, z, c] to the list
			else if (row.size() == 4)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3])));

			//Add point P = [x, y, z, fc, r, g, b] to the list
			else if (row.size() == 6)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), fc, std::stod(row[3]), std::stoi(row[4]), std::stoi(row[5])));

			//Add point P = [x, y, z, fc, r, g, b] to the list
			else if (row.size() == 7)
				U.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]), std::stoi(row[4]), std::stoi(row[5]), std::stoi(row[6])));


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