// Description: Export output data into DXF fike

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


#include "dxfexport.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "point3d.h"

void DXFExport::exportClustersToDXF(const std::string& file_name, const TVector <Facility>& F, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP)
{
	//Export clusters represented by facilities to DXF file
	const unsigned int color_points = 5, color_facilities = 3, color_facilities_lines = 3, color_normals = 1, color_text = 3;
	const double font_height = 0.001;
	const std::string level_points = "Source_points", level_facilities = "Facilities", level_facilities_lines = "Facilities_lines", level_normals = "Normals", level_text = "Text";

	std::ofstream file;

	try
	{
		file.open(file_name);

		//Create header section
		createHeaderSection(file);

		//Create table section
		createTableSection(file);

		//Create layer for source points
		createLayerSection(file, level_points, color_points);

		//Create layer for facilities
		createLayerSection(file, level_facilities, color_facilities);

		//Create layer for lines conecting facilities
		createLayerSection(file, level_facilities_lines, color_facilities_lines);

		//Create layer for normals
		createLayerSection(file, level_normals, color_normals);

		//Create layer for text
		createLayerSection(file, level_text, color_text);

		//End table header
		endTableSection(file);

		//Create entity section
		createEntitySection(file);

		//Process all facilities
		for (Facility f : F)
		{
			//Get point id
			const int p_idx2 = abs(f.p_idx) - 1;

			//Create facility
			createPoint(file, level_facilities, points[p_idx2].x, points[p_idx2].y, points[p_idx2].z, color_facilities);

			//Create facility label
			std::string facility_label = std::to_string(p_idx2);
			createText(file, level_text, facility_label, points[p_idx2].x + 0.5 * font_height, points[p_idx2].y - 0.5 * font_height, points[p_idx2].z, 0.0, font_height, color_text);

			//Create connected points and connecting lines
			for (auto u_idx : f.U_idxs)
			{
				//Get client ID
				const int u_idx2 = abs(u_idx) - 1;

				//Create connected point
				createPoint(file, level_points, points[u_idx2].x, points[u_idx2].y, points[u_idx2].z, color_points);

				//Create connecting line
				createLine(file, level_facilities_lines, points[p_idx2].x, points[p_idx2].y, points[p_idx2].z, points[u_idx2].x, points[u_idx2].y, points[u_idx2].z, color_facilities_lines);

				//if (points[p_idx2].x == points[u_idx2].x && points[p_idx2].y == points[u_idx2].y)
				//	std::cout << "identical";

				//Create point label
				std::string point_label = std::to_string(u_idx2);
				createText(file, level_text, point_label, points[u_idx2].x + 0.5 * font_height, points[u_idx2].y - 0.5 * font_height, points[u_idx2].z, 0.0, font_height, color_text);

				//Create normal at a point
				createLine(file, level_normals, points[u_idx2].x, points[u_idx2].y, points[u_idx2].z, points[u_idx2].x + 0.1 * RP[u_idx2].a, points[u_idx2].y + 0.1 * RP[u_idx2].b, points[u_idx2].z + 0.1 * RP[u_idx2].c, color_normals);
			}

			//Create normal at facility
			createLine(file, level_normals, points[p_idx2].x, points[p_idx2].y, points[p_idx2].z, points[p_idx2].x + 0.1 * RP[p_idx2].a, points[p_idx2].y + 0.1 * RP[p_idx2].b, points[p_idx2].z + 0.1 * RP[p_idx2].c, color_normals);
		}

		//End header section
		endHeaderSection(file);

		//Close file
		file.close();
	}

	//Any error has appeared
	catch (std::ifstream::failure e)
	{
		std::cerr << "Exception writing DXF file\n";
	}
}


void DXFExport::createHeaderSection(std::ofstream& file)
{
	//Create header section
	const std::string object_type_id = "0\n";
	const std::string object_type = "SECTION\n";
	const std::string header_id = "2\n";
	const std::string header_name = "HEADER\n";
	const std::string variable_name = "9\n";
	const std::string acad_version = "$ACADVER\n";
	const std::string acad_version_id = "1\n";
	const std::string acad_version_val = "AC1006\n";
	const std::string section_terminator = "ENDSEC\n";

	/* Add to file */
	file << object_type_id;
	file << object_type;
	file << header_id;
	file << header_name;
	file << variable_name;
	file << acad_version;
	file << acad_version_id;
	file << acad_version_val;
	file << object_type_id;
	file << section_terminator;
}


void DXFExport::endHeaderSection(std::ofstream& file)
{
	//Create end of the header section
	const std::string entity_id = "0\n";
	const std::string section_terminator = "ENDSEC\n";
	const std::string file_terminator = "EOF\n";

	/* Add to file */
	file << entity_id;
	file << section_terminator;
	file << entity_id;
	file << file_terminator;
}


void DXFExport::createTableSection(std::ofstream& file)
{
	//Create table section
	const std::string object_type_id = "0\n";
	const std::string object_name = "SECTION\n";
	const std::string table_id = "2\n";
	const std::string table = "TABLES\n";
	const std::string table_name = "TABLE\n";
	const std::string layer_name = "LAYER\n";
	const std::string max_number_of_entries_id = "70\n";
	const std::string max_number_of_entries = "0\n";

	file << object_type_id;
	file << object_name;
	file << table_id;
	file << table;
	file << object_type_id;
	file << table_name;
	file << table_id;
	file << layer_name;
	file << max_number_of_entries_id;
	file << max_number_of_entries;
}


void DXFExport::endTableSection(std::ofstream& file)
{
	//Write end of the table section
	const std::string object_type_id = "0\n";
	const std::string table_end_name = "ENDTAB\n";
	const std::string section_end_name = "ENDSEC\n";

	/* Add to file */
	file << object_type_id;
	file << table_end_name;
	file << object_type_id;
	file << section_end_name;
}


void DXFExport::createLayerSection(std::ofstream& file, const std::string& layer_name, const unsigned int color)
{
	//Add section for one layer
	const std::string object_type_id = "0\n";
	const std::string object_name = "LAYER\n";
	const std::string layer_name_id = "2\n";
	const std::string layer_flag_id = "70\n";
	const std::string layer_flag = "0\n";
	const std::string color_number_id = "62\n";
	const std::string line_type_id = "6\n";
	const std::string line_type_name = "CONTINUOUS\n";

	/* Add to file */
	file << object_type_id;
	file << object_name;
	file << layer_name_id;
	file << layer_name << '\n';
	file << layer_flag_id;
	file << layer_flag;
	file << color_number_id;
	file << color << '\n';
	file << line_type_id;
	file << line_type_name;
}


void DXFExport::createEntitySection(std::ofstream& file)
{
	//Create section for entities
	const std::string object_type_id = "0\n";
	const std::string object_name = "SECTION\n";
	const std::string entity_name_id = "2\n";
	const std::string entity_name = "ENTITIES\n";

	/* Add to file */
	file << object_type_id;
	file << object_name;
	file << entity_name_id;
	file << entity_name;
}


void DXFExport::createLine(std::ofstream& file, const std::string& layer_name, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const int color)
{
	//Write line to DXF file
	const std::string entity_id = "0\n";
	const std::string entity_name = "LINE\n";
	const std::string level_id = "8\n";
	const std::string color_id = "62\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";
	const std::string xii_id = "11\n";
	const std::string yii_id = "21\n";
	const std::string zii_id = "31\n";

	//Set format and precision
	file << std::fixed << std::setprecision(3);

	/* Add to file */
	file << entity_id;
	file << entity_name;
	file << level_id;
	file << layer_name << '\n';
	file << color_id;
	file << color << '\n';
	file << xi_id;
	file << x1 << '\n';
	file << yi_id;
	file << y1 << '\n';
	file << zi_id;
	file << z1 << '\n';
	file << xii_id;
	file << x2 << '\n';
	file << yii_id;
	file << y2 << '\n';
	file << zii_id;
	file << z2 << '\n';
}


void DXFExport::createPoint(std::ofstream& file, const std::string& layer_name, const double x, const double y, const double z, const int color)
{
	//Write point to DXF file
	const std::string entity_id = "0\n";
	const std::string entity_name = "POINT\n";
	const std::string level_id = "8\n";
	const std::string color_id = "62\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";

	//Set format and precision
	file << std::fixed << std::setprecision(3);

	/* Add to file */
	file << entity_id;
	file << entity_name;
	file << level_id;
	file << layer_name << '\n';
	file << color_id;
	file << color << '\n';
	file << xi_id;
	file << x << '\n';
	file << yi_id;
	file << y << '\n';
	file << zi_id;
	file << z << '\n';
}


void DXFExport::createText(std::ofstream& file, const std::string& layer_name, const std::string& text, const double x, const double y, const double z, const double rotation, const double height, const unsigned int color)
{
	//Create text
	const std::string entity_id = "0\n";
	const std::string entity_name = "TEXT\n";
	const std::string style_id = "7\n";
	const std::string text_style = "PNTNUM\n";
	const std::string rotation_id = "50\n";
	const std::string level_id = "8\n";
	const std::string color_id = "62\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";
	const std::string height_id = "40\n";
	const std::string text_id = "1\n";

	//Set format and precision
	file << std::fixed << std::setprecision(3);

	/* Add to file */
	file << entity_id;
	file << entity_name;
	file << style_id;
	file << text_style;
	file << rotation_id;
	file << rotation << '\n';
	file << level_id;
	file << layer_name << '\n';
	file << color_id;
	file << color << '\n';
	file << xi_id;
	file << x << '\n';
	file << yi_id;
	file << y << '\n';
	file << zi_id;
	file << z << '\n';
	file << height_id;
	file << height << '\n';
	file << text_id;
	file << text << '\n';
}


std::string DXFExport::to_string(const double value, const unsigned short dec_places)
{
	//Convert to string,  set the decimal places
	std::ostringstream out;
	out << std::setprecision(dec_places) << value;
	return out.str();
}