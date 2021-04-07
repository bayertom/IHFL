// Description: Export lines, points, polygons to 2D/3D DXF file

// Copyright (c) 2015 - 2016
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



#ifndef DXFExport_HPP
#define DXFExport_HPP

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>


#include "libalgo/source/const2/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/exceptions/FileWriteException.h"

template <typename T>
void DXFExport::exportDTToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &nl, const TVector <std::shared_ptr<HalfEdge <T> > > &dt, const T font_height)
{
	//Export DT to DXF file
	const unsigned int color_dt = 1, color_points = 5, color_points_labels = 5;
	const std::string level_dt = "dt", level_points = "dt_points", level_points_labels = "dt_points_labels";

	std::ofstream file;
	//std::setprecision(5) << std::fixed;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for dt
			createLayerSection(file, level_dt, color_dt);

			//Create layer for dt points
			createLayerSection(file, level_points, color_points);

			//Create layer for dt point labels
			createLayerSection(file, level_points_labels, color_points_labels);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all halfedges
			processHalfEdges(file, dt, level_dt, color_dt);
			
			//Process all test points
			/*
			for (const std::shared_ptr <Node3DCartesian <T> > node : nl)
			{
				//Create test point
				createPoint(file, level_points, node->getX(), node->getY(), node->getZ(), color_points);

				//Create test point label
				//std::string point_label = std::to_string(node->getPointID());
				createText(file, level_points_labels, point_label, node->getX() + 0.5 * font_height, node->getY() - 0.5 * font_height, node->getZ(), 0.0, font_height, color_points_labels);
			}
			*/
			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportContourLinesToDXF(const std::string &file_name, const TVector <std::shared_ptr <Edge <Node3DCartesian <T> > > > &contours, const T font_height)
{
	//Export contour lines to DXF file
	const unsigned int color_cont = 1, color_cont_points = 5;
	const std::string level_cont = "contour_lines", level_cont_points = "contour_lines_points", level_points_labels = "dt_points_labels";

	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for contour lines
			createLayerSection(file, level_cont, color_cont);

			//Create layer for dt contour lines labels
			createLayerSection(file, level_cont_points, color_cont_points);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all edges
			processEdges(file, contours, level_cont, color_cont);
			/*
			//Process all test points
			for (const std::shared_ptr <Node3DCartesian <T> > node : nl)
			{
			//Create test point
			createPoint(file, level_points, node->getX(), node->getY(), node->getZ(), color_points);

			//Create test point label
			std::string point_label = std::to_string(node->getPointID());
			createText(file, level_points_labels, point_label, node->getX() + 0.5 * font_height, node->getY() - 0.5 * font_height, node->getZ(), 0.0, font_height, color_points_labels);
			}
			*/
			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}

template <typename T>
void DXFExport::exportContourLinesToDXF(const std::string &file_name, TVector2D <std::shared_ptr <Node3DCartesian <T> > > contours_polylines, const T font_height)
{
	//Export contour lines given by polylines to DXF file
	const unsigned int color_cont = 1, color_cont_points = 5;
	const std::string level_cont = "contour_lines", level_cont_points = "contour_lines_points", level_points_labels = "dt_points_labels";

	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for contour lines
			createLayerSection(file, level_cont, color_cont);

			//Create layer for dt contour lines labels
			createLayerSection(file, level_cont_points, color_cont_points);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process contour lines
			for (TVector <std::shared_ptr<Node3DCartesian <T> > > polyline : contours_polylines)
			{
				processPolyline(file, polyline, level_cont, color_cont);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}



template <typename T>
void DXFExport::exportGraticuleToDXF(const std::string &file_name, const TVector <Meridian <T> > &meridians, const TVector2D <Point3DCartesian<T> > & meridians_proj, const TVector <Parallel <T> > &parallels, const TVector2D <Point3DCartesian<T> > & parallels_proj, const TVector <Point3DCartesian<T> > & test_points, const TVector <Point3DCartesian<T> > & reference_points_proj, const T font_height, const T step)
{
	//Export the graticule formed by meridians/parallels into DXF file
	//Export test points, reference points, and projected points
	const unsigned int color_graticule = 8, color_test_points = 5, color_reference_points = 1;

	const std::string  level_meridians = "Meridians", level_parallels = "Parallels", level_meridian_labels = "Meridian_labels",
		level_parallel_labels = "Parallel_labels", level_test_points = "Test_points", level_test_point_labels = "Test_point_labels",
		level_proj_reference_points = "Reference_points_proj", level_proj_reference_point_labels = "Reference_points_proj_labels";
	
	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for meridians
			createLayerSection(file, level_meridians, color_graticule);

			//Create layer for parallels
			createLayerSection(file, level_parallels, color_graticule);

			//Create layer for meridian labels
			createLayerSection(file, level_meridian_labels, color_graticule);

			//Create layer for parallel labels
			createLayerSection(file, level_parallel_labels, color_graticule);

			//Create layer for test points
			createLayerSection(file, level_test_points, color_test_points);

			//Create layer for test points
			createLayerSection(file, level_test_point_labels, color_test_points);

			//Create layer for projected reference points
			createLayerSection(file, level_proj_reference_points, color_reference_points);

			//Create layer for projected reference points
			createLayerSection(file, level_proj_reference_point_labels, color_reference_points);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all meridians
			const unsigned int nm = meridians_proj.size();
			for (unsigned int i = 0; i < nm; i++)
			{
				processGraticuleElements(file, meridians_proj[i], meridians[i].getLon(), level_meridians, level_meridian_labels, font_height, step, color_graticule);
			}

			//Process all parallels
			const unsigned int np = parallels_proj.size();
			for (unsigned int i = 0; i < np; i++)
			{
				processGraticuleElements(file, parallels_proj[i], parallels[i].getLat(), level_parallels, level_parallel_labels, font_height, step, color_graticule);
			}
			
			//Process all test points
			unsigned int index_test = 0;
			for (const Point3DCartesian <T> point : test_points)
			{
				//Create test point
				createPoint(file, level_test_points, point.getX(), point.getY(), point.getZ(), color_test_points);

				//Create test point label
				std::string point_label = std::to_string(++index_test);
				createText(file, level_test_point_labels, point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0.0, font_height, color_test_points);
			}

			//Process all projected reference points
			unsigned int index_reference = 0;
			for (const Point3DCartesian <T> point : reference_points_proj)
			{
				//Create reference point
				createPoint(file, level_proj_reference_points, point.getX(), point.getY(), point.getZ(), color_reference_points);

				//Create reference point label
				std::string point_label = std::to_string( ++index_reference);
				createText(file, level_proj_reference_point_labels, point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0.0, font_height, color_reference_points);
			}
			
			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportStraightWalkToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &nl_dt, const TVector <std::shared_ptr<HalfEdge <T> > > &dt, const TVector <std::shared_ptr<Node3DCartesian <T> > > &walk, const TVector <std::shared_ptr<Node3DCartesian <T> > > &se)
{
	//Export straight walk, start/end points and the triangulation
	//Export also DT to DXF file
	const unsigned int color_dt = 1, color_walk = 5, color_se = 5;
	const std::string level_dt = "dt", level_walk = "Straigh_walk", level_se = "endpoints";

	std::ofstream file;
	//std::setprecision(5) << std::fixed;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for dt
			createLayerSection(file, level_dt, color_dt);

			//Create layer for walk
			createLayerSection(file, level_walk, color_walk);

			//Create layer for walk end points
			createLayerSection(file, level_se, color_se);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all halfedges
			processHalfEdges(file, dt, level_dt, color_dt);

			//Process the walk
			processPolyline(file, walk, level_walk, color_walk);

			//Process its endpoints
			for (const std::shared_ptr <Node3DCartesian <T> > node : se)
			{
				//Create test point
				createPoint(file, level_se, node->getX(), node->getY(), node->getZ(), color_se);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportStraightWalkToDXF2(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &walk, const TVector <std::shared_ptr<Node3DCartesian <T> > > &se)
{
	//Export straight walk, start/end points and the triangulation
	//Export DT to DXF file
	const unsigned int color_walk = 3, color_se = 5;
	const std::string level_walk = "Straigh_walk", level_se = "endpoints";

	std::ofstream file;
	//std::setprecision(5) << std::fixed;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for walk
			createLayerSection(file, level_walk, color_walk);

			//Create layer for walk end points
			createLayerSection(file, level_se, color_se);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process the walk
			processPolyline(file, walk, level_walk, color_walk);

			//Process its endpoints
			createLine(file, level_se, se[0]->getX(), se[0]->getY(), se[0]->getZ(), se[1]->getX(), se[1]->getY(), se[1]->getZ(), color_se);

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportClustersToDXF(const std::string& file_name, const TVector <Facility>& F, const TVector2D <double>& ANN)
{
	//Export clusters represented by facilities to DXF file
	const unsigned int color_points = 5, color_facilities = 3, color_facilities_lines = 3, color_normals = 1, color_text = 3;
	const double font_height = 0.2;
	const std::string level_points = "Source_points", level_facilities = "Facilities", level_facilities_lines = "Facilities_lines", level_normals = "Normals", level_text = "Text";

	std::ofstream file;

	try
	{
		file.open(file_name);

		if (file.is_open())
		{
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
			int i = 0;
			for (Facility f : F)
			{
				//Create facility
				createPoint(file, level_facilities, f.u.x, f.u.y, f.u.z, color_facilities);

				//Create facility label
				std::string facility_label = std::to_string(f.u.id);
				createText(file, level_text, facility_label, f.u.x + 0.5 * font_height, f.u.y - 0.5 * font_height, f.u.z, 0.0, font_height, color_text);

				//Create connected points and connecting lines

				for (const Point3D p : f.C)
				{
					//Create connected point
					createPoint(file, level_points, p.x, p.y, p.z, color_points);

					//Create connecting line
					createLine(file, level_facilities_lines, f.u.x, f.u.y, f.u.z, p.x, p.y, p.z, color_facilities_lines);

					//Create point label
					std::string point_label = std::to_string(p.id);
					createText(file, level_text, point_label, p.x + 0.5 * font_height, p.y - 0.5 * font_height, p.z, 0.0, font_height, color_text);

					//Create normal at a point
					createLine(file, level_normals, p.x, p.y, p.z, p.x + 0.1 * ANN[p.id][0], p.y + 0.1 * ANN[p.id][1], p.z + 0.1 * ANN[p.id][2], color_normals);
				}

				//Create normal at facility
				createLine(file, level_normals, f.u.x, f.u.y, f.u.z, f.u.x + 0.1 * ANN[f.u.id][0], f.u.y + 0.1 * ANN[f.u.id][1], f.u.z + 0.1 * ANN[f.u.id][2], color_normals);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure&)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}



template <typename T>
void DXFExport::exportPolylinesToDXF(const std::string& file_name, TVector2D <Point3DCartesian <T> > polylines, const T font_height)
{
	//Export contour lines given by polylines to DXF file
	const unsigned int color_poly = 1, color_poly_points = 5;
	const std::string level_poly = "polylines", level_poly_points = "polylines_points";

	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for contour lines
			createLayerSection(file, level_poly, color_poly);

			//Create layer for dt contour lines labels
			createLayerSection(file, level_poly_points, color_poly_points);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process contour lines
			for (TVector <Point3DCartesian <T> > polyline : polylines)
			{
				processPolyline(file, polyline, level_poly, color_poly);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure&)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportPolygonToDXF(const std::string& file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > >& polygon, const unsigned int color)
{
	//Export polygon to DXF
	TVector <std::shared_ptr<Node3DCartesian <T> > > polygon_temp = polygon;

	//Add first point
	polygon_temp.push_back(polygon_temp[0]);

	//Export
	exportPolylineToDXF(file_name, polygon_temp, color);
}

template <typename T>
void DXFExport::exportPolylineToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &polyline, const unsigned int color)
{
	//Export polyline to DXF
	const std::string  level_pol = "Polyline";

	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer
			createLayerSection(file, level_pol, color);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process contour lines
			processPolyline(file, polyline, level_pol, color);
			
			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename T>
void DXFExport::exportPolylineToDXF(const std::string& file_name, const TVector <Point3DCartesian <T> > & polyline, const unsigned int color)
{
	//Export polyline to DXF
	const std::string  level_pol = "Polyline";

	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer
			createLayerSection(file, level_pol, color);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process contour lines
			processPolyline(file, polyline, level_pol, color);

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure&)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename Point, typename T>
void DXFExport::exportLinesToDXF(const std::string& file_name, const TVector <Point> &start, const TVector <Point>& end, const T font_height, const unsigned int color)
{
	//Export list of points to DXF file
	std::ofstream file;

	//Export lines to DXF
	try
	{
		file.open(file_name);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for points
			createLayerSection(file, "Lines", color);

			//Create layer for point labels
			//createLayerSection(file, "Point_labels", color);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all points
			for (int i = 0; i < start.size(); i++)
			{
				//Create point
				createLine(file, "Lines", start[i].x, start[i].y, start[i].z, end[i].x, end[i].y, end[i].z, color);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure&)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}

}



template <typename Point, typename T>
void DXFExport::exportPointsToDXF(const std::string &file_name, const TVector <Point> &points, const T font_height, const unsigned int color)
{
	//Export list of points to DXF file
	std::ofstream file;

	try
	{
		file.open(file_name);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for points
			createLayerSection(file, "Points", color);

			//Create layer for point labels
			//createLayerSection(file, "Point_labels", color);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all points
			for (const Point point:points)
			{
				//Create point
				createPoint(file, "Points", point.x, point.y, point.z, color);
				//createPoint(file, "Points", point.getX(), point.getY(), point.getZ(), color);
				
				/*
				//Create point label
				char point_id_text[255];
				sprintf(point_id_text, "%d", point.id);
				//sprintf(point_id_text, "%d", point.getPointID());
				std::string point_label = point_id_text;
				createText(file, "Point_labels", point_label, point.x + 0.5 * font_height, point.y - 0.5 * font_height, point.z, 0, font_height, color);
				//createText(file, "Point_labels", point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0, font_height, color);
				*/
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


inline void DXFExport::createHeaderSection ( std::ofstream & file )
{
        //Create header section
	const std::string object_type_id = "0\n";
	const std::string object_type = "SECTION\n";
	const std::string header_id = "2\n";
	const std::string header_name = "HEADER\n";
	const std::string variable_name = "9\n";
	const std::string acad_version = "$ACADVER\n";
	const std::string acad_version_id =  "1\n";
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


inline void DXFExport::endHeaderSection (std::ofstream & file)
{
        //Create end of the header section
	const std::string entity_id = "0\n";
	const std::string section_terminator ="ENDSEC\n";
	const std::string file_terminator = "EOF\n";
	
	/* Add to file */
	file << entity_id;
        file << section_terminator;
        file << entity_id;
        file << file_terminator;
}


inline void DXFExport::createTableSection (std::ofstream & file )
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


inline void DXFExport::endTableSection ( std::ofstream & file )
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


inline void DXFExport::createLayerSection (std::ofstream & file, const std::string &layer_name, const unsigned int color )
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


inline void DXFExport::createEntitySection ( std::ofstream & file )
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


template <typename T>
void DXFExport::createLine (std::ofstream & file, const std::string &layer_name, const T x1, const T y1, const T z1, const T x2, const T y2, const T z2, const int color)
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


template <typename T>
void DXFExport::createPoint (std::ofstream & file, const std::string & layer_name, const T x, const T y, const T z, const int color )
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


template <typename T>
void DXFExport::createText (std::ofstream & file, const std::string & layer_name, const std::string & text, const T x, const T y, const T z, const T rotation, const T height, const unsigned int color )
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


template <typename T>
void DXFExport::processHalfEdges(std::ofstream & file, const TVector <std::shared_ptr<HalfEdge <T> > > &hl, const std::string &layer_name, const unsigned int color)
{
	//Process all half edges
	const unsigned int n = hl.size();

	//Process halfedges
	for (unsigned int i = 0; i < n; i++)
	{
		// Get edge
		std::shared_ptr <HalfEdge <T> > e = hl[i];

		//Use simplex indentificator to eliminate T processing of the edge
		if (!e->isSimplexEdge())
		{
			// Set twin edge as simplex, i. e. processed
			if (e->getTwinEdge()) e->getTwinEdge()->setSimplexEdge(true);

			// Get start point
			const T x1 = e->getPoint()->getX();
			const T y1 = e->getPoint()->getY();
			const T z1 = e->getPoint()->getZ();

			// Get end point
			const T x2 = e->getNextEdge()->getPoint()->getX();
			const T y2 = e->getNextEdge()->getPoint()->getY();
			const T z2 = e->getNextEdge()->getPoint()->getZ();

			//Create line
			createLine(file, layer_name, x1, y1, z1, x2, y2, z2, color);

			//std::string point_id_text = std::to_string(i);
			//createText(file, layer_name, point_id_text, 0.5 *(x1 + x2), 0.5 * (y1 + y2), 0.5 * (z1 + z2), 0.0, 0.5, 0);
		}
	}

	// Set all edges as non simplex
	for (unsigned int i = 0; i < n; i++)
	{
		hl[i]->setSimplexEdge(false);
	}
}


template <typename T>
void DXFExport::processEdges(std::ofstream & file, const TVector <std::shared_ptr <Edge <Node3DCartesian<T> > > > &hl, const std::string &layer_name, const unsigned int color)
{
	//Process all edges
	const unsigned int n = hl.size();

	//Process halfedges
	for (unsigned int i = 0; i < n; i++)
	{
		// Get edge
		std::shared_ptr <Edge <Node3DCartesian<T> > > e = hl[i];

		// Get start point
		const T x1 = e->getP1()->getX();
		const T y1 = e->getP1()->getY();
		const T z1 = e->getP1()->getZ();

		// Get end point
		const T x2 = e->getP2()->getX();
		const T y2 = e->getP2()->getY();
		const T z2 = e->getP2()->getZ();

		//Create line
		createLine(file, layer_name, x1, y1, z1, x2, y2, z2, color);
	}
}


template <typename T>
void DXFExport::processPolyline(std::ofstream & file, TVector <std::shared_ptr<Node3DCartesian <T> > > polyline, const std::string &layer_name, const unsigned int color)
{
	const unsigned int n = polyline.size();

	//Process halfedges
	for (unsigned int i = 0; i < n - 1; i++)
	{
		// Get start point
		const T x1 = polyline[i]->getX();
		const T y1 = polyline[i]->getY();
		const T z1 = polyline[i]->getZ();

		// Get end point
		const T x2 = polyline[i+1]->getX();
		const T y2 = polyline[i+1]->getY();
		const T z2 = polyline[i+1]->getZ();

		//Create line
		createLine(file, layer_name, x1, y1, z1, x2, y2, z2, color);
	}
}


template <typename T>
void DXFExport::processPolyline(std::ofstream& file, TVector <Point3DCartesian <T> > polyline, const std::string& layer_name, const unsigned int color)
{
	const unsigned int n = polyline.size();

	//Process halfedges
	for (unsigned int i = 0; i < n - 1; i++)
	{
		// Get start point
		const T x1 = polyline[i].getX();
		const T y1 = polyline[i].getY();
		const T z1 = polyline[i].getZ();

		// Get end point
		const T x2 = polyline[i + 1].getX();
		const T y2 = polyline[i + 1].getY();
		const T z2 = polyline[i + 1].getZ();

		//Create line
		createLine(file, layer_name, x1, y1, z1, x2, y2, z2, color);
	}
}



template <typename GraticulePart, typename T>
void DXFExport::processGraticuleElements( std::ofstream & file, GraticulePart & part, const T val, const std::string &layer_graticule_name, const std::string &layer_labels_name, const T font_height, const T step, const unsigned int color)
{
	//Export meridian or parallel (defined as a template parameter GraticulePart)

	//Process all points
	const unsigned int n = part.size();

	for (unsigned int i = 1; i < n; i++)
	{
		//Get actual point
		Point3DCartesian<T> p = &part[i];

		//Get previous point
		Point3DCartesian<T>  p_previous = &part[i - 1];

		//Create a line> part of the meridian or parallel
		createLine(file, layer_graticule_name, p_previous.getX(), p_previous.getY(), p_previous.getZ(), p.getX(), p.getY(), p.getZ(), color);
	}

	//Create label
	char point_id_char[255];

	if (n > 1)
	{
		//Set accuracy depending on a step
		if (step > 1.0) sprintf(point_id_char, "%3.1f", val);
		else if (step > 0.1) sprintf(point_id_char, "%3.2f", val);
		else sprintf(point_id_char, "%3.3f", val);

		std::string point_id_text = point_id_char;

		//Compute bearing
		Point3DCartesian<T>  p1 = &part[0.5 * n - 1], p2 =  &part[0.5 * n];
		const T bearing = atan2(p2.getY() - p1.getY(), p2.getX() - p1.getX()) * RO;

		//Create label for meridian/parallel
		createText(file, layer_labels_name, point_id_text, part[0.5 * n - 1].getX() + 0.5 * font_height * cos(bearing * M_PI / 180), 
			part[0.5 * n - 1].getY() + 0.5 * font_height * sin(bearing * M_PI / 180), part[0].getZ(), bearing, font_height, color);
	}
}


template <typename T>
std::string DXFExport::to_string(const T value, const unsigned short dec_places)
{
	//Convert to string,  set the decimal places
	std::ostringstream out;
	out << std::setprecision(dec_places) << value;
	return out.str();
}

#endif
