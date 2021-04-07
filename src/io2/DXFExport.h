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


#ifndef DXFExport_H
#define DXFExport_H

#include <memory>

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TVector2D.h"

#include "libalgo/source/structures/edge/Edge.h"
#include "libalgo/source/structures/edge/HalfEdge.h"
#include "libalgo/source/structures/graticule/Meridian.h"
#include "libalgo/source/structures/graticule/Parallel.h"

struct Facility;
struct Point3D;

//Export to DXF file
class DXFExport
{
        public:
              
		template <typename T>
		static void exportDTToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &nl, const TVector <std::shared_ptr<HalfEdge <T> > > &dt, const T font_height);

		template <typename T>
		static void exportContourLinesToDXF(const std::string &file_name, const TVector <std::shared_ptr< Edge <Node3DCartesian <T> > > > &contours, const T font_height);

		template <typename T>
		static void exportContourLinesToDXF(const std::string &file_name, TVector2D <std::shared_ptr <Node3DCartesian <T> > > contours_polylines, const T font_height);

		template <typename T>
		static void exportGraticuleToDXF(const std::string &file_name, const TVector <Meridian <T> > &meridians, const TVector2D <Point3DCartesian<T> > & meridians_proj, const TVector <Parallel <T> > &parallels, const TVector2D <Point3DCartesian<T> > & parallels_proj, const TVector <Point3DCartesian<T> > & test_points, const TVector <Point3DCartesian<T> > & reference_points_proj, const T font_height, const T step);

		template <typename T>
		static void exportStraightWalkToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &nl_dt, const TVector <std::shared_ptr<HalfEdge <T> > > &dt, const TVector <std::shared_ptr<Node3DCartesian <T> > > &walk, const TVector <std::shared_ptr<Node3DCartesian <T> > > &se);

		template <typename T>
		static void exportStraightWalkToDXF2(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &walk, const TVector <std::shared_ptr<Node3DCartesian <T> > > &se);

		template <typename T>
		static void exportClustersToDXF(const std::string& file_name, const TVector <Facility> &F, const TVector2D <double>& ANN);

		template <typename T>
		static void exportPolylinesToDXF(const std::string& file_name, TVector2D <Point3DCartesian <T> > polylines, const T font_height);

		template <typename T>
		static void exportPolygonToDXF(const std::string& file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > >& polylgon, const unsigned int color);

		template <typename T>
		static void exportPolylineToDXF(const std::string &file_name, const TVector <std::shared_ptr<Node3DCartesian <T> > > &polyline, const unsigned int color);

		template <typename T>
		static void exportPolylineToDXF(const std::string& file_name, const TVector <Point3DCartesian <T>>& polyline, const unsigned int color);

		template <typename Point, typename T>
		static void exportPointsToDXF(const std::string &file_name, const TVector <Point> &points, const T font_height, const unsigned int color);

		template <typename Point, typename T>
		static void exportLinesToDXF(const std::string& file_name, const TVector <Point>& start, const TVector <Point>& end, const T font_height, const unsigned int color);


        private:

                static void createHeaderSection (std::ofstream & file);
                static void createTableSection (std::ofstream & file);
                static void endTableSection (std::ofstream & file);
                static void createLayerSection (std::ofstream & file, const std::string &layer_name, const unsigned int color );
                static void createEntitySection (std::ofstream & file);
		static void endHeaderSection(std::ofstream & file);

                template <typename T>
                static void createLine (std::ofstream & file, const std::string &layer_name, const T x1, const T y1, const T z1, const T x2, const T y2, const T z2, const int color);

                template <typename T>
                static void createPoint (std::ofstream & file, const std::string &layer_name, const T x, const T y, const T z, const int color );

                template <typename T>
                static void createText (std::ofstream & file, const std::string &layer_name, const std::string & text, const T x, const T y, const T z, const T rotation, const T height, const unsigned int color );

		template <typename T>
		static void processHalfEdges(std::ofstream & file, const TVector <std::shared_ptr<HalfEdge <T> > > &hl, const std::string &layer_name,  const unsigned int color);

		template <typename T>
		static void processEdges(std::ofstream & file, const TVector <std::shared_ptr <Edge <Node3DCartesian <T> > > > &hl, const std::string &layer_name, const unsigned int color);
		
		template <typename T>
		static void processPolyline(std::ofstream & file, TVector <std::shared_ptr<Node3DCartesian <T> > > polyline, const std::string &layer_name, const unsigned int color);

		template <typename T>
		static void processPolyline(std::ofstream& file, TVector <Point3DCartesian <T> > polyline, const std::string& layer_name, const unsigned int color);


		template <typename GraticulePart, typename T>
		static void processGraticuleElements(std::ofstream & file, GraticulePart & part, const T val, const std::string &layer_graticule_name, const std::string &layer_labels_name, const T font_height, const T step, const unsigned int color);

		template <typename T>
		static std::string to_string(const T value, const unsigned short dec_places);

};

#include "DXFExport.hpp"

#endif
