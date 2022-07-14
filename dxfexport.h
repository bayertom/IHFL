#ifndef DXFExport_H
#define DXFExport_H

#include <memory>
#include <string>

#include "tvector.h"
#include "facility.h"
#include "regressionplane.h"


//Export to DXF file
class DXFExport
{
public:

        static void exportClustersToDXF(const std::string& file_name, const TVector <Facility>& F, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP);

private:

        static void createHeaderSection(std::ofstream& file);
        static void createTableSection(std::ofstream& file);
        static void endTableSection(std::ofstream& file);
        static void createLayerSection(std::ofstream& file, const std::string& layer_name, const unsigned int color);
        static void createEntitySection(std::ofstream& file);
        static void endHeaderSection(std::ofstream& file);

        static void createLine(std::ofstream& file, const std::string& layer_name, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const int color);
        static void createPoint(std::ofstream& file, const std::string& layer_name, const double x, const double y, const double z, const int color);
        static void createText(std::ofstream& file, const std::string& layer_name, const std::string& text, const double x, const double y, const double z, const double rotation, const double height, const unsigned int color);
        static std::string to_string(const double value, const unsigned short dec_places);
};

#endif
