#ifndef IO_H
#define IO_H

#include <string>
#include "point3d.h"
#include "tvector.h"

//Input - output functions
class IO
{
	public:
		static void loadPointCloud(const std::string& file_name, TVector <Point3D>& U, const double fc = 1.0); 
		static void savePointCloud(const std::string& file_name, const TVector <Point3D>& U);
		static void saveClientsToFacilites(const std::string& file_name, const TVector <int>& CL);
		static void savePointCloudAndStatistics(const std::string& file_name, const TVector <Point3D>& U, const TVector <int>& NC, const TVector <double>& RAD, const TVector <double>& ABN, const TVector <double>& DFP, const TVector <double>& ASP, const TVector <int>& DIM, const TVector <int>& OVER);

};

#endif
