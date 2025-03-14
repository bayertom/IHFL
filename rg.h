// Description: Region growing algorithm
// Uniform and non-uniform version

// Copyright (c) 2024 - 2025
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


#ifndef RG_H
#define RG_H

#include <queue>

#include "pfregionmetric.h"

#include "point3d.h"
#include "tvector2D.h"
#include "region.h"

class RG
{
	public:
		static void regionGrow( TVector <Point3D>& points, const pfRegionMetric& pf_edge_metric, const TVector<double>& edge_thresholds, const int k, TVector2D <Point3D>& regions_points);
		static void regionGrowDyn(TVector <Point3D>& points, const pfRegionMetric& pf_edge_metric, const TVector<double>& edge_thresholds, const int k, TVector2D <Point3D>& regions_points);
		static bool inSameRegion(const TVector <double>& p1_vals, const TVector <double>& p2_vals, const TVector <double>& edge_thresholds);
	
private:
		static TVector2D <double> computeFeatures(const TVector <Point3D>& points, const TVector2D <size_t>& knn_idxs);
		static TVector <bool> findEdgePointsCurv(const TVector2D <double>& features, const double threshold);
		static double get2VectorsAngle(const double ux, const double uy, const double uz, const double vx, const double vy, const double vz);
		static Region mergeRegions(Region& r1, const Region& r2);
		static std::tuple<double, double, double> measureSimilarity(const double curv1, const double tx1, const double ty1, const double tz1, const double nx1, const double ny1, const double nz1, const double curv2, const double tx2, const double ty2, const double tz2, const double nx2, const double ny2, const double nz2);
		static double getHausdorffDist(const Region& r1, const Region& r2, const TVector <Point3D>& points);
		static double hausdorffDist(const TVector <int>& idxs1, const TVector <int>& idxs2, const TVector <Point3D>& points);

};

#endif
