// Description: Incremental heuristic facility location clustering
// Uniform and non-uniform version

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


#ifndef IHFL_H
#define IHFL_H

#include <tuple>

#include "point3d.h"
#include "pfnorm.h"
#include "regressionplane.h"
#include "tvector2D.h"
#include "facility.h"
#include "gridindexing.h"


#include <iostream>
#include <fstream>
#include <iomanip>

//Incremental heuristic facility location algorithm
//Clustering according to the given (pseudo) norm
class IHFL
{
	private:
		bool non_uniform_cl;		//Non-uniform clusterization (recompute cost of points according to normal vectors)
		int k;				//Amount of nearest neighbors, estimation of normal vectors
		double  lambda,			//Radius of the ball
			bin,			//Size of the bin of the grid (performance of the spatial indexing)
			mju,			//Isotropic factor
			l,			//Penalization outside the ball, power of distance differences
			x_scale,		//Scale factor in X direction
			y_scale,		//Scale factor in Y direction
			z_scale;		//Scale factor in Z direction

		const pfnorm &fnorm;		//Reference to the (pseudo) norm


	public:
		 IHFL(const bool non_uniform_cl_, const int k_, const double lambda_, const double bin_, const double mju_, const double l_, const double x_scale_, 
			const double y_scale_, const double z_scale_, const pfnorm &fnorm_) : non_uniform_cl(non_uniform_cl_), k(k_), lambda(lambda_), bin(bin_), 
			 mju(mju_), l(l_), x_scale(x_scale_), y_scale(y_scale_), z_scale(z_scale_), fnorm(fnorm_) {}
		
		 void generateClusters(const double w, const double h, const double rad, const int nc, const int n, TVector <Point3D>& U);
		 void generateCone(const double a, const double b, const int n, TVector <Point3D>& U);
		 void generateCylinder(const double a, const double b, const int n, TVector <Point3D>& U);
		 void generateCube(const double a, const int n, TVector <Point3D>& U);
		 
		 void clusterizeIHFL(TVector <Point3D>& U, const double fc, const GridIndexing &gi, TVector2D <Facility> &FG, TVector <RegressionPlane> & RP);
		 void clusterStatistics(const TVector <Point3D>& points, const TVector2D <Facility> &FG, const GridIndexing& gi, const TVector <RegressionPlane>& RP, 
		     TVector <int>& NC, TVector <double>& RAD, TVector <double>& ABN, TVector <double>& DFP, TVector <double>& ASP, TVector <int>& DIM, 
		     TVector <int>& OVER, TVector <double>& SLO);

		 double nL2(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nL1(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nEll(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nDIS(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nABN(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nABLP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nDFP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double computeCost(const TVector <Facility>& F, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP);

private:
		 double distL2(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
		 double distL1(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
		 double distEll(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
		 double pointPlaneDist(const double qx, const double qy, const double qz, const double a, const double b, const double c, const double d);
		 double abn(const RegressionPlane & pa, const RegressionPlane & pb);
		 double ablp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double dfp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		
		 void sideConePoint(const double a, const double b, double& h, double& r, double& t);
		 void baseConePoint(const double a, const double b, double& h, double& r, double& t);
		 void sideCylinderPoint(const double a, const double b, double& h, double& r, double& t);
		 void baseCylinderPoint(const double a, const double b, double& h, double& r, double& t);

		 void updateClusters(const int i, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP, const GridIndexing & gi, TVector2D <Facility >& FG);
		 void recomputeFacilityCosts (const double fc, double rat, const TVector <RegressionPlane>& RP, const pfnorm& fnorm, TVector <Point3D>& U);
		 void getAveragePointNormal(TVector <Point3D>& U, const TVector2D <size_t>& knn_id, TVector <RegressionPlane>& RP);
};

#endif