#ifndef IHFL2_H
#define IHFL2_H

#include <tuple>
#include <memory>

#include "point3d.h"
#include "pfnorm.h"
#include "regressionplane.h"
#include "tvector2D.h"
#include "facility.h"
#include "gridindexing.h"

//Incremental heuristic facility location algorithm
//Clustering according to the given (pseudo) norm
class IHFL2
{
	private:
		bool non_uniform_cl;		//Non-uniform clusterization (recompute cost of points according to normal vectors)
		int k;				//Amount of nearest neighbors, estimation of normal vectors
		double  lambda,			//Radius of the ball
			mju,			//Isotropic factor
			l;			//Penalization outside the ball, power of distance differences
		const pfnorm &fnorm;		//Reference to the (pseudo) norm


	public:
		 IHFL2(const bool non_uniform_cl_, const int k_, const double lambda_, const double mju_, const double l_, const pfnorm &fnorm_) :
			non_uniform_cl(non_uniform_cl_), k(k_), lambda(lambda_), mju(mju_), l(l_), fnorm(fnorm_) {}
		
		 void generateClusters(const double w, const double h, const double rad, const int nc, const int n, TVector <Point3D>& U);
		 void generateCone(const double a, const double b, const int n, TVector <Point3D>& U);
		 void generateCylinder(const double a, const double b, const int n, TVector <Point3D>& U);
		 void generateCube(const double a, const int n, TVector <Point3D>& U);
		 void clusterizeIHFL(TVector <Point3D>& U, const double fc, TVector<Facility> &F, TVector <RegressionPlane>& RP);
		 void clusterStatistics(const TVector <Point3D>& points, const TVector <Facility>& F, const TVector <RegressionPlane>& RP, 
		      TVector <int>& NC, TVector <double>& RAD, TVector <double>& ABN, TVector <double>& DFP, TVector <double>& ASP, 
		      TVector <int>& DIM, TVector <int>& OVER, TVector <double>& SLO);

		 double nL2(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nDIS(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nABN(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nABLP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double nDFP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double computeCost(const TVector <Facility>& F, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP);

private:
		 double dist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
		 double pointPlaneDist(const double qx, const double qy, const double qz, const double a, const double b, const double c, const double d);
		 double abn(const RegressionPlane & pa, const RegressionPlane & pb);
		 double ablp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		 double dfp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb);
		
		 void sideConePoint(const double a, const double b, double& h, double& r, double& t);
		 void baseConePoint(const double a, const double b, double& h, double& r, double& t);
		 void sideCylinderPoint(const double a, const double b, double& h, double& r, double& t);
		 void baseCylinderPoint(const double a, const double b, double& h, double& r, double& t);

		 void updateClusters(const int pi, const TVector <Point3D>& points, TVector <Facility>& F, const TVector<RegressionPlane>& RP);
		 void updateClusters2(const int i, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP, GridIndexing& gi, TVector2D <std::shared_ptr <Facility> >& F);
		 void updateClusters3(const int i, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP, GridIndexing& gi, TVector2D <Facility>& F);
		 void recomputeFacilityCosts (const double fc, double rat, const TVector <RegressionPlane>& RP, const pfnorm& fnorm, TVector <Point3D>& U);
		 void getAveragePointNormal(TVector <Point3D>& U, const TVector2D <size_t>& knn_id, TVector <RegressionPlane>& RP);
};

#endif