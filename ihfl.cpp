// Description: Incremental heuristic facility location clustering

// Copyright (c) 2021 - 2025
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


#include "ihfl.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>

#include "knnsearch.h"
#include "regressionplane.h"
#include "pca.h"
#include "isfacilitysefordeletion.h"
#include "sortpoints3dbyz.h"

#include "Eigen/Dense"                               
#include "Eigen/Sparse"                              
#include "Eigen/Core"

double IHFL::distL2(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute L2 (Euclidean) distance
	const double dx = x2 - x1;
	const double dy = y2 - y1;
	const double dz = z2 - z1;

	return sqrt(dx * dx + dy * dy + dz * dz);
}


double IHFL::distL1(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute L1 distance
	const double dx = x2 - x1;
	const double dy = y2 - y1;
	const double dz = z2 - z1;

	return fabs(dx) + fabs(dy) + fabs(dz);
}


double IHFL::distL22(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute square of L2 (Euclidean) distance
	const double dx = x2 - x1;
	const double dy = y2 - y1;
	const double dz = z2 - z1;

	return dx * dx + dy * dy + dz * dz;
}


double IHFL::distEll(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute elliptic distance
	const double dx = x_scale * (x2 - x1);
	const double dy = y_scale * (y2 - y1);
	const double dz = z_scale * (z2 - z1);

	return sqrt(dx * dx + dy * dy + dz * dz);
}


double IHFL::pointPlaneDist(const double qx, const double qy, const double qz, const double a, const double b, const double c, const double d)
{
	//Distance of the point q = [qz, qy, qz] from the plane given by (a, b, c, d)
	return fabs(a * qx + b * qy + c * qz + d) / sqrt(a * a + b * b + c * c);
}


double IHFL::abn(const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Compute ABN (Angle Between Planes)
	const double norma = sqrt(pa.a * pa.a + pa.b * pa.b + pa.c * pa.c);
	const double normb = sqrt(pb.a * pb.a + pb.b * pb.b + pb.c * pb.c);
	const double dot_ab = fabs(pa.a * pb.a + pa.b * pb.b + pa.c * pb.c);

	return acos(std::max(-1.0, std::min(1.0, dot_ab / (norma * normb))));
}


double IHFL::ablp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Compute ABLP (Angle Between Line and RegressionPlane)
	const double ux = b.x - a.x;
	const double uy = b.y - a.y;
	const double uz = b.z - a.z;

	//Norms
	const double norma = sqrt(pa.a * pa.a + pa.b * pa.b + pa.c * pa.c);
	const double normb = sqrt(pb.a * pb.a + pb.b * pb.b + pb.c * pb.c);
	const double normu = sqrt(ux * ux + uy * uy + uz * uz);

	//Dot products
	const double dot_pau = fabs(ux * pa.a + uy * pa.b + uz * pa.c);
	const double dot_pbu = fabs(-ux * pb.a - uy * pb.b - uz * pb.c);

	//Compute angles
	const double ablp1 = asin(std::min(1.0, dot_pau / (norma * normu)));
	const double ablp2 = asin(std::min(1.0, dot_pbu / (normb * normu)));

	return ablp1 + ablp2;
}


double IHFL::dfp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Compute DFP (Distance from plane)
	const double dfa = pointPlaneDist(a.x, a.y, a.z, pb.a, pb.b, pb.c, pb.d);
	const double dfb = pointPlaneDist(b.x, b.y, b.z, pa.a, pa.b, pa.c, pa.d);

	//Norm
	return sqrt(dfa * dfa + dfb * dfb);
}


double IHFL::lin(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute linearity
	const double lina = (l1a - l2a) / l1a;
	const double linb = (l1b - l2b) / l1b;

	return fabs(lina - linb);
	//return 0.5 * (lina + linb);
}


double IHFL::pla(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute planarity
	const double plaa = (l2a - l3a) / l1a;
	const double plab = (l2b - l3b) / l1b;

	return 0.5 * (plaa + plab);
}


double IHFL::sph(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute sphericity
	const double spha = l3a / l1a;
	const double sphb = l3b / l1b;

	return fabs(spha - sphb);
}


double IHFL::omn(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute omnivariance
	const double omna = pow(l1a * l2a * l3a, 1 / 3.0);
	const double omnb = pow(l1b * l2b * l3b, 1 / 3.0);

	return fabs(omna - omnb);
}


double IHFL::ani(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute anisotropy
	const double ania = (l1a - l3a) / l1a;
	const double anib = (l1b - l3b) / l1b;

	return fabs(ania - anib);
}


double IHFL::cur(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Compute curvature change
	const double curva = l3a / (l1a + l2a + l3a);
	const double curvb = l3b / (l1b + l2b + l3b);

	return fabs(curva - curvb);
}


double IHFL::ent(const double l1a, const double l2a, const double l3a, const double l1b, const double l2b, const double l3b)
{
	//Eigen entropy
	const double enta = - (l1a * log(l1a) + l2a * log(l2a) + l3a * log(l3a));
	const double entb = - (l1b * log(l1b) + l2b * log(l2b) + l3b * log(l3b));

	return fabs(enta - entb);
}


double IHFL::ver(const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Verticality, measured by the slope
	const double norm1 = sqrt(pa.a * pa.a + pa.b * pa.b + pa.c * pa.c);
	const double f_slope1 = acos(pa.c / norm1);

	const double norm2 = sqrt(pb.a * pb.a + pb.b * pb.b + pb.c * pb.c);
	const double f_slope2 = acos(pb.c / norm2);

	return 0.5 * fabs(log(f_slope1) + log(f_slope2));
}


double IHFL::hor(const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Horizontality, measured by the slope
	const double norm1 = sqrt(pa.a * pa.a + pa.b * pa.b + pa.c * pa.c);
	const double f_slope1 = acos(pa.c / norm1);

	const double norm2 = sqrt(pb.a * pb.a + pb.b * pb.b + pb.c * pb.c);
	const double f_slope2 = acos(pb.c / norm2);

	return 0.5 * (log(f_slope1) + log(f_slope2));
}


double IHFL::mL2(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//L2 metric
	return distL2(a.x, a.y, a.z, b.x, b.y, b.z);
}


double IHFL::mL1(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//L1 metric
	return distL1(a.x, a.y, a.z, b.x, b.y, b.z);
}


double IHFL::mL22(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//L22 metric
	return distL22(a.x, a.y, a.z, b.x, b.y, b.z);
}


double IHFL::mEll(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Elliptic distance
	return distEll(a.x, a.y, a.z,  b.x, b.y, a.z);
}


double IHFL::pmGeo(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Geographic pseudometric by prof Sykora
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Constrained pseudometric
	if (dab < lambda)
		return 0.0;
	else
		return dab - lambda;
}


double IHFL::pmABN(const Point3D& a, const Point3D& b, const RegressionPlane & pa, const RegressionPlane & pb)
{
	//Pseudometric: Angle between normals (ABN), tangent model, g1
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute abn criterion
	const double dabn = abn(pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dabn + (1 - mju) * dab);
	else
		return (mju * dabn + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmDIS(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Euclidean distance between p and tangent plane, tangent model (g2)
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute abn criterion
	const double dabn = abn(pa, pb);
	const double dh = sin(0.5 * dabn) * dab;

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmABLP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Angle between line and planes, secant model (g3)
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute ablp criterion
	const double dablp = ablp(a, b, pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dablp  + (1 - mju) * dab);
	else
		return (mju * dablp + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmDFP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Distance from plane (DFP), secant model (g4)
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute dfp criterion
	const double dh = dfp(a, b, pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmLIN(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Linearity
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute linearity criterion
	const double dh = lin(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmPLA(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Planarity
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute panarity criterion
	const double dh = pla(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmSPH(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Sphericity
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute sphericity criterion
	const double dh = sph(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmOMN(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Omnivariance
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute omnivariance criterion
	const double dh = omn(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmANI(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Anisotropy
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute anisotropy criterion
	const double dh = ani(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmCUR(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Curvature change
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute curvature criterion
	const double dh = cur(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmENT(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Eigen entropy
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute eigen entropy
	const double dh = ent(pa.lambda1, pa.lambda2, pa.lambda3, pb.lambda1, pb.lambda2, pb.lambda3);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmVER(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: verticality
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute verticality
	const double dh = ver(pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


double IHFL::pmHOR(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: horizontality
	const double dab = mL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute horizontality
	const double dh = hor(pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dh + (1 - mju) * dab);
	else
		return (mju * dh + (1 - mju) * dab) + mju * pow(dab - lambda, l);
}


void IHFL::generateClusters(const double w, const double h, const double rad, const int nc, const int n, TVector <Point3D>& U)
{
	//Generate nc clusters, each is represented by n points
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu1(0.0, 1.0);

	//Generate points
	for (int i = 0; i < nc; i++)
	{
		//Generate random center of the cluster
		const double xc = disu1(gen) * w - rad;
		const double yc = disu1(gen) * h - rad;

		//Generate points of the cluster
		std::uniform_real_distribution<double> disu2(xc, xc + rad);
		std::uniform_real_distribution<double> disu3(yc, rad);

		for (int j = 0; j < n; j++)
			U.push_back(Point3D(0, disu2(gen), disu3(gen), 0));
	}
}


void IHFL::sideConePoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the side of cone
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a * sqrt(disu(gen));
	r = (b / a) * h;
	t = 2.0 * std::numbers::pi * disu(gen);
}


void IHFL::baseConePoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the base of cone
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a;
	r = b * sqrt(disu(gen));
	t = 2.0 * std::numbers::pi * disu(gen);
}


void IHFL::sideCylinderPoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the side of cylinder
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a * disu(gen);
	r = b;
	t = 2.0 * std::numbers::pi * disu(gen);
}


void IHFL::baseCylinderPoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the base of cylinder
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> disi(0, 1);
	std::uniform_real_distribution<double> disu(0, 1);

	h = disi(gen) * a;
	r = b * sqrt(disu(gen));
	t = 2.0 * std::numbers::pi * disu(gen);
}


void IHFL::generateCone(const double a, const double b, const int n, TVector <Point3D>& U)
{
	//Generate points on the cone surface
	double t = 0;

	//Initialize random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		double h, r, t;

		//Point on the base
		if (disu(gen) < b / (b + sqrt(a * a + b * b)))
			baseConePoint(a, b, h, r, t);

		//Point on the side
		else
			sideConePoint(a, b, h, r, t);

		//Cone coordipates
		const double x = r * cos(t);
		const double y = h;
		const double z = r * sin(t);

		//Apend to the list
		U.push_back(Point3D(0, x, z, -y));
	}
}


void IHFL::generateCylinder(const double a, const double b, const int n, TVector <Point3D>& U)
{
	//Generate cylinder
	double t = 0;

	//Initialize random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		double h, r, t;

		//Point on the base
		if (disu(gen) < b / (2 * (a + b)))
			baseCylinderPoint(a, b, h, r, t);

		//Point on the side
		else
			sideCylinderPoint(a, b, h, r, t);

		//Cone coordipates
		const double x = r * cos(t);
		const double y = h;
		const double z = r * sin(t);

		//Apend to the list
		U.push_back(Point3D(0, x, z, -y));
	}
}


void IHFL::generateCube(const double a, const int n, TVector <Point3D>& U)
{
	//Generate points on the surface of cube
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_int_distribution <int> disi(0, 5);
	std::uniform_real_distribution<double> disu(0, 1);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		//Select random side of the cube
		int side = disi(gen);

		//Create empty row
		double xyz[] = { 0, 0, 0 };

		//Coordipate index
		int c = side % 3;

		//Create a point
		if (side > 2)
			xyz[c] = a;
		else
			xyz[c] = 0;

		xyz[(c + 1) % 3] = a * disu(gen);
		xyz[(c + 2) % 3] = a * disu(gen);

		//Apend to the list
		U.push_back(Point3D(0, xyz[0], xyz[1], xyz[2]));
	}
}


void IHFL::recomputeFacilityCostsByPCA(const TVector <RegressionPlane>& RP, TVector <Point3D> &U)
{
	//Recompute facility cost according to the normals for non-uniform clustering given by PCA
	const double eps = 0.0001, max_fc = 10000;
	
	//Divide by the height of the PCA bounding box 
	for (int i = 0; i < U.size(); i++)
	{
		//Compute ratio
		const double ratio = U[i].fc / (RP[i].height + eps);

		//Use square ratio
		U[i].fc = std::min(ratio * ratio, max_fc);
	}
}


void IHFL::getAveragePointNormal(const TVector <Point3D>& P, const TVector2D <size_t>& knn_idxs, TVector <RegressionPlane>& RP)
{
	//Compute average point normal using SVD decomposition, regression error and ABN
	const int n = P.size();

	std::cout << ">> SVD decomposition: ";
	for (int i = 0; i < n; i++)
	{
		//Convert nearest neighbors to matrix A
		const int m = knn_idxs[i].size();
		Eigen::MatrixXd A(m, 3);

		for (int j = 0; j < m; j++)
		{
			const int knn_idx = knn_idxs[i][j];
			A(j, 0) = P[knn_idx].x;
			A(j, 1) = P[knn_idx].y;
			A(j, 2) = P[knn_idx].z;
		}

		//Compute PCA
		auto [U, S, V] = PCA::pca(A);

		//Compute regression plane using SVD
		RegressionPlane plane;
		plane.computeRegressionPlane(A, S, V);

		//Add average normal vector and regression error to the list
		RP.push_back(plane);
	}

	std::cout << "OK \n";
}


void IHFL::clusterizeIHFL(TVector <Point3D>&U,const GridIndexing& gi, TVector2D <Facility> &FG, TVector <RegressionPlane> &RP)
{
	//Perform incremental heuristic facility clustering by (IHFL) according to a given metric
	//Proposed incremental algorithms selecting one of four strategies
	const int n = U.size();

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_idxs;
	TVector2D <float> knn_dists;
	TVector <int> K(1, k);
	KNNSearch search (U);
	search.findAllKNN(U, K, knn_idxs, knn_dists);

	//Compute average normals, recompute facility cost
	getAveragePointNormal(U, knn_idxs, RP);
	
	//Recompute facility costs according to normals (replace old values)]
	const double multiplier = 10.0;

	if (non_uniform_cl && recompute_cost)
	{
		recomputeFacilityCostsByPCA( RP, U);
	}

	//Process all points
	std::cout << ">> Clusterization: ";

	//First point becomes the facility
	Facility f0(1, U[0].fc);

	//Compute its grid index (position in the 3D grid)
	int idx_p = gi.get1DIndex(U[0]);

	//Add facility to the grid index structure
	FG[idx_p].push_back(f0);

	//Proces remaining points one by one
	for (int i = 1; i < n; i++)
	{
		//Print status
		if (i % (25000) == 0)
		{
			std::cout << i << " ";
		}

		//Update clusters using the IHFL method, heuristic approach
		updateClusters(i, U, RP, gi, FG);
	}
}


void IHFL::updateClusters(const int i, const TVector <Point3D>& points, const TVector <RegressionPlane>& RP, const GridIndexing& gi, TVector2D <Facility >& FG)
{
	//Incremental method, update of the clusterization takes one of four strategies:
	//   1) Create new facility at p: S1
	//   2) Connect p to the (pseudo) nearest facility: S2
	//   3) Reallocate all clusters (pseudo) nearest to p. Create facility at p + multiple full reallocations: S3
	//   4) Reallocate parts of clusters (pseudo) nearest to p. Create facility at p + multiple partial reallocations: S4

	//Initialize index of the nearest cluster
	std::pair<int, int> idx_fac_nearest = { -1, -1 };

	//Initial costs for different strategies
	double c_nearest = std::numeric_limits<float>::max();									//Cost of the nearest facility, strategy S1
	double c_new = points[i].fc;												//Cost for creation of the new facility, strategy S2
	double c_reallocate_facilities = points[i].fc;										//Cost for the cluster modification (full or partial reassignment to another facility), strategy S3 a) b)

	//Get grid index of the added point
	int idx_p = gi.get1DIndex(points[i]);

	//Get indices of all facilities in the adjacent cells of the grid
	TVector <int> idxs_fac = gi.getAdjacentIndices(points[i]);

	//Process all facilities in the adjacent cells of the grid
	TVector <std::pair <int, int> > reallocate_facilities;
	for (int &idx_fac : idxs_fac)
	{
		//Process all facilities in the same bin
		for (int j = 0; j < FG[idx_fac].size(); j++)
		{
			//Get facility
			Facility & fac = FG[idx_fac][j];

			//Reset sign and shift
			const int p_idx2 = abs(fac.p_idx) - 1;

			//Distance and pseudodistance between point p[i] and facility
			const double dist_pf = distL2(points[i].x, points[i].y, points[i].z, points[p_idx2].x, points[p_idx2].y, points[p_idx2].z);    
			const double p_dist_pf = (this->*pf_clust_metric)(points[i], points[p_idx2], RP[i], RP[p_idx2]);							

			//const double dist_pf = distL2( points[p_idx2].x, points[p_idx2].y, points[p_idx2].z, points[i].x, points[i].y, points[i].z);
			//const double p_dist_pf = (this->*pf_clust_metric)(points[p_idx2], points[i], RP[p_idx2], RP[i]);			

			//Cummulated values
			double dc_all = points[i].fc - fac.fc;			//Cost difference: reassignment of all cluster to  p - cost for the old facility deletition
			double dc_closer = points[i].fc - p_dist_pf;		//Cost difference: reassignment of cluster points closer to p

			//Reallocate only according to near facilities (closer than 2 * lambda)
			if (dist_pf <= 1.0 * lambda) //Was 3.0 * lambda
			{
				//Distance of a point and the current facility: Strategy S1
				if (p_dist_pf < c_nearest && p_dist_pf > 0)	//Actualize distance to the nearest facility representing cost
				{
					//Actualize cost
					c_nearest = p_dist_pf;
					
					//Store facility ID
					idx_fac_nearest.first = idx_fac;
					idx_fac_nearest.second = j;
				}

				//Browse all points ui of the facility
				for (int& u_idx : fac.U_idxs)
				{
					//Reset sign and shift of the index
					const int u_idx2 = abs(u_idx) - 1;

					//const double p_dist_up = (this->*pf_clust_metric)(points[u_idx2], points[i], RP[u_idx2], RP[i]);		//Pseudodistance of the cluster point u to the proposed new center p[i]
					//const double p_dist_uf = (this->*pf_clust_metric)(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);	//Pseudodistance of the cluster point u to its cluster center f.u

					//Pseudodistance of the cluster point u to the proposed new center p[i]
					const double p_dist_up = (this->*pf_clust_metric)(points[i], points[u_idx2], RP[i], RP[u_idx2]);
					
					//Pseudodistance of the cluster point u to its facility (cluster) center f.u
					const double p_dist_uf = (this->*pf_clust_metric)(points[p_idx2], points[u_idx2], RP[p_idx2], RP[u_idx2]);	//Pseudodistance of the cluster point u to its cluster center f.u

					//Current facility point u is closer to p: cost for the reassignment u to the new center p
					//Strategy S4, compute new cost increment
					if (p_dist_up < p_dist_uf)
					{
						//Cost difference
						dc_closer += p_dist_up - p_dist_uf;

						//Change sign to plus (reallocate to p)
						u_idx = abs(u_idx2) + 1;
					}

					//Change sign to minus: (remain connected to old facility)
					else
						u_idx = -abs(u_idx2) - 1;

					//Cost for the reassignment of any facility point to p
					//Strategy S3
					dc_all += p_dist_up - p_dist_uf;
				}

				//Entire facility is reallocated to p: does the newly created cluster at p decrease the cost?
				//Strategy S3
				if (dc_all < dc_closer)
				{
					//Cost difference
					c_reallocate_facilities = c_reallocate_facilities + dc_all - points[i].fc + p_dist_pf;			//Expected cost increment

					//Change facility ID to plus
					fac.p_idx = abs(fac.p_idx);

					//Set ID of the reallocated cluster
					reallocate_facilities.push_back(std::pair<int, int>(idx_fac, j));
				}

				//Reallocate part of the facility to p: is there any improvement?
				//Strategy S4
				else if (dc_closer < dc_all)
				{
					//Cost difference
					c_reallocate_facilities = c_reallocate_facilities + dc_closer - points[i].fc + p_dist_pf;			//Expected cost increment

					//Change facility ID to minus
					fac.p_idx = -abs(fac.p_idx);

					//Set ID of the reallocated facility
					reallocate_facilities.push_back(std::pair<int, int>(idx_fac, j));
				}
			}
		}
	}

	//Heuristic decision: find minimum cost increment and choose the optimal strategy (S1-S3 a) b))
	const double c_min = std::min(c_new, std::min(c_nearest, c_reallocate_facilities));

	//Strategy S1: Create new facility at p
	if (c_min == c_new)
	{
		Facility fac_new(i + 1, points[i].fc);
		FG[idx_p].push_back(fac_new);
	}

	//Strategy S2: Connect p to the nearest facility
	else if (c_min == c_nearest)
	{
		FG[idx_fac_nearest.first][idx_fac_nearest.second].U_idxs.push_back(i + 1);
	}

	//Strategy S3 or S4
	else
	{
		//Create new facility fp at p
		Facility fac_new(i + 1, points[i].fc);

		//Process all facilities affected by the reallocation one by one
		for (const auto& k : reallocate_facilities)
		{
			//Partial realloation: Strategy S4
			if (FG[k.first][k.second].p_idx < 0)
			{
				//Reconnect all clients
				TVector <int> U_old_idxs;
				for (const int& u_idx : FG[k.first][k.second].U_idxs)
				{
					//Reconnect client to the new facility
					if (u_idx > 0)
						fac_new.U_idxs.push_back(abs(u_idx));

					//Client remains to the old facility
					else
						U_old_idxs.push_back(abs(u_idx));
				}

				//Connect remaining points to the old facility (faster then delete)
				FG[k.first][k.second].U_idxs = std::move(U_old_idxs);
			}

			//Full reallocation: Strategy S3
			else
			{
				//Connect all clients of the old facility to the new facility fp at p
				fac_new.U_idxs.insert(fac_new.U_idxs.end(), std::make_move_iterator(FG[k.first][k.second].U_idxs.begin()), std::make_move_iterator(FG[k.first][k.second].U_idxs.end()));

				//Connect old facility to the new facility fp at p
				fac_new.U_idxs.push_back(abs(FG[k.first][k.second].p_idx));

				//Mark old facility for the deletetion
				FG[k.first][k.second].del = true;
			}
		}

		//Add new facility at p to the list
		FG[idx_p].push_back(fac_new);
	}

	//Delete marked old facilities (now they are fully reallocated to p, strategy S3)
	for (const int & idx_fac : idxs_fac)
		FG[idx_fac].erase(std::remove_if(FG[idx_fac].begin(), FG[idx_fac].end(), IsFacilitySetForDeletion()), FG[idx_fac].end());
}


double IHFL::computeCost(const TVector <Facility>& F, const TVector <Point3D> &points, const TVector <RegressionPlane>& AN)
{
	//Compute cost of the clusterization
	double total_cost = 0;

	for (const auto &f : F)
	{
		//Add the facility cost
		total_cost += f.fc;

		//Process all connected points
		const int p_idx2 = abs(f.p_idx) - 1;
		for (const auto c_id : f.U_idxs)
		{
			const int c_id2 = fabs(c_id) - 1;
			total_cost += (this->*pf_clust_metric)(points[p_idx2], points[c_id2], AN[p_idx2], AN[c_id2]);
		}
	}

	return total_cost;
}


void IHFL::clusterStatistics(const TVector <Point3D>& points, const TVector2D <Facility> &FG, const GridIndexing& gi, const TVector <RegressionPlane>& RP, TVector <int>& NC, TVector <double> &RAD, TVector <double>& ABN, TVector <double>& DFP, TVector <double>& ASP, TVector <int>& DIM, TVector <int>& OVER, TVector <double>& SLO)
{
	//Compute parameters of the cluster
	for (const auto &fg : FG)
	{ 
		//Process all facilities inside the grid bin
		for (const Facility& f : fg)
		{
			//Browse all points ui of the facility
			int i = 0;
			double f_radius = 0, f_abn = 0, f_dfp = 0, f_aspect = 1, f_dim = 0;

			//Reset sign and shift of the index
			const int p_idx2 = abs(f.p_idx) - 1;

			//Process all clients
			Eigen::MatrixXd M(f.U_idxs.size(), 3);
			for (const int u_idx : f.U_idxs)
			{
				//Reset sign and shift of the index
				const int u_idx2 = abs(u_idx) - 1;

				//Convert clients to matrix
				M(i, 0) = points[u_idx2].x;
				M(i, 1) = points[u_idx2].y;
				M(i, 2) = points[u_idx2].z;

				//Compute pseudometrics
				f_radius = std::max(mL2(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]), f_radius);
				f_abn += pmABN(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);
				f_dfp += pmDFP(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);

				i++;
			}

			//PCA Analysis
			if (f.U_idxs.size() > 0)
			{
				//Perform PCA
				auto [U, S, V] = PCA::pca(M);

				//Store singular values: descendent sort
				double lambda1 = S(0, 0);
				double lambda2 = S(1, 0);
				double lambda3 = S(2, 0);

				//Compute cluster aspect
				f_aspect = (fabs(lambda1) < 1.0e-6 ? 1 : lambda2 / lambda1);

				//Sum of singular values
				const double lambda_sum = lambda1 + lambda2 + lambda3;

				//Cluster dimension 0
				if (fabs(lambda_sum) < 1.0e-6)
					f_dim = 0;

				//Cluster dimension: 0 - 3
				else
				{
					//Fractions of singular values
					const double lambda1f = lambda1 / lambda_sum;
					const double lambda2f = lambda2 / lambda_sum;
					const double lambda3f = lambda3 / lambda_sum;

					//Compute sum
					const double lambda_sum = lambda1f + lambda2f + lambda3f;

					//Set cluster dimensions
					//Point
					if (lambda1f < 0.01)	
						f_dim = 0;

					//Line	
					else if (lambda2f / lambda_sum < 0.01)	
						f_dim = 1;

					//Circle
					else if (lambda3f / lambda_sum < 0.01)
						f_dim = 2;

					//Sphere
					else				
						f_dim = 3;
				}
			}

			//Compute overlap ratio
			double f_over = 0;

			//Get indices of all facilities in the adjacent cells of the grid
			TVector <int> idxs_fac = gi.getAdjacentIndices(points[p_idx2]);

			//Process all facilities in the adjacent cells of the grid
			for (const int & idx_fac : idxs_fac)
			{
				//Process all facilities in the same bin
				for (const Facility &f2 : FG[idx_fac])
				{
					//Different facilities
					if (f.p_idx != f2.p_idx)
					{
						//Take all connected points
						for (int k = 0; k < f2.U_idxs.size(); k++)
						{
							//Reset sign and shift of the index
							const int u_idx2 = abs(f2.U_idxs[k]) - 1;

							//Measure distance between the facility and a client
							const double d = mL2(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);

							//Point closer to centroid than radius, increment overlap ratio
							if (d < f_radius)
								f_over += 1;
						}
					}
				}
			}

			//Compute slope
			const double norm = sqrt(RP[p_idx2].a * RP[p_idx2].a + RP[p_idx2].b * RP[p_idx2].b + RP[p_idx2].c * RP[p_idx2].c);
			const double f_slope = acos(RP[p_idx2].c / norm) * 180 / std::numbers::pi;
			
			//Store statistical values
			NC.push_back(f.U_idxs.size());
			RAD.push_back(f_radius);
			ABN.push_back(f_abn / std::max(1, (int)f.U_idxs.size()));
			DFP.push_back(f_dfp / std::max(1, (int)f.U_idxs.size()));
			ASP.push_back(f_aspect);
			DIM.push_back(f_dim);
			SLO.push_back(f_slope);
			OVER.push_back(f_over);
		}
	}
}


TVector <int> IHFL::outputFaciliesToIDXs(const TVector <Point3D>& points, const TVector2D <Facility>& FG)
{
	//Convert output faciities to indices
	TVector <int> ID;

	for (const auto& fg : FG)
	{
		//Process all facilities inside the grid bin
		for (const Facility& f : fg)
		{
			//Reset sign and shift of the index
			const int p_idx2 = abs(f.p_idx) - 1;

			//Add to the list of point indices within the point cloud
			ID.push_back(points[p_idx2].id);
		}
	}

	return ID;
}


TVector2D <int> IHFL::facilitiesToIDXs(const TVector <Point3D>& points, const TVector <Facility>& F)
{
	//Convert clusters to lists
	TVector2D <int> cloud_idxs;

	for (auto f : F)
	{
		//Convert to the list
		TVector <int> cluster_idxs = f.U_idxs;

		//Add to the first position
		cluster_idxs.emplace(cluster_idxs.begin(), f.p_idx);

		//Process all points
		for (int i = 0; i < cluster_idxs.size(); i++)
		{
			//Convert to point tile index
			const int u_idx2 = abs(cluster_idxs[i]) - 1;

			//Get point cloud tile index
			const int ucl_id = points[u_idx2].id;

			//Add to the list
			cluster_idxs[i] = ucl_id;
		}

		//Add to the list
		cloud_idxs.push_back(cluster_idxs);
	}

	return cloud_idxs;
}


std::map <int, int> IHFL::clientsToFacilitiesIDXs(const TVector <Point3D>& kd_points_tile, const TVector <Facility>& F)
{
	//Assign clients to facilities within a point cloud
	std::map <int, int> clients_to_facilities_cloud;

	for (int i = 0; i < F.size(); i++)
	{
		//Get point id
		const int p_idx2 = abs(F[i].p_idx) - 1;

		//Asign facility to itself
		clients_to_facilities_cloud[kd_points_tile[p_idx2].id] = kd_points_tile[p_idx2].id;

		//Browse all clients connected to the facility
		for (int u_idx : F[i].U_idxs)
		{
			//Get client ID
			const int u_idx2 = abs(u_idx) - 1;

			//Asign client to the facility
			clients_to_facilities_cloud[kd_points_tile[u_idx2].id] = kd_points_tile[p_idx2].id;
		}
	}

	return clients_to_facilities_cloud;
}
