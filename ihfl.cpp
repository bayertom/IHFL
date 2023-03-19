#include "ihfl.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>

#include "knnsearch.h"
#include "regressionplane.h"
#include "isfacilitysefordeletion.h"

#include <Eigen/Dense>                               
#include <Eigen/Sparse>                              
#include <Eigen/Core>

double IHFL::dist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute Euclidean distance
	const double dx = x2 - x1;
	const double dy = y2 - y1;
	const double dz = z2 - z1;

	return sqrt(dx * dx + dy * dy + dz * dz);
}


double IHFL::pointPlaneDist(const double qx, const double qy, const double qz, const double a, const double b, const double c, const double d)
{
	//Distance of the point q = [qz, qy, qz] from the plane given by (a, b, c, d)
	return fabs(a * qx + b * qy + c * qz + d) / sqrt(a * a + b * b + c * c);
}


double IHFL::abn(const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Compute ABN criterion (Angle Between Planes)
	const double norma = sqrt(pa.a * pa.a + pa.b * pa.b + pa.c * pa.c);
	const double normb = sqrt(pb.a * pb.a + pb.b * pb.b + pb.c * pb.c);
	const double dot_ab = fabs(pa.a * pb.a + pa.b * pb.b + pa.c * pb.c);

	return acos(std::min(1.0, dot_ab / (norma * normb)));
}


double IHFL::ablp(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Compute ABLP criterion (Angle Between Line and RegressionPlane)
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
	//Compute DFP criterion (Distance from plane)
	const double dfa = pointPlaneDist(a.x, a.y, a.z, pb.a, pb.b, pb.c, pb.d);
	const double dfb = pointPlaneDist(b.x, b.y, b.z, pa.a, pa.b, pa.c, pa.d);

	//Norm
	return sqrt(dfa * dfa + dfb * dfb);
}


double IHFL::nL2(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//L2 norm
	return dist(a.x, a.y, a.z, b.x, b.y, b.z);
}


double IHFL::nABN(const Point3D& a, const Point3D& b, const RegressionPlane & pa, const RegressionPlane & pb)
{
	//Pseudometric: Angle between normals (ABN), tangent model, g1
	const double dab = nL2(a, b, pa, pb);

	//Identical points
	if (dab < std::numeric_limits<float>::min())
		return 0;

	//Compute abn criterion
	const double dabn = abn(pa, pb);

	//Constrained hybrid pseudometric
	if (dab < lambda)
		return (mju * dabn + (1 - mju) * dab);
	else
		return (mju * dabn + (1 - mju) * dab) + mju * (dab - lambda, l);
}


double IHFL::nDIS(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Euclidean distance between p and tangent plane, tangent model (g2)
	const double dab = nL2(a, b, pa, pb);

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


double IHFL::nABLP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Angle between line and planes, secant model (g3)
	const double dab = nL2(a, b, pa, pb);

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


double IHFL::nDFP(const Point3D& a, const Point3D& b, const RegressionPlane& pa, const RegressionPlane& pb)
{
	//Pseudometric: Distance from plane (DFP), secant model (g4)
	const double dab = nL2(a, b, pa, pb);

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
	//Points on the base of cone
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
	//Generate cone
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


void IHFL::recomputeFacilityCosts(const double fc, double rat, const TVector <RegressionPlane>& RP, const pfnorm& fnorm, TVector <Point3D> &U)
{
	//Recompute facility cost according to the normals for non-uniform clustering
	const double eps = 0.0001;

	for (int i = 0; i < U.size(); i++)
	{
		//ABN
		if (fnorm == &IHFL::nDIS || fnorm == &IHFL::nABN || fnorm == &IHFL::nABLP)
		{
			//Compute new facility cost
			U[i].fc = std::max(std::min(fc * fc / (RP[i].abn + eps), rat * fc), fc / rat);
		}
		
		//DFP + L2
		else if (fnorm == &IHFL::nDFP || fnorm == &IHFL::nL2)
		{
			U[i].fc = std::max(std::min(fc * fc / (RP[i].sigma + eps), rat * fc), fc / rat);
		}
	}
}


void IHFL::getAveragePointNormal(TVector <Point3D>& U, const TVector2D <size_t>& knn_id, TVector <RegressionPlane>& RP)
{
	//Compute average point normal using SVD decomposition, regression error and ABN
	const int n = U.size();

	std::cout << ">> SVD decomposition: ";
	for (int i = 0; i < n; i++)
	{
		//Add all neighbors to the list
		TVector <Point3D> KNN;
		for (size_t index : knn_id[i])
			KNN.push_back(U[index]);

		//Compute regression plane using SVD
		RegressionPlane plane;
		plane.computeRegressionPlane(KNN);

		//Add average normal vector and regression error to the list
		//ABN will be computed later
		RP.push_back(plane);
	}

	//Compute max ABN for points and their neighbors
	for (int i = 0; i < n; i++)
	{
		//Find max ABN
		const size_t index0 = knn_id[i][0];
		double abn_max = 0;
		for (size_t index : knn_id[i])
		{
			//const double abn_j = abn(AN[index0], AN[index]);
			const double abn_j = nL2(U[index0], U[index], RP[index0], RP[index]) * sin(abn(RP[index0], RP[index]) / 2);

			//Actualize maximum
			if (abn_j > abn_max)
				abn_max = abn_j;
		}

		//Store ABN max value to the list
		RP[i].abn = abn_max;
	}

	std::cout << "OK \n";
}


void IHFL::clusterizeIHFL(TVector <Point3D>&U, const double fc, TVector <Facility>&F, TVector <RegressionPlane> &RP)
{
	//Perform incremental heuristic facility clustering by (IHFL) according to a given norm
	//Proposed incremental algorithms selecting one of four strategies
	const int n = U.size();

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_id;
	TVector2D <float> knn_dist;
	TVector <int> K(1, k);
	KNNSearch search (U);
	search.findAllKNN(U, K, knn_id, knn_dist);

	//Compute average normals, regression errors
	getAveragePointNormal(U, knn_id, RP);

	//Recompute facility costs according to normals (replace old values)]
	const double multiplier = 10.0;
	if (non_uniform_cl)
	{
		//if (fnorm != &IHFL::nL2)
			recomputeFacilityCosts(fc, multiplier, RP, fnorm, U);
		//else
			//recomputeFacilityCosts(fc, multiplier * multiplier * multiplier, RP, fnorm, U);
			//recomputeFacilityCosts(lambda, fc, 5.0, RP, &IHFL::nL2, U);
	}

	//Process all points
	std::cout << ">> Clusterization: ";

	//First point becomes the facility
	Facility f0(U[0].id + 1, U[0].fc);
	F.push_back(f0);

	//Proces remaining points one by one
	for (int i = 1; i < n; i++)
	{
		//Print 
		if (i % (25000) == 0)
			std::cout << i << " ";

		//Update clusters using the IHFL method, heuristic approach
		updateClusters(i, U, F, RP);
	}
}


void IHFL::updateClusters(const int i, const TVector <Point3D>& points, TVector <Facility>& F, TVector <RegressionPlane>& RP)
{
	//Incremental method, update of the clusterization takes one of four strategies:
	//   1) Create new facility at p: S1
	//   2) Connect p to the (pseudo) nearest facility: S2
	//   3) Reallocate all clusters (pseudo) nearest to p. Create facility at p + multiple full reallocations: S3
	//   4) Reallocate parts of clusters (pseudo) nearest to p. Create facility at p + multiple partial reallocations: S4

	//Initialize index
	int j = 0, idx_nearest = -1;

	//Initial costs for different strategies
	double c_nearest = std::numeric_limits<float>::max();   //Cost of the nearest facility, strategy S1
	double c_new = points[i].fc;					//Cost for creation of the new facility, strategy S2
	double c_reallocate_clusters = points[i].fc;			//Cost for the cluster modification (full or partial reassignment to another facility), strategy S3 a) b)

	//Browse all facilities
	TVector <int> reallocate_oper;				//Changed cluster id and type of reassignment operation (false = partial, true = full)
	for (Facility& f : F)
	{
		//Reset sign and shift
		const int p_idx2 = abs(f.p_idx) - 1;

		//Norm and pseudonorm
		const double dist_pf = dist(points[i].x, points[i].y, points[i].z, points[p_idx2].x, points[p_idx2].y, points[p_idx2].z);       //Distance between point p and facility
		const double dpf = (this->*fnorm)(points[i], points[p_idx2], RP[i], RP[p_idx2]);		//Pseudonorm between point p and facility

		//Cummulated values
		double dc_all = points[i].fc - f.fc;					//Cost diifference: reassignment of all cluster to  p - cost for the old facility deletition
		double dc_closer = points[i].fc - dpf;					//Cost difference: reassignment of cluster points closer to p

		//Reallocate only according to near facilities
		if (dist_pf < 3.0 * lambda)
		{
			//Distance point and the current facility: Strategy S1
			if (dpf < c_nearest && dpf > 0)  //Actualize distance to the nearest facility
			{
				c_nearest = dpf;
				idx_nearest = j;
			}

			//Browse all points ui of the facility
			for (int& u_idx : f.U_idxs)
			{
				//Reset sign and shift of the index
				const int u_idx2 = abs(u_idx) - 1;

				//Compute pseudonorms using pointers to member functions
				const double dup = (this->*fnorm)(points[u_idx2], points[i], RP[u_idx2], RP[i]);           //Pseudonorm of the cluster point u to the proposed new center p
				const double duf = (this->*fnorm)(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);       //Pseudonorm of the cluster point u to its cluster center f.u

				//Current cluster point u is closer to p: cost for the reassignment u to the new center p
				//Strategy S4, compute new cost increment
				if (dup < duf)
				{
					dc_closer += dup - duf;

					//Change sign to plus (reallocate to p)
					u_idx = u_idx2 + 1;
				}

				//Change sign to minus: (connected to old facility)
				else
					u_idx = -u_idx2 - 1;

				//Cost for the reassignment of any cluster point to p
				//Strategy S3
				dc_all += dup - duf;
			}

			//Entire cluster is reallocated to p: does the newly created cluster at p decrese the cost?
			//Strategy S3)
			if ((dc_all < dc_closer))
			{
				c_reallocate_clusters = c_reallocate_clusters + dc_all - points[i].fc + dpf;	     //Expected cost increment

				//Set +, S3 strategy
				reallocate_oper.push_back(j + 1);
			}

			//Reallocate part of the cluster to p: is there any improvement?
			//Strategy S4)
			else if (dc_closer < dc_all)
			{
				c_reallocate_clusters = c_reallocate_clusters + dc_closer - points[i].fc + dpf;       //Expected increment

				//Set -, S4 strategy
				reallocate_oper.push_back(-j - 1);
			}
		}

		//Increment j
		j++;
	}

	//Heuristic decision: find minimum cost increment and choose the optimal strategy (S1-S3 a) b))
	const double c_min = std::min(c_new, std::min(c_nearest, c_reallocate_clusters));

	//S1: Create new facility at p is the cheapest
	if (c_min == c_new)
	{
		Facility fac_new(i + 1, points[i].fc);
		F.push_back(fac_new);
	}

	//S2: Connect p to the nearest facility
	else if (c_min == c_nearest)
	{
		F[idx_nearest].U_idxs.push_back(i + 1);
	}

	//S3 or S4
	else
	{
		//Create new facility fp at p
		Facility fac_new(i + 1, points[i].fc);

		//Process all clusters affected by the reallocation one by one
		for (const int& k : reallocate_oper)
		{
			//Reset sign and shift
			const int k2 = abs(k) - 1;

			//Partial realloation: S4
			if (k < 0)
			{
				//Reconnect all clients
				TVector <int> U_old_idx;
				for (const int & u_idx : F[k2].U_idxs)
				{
					//Reconnect to the new facility
					if (u_idx > 0)
						fac_new.U_idxs.push_back(abs(u_idx));

					//Remain to the old facility
					else
						U_old_idx.push_back(abs(u_idx));
				}

				//Connect points to the old facility
				F[k2].U_idxs = std::move(U_old_idx);
			}

			//Full reallocation: S3
			else
			{
				//Insert all clients to the new facility fp at p
				fac_new.U_idxs.insert(fac_new.U_idxs.end(), std::make_move_iterator(F[k2].U_idxs.begin()), std::make_move_iterator(F[k2].U_idxs.end()));

				//Connect old facility to the new facility fp at p
				fac_new.U_idxs.push_back(F[k2].p_idx);

				//Mark old facilty for the deletetion
				F[k2].del = true;
			}
		}

		//Add new facility at p to the list
		F.push_back(fac_new);
	}

	//Delete marked facilities of the old cluster (now it is rerallocated)
	F.erase(std::remove_if(F.begin(), F.end(), IsFacilitySetForDeletion()), F.end());
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
			total_cost += (this->*fnorm)(points[p_idx2], points[c_id2], AN[p_idx2], AN[c_id2]);
		}
	}

	return total_cost;
}

void IHFL::clusterStatistics(const TVector <Point3D>& points, const TVector <Facility>& F, const TVector <RegressionPlane>& RP, TVector <int>& NC, TVector <double> &RAD, TVector <double>& ABN, TVector <double>& DFP, TVector <double>& ASP, TVector <int>& DIM, TVector <int>& OVER)
{
	//Compute parameters of the cluster
	for (const auto f : F)
	{
		//Browse all points ui of the facility
		int i = 0;
		double f_radius = 0, f_abn = 0, f_dfp = 0, f_aspect = -1, f_dim = 0;
		Eigen::MatrixXd M(f.U_idxs.size(), 3);
		for (int u_idx : f.U_idxs)
		{
			//Reset sign and shift of the index
			const int u_idx2 = abs(u_idx) - 1;
			const int p_idx2 = abs(f.p_idx) - 1;

			//Convert clients to matrix
			M(i, 0) = points[u_idx2].x;
			M(i, 1) = points[u_idx2].y;
			M(i, 2) = points[u_idx2].z;

			//Compute pseudonorms
			f_radius = std::max(nL2(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]), f_radius);
			f_abn += nABN(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);
			f_dfp += nDFP(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);

			i++;
		}

		//PCA decomposition
		if (f.U_idxs.size() > 0)
		{
			//PCA analysis
			Eigen::MatrixXd centered = M.rowwise() - M.colwise().mean();
			Eigen::MatrixXd cov = centered.adjoint() * centered;
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

			//Eigen vectors and eigen values
			Eigen::MatrixXd eivect = eig.eigenvectors();
			Eigen::MatrixXd eival = eig.eigenvalues();

			//Compute cluster aspect
			f_aspect = (eival(2,0) == 0 ? -1 : eival(1, 0) / eival(2, 0));

			//Sum of eigen values
			const double lambda_sum = eival(0, 0) + eival(1, 0) + eival(2, 0);

			//Cluster dimension 0
			if (lambda_sum == 0)
				f_dim = 0;

			//Cluster dimension: 0 - 3
			else
			{
				//Fractions of eigen values
				const double lambda1f = eival(0, 0) / lambda_sum;
				const double lambda2f = eival(1, 0) / lambda_sum;
				const double lambda3f = eival(2, 0) / lambda_sum;

				//Set cluster dimensions
				if (lambda3f < 0.01)
					f_dim = 0;
				else if (lambda2f < 0.01)
					f_dim = 1;
				else if (lambda1f < 0.01)
					f_dim = 2;
				else
					f_dim = 3;
			}
		}

		//Store values
		NC.push_back(f.U_idxs.size());
		RAD.push_back(f_radius);
		ABN.push_back(f_abn / std::max(1, (int)f.U_idxs.size()));
		DFP.push_back(f_dfp / std::max(1, (int)f.U_idxs.size()));
		ASP.push_back(f_aspect);
		DIM.push_back(f_dim);

		/*
		Eigen::MatrixXd mat(points.size(), 3);

		for (int i = 0; i < points.size(); i++)
		{
			mat(i, 0) = points[i].x;
			mat(i, 1) = points[i].y;
			mat(i, 2) = points[i].z;
		}
		*/
		/*
		double l = 1, w = 2, h = 3;
		Eigen::MatrixXd mat(1, 3);
		mat(0, 0) = 0;
		/*
		mat(0, 1) = 0;
		mat(0, 2) = 0;

		mat(1, 0) = l;
		mat(1, 1) = 0;
		mat(1, 2) = 0;

		mat(2, 0) = 0;
		mat(2, 1) = w;
		mat(2, 2) = 0;

		mat(3, 0) = l;
		mat(3, 1) = w;
		mat(3, 2) = 0;

		mat(4, 0) = 0;
		mat(4, 1) = 0;
		mat(4, 2) = h;

		mat(5, 0) = l;
		mat(5, 1) = 0;
		mat(5, 2) = h;

		mat(6, 0) = 0;
		mat(6, 1) = w;
		mat(6, 2) = h;

		mat(7, 0) = l;
		mat(7, 1) = w;
		mat(7, 2) = h;
		

		Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
		Eigen::MatrixXd cov = centered.adjoint() * centered;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
		Eigen::MatrixXd eivect = eig.eigenvectors();
		Eigen::MatrixXd eival = eig.eigenvalues();

		std::cout << eivect;
		std::cout << eival;
		*/
	}

	//Compute overlap ratio
	OVER.resize(F.size());
	std::fill(OVER.begin(), OVER.end(), 0);

	//Process facilities one by one
	for (int i = 0; i < F.size(); i++)
	{
		//Reset sign and shift of the index
		const int p_idx2 = abs(F[i].p_idx) - 1;

		//Check all facilities
		for (int j = 0; j < F.size(); j++)
		{
			if (i != j)
			{
				//Take all connected points
				for (int k = 0; k < F[j].U_idxs.size(); k++)
				{
					const int u_idx2 = abs(F[j].U_idxs[k]) - 1;

					//Increment overlap ratio
					const double d = nL2(points[u_idx2], points[p_idx2], RP[u_idx2], RP[p_idx2]);

					if (d < RAD[i])
						OVER[i] = OVER[i] + 1;

				}
			}
		}
	}
	
}
