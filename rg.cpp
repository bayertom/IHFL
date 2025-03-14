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


#include "rg.h"

#include <queue>
#include <iostream>
#include <tuple>
#include <numbers>

#include "knnsearch.h"
#include "pca.h"
#include "regressionplane.h"
#include "dxfexport.h"
#include "facility.h"

#include "Eigen/Dense"                               
#include "Eigen/Sparse"                              
#include "Eigen/Core"

void RG::regionGrow(TVector <Point3D>& points, const pfRegionMetric& pf_region_metric, const TVector<double>& edge_thresholds, const int k, TVector2D <Point3D>& regions_points)
{
	//Region-growing strategy according to the computed values of point features
	const int n = points.size();

	//Auxilliary data structures
	TVector<bool> visited(n, false);
	TVector<bool> open(n, false);
	TVector <Region> regions;

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_idxs;
	TVector2D <float> knn_dists;
	TVector <int> K(1, k);
	KNNSearch search(points);
	search.findAllKNN(points, K, knn_idxs, knn_dists);

	//Convert to facilities
	TVector <Facility> facilities;
	/*
	Facility f(1, 1);
	TVector <RegressionPlane> planes(n);
	for (int i = 0; i < n; i++)
	{
		Facility f(i+1, 1);
		for (int j = 0; j < k; j++)
		{
			f.U_idxs.push_back(knn_idxs[i][j]+1);
		}

		facilities.push_back(f);
	}

	DXFExport::exportClustersToDXF("nn_clusters.dxf", facilities, points, planes);
	*/
	
	//Compute point features using PCA
	TVector2D <double> point_features =  computeFeatures(points, knn_idxs);

	//Find edge points (use a curvature)
	TVector <bool> edge_points = findEdgePointsCurv(point_features, edge_thresholds[0]);

	//Process the points of the cloud one by one
	for (int i = 0; i < n; i++)
	{
		//Proces unvisited point
		if (!visited[i])
		{
			//Initialize new region
			Region new_region;

			//Create queue
			std::queue<int> Q;

			//Add first point
			Q.push(i);

			//Process, until Q is empty
			while (!Q.empty())
			{
				//Get and pop item
				int idx = Q.front();
				Q.pop();

				// Create facility
				// Facility f(idx + 1, 1);

				//Point has not been visited
				if (!visited[idx])
				{
					//Set as visited
					visited[idx] = true;
					
					//Add point to the region and update its properties
					new_region.add(idx, point_features[idx][0], point_features[idx][1], point_features[idx][2], point_features[idx][3], point_features[idx][4], point_features[idx][5], point_features[idx][6]);

					//Process all nearest neighbors
					for (size_t k_idx : knn_idxs[idx])
					{
						//Point not visited and in edge points
						if ((!visited[k_idx]) && (edge_points[k_idx]) /* && (!open[k_idx])*/)
						{
							//Add adjacent point given by the index k_idx to Q, has same properties as current region
							TVector new_reg_values = { new_region.curv_aver, new_region.tx_aver, new_region.ty_aver, new_region.tz_aver, new_region.nx_aver, new_region.ny_aver, new_region.nz_aver };
							if (pf_region_metric(new_reg_values, point_features[k_idx], edge_thresholds))
							{
								// Add nodes to the facility
								// f.U_idxs.push_back(k_idx + 1);

								open[k_idx] = true;
								Q.push(k_idx);
							}
						}
					}
					
					//Add facilities to the list
					/*
						if (f.U_idxs.size() > 0)
					{
						facilities.push_back(f);
					}
					*/
				}
			}

			//Ignore small regions containing less than 10 points
			if (new_region.points_idx.size() > 10)
			{
				
				//Find regions that need to be merged
				TVector new_reg_values = { new_region.curv_aver, new_region.tx_aver, new_region.ty_aver, new_region.tz_aver, new_region.nx_aver, new_region.ny_aver, new_region.nz_aver };
				TVector <int> mergeable_regions;
				for (int j = 0; j < regions.size(); j++)
				{
					//Get features of the j-th region
					TVector regj_values = { regions[j].curv_aver, regions[j].tx_aver, regions[j].ty_aver, regions[j].tz_aver, regions[j].nx_aver, regions[j].ny_aver, regions[j].nz_aver };

					//New region and j-th region could be merged
					if (pf_region_metric(new_reg_values, regj_values, edge_thresholds))
					{
						double hd =  getHausdorffDist(new_region, regions[j], points);
						std::cout << hd << ' ';
						if (hd < 0.1)
							mergeable_regions.push_back(j);
					}
				}
				
				//Some mergeable regions were found
				if (mergeable_regions.size() > 0)
				{
					//Create new region
					Region merged_region = new_region;

					for (int j = mergeable_regions.size() - 1; j >= 0; --j) 
					{
						//Merge two regions
						merged_region = mergeRegions(merged_region, regions[mergeable_regions[j]]);
						
						//Erase old region
						regions.erase(regions.begin() + mergeable_regions[j]);
					}

					//Add merged region to the list
					regions.push_back(merged_region);

					break;
				}
				
				//No regions to merge, add new region
				else
					regions.push_back(new_region);
			}
		}		
	}

	//Convert region indices to points
	for (int i = 0; i < regions.size(); i++)
	{
		//Create list of points representing a region
		TVector <Point3D> region_points;
		for (int j = 0; j < regions[i].points_idx.size(); j++)
			region_points.push_back(points[regions[i].points_idx[j]]);

		//Add region to the list of regions
		regions_points.push_back(region_points);
	}

	//Export facilities
	//TVector <RegressionPlane> planes(n);
	//DXFExport::exportClustersToDXF("nn_clusters_edges.dxf", facilities, points, planes);
}


void RG::regionGrowDyn(TVector <Point3D>& points, const pfRegionMetric& pf_region_metric, const TVector<double>& edge_thresholds, const int k, TVector2D <Point3D>& regions_points)
{
	//Region-growing strategy according to the computed values of point features
	const int n = points.size();

	//Auxilliary data structures
	TVector<bool> visited(n, false);
	TVector<bool> open(n, false);
	TVector <Region> regions;

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_idxs;
	TVector2D <float> knn_dists;
	TVector <int> K(1, k);
	KNNSearch search(points);
	search.findAllKNN(points, K, knn_idxs, knn_dists);

	//Compute point features using PCA
	TVector2D <double> point_features = computeFeatures(points, knn_idxs);

	//Find edge points (use a curvature)
	TVector <bool> edge_points = findEdgePointsCurv(point_features, edge_thresholds[0]);

	//Process the points of the cloud one by one
	std::unordered_map<int, int> point_to_region;
	for (int i = 0; i < n; i++)
	{
		//Proces unvisited point
		if (!visited[i])
		{
			//Create queue
			std::queue<int> Q;

			//Add first point
			Q.push(i);

			//Set point as visited
			visited[i] = true;

			//Create new region
			Region new_region;

			//Add point to the region and update its properties
			new_region.add(i, point_features[i][0], point_features[i][1], point_features[i][2], point_features[i][3], point_features[i][4], point_features[i][5], point_features[i][6]);

			//Process, until Q is empty
			while (!Q.empty())
			{
				//Get and pop item
				int idx = Q.front();
				Q.pop();

				//Process all nearest neighbors
				for (size_t k_idx : knn_idxs[idx])
				{
					//Point not visited and in edge points
					if ((!visited[k_idx]) && (edge_points[k_idx]) /* && (!open[k_idx])*/)
					{
						//Add adjacent point given by the index k_idx to Q
						//It has same properties as current region
						TVector new_reg_values = { new_region.curv_aver, new_region.tx_aver, new_region.ty_aver, new_region.tz_aver, new_region.nx_aver, new_region.ny_aver, new_region.nz_aver };
						if (pf_region_metric(new_reg_values, point_features[k_idx], edge_thresholds))
						//if (pf_region_metric(point_features[idx], point_features[k_idx], edge_thresholds))
						{
							//Set as visited
							visited[k_idx] = true;

							//Add to the queue
							Q.push(k_idx);

							//Check if the neighbor already belongs to an existing region
							auto it_region = point_to_region.find(k_idx);

							//Point belongs to another region
							if (it_region != point_to_region.end())
							{
								//Get region index
								int region_idx = point_to_region[k_idx];

								//Merge both regions
								new_region = mergeRegions(new_region, regions[region_idx]);
							}

							//Add point to the new region
							else
							{
								//Add point to the new region
								new_region.add(k_idx, point_features[k_idx][0], point_features[k_idx][1], point_features[k_idx][2], point_features[k_idx][3], point_features[k_idx][4], point_features[k_idx][5], point_features[k_idx][6]);

								//Update index
								point_to_region[k_idx] = regions.size();
							}
						}
					}
				}
			}

			//Ignore small regions containing less than 10 points
			if (new_region.points_idx.size() > 10)
			{
				regions.push_back(new_region);
			
			}
		}				
	}

	//Convert region indices to points
	for (int i = 0; i < regions.size(); i++)
	{
		//Create list of points representing a region
		TVector <Point3D> region_points;
		for (int j = 0; j < regions[i].points_idx.size(); j++)
			region_points.push_back(points[regions[i].points_idx[j]]);

		//Add region to the list of regions
		regions_points.push_back(region_points);
	}

	//Export facilities
	//TVector <RegressionPlane> planes(n);
	//DXFExport::exportClustersToDXF("nn_clusters_edges.dxf", facilities, points, planes);
}


Region RG::mergeRegions(Region& r1, const Region& r2)
{
	//Merge two regions into new one
	Region merged = r1;

	//Compute new average values
	merged.curv_aver = 0.5 * (r1.curv_aver + r2.curv_aver);
	merged.tx_aver = 0.5 * (r1.tx_aver + r2.tx_aver);
	merged.ty_aver = 0.5 * (r1.ty_aver + r2.ty_aver);
	merged.tz_aver = 0.5 * (r1.tz_aver + r2.tz_aver);
	merged.nx_aver = 0.5 * (r1.nx_aver + r2.nx_aver);
	merged.ny_aver = 0.5 * (r1.ny_aver + r2.ny_aver);
	merged.nz_aver = 0.5 * (r1.nz_aver + r2.nz_aver);

	//Copy items
	copy(r2.points_idx.begin(), r2.points_idx.end(), back_inserter(merged.points_idx));

	return merged;
}


TVector2D <double> RG::computeFeatures(const TVector <Point3D>& points, const TVector2D <size_t>& knn_idxs)
{
	//Compute features using PCA at given point
	//Curvature, (tx, ty, tz), (nx, ny, nz)
	const double n = points.size();
	TVector2D <double> point_features;

	//Process points one by one
	for (int i = 0; i < n; i++)
	{
		//Convert nearest neighbors to matrix A
		const int m = knn_idxs[i].size();
		Eigen::MatrixXd A(m, 3);

		for (int j = 0; j < m; j++)
		{
			const int knn_idx = knn_idxs[i][j];
			A(j, 0) = points[knn_idx].x;
			A(j, 1) = points[knn_idx].y;
			A(j, 2) = points[knn_idx].z;
		}

		//Compute PCA
		auto [U, S, V] = PCA::pca(A);

		//Get eigen-values: descendent sort
		const double lambda1 = S(0, 0);
		const double lambda2 = S(1, 0);
		const double lambda3 = S(2, 0);

		//Compute curvature
		const double curv = lambda3 / (lambda1 + lambda2 + lambda3);

		//Add principal components to the list
		point_features.push_back({curv, V(0, 0), V(1, 0), V(2, 0), V(0, 2), V(1, 2), V(2, 2) });
	}

	return point_features;
}


TVector <bool> RG::findEdgePointsCurv(const TVector2D <double>& point_features, const double threshold)
{
	//Detect edge points using the curvature criteria
	TVector <bool> edge_points(point_features.size(), false);

	//Process points one by one
	for (int i = 0; i < point_features.size(); i++)
	{
		//Point exceeds the threshold -> edge point
		if (point_features[i][0] > threshold)
			edge_points[i] = true;
	}

	return edge_points;
}


bool RG::inSameRegion(const TVector <double>& p1_vals, const TVector <double>& p2_vals, const TVector <double>& edge_thresholds)
{
	//Returns, if two points belong to the same region according to:
	//   a) curvature differences
	//   b) direction vectors (abt) differences
	//   b) normal vectors (abt) differences

	auto [dcurv, abt, abn] = measureSimilarity(p1_vals[0], p1_vals[1], p1_vals[2], p1_vals[3], p1_vals[4], p1_vals[5], p1_vals[6],
		                               p2_vals[0], p2_vals[1], p2_vals[2], p2_vals[3], p2_vals[4], p2_vals[5], p2_vals[6]);

	//Compare thresholds
	bool res = /*dcurv < edge_thresholds[1] &&*/ abt < edge_thresholds[2] /* && abn < edge_thresholds[3]*/;
	return res;
}


std::tuple<double, double, double> RG::measureSimilarity(const double curv1, const double tx1, const double ty1, const double tz1, const double nx1, const double ny1, const double nz1, const double curv2, const double tx2, const double ty2, const double tz2, const double nx2, const double ny2, const double nz2)
{
	//Measure similarity between two points/regions
	
	//Compute curvature differences
	const double dcurv = fabs(curv1 - curv2);

	//Compute angle between tangents
	const double abt = get2VectorsAngle(tx1, ty1, tz1, tx2, ty2, tz2);
	const double abt2 = std::min(abt, fabs(abt - std::numbers::pi));

	//Compute angle between normals
	const double abn = get2VectorsAngle(nx1, ny1, nz1, nx2, ny2, nz2);
	const double abn2 = std::min(abn, fabs(abn - std::numbers::pi));

	return { dcurv, abt2, abn2 };
}


double RG::get2VectorsAngle(const double ux, const double uy, const double uz, const double vx, const double vy, const double vz)
{
	//Get angle between two vectors
	const double nu = sqrt(ux * ux + uy * uy + uz * uz);
	const double nv = sqrt(vx * vx + vy * vy + vz * vz);
	const double dot = fabs(ux * vx + uy * vy + uz * vz);

	return acos(std::max(-1.0, std::min(1.0, dot / (nu * nv))));
}


double RG::getHausdorffDist(const Region& r1, const Region& r2, const TVector <Point3D>& points)
{
	//Hausdorff distance from r1 to r2
	const double dmax1 = hausdorffDist(r1.points_idx, r2.points_idx, points);

	//Hausdorff distance from r2 to r1
	const double dmax2 = hausdorffDist(r2.points_idx, r1.points_idx, points);

	//Return maximum of both distances
	return std::max(dmax1, dmax2);
}


double RG::hausdorffDist(const TVector <int>& idxs1, const TVector <int>& idxs2, const TVector <Point3D>& points)
{
	//Compute Hausdorff distance from r1 to r2
	double dmax = 0;
	int ttt = idxs1.size();
	int uuu = idxs2.size();

	for (unsigned int i = 0; i < idxs1.size(); i++)
	{
		double dmin = DBL_MAX;
		for (unsigned int j = 0; j < idxs2.size(); j++)
		{
			//Compute distance
			const double dx = points[idxs1[i]].x - points[idxs2[j]].x;
			const double dy = points[idxs1[i]].y - points[idxs2[j]].y;
			const double dz = points[idxs1[i]].z - points[idxs2[j]].z;
			const double d = sqrt(dx * dx + dy * dy + dz * dz);

			//if (ttt == 109)
			//	std::cout << points[idxs1[i]].x << ' ' << points[idxs1[i]].y <<  ' ' << points[idxs2[j]].x << ' ' <<points[idxs2[j]].y << ' ' << d << '\n';

			//Update minimum
			if (d < dmin)
			{
				dmin = d;

				//if (ttt == 109 && uuu == 154)
				//	std::cout << points[idxs1[i]].x << ' ' << points[idxs1[i]].y << ' ' <<  points[idxs1[i]].z << ' ' << points[idxs2[j]].x << ' ' << points[idxs2[j]].y << points[idxs2[j]].z << ' ' << ' ' << d << '\n';

			}
		}

		//Update dmax
		if (dmin > dmax)
			dmax = dmin;
	}

	return dmax;
}


