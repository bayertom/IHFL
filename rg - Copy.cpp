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

#include "knnsearch.h"
#include "pca.h"
#include "regressionplane.h"
#include "dxfexport.h"
#include "facility.h"

#include "Eigen/Dense"                               
#include "Eigen/Sparse"                              
#include "Eigen/Core"

void RG::regionGrow(TVector <Point3D>& points, const pfEdgeMetric& pf_edge_metric, const TVector<double>& edge_thresholds, const int k, TVector2D <Point3D>& regions_points)
{
	//Region-growing strategy
	const int n = points.size();

	//Auxilliary data structures
	TVector<bool> visited(n, false);
	TVector<bool> open(n, false);
	TVector2D <int> regions;

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

	//Find edge points
	TVector <bool> edge_points(n, false);
	TVector2D <double> edge_values;
	findEdgePointsPCA(points, knn_idxs, edge_thresholds[0], edge_points, edge_values);
	/*
	TVector <bool> edge_points(n, false);
	std::priority_queue<std::pair<double, int>> edge_points_pq;
	TVector2D <double> edge_values;
	findEdgePointsPCA2(points, knn_idxs, edge_thresholds[0], edge_points, edge_points_pq, edge_values);
	//sortAAccordingToB(points, edge_values, 0);
	//sortAAccordingToB(edge_points, edge_values, 0);
	//sortAAccordingToB(knn_idxs, edge_values, 0);
	*/

	//Process the point cloud
	for (int i = 0; i < n; i++)
	//while (!edge_points_pq.empty()) 
	{
		//int i = edge_points_pq.top().second;
		//edge_points_pq.pop();

		//Proces unvisited point
		if (!visited[i])
		{
			//Initialize region
			TVector <int> region;

			//Create query
			std::queue<int> Q;

			//Add first point
			Q.push(i);

			//Process, until Q is empty
			while (!Q.empty())
			{
				//Get and pop item
				int idx = Q.front();
				Q.pop();


				Facility f(idx + 1, 1);

				//Point has not been visited
				if (!visited[idx])
				{
					//Set as visited
					visited[idx] = true;
					/*
					if (region.size() > 0)
					{
						//Test distance
						double dx = points[idx].x - points[region.back()].x;
						double dy = points[idx].y - points[region.back()].y;
						double dz = points[idx].z - points[region.back()].z;

						double d = sqrt(dx * dx + dy * dy + dz * dz);
						std::cout << points[idx].x << ' ' << points[idx].y << ' ' << points[idx].z << ' ' << points[region.back()].x << ' ' << points[region.back()].y << ' ' << points[region.back()].z << ' ' << d << '\n';
						if (d > 0.2)
							std::cout << "\nError\n";

					}
					*/
					//Add to the region
					region.push_back(idx);

					//Process all nearest neighbors
					for (size_t k_idx : knn_idxs[idx])
					{
						//Point not visited and in edge points
						if ((!visited[k_idx]) && (edge_points[k_idx]) /* && (!open[k_idx])*/)
						{
							
							//Add point to Q
							if (pf_edge_metric(edge_values[idx], edge_values[k_idx], edge_thresholds))
							{
								f.U_idxs.push_back(k_idx + 1);

								open[k_idx] = true;
								Q.push(k_idx);
								//std::cout  << points[k_idx].x << ' ' << points[k_idx].y << ' ' << points[k_idx].z  << '\n';
								/*
								double dx = points[k_idx].x - points[idx].x;
								double dy = points[k_idx].y - points[idx].y;
								double dz = points[k_idx].z - points[idx].z;

								double d = sqrt(dx * dx + dy * dy + dz * dz);
								std::cout << points[idx].x << ' ' << points[idx].y << ' ' << points[idx].z<< ' ' << points[k_idx].x << ' ' << points[k_idx].y << ' ' << points[k_idx].z << ' ' << d << '\n';
								if (d > 0.2)
									std::cout << "\nError\n";
									*/
							}
						}
					}
					double xxx = 37;
					if (f.U_idxs.size() > 0)
					{
						facilities.push_back(f);
					}
				}
			}
			//Ignore small regins containing less than 5 points
			if (region.size() > 10)
			{
				regions.push_back(region);
				break;
			}
		}		
	}

	//Convert region indices to points
	for (int i = 0; i < regions.size(); i++)
	{
		//Create list of points representing a region
		TVector <Point3D> region_points;
		for (int j = 0; j < regions[i].size(); j++)
			region_points.push_back(points[regions[i][j]]);

		//Add region to the list of regions
		regions_points.push_back(region_points);
	}

	//Export facilities
	TVector <RegressionPlane> planes(n);
	DXFExport::exportClustersToDXF("nn_clusters_edges.dxf", facilities, points, planes);
}

void RG::findEdgePointsPCA(TVector <Point3D>& points, const TVector2D <size_t>& knn_idxs, const double threshold, TVector <bool>& edge_points, TVector2D <double> & edge_values)
{
	//Find edge points using PCA: l3/(l1 + l2 + l3) > threshold
	const double n = points.size();

	std::cout << ">> SVD decomposition: ";
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

		//Add to the list: curvature + normal vector
		edge_values.push_back({curv, V(0, 2), V(1, 2), V(2, 2) });

		//std::cout << curv << '\n';

		//Add point to the edge points
		if (curv > threshold)
			edge_points[i] = true;
	}

	//Sort according to the first column in ascending order
	
	//std::sort(edge_values.begin(), edge_values.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
	//	return a[0] < b[0]; });

	std::cout << "OK \n";
}


void RG::findEdgePointsPCA2(TVector <Point3D>& points, const TVector2D <size_t>& knn_idxs, const double threshold, TVector <bool>& edge_points, std::priority_queue<std::pair<double, int>> &edge_points_pq, TVector2D <double>& edge_values)
{
	//Find edge points using PCA: l3/(l1 + l2 + l3) > threshold
	const double n = points.size();

	std::cout << ">> SVD decomposition: ";
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

		//Add to the list: curvature + normal vector
		edge_values.push_back({ curv, V(0, 2), V(1, 2), V(2, 2) });

		//std::cout << curv << '\n';

		//Add point to the edge points
		if (curv > threshold)
		{
			edge_points_pq.push({ curv, i });
			edge_points[i] = true;
		}
	}

	//Sort according to the first column in ascending order

	//std::sort(edge_values.begin(), edge_values.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
	//	return a[0] < b[0]; });

	std::cout << "OK \n";
}



bool RG::curvABNEdge(const TVector <double>& p1_vals, const TVector <double>& p2_vals, const TVector <double>& edge_thresholds)
{
	//Returns, if two points belong to the same edge according to:
	//   a) curvature differences
	//   b) abn criteria differences

	//Compute curvature
	const double dcurv = fabs(p1_vals[0] - p2_vals[0]);

	//Compute ABN
	const double a1 = p1_vals[1], b1 = p1_vals[2], c1 = p1_vals[3];
	const double a2 = p2_vals[1], b2 = p2_vals[2], c2 = p2_vals[3];
	
	const double norma = sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	const double normb = sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	const double dot_ab = fabs(a1 * a2 + b1 * b2 + c1 * c2);
	const double abn = acos(std::min(1.0, dot_ab / (norma * normb)));
	//std::cout << abn << ' ';

	//Compare thresholds
	bool res = dcurv < edge_thresholds[1] && abn < edge_thresholds[2];
	if (res)
		std::cout << abn << '\n';
	return res;
}

template <typename T1, typename T2>
void RG::sortAAccordingToB(TVector<T1>& a, TVector2D<T2>& b, const unsigned int k) 
{
	//Sort vector a according to the k-th column of b, written for a generic type
	const int n = b.size();

	//Create an index vector {0, 1, 2, ...}
	TVector<int> idxs(n);
	std::iota(idxs.begin(), idxs.end(), 0);

	//Sort indices based on the k-t column of b
	std::sort(idxs.begin(), idxs.end(), [&](int i, int j) {
		return b[i][k] > b[j][k];
		});

	//Apply the same sorting to a and b
	TVector<T1> a_s(n);
	TVector2D<T2> b_s(n);

	for (size_t i = 0; i < n; ++i) 
	{
		a_s[i] = a[idxs[i]];
		b_s[i] = b[idxs[i]];
	}

	//Move sorted data back to the original vectors
	a = std::move(a_s);
	b = std::move(b_s);
}
