// Description: KNN search using the nanoflann library

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


#include "knnsearch.h"

#include <iostream>

#include "nanoflann.hpp"
#include "utils.h"

using namespace nanoflann;


void KNNSearch::findAllKNN(const TVector <Point3D>&qpoints, const TVector<int>&K, TVector2D <size_t>&knn_id, TVector2D <float>&knn_dist)
{
	//Find all k nearest neighbors to any qpoint in points
	const int n = points.size(), nq = qpoints.size(), nk = K.size();

	std::cout << ">> Convert points to the cloud: ";

	//Convert points to the cloud
	PointCloud <double> pointsc;
	pointsc.pts.resize(n);
	for (int i = 0; i < n; i++)
	{
		pointsc.pts[i].x = points[i].x;
		pointsc.pts[i].y = points[i].y;
		pointsc.pts[i].z = points[i].z;
	}

	std::cout << "OK \n";

	//Create kd tree
	std::cout << ">> Building KD-tree: ";
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3> my_kd_tree_t;
	my_kd_tree_t  index(3, pointsc, KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();
	std::cout << "OK \n";

	//Perform knn search
	std::cout << ">> Perform knn-search: ";
	for (int i = 0; i < nq; i++)
	{
		//Fixed or variable amount of neighbors
		const int k = std::max((nk == 1 ? K[0] : K[i]), 1);

		//Auxilliary data structures
		TVector<size_t> ret_indices(k);
		TVector<float> out_dists_sqr(k);

		//Initialize result set
		nanoflann::KNNResultSet<float> resultSet(k);
		resultSet.init(&ret_indices[0], &out_dists_sqr[0]);

		//Perform knn search
		double query_pt[3] = { qpoints[i].x, qpoints[i].y, qpoints[i].z };
		index.findNeighbors(resultSet, &query_pt[0]);

		//Add to the list
		knn_id.push_back(ret_indices);
		knn_dist.push_back(out_dists_sqr);
	}

	std::cout << "OK \n";
}