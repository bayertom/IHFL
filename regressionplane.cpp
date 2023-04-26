// Description: Compute the regression plane using SVD

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


#include "regressionplane.h"

#include "Eigen/Core"
#include "Eigen/SVD"

void RegressionPlane::computeRegressionPlane(const TVector <Point3D>& P)
{
	//Compute regression plane using SVD
	const int n = P.size();

	//Not enough points for the regression plane
	if (n < 3)
		return;

	//Create matrix of points
	Eigen::MatrixXd X(n, 3);

	//Convert points to matrix
	for (int i = 0; i < n; i++)
	{
		X(i, 0) = P[i].x;
		X(i, 1) = P[i].y;
		X(i, 2) = P[i].z;
	}

	//Compute centroid
	Eigen::RowVector3d C = Eigen::RowVector3d(X.col(0).mean(), X.col(1).mean(), X.col(2).mean()).transpose();

	//Subtract centroid
	Eigen::MatrixXd CX = X.rowwise() - C;

	//Compute SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CX, Eigen::ComputeFullV | Eigen::ComputeFullU/*Eigen::ComputeThinU | Eigen::ComputeThinV*/);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd S = svd.singularValues();
	Eigen::MatrixXd V = svd.matrixV();

	//Compute sigma: S[min] cummulated distance of points from the plane
	sigma = S(2, 0);

	//Compute variability of normals S[min] / (S[min] + S[mid] + S[max])
	var = S(2, 0) / (S(0, 0) + S(1, 0) + S(2, 0));

	//Compute parameters a, b, c of the regression plane (last column of V)
	a = V(0, 2);
	b = V(1, 2);
	c = V(2, 2);

	//Compute parameter d of the regression plane
	d = -(a * C(0, 0) + b * C(0, 1) + c * C(0, 2));
}

