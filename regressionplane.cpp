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
#include <iostream>

void RegressionPlane::computeRegressionPlane(const Eigen::MatrixXd& A, const Eigen::MatrixXd &S, Eigen::MatrixXd &V)
{
	//Compute regression plane using PCA

	//Compute points centroid (mean of their coordinates represented by columns)
	Eigen::RowVectorXd M = A.colwise().mean();

	//Store eigen values: descendent sort
	lambda1 = S(0, 0);
	lambda2 = S(1, 0);
	lambda3 = S(2, 0);

	//Rotate the input set parallel with the coordinate axes
	Eigen::MatrixXd AT = A * V;

	//Compute extreme values of rotated coordinates
	Eigen::MatrixXd MAX = AT.colwise().maxCoeff();
	Eigen::MatrixXd MIN = AT.colwise().minCoeff();

	//Compute height of the min-max box (difference of Z-coordinates)
	height = MAX(0, 2) - MIN(0, 2);
	
	//Compute parameters a, b, c of the regression plane (last column of V)
	a = V(0, 2);
	b = V(1, 2);
	c = V(2, 2);

	//Compute parameter d of the regression plane passing through M
	d = -(a * M(0, 0) + b * M(0, 1) + c * M(0, 2));
}