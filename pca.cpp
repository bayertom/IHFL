// Description: Principal component analysis

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

#include "pca.h"
#include <iostream>

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> PCA::pca(const Eigen::MatrixXd& A)
{
	//Compute PCA using SVD

	//Compute means of coordinates over columns
	Eigen::RowVectorXd M = A.colwise().mean();

	//Subtract mean: B = A - M
	Eigen::MatrixXd B = A.rowwise() - M;

	//Covariance matrix: C = B' * B / (m - 1)
	Eigen::MatrixXd C = (B.adjoint() * B) / double(A.rows() - 1);

	//Compute SVD, full version: [U, S, V] = svd(C)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullV | Eigen::ComputeFullU);

	//Get matrices
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd S = svd.singularValues();
	Eigen::MatrixXd V = svd.matrixV();

	return {U, S, V};
}