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


#ifndef RegressionPlane_H
#define RegressionPlane_H

#include "Eigen/Core"
#include "Eigen/SVD"

#include "tvector.h"
#include "point3d.h"

//Parameters of the regression plane given by P
struct RegressionPlane
{
	double  a,	 //Plane equation, a
		b,	 //Plane equation, b
		c,	 //Plane equation, c
		d,	 //Plane equation, d
		height,	 //Height of the PCA min-max box over the points
		lambda1, //Eigenvalue 1
		lambda2, //Eigenvalue 2
		lambda3; //Eigenvalue 3

	RegressionPlane() : a(1.0), b(0.0), c(0.0), d(0.0), height(0.0), lambda1(1.0), lambda2(0.0), lambda3(0.0) {}

	void computeRegressionPlane(const Eigen::MatrixXd& A, const Eigen::MatrixXd& S, Eigen::MatrixXd& V);
};

#endif