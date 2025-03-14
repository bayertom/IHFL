// Description: Pointer to the norm function 

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

#ifndef Ppf_clust_metric_H
#define Ppf_clust_metric_H

#include "point3d.h"
#include "regressionplane.h"

//Forward declaration
class IHFL;

//Pointer to the (pseudo) metric function used for the clusterization
typedef double (IHFL:: * pfClustMetric) (const Point3D&, const Point3D&, const RegressionPlane&, const RegressionPlane&);

#endif