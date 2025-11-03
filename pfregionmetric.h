// Description: Pointer to the metric (static function) measuring if two points belong to the same region
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

#ifndef PFRegionMetric_H
#define PFRegionMetric_H

#include "tvector.h"

//Forward declaration
class RG;

//Pointer to the (pseudo) metric represented by a static function measuring if two points belong to the same region (edge)
typedef bool (*pfRegionMetric) (const TVector <double> &, const TVector <double>&, const TVector <double>&);

#endif