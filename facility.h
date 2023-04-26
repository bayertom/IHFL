// Description: Facility location clustering

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


#ifndef Facility_H
#define Facility_H

//Facility definition
struct Facility {
	int p_idx;				//Facility center p, index of point
	double fc;				//Facility cost, non-uniform clusterization
	TVector <int> U_idxs;			//Connected clients ui, indexes of points
	bool del;				//Facility marked as deleted
	Facility(const unsigned int& p_idx_, const double& fc_) : p_idx(p_idx_), fc(fc_), del(false) {}
};

#endif