// Description: Mark faciity for the deletiotion (S3 strategy)

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


#ifndef IsFacilitySetForDeletion_H
#define IsFacilitySetForDeletion_H

#include "facility.h"

//Mark facility for deletion
class IsFacilitySetForDeletion
{
	public:
		bool operator() (const Facility & f) const
		{
			return f.del;
		}
};


#endif
