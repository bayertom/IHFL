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
