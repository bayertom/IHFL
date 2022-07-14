
#ifndef PFNorm_H
#define PFNorm_H

#include "point3d.h"
#include "regressionplane.h"

//Forward ddeclaration
class IHFL;

//Pointer to the (pseudo) norm function
typedef double (IHFL:: * pfnorm) (const Point3D&, const Point3D&, const RegressionPlane&, const RegressionPlane&);

#endif